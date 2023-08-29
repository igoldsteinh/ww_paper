# Turing model for SEIRR-ww model
using LinearAlgebra
using ForwardDiff
prob = ODEProblem(seirr_ode_log!,
zeros(6),
(0.0, obstimes[end]),
ones(4))


@model function bayes_seirr_student(data_log_copies, obstimes, param_change_times, extra_ode_precision, prob, abs_tol, rel_tol)
    # Calculate number of observed datapoints timepoints
    l_copies = length(obstimes)
    l_param_change_times = length(param_change_times)
  
    # Priors
    R0_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), Diagonal(ones(l_param_change_times + 2))) # +2, 1 for var, 1 for init
    S_SEIR1_non_centered ~ Normal()
    I_EIR1_non_centered ~ Normal()
    R1_ER1_non_centered ~ Normal()
    gamma_non_centered ~ Normal() # rate to I
    nu_non_centered ~ Normal() # rate to R1
    eta_non_centered ~ Normal() # rate to R2
    rho_gene_non_centered ~ Normal() # gene detection rate
    τ_non_centered ~ Normal() # standard deviation for log scale data
    lambda_non_centered ~ Normal() # percentage of emissions from the I vs R1 compartment
    df ~ Gamma(df_shape, df_scale)
  
    # Transformations
    gamma = exp(gamma_non_centered * gamma_sd + gamma_mean)
    nu = exp(nu_non_centered * nu_sd + nu_mean)
    eta = exp(eta_non_centered * eta_sd + eta_mean)
    rho_gene = exp(rho_gene_non_centered * rho_gene_sd + rho_gene_mean)
    τ = exp(τ_non_centered * tau_sd + tau_mean)
    lambda = logistic(lambda_non_centered * lambda_sd + lambda_mean)
    sigma_R0_non_centered = R0_params_non_centered[1]
    sigma_R0 = exp(sigma_R0_non_centered * sigma_R0_sd + sigma_R0_mean)
    r0_init_non_centered = R0_params_non_centered[2]
    r0_init = exp(r0_init_non_centered * r0_init_sd + r0_init_mean)
    beta_init = r0_init * nu
    log_R0_steps_non_centered = R0_params_non_centered[3:end]
    S_SEIR1 = logistic(S_SEIR1_non_centered * S_SEIR1_sd + S_SEIR1_mean)
    I_EIR1 = logistic(I_EIR1_non_centered * I_EIR1_sd + I_EIR1_mean)
    R1_ER1 = logistic(R1_ER1_non_centered * R1_ER1_sd + R1_ER1_mean)
    S_init = S_SEIR1 * active_pop
    I_init = max(I_EIR1 * (active_pop - S_init), 1) # Make sure at least 1 Infectious
    R1_init = max(R1_ER1 * (active_pop - S_init - I_init), 1) # Make sure at least 1 Infection
    E_init = max(active_pop - (S_init + I_init + R1_init), 1) # Make sure at least 1 Exposed
    u0 = [S_init, E_init, I_init, R1_init, 1.0, I_init] # Intialize with 1 in R2 so there are no problems when we log for ODE
    log_u0 = log.(u0)
    p0 = [beta_init, gamma, nu, eta]

    # Time-varying parameters
    beta_t_values_no_init = exp.(log(r0_init) .+ cumsum(vec(log_R0_steps_non_centered) * sigma_R0)) * nu
    beta_t_values_with_init = vcat(beta_init, beta_t_values_no_init)
  
  
    function param_affect_beta_IFR!(integrator)
      ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
      integrator.p[1] = beta_t_values_no_init[ind_t] # Replace beta with a new value from beta_t_values
    end
  
    param_callback = PresetTimeCallback(param_change_times, param_affect_beta_IFR!, save_positions = (false, false))
  
    # Solve the ODE  at obstimes
    if extra_ode_precision
      sol = solve(prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false, abstol = abs_tol, reltol = rel_tol,
      u0=log_u0, 
      p=p0, 
      tspan=(0.0, obstimes[end]))   
    else
      sol = solve(prob, Tsit5(), callback = param_callback, saveat = 1.0, save_start = true, verbose = false,
      u0=log_u0, 
      p=p0, 
      tspan=(0.0, obstimes[end])) 
    end
      
    # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    if sol.retcode != :Success
      Turing.@addlogprob! -Inf
      return
    end
  
    sol_reg_scale_array = exp.(Array(sol))
  
  
  
    log_genes_mean = log.(sol_reg_scale_array[3,2:end] .* lambda + (1 - lambda) .* sol_reg_scale_array[4, 2:end]) .+ log(rho_gene) # first entry is the initial conditions, we want 2:end
    new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]
    
    # Likelihood
    for i in 1:l_copies
      index = obstimes[i] # what time in the ode matches the obs time?
      data_log_copies[i] ~ GeneralizedTDist(log_genes_mean[round(Int64,index)], τ, df) 
    end
  
    # Generated quantities
    S = sol_reg_scale_array[1, :]
    r0_t_values_with_init = beta_t_values_with_init / nu
    R0_full_values = zeros(Real, round(Int64,obstimes[end]))
    # ok here's the idea 
    # Rt is actually a function of the S compartment, it changes subtly as S changes, even though R0 is flat for a particular week
    # so for reach flat week of R0, we should still get 7 different values of Rt because of the changes in S
    # r0_t_values_with_init = ones(length(param_change_times) + 1)
    for i in 1:(round(Int64,obstimes[end]))
        R0_full_values[i] = r0_t_values_with_init[floor(Int64, (i-1)/7) + 1]
    end 
    
    # there are too many Rt values if we use all of S b/c S[1] is 0
    # and what we need is S[133] = 132 for the time period (132,133]
    Rt_t_values = R0_full_values .* S[1:end-1] / popsize
  
    return (
      gamma = gamma,
      nu = nu,
      eta = eta,
      rho_gene = rho_gene,
      r0_init = r0_init,
      sigma_R0 = sigma_R0, 
      τ = τ,
      lambda = lambda, 
      df = df,
      beta_t_values = beta_t_values_with_init,
      r0_t_values = r0_t_values_with_init,
      R0_full_values,
      rt_t_values = Rt_t_values,
      S_SEIR1,
      I_EIR1,
      R1_ER1,
      S_init,
      E_init,
      I_init,
      R1_init,
      S = sol_reg_scale_array[1, :],
      E = sol_reg_scale_array[2, :],
      I = sol_reg_scale_array[3, :],
      R1 = sol_reg_scale_array[4, :],
      R2 = sol_reg_scale_array[5, :],
      new_cases = new_cases,
      log_genes_mean = log_genes_mean
    )
  end
  
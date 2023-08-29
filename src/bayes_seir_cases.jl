# Turing model for SEIR-cases model
using LinearAlgebra
prob = ODEProblem{true}(seir_ode_log!,
zeros(5),
(0.0, obstimes[end]),
ones(3))


@model function bayes_seir_cases(data_cases, obstimes, param_change_times, extra_ode_precision, prob, abs_tol, rel_tol)
    # Calculate number of observed datapoints timepoints
    l_copies = length(obstimes)
    l_param_change_times = length(param_change_times)
  
    # Priors
    R0_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), Diagonal(ones(l_param_change_times + 2))) # +2, 1 for init, 1 for sigma 
    S_SEI_non_centered ~ Normal()
    I_EI_non_centered ~ Normal()
    gamma_non_centered ~ Normal() # rate to I
    nu_non_centered ~ Normal() # rate to R1
    rho_case_non_centered ~ Normal() # gene detection rate
    phi_non_centered ~ Normal()
    sigma_R0_non_centered = R0_params_non_centered[2]

  
    # Transformations
    gamma = exp(gamma_non_centered * gamma_sd + gamma_mean)
    nu = exp(nu_non_centered * nu_sd + nu_mean)
    rho_case = logistic(rho_case_non_centered * rho_case_sd + rho_case_mean)
    phi_cases = exp(phi_non_centered * phi_sd + phi_mean)
    sigma_R0 = exp(sigma_R0_non_centered * sigma_R0_sd + sigma_R0_mean)  
    r0_init_non_centered = R0_params_non_centered[1]
    r0_init = exp(r0_init_non_centered * r0_init_sd + r0_init_mean)
    beta_init = r0_init * nu
  
  
    S_SEI = logistic(S_SEI_non_centered * S_SEI_sd + S_SEI_mean)
    I_EI = logistic(I_EI_non_centered * I_EI_sd + I_EI_mean)
  
    log_R0_steps_non_centered = R0_params_non_centered[3:end]
    S_init = S_SEI * active_pop
    I_init = max(I_EI * (active_pop - S_init), 1) # Make sure at least 1 Infectious
    E_init = max(active_pop - (S_init + I_init), 1) # Make sure at least 1 Exposed
    u0 = [S_init, E_init, I_init, 1.0, I_init] # Intialize with 1 in R so there are no problems when we log for ODE
    u0 = [S_init, E_init, I_init, 1.0, I_init] # Intialize with 1 in R so there are no problems when we log for ODE
    log_u0 = log.(u0) 
    p0 = [beta_init, gamma, nu]

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
      sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false, abstol = abs_tol, reltol = rel_tol,
                    u0=log_u0, 
                    p=p0, 
                    tspan=(0.0, obstimes[end]))  
    else
      sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false,
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
    
    total_E2I = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]
    cases_mean = total_E2I .* rho_case
    # Likelihood
    for i in 1:l_copies
      data_cases[i] ~ NegativeBinomial2(max(cases_mean[i], 0.0), phi_cases)
    end
  
    # Generated quantities
    S = sol_reg_scale_array[1, :]
    r0_t_values_with_init = beta_t_values_with_init / nu    
    rt_t_values = r0_t_values_with_init .* S[1:(end-1)] / popsize
  
    return (
      gamma = gamma,
      nu = nu,
      rho_case = rho_case,
      r0_init = r0_init,
      sigma_R0 = sigma_R0, 
      phi_cases = phi_cases,
      beta_t_values = beta_t_values_with_init,
      r0_t_values = r0_t_values_with_init,
      rt_t_values = rt_t_values,
      S_SEI,
      I_EI,
      S_init,
      E_init,
      I_init,
      S = sol_reg_scale_array[1, :],
      E = sol_reg_scale_array[2, :],
      I = sol_reg_scale_array[3, :],
      R = sol_reg_scale_array[4, :],
      C = sol_reg_scale_array[5, :],
      total_E2I = total_E2I,
      cases_mean = cases_mean
    )
  end
  
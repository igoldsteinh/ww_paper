# Turing model for EIRR-ww with ODE solver
using LinearAlgebra
prob = ODEProblem{true}(eirr_ode_log!,
zeros(5),
(0.0, obstimes[end]),
ones(4))

@model function bayes_eirr_student(data_log_copies, obstimes, param_change_times, extra_ode_precision, prob, abs_tol, rel_tol, savetimes, index)
  # Calculate number of observed datapoints timepoints
  l_copies = length(obstimes)
  l_param_change_times = length(param_change_times)

  # Priors
  rt_params_non_centered ~ MvNormal(l_param_change_times + 2, 1) # +2, 1 for var, 1 for init
  I_init_non_centered ~ Normal()
  R1_init_non_centered ~ Normal()
  E_init_non_centered ~ Normal()
  γ_non_centered ~ Normal() # rate to I
  ν_non_centered ~ Normal() # rate to R1
  η_non_centered ~ Normal() # rate to R2
  ρ_gene_non_centered ~ Normal() # scales concentrations
  τ_non_centered ~ Normal() # standard deviation for log scale data
  lambda_non_centered ~ Normal() # proportion of concentration from the I vs R1 compartment
  df ~ Gamma(df_shape, df_scale)

  # Transformations
  γ = exp(γ_non_centered * gamma_sd + gamma_mean)
  ν = exp(ν_non_centered * nu_sd + nu_mean)
  η = exp(η_non_centered * eta_sd + eta_mean)
  ρ_gene = exp(ρ_gene_non_centered * rho_gene_sd + rho_gene_mean)
  τ = exp(τ_non_centered * tau_sd + tau_mean)
  lambda = logistic(lambda_non_centered * lambda_sd + lambda_mean)
  σ_rt_non_centered = rt_params_non_centered[1]
  σ_rt = exp(σ_rt_non_centered * sigma_rt_sd + sigma_rt_mean)
  rt_init_non_centered = rt_params_non_centered[2]
  rt_init = exp(rt_init_non_centered * rt_init_sd + rt_init_mean)
  alpha_init = rt_init * ν

  log_rt_steps_non_centered = rt_params_non_centered[3:end]

  I_init = I_init_non_centered * I_init_sd + I_init_mean
  R1_init = R1_init_non_centered * R1_init_sd + R1_init_mean
  E_init = E_init_non_centered * E_init_sd + E_init_mean
  C_init = I_init 
  u0 = [E_init, I_init, R1_init, 1.0, C_init] # Intialize with 1 in R2 so there are no problems when we log for ODE
  log_u0 = log.(u0) 
  p0 = [alpha_init, γ, ν, η]

  # Time-varying parameters
  alpha_t_values_no_init = exp.(log(rt_init) .+ cumsum(vec(log_rt_steps_non_centered) * σ_rt)) * ν
  alpha_t_values_with_init = vcat(alpha_init, alpha_t_values_no_init)

  function param_affect_β_IFR!(integrator)
    ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
    integrator.p[1] = alpha_t_values_no_init[ind_t] # Replace β with a new value from β_t_values
  end

  param_callback = PresetTimeCallback(param_change_times, param_affect_β_IFR!, save_positions = (false, false))

  # Solve the ODE at intervals of 1.0, could also solve at obstimes
  if extra_ode_precision
    sol = solve(prob, Tsit5(), callback = param_callback, saveat = savetimes, save_start = true, verbose = false, abstol = abs_tol, reltol = rel_tol,
                  u0=log_u0, 
                  p=p0, 
                  tspan=(0.0, obstimes[end]))  
  else
    sol = solve(prob, Tsit5(), callback = param_callback, saveat = savetimes, save_start = true, verbose = false,
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

  log_genes_mean = log.(sol_reg_scale_array[2,2:end] .* lambda + (1 - lambda) .* sol_reg_scale_array[3, 2:end]) .+ log(ρ_gene) # first entry is the initial conditions, we want 2:end
  new_cases = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]

  # Likelihood
  for i in 1:l_copies
    data_log_copies[i] ~ GeneralizedTDist(log_genes_mean[round(Int64,index[i])], τ, df) 
  end

  # Generated quantities
  Rₜ_t_values_with_init = alpha_t_values_with_init/ ν

  return (
    γ = γ,
    ν = ν,
    η = η,
    ρ_gene = ρ_gene,
    rt_init = rt_init,
    σ_rt = σ_rt,
    τ = τ,
    lambda = lambda,
    df = df,
    alpha_t_values = alpha_t_values_with_init,
    Rₜ_t_values = Rₜ_t_values_with_init,
    E_init,
    I_init,
    R1_init,
    E = sol_reg_scale_array[1, :],
    I = sol_reg_scale_array[2, :],
    R1 = sol_reg_scale_array[3, :],
    R2 = sol_reg_scale_array[4, :],
    new_cases = new_cases,
    log_genes_mean = log_genes_mean
  )
end

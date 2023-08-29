# Turing model for EIR-cases model
using LinearAlgebra
using ForwardDiff
@model function bayes_eir_cases!(outs_tmp, data_cases, obstimes, param_change_times)
  # Calculate number of observed datapoints timepoints
  l_copies = length(obstimes)
  l_param_change_times = length(param_change_times)

  # Priors
  rt_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), Diagonal(ones(l_param_change_times + 2))) # +2, 1 for var, 1 for init
  gamma_non_centered ~ Normal() # rate to I
  nu_non_centered ~ Normal() # rate to Re
  I_init_non_centered ~ Normal()
  E_init_non_centered ~ Normal()
  rho_case_non_centered ~ Normal() # case detection rate
  phi_non_centered ~ Normal() # NB overdispersion

  # Transformations
  gamma = exp(gamma_non_centered * gamma_sd + gamma_mean)
  nu = exp(nu_non_centered * nu_sd + nu_mean)
  rho_case = logistic(rho_case_non_centered * rho_case_sd + rho_case_mean)
  phi_cases = exp(phi_non_centered * phi_sd + phi_mean)
  sigma_rt_non_centered = rt_params_non_centered[1]
  sigma_rt = exp(sigma_rt_non_centered * sigma_rt_sd + sigma_rt_mean)
  rt_init_non_centered = rt_params_non_centered[2]
  rt_init = exp(rt_init_non_centered * rt_init_sd + rt_init_mean)
  alpha_init = rt_init * nu
  log_rt_steps_non_centered = rt_params_non_centered[3:end]
  I_init = I_init_non_centered * I_init_sd + I_init_mean
  E_init = E_init_non_centered * E_init_sd + E_init_mean
  C_init = I_init 
  u0 = [E_init, I_init, 1.0, C_init] # Intialize with 1 in R2 so there are no problems when we log for ODE

  # Time-varying parameters
  alpha_t_values_no_init = exp.(log(rt_init) .+ cumsum(vec(log_rt_steps_non_centered) * sigma_rt)) * nu
  alpha_t_values_with_init = vcat(alpha_init, alpha_t_values_no_init)
  
  #Solve the ODE
  sol_reg_scale_array = eir_closed_solution!(outs_tmp, 1:obstimes[end], param_change_times, 0.0, alpha_t_values_with_init, u0, gamma, nu)
  
  # If something goes wrong, discard the sample
  if any(!isfinite, sol_reg_scale_array)
    Turing.@addlogprob! -Inf
    return
  end

  # Create case mean
  new_cases = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]
  cases_mean = new_cases .* rho_case
  
  # Likelihood
  for i in 1:round(Int64, l_copies)
    data_cases[i] ~ NegativeBinomial2(max(cases_mean[i], 0.0), phi_cases)
  end
  
  # Generated quantities
  rt_t_values_with_init = alpha_t_values_with_init/ nu

  return (
    gamma = gamma,
    nu = nu,
    rho_cases = rho_case,
    rt_init = rt_init,
    sigma_rt = sigma_rt,
    phi_cases = phi_cases,
    alpha_t_values = alpha_t_values_with_init,
    rt_t_values = rt_t_values_with_init,
    I_init,
    E_init,
    E = sol_reg_scale_array[2, :],
    I = sol_reg_scale_array[3, :],
    R = sol_reg_scale_array[4, :],
    new_cases = new_cases,
    cases_mean = cases_mean
  )
end

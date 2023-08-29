# Turing model for EIRR-ww model
using LinearAlgebra
using ForwardDiff
using NaNMath
@model function bayes_eirrc_closed!(outs_tmp,data_log_copies, obstimes, param_change_times)
  # Calculate number of observed datapoints timepoints
  l_copies = length(obstimes)
  l_param_change_times = length(param_change_times)

  # Priors
  rt_params_non_centered ~ MvNormal(zeros(l_param_change_times + 2), Diagonal(ones(l_param_change_times + 2))) # +2, 1 for var, 1 for init
  I_init_non_centered ~ Normal()
  R1_init_non_centered ~ Normal()
  E_init_non_centered ~ Normal()
  gamma_non_centered ~ Normal() # rate to I
  nu_non_centered ~ Normal() # rate to Re
  eta_non_centered ~ Normal() # rate to Rd
  rho_gene_non_centered ~ Normal() # gene detection rate
  tau_non_centered ~ Normal() # standard deviation for log scale data
  lambda_non_centered ~ Normal() # percentage of emissions from the I vs R1 compartment
  df ~ Gamma(df_shape, df_scale)

  # Transformations
  gamma = exp(gamma_non_centered * gamma_sd + gamma_mean)
  nu = exp(nu_non_centered * nu_sd + nu_mean)
  eta = exp(eta_non_centered * eta_sd + eta_mean)
  rho_gene = exp(rho_gene_non_centered * rho_gene_sd + rho_gene_mean)
  tau = exp(tau_non_centered * tau_sd + tau_mean)
  lambda = logistic(lambda_non_centered * lambda_sd + lambda_mean)
  sigma_rt_non_centered = rt_params_non_centered[1]
  sigma_rt = exp(sigma_rt_non_centered * sigma_rt_sd + sigma_rt_mean)
  rt_init_non_centered = rt_params_non_centered[2]
  rt_init = exp(rt_init_non_centered * rt_init_sd + rt_init_mean)
  alpha_init = rt_init * nu
  log_rt_steps_non_centered = rt_params_non_centered[3:end]
  I_init = I_init_non_centered * I_init_sd + I_init_mean
  R1_init = R1_init_non_centered * R1_init_sd + R1_init_mean
  E_init = E_init_non_centered * E_init_sd + E_init_mean
  u0 = [E_init, I_init, R1_init, 1.0, 1.0] # Intialize with 1 in R2 b/c it doesn't matter

  # Time-varying parameters
  alpha_t_values_no_init = exp.(log(rt_init) .+ cumsum(vec(log_rt_steps_non_centered) * sigma_rt)) * nu
  alpha_t_values_with_init = vcat(alpha_init, alpha_t_values_no_init)

  # Solve the ODE 
  sol_reg_scale_array = new_eirrc_closed_solution!(outs_tmp, 1:maximum(obstimes), param_change_times, 0.0, alpha_t_values_with_init, u0, gamma, nu, eta)
  
  log_genes_mean = NaNMath.log.(sol_reg_scale_array[3, 2:end] .* lambda + (1 - lambda) .* sol_reg_scale_array[4, 2:end]) .+ log(rho_gene) # first entry is the initial conditions, we want 2:end
  incid = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]
  # Likelihood 
  for i in 1:l_copies
    index = obstimes[i] # what time in the ode matches the obs time?
    data_log_copies[i] ~ GeneralizedTDist(log_genes_mean[round(Int64,index)], tau, df) 
  end

  # Generated quantities
  rt_t_values_with_init = alpha_t_values_with_init/ nu

  return (
    gamma = gamma,
    nu = nu,
    eta = eta,
    rho_gene = rho_gene,
    rt_init = rt_init,
    sigma_rt = sigma_rt,
    tau = tau,
    lambda = lambda,
    df = df,
    alpha_t_values = alpha_t_values_with_init,
    rt_t_values = rt_t_values_with_init,
    E_init,
    I_init,
    R1_init,
    E = sol_reg_scale_array[2, :],
    I = sol_reg_scale_array[3, :],
    R1 = sol_reg_scale_array[4, :],
    R2 = sol_reg_scale_array[5, :],
    C = sol_reg_scale_array[6, :],
    incid = incid,
    log_genes_mean = log_genes_mean
  )
end

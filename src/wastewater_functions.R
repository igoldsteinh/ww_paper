# wasterwater functions
library(stemr)
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)
library(EpiEstim)
library(zoo)
library(sdprisk)
library(brms)
library(rstan)
library(truncnorm)
library(tidybayes)
set.seed(1234)
# source(here::here("src", "huisman_functions.R"))

# first function seirr using stemr ----------------------------------------
sim_ww_data <- function(r0_trajectory,
                        tests,
                          num_sim,
                          seed = 225,
                          pop_size = 300000,
                          init_incid = 10,
                          init_latent = 0,
                          init_recover1 = 0,
                          init_recover2 = 0,
                          mu_param = 1/7,
                          omega_param = 0,
                          ksi_param = 1/18,
                          nu_param = 1/4,
                          gene_rho_param = 0.01,
                          sd_param = 0.69,
                          frac_infec_param= 0.64,
                          case_rho = 9E-5, 
                          case_kappa = 5,
                          obs_dist = "normal",
                          ODE = TRUE){
  # set.seed(1)
  # r0_trajectory = r0
  # num_sim=1
  # seed = 225
  # pop_size = 300000
  # init_incid = 10
  # init_latent = 0
  # init_recover1 = 0
  # init_recover2 = 0
  # mu_param = 1/10
  # omega_param = 0
  # ksi_param = 1/24
  # nu_param = 1/4
  # rhon1_param = 0.003
  # rhon2_param = 0.0006
  # frac_infec_param= 0.64
  # kappa_param = 1.31
  # case_rho = 9E-5
  # case_kappa = 5

  pop_size <- pop_size
  
  R0 <- c(r0_trajectory[1], r0_trajectory) # we need to have an initial r0_trajectory value
  strata <- NULL
  # in R1, we still have emission, in R2, there is no emission
  compartments <- c("S", "E",  "I", "R1", "R2")
  rates <- list(rate("beta_t * I", "S", "E", incidence = T),
                rate("nu", "E", "I", incidence = T),
                rate("mu", "I", "R1", incidence = T),
                rate("ksi", "R1", "R2", incidence = T),
                rate("omega", "R1", "S", incidence = T))
  
  state_initializer <- list(stem_initializer(c(S = pop_size - init_latent - init_incid  - init_recover1 - init_recover2, 
                                               E = init_latent, 
                                               I = init_incid, 
                                               R1 = init_recover1,
                                               R2 = init_recover2),
                                             fixed = T, 
                                             prior = c(pop_size - init_latent - init_incid - init_recover1 - init_recover2,
                                                       init_latent, 
                                                       init_incid,
                                                       init_recover1,
                                                       init_recover2)))
  adjacency <- NULL
  
  
  # setting up time varying R0, translate through Beta which is the actual parameter in the model
  
  parameters = c(mu = mu_param, 
                 omega = omega_param, 
                 nu = nu_param, 
                 gene_rho = gene_rho_param,
                 sd = sd_param,
                 ksi = ksi_param,
                 pop_size = pop_size,
                 frac_infec = frac_infec_param,
                 case_rho = case_rho, 
                 case_kappa = case_kappa)
  
  time <- 0:(length(r0_trajectory))
  tcovar <- cbind(time = time,
                  beta_t = R0 * parameters[["mu"]] / pop_size,
                  tests = c(0,tests))
  
  constants <- c(t0 = 0)
  
  t0 <- 0; 
  tmax <- length(r0_trajectory);
  
  dynamics <-
    stem_dynamics(
      rates = rates,
      tmax = tmax,
      parameters = parameters,
      state_initializer = state_initializer,
      compartments = compartments,
      constants = constants,
      strata = strata,
      adjacency = adjacency,
      tcovar = tcovar,
      messages = T,
      compile_ode = T,
      compile_rates = T,
      compile_lna = F,
      rtol = 1e-6,
      atol = 1e-6,
      step_size = 1e-6
    )
  
  if (obs_dist == "normal"){
    emissions <-
      list(emission(meas_var = "log_gene_copies1", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies2", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies3", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies4", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies5", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies6", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies7", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies8", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies9", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies10", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "cases", 
                    distribution = "negbinomial", 
                    emission_params = c("case_kappa", "E2I * case_rho * tests"),
                    incidence = TRUE,
                    obstimes = seq(1, tmax, by = 1)))
    
  }
  
  if (obs_dist == "student") {
    emissions <-
      list(emission(meas_var = "cases", 
                    distribution = "negbinomial", 
                    emission_params = c("case_kappa", "E2I * case_rho * tests"),
                    incidence = TRUE,
                    obstimes = seq(1, tmax, by = 1)))
    
  }
  measurement_process <- stem_measure(emissions = emissions, 
                                      dynamics = dynamics, 
                                      messages = T)
  
  stem_object <- make_stem(dynamics = dynamics, 
                           measurement_process = measurement_process)
  
  if (ODE == TRUE) {
    stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                          method       = "ode",
                                          paths        = TRUE,
                                          observations = T,
                                          nsim         = num_sim,
                                          census_times = unique(c(0:tmax)))
    
  } else {
    stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                          method       = "gillespie",
                                          paths        = TRUE,
                                          observations = T,
                                          nsim         = num_sim,
                                          census_times = unique(c(0:tmax)))
    
  }
  
  outs <- list(stem_data_stochastic, tcovar)
  
}


#  seiirr using stemr ----------------------------------------
sim_ww_data_seiirr <- function(r0_trajectory,
                        tests,
                        num_sim,
                        seed = 225,
                        pop_size = 300000,
                        init_infec1 = 10,
                        init_infec2 = 10,
                        init_latent = 0,
                        init_recover1 = 0,
                        init_recover2 = 0,
                        nu_param = 1/7,
                        omega_param = 0,
                        eta_param = 1/18,
                        gamma_param = 1/4,
                        case_rho = 9E-5, 
                        case_kappa = 5,
                        ODE = TRUE){
  # set.seed(1)
  # r0_trajectory = r0
  # num_sim=1
  # seed = 225
  # pop_size = 300000
  # init_incid = 10
  # init_latent = 0
  # init_recover1 = 0
  # init_recover2 = 0
  # mu_param = 1/10
  # omega_param = 0
  # ksi_param = 1/24
  # nu_param = 1/4
  # rhon1_param = 0.003
  # rhon2_param = 0.0006
  # frac_infec_param= 0.64
  # kappa_param = 1.31
  # case_rho = 9E-5
  # case_kappa = 5
  
  pop_size <- pop_size
  
  R0 <- c(r0_trajectory[1], r0_trajectory) # we need to have an initial r0_trajectory value
  strata <- NULL
  # in R1, we still have emission, in R2, there is no emission
  # I1 and I2 will have the same infectiousness, same durations, but different emissions
  compartments <- c("S", "E",  "I1", "I2", "R1", "R2")
  rates <- list(rate("beta_t * (I1 + I2)", "S", "E", incidence = T),
                rate("gamma", "E", "I1", incidence = T),
                rate("nu", "I1", "I2", incidence = T),
                rate("nu", "I2", "R1", incidence = T),
                rate("eta", "R1", "R2", incidence = T),
                rate("omega", "R2", "S", incidence = T))
  
  state_initializer <- list(stem_initializer(c(S = pop_size - init_latent - init_infec1 - init_infec2  - init_recover1 - init_recover2, 
                                               E = init_latent, 
                                               I1 = init_infec1,
                                               I2 = init_infec2,
                                               R1 = init_recover1,
                                               R2 = init_recover2),
                                             fixed = T, 
                                             prior = c(pop_size - init_latent - init_infec1 - init_infec2 - init_recover1 - init_recover2,
                                                       init_latent, 
                                                       init_infec1,
                                                       init_infec2,
                                                       init_recover1,
                                                       init_recover2)))
  adjacency <- NULL
  
  
  # setting up time varying R0, translate through Beta which is the actual parameter in the model
  
  parameters = c(gamma = gamma_param, 
                 omega = omega_param, 
                 nu = nu_param, 
                 eta = eta_param,
                 pop_size = pop_size,
                 case_rho = case_rho, 
                 case_kappa = case_kappa)
  
  time <- 0:(length(r0_trajectory))
  tcovar <- cbind(time = time,
                  beta_t = R0 * 1/(2 * 1/parameters[["nu"]]) / pop_size,
                  tests = c(0,tests))
  
  constants <- c(t0 = 0)
  
  t0 <- 0; 
  tmax <- length(r0_trajectory);
  
  dynamics <-
    stem_dynamics(
      rates = rates,
      tmax = tmax,
      parameters = parameters,
      state_initializer = state_initializer,
      compartments = compartments,
      constants = constants,
      strata = strata,
      adjacency = adjacency,
      tcovar = tcovar,
      messages = T,
      compile_ode = T,
      compile_rates = T,
      compile_lna = F,
      rtol = 1e-6,
      atol = 1e-6,
      step_size = 1e-6
    )
  

    emissions <-
      list(emission(meas_var = "cases", 
                    distribution = "negbinomial", 
                    emission_params = c("case_kappa", "E2I1 * case_rho * tests"),
                    incidence = TRUE,
                    obstimes = seq(1, tmax, by = 1)))
    
  
  measurement_process <- stem_measure(emissions = emissions, 
                                      dynamics = dynamics, 
                                      messages = T)
  
  stem_object <- make_stem(dynamics = dynamics, 
                           measurement_process = measurement_process)
  
  if (ODE == TRUE) {
    stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                          method       = "ode",
                                          paths        = TRUE,
                                          observations = T,
                                          nsim         = num_sim,
                                          census_times = unique(c(0:tmax)))
    
  } else {
    stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                          method       = "gillespie",
                                          paths        = TRUE,
                                          observations = T,
                                          nsim         = num_sim,
                                          census_times = unique(c(0:tmax)))
    
  }
  
  outs <- list(stem_data_stochastic, tcovar)
  
}


# sim seirrs model --------------------------------------------------------

sim_seirrs_data <- function(r0_trajectory,
                        tests,
                        num_sim,
                        seed = 225,
                        pop_size = 300000,
                        init_incid = 10,
                        init_latent = 0,
                        init_recover1 = 0,
                        init_recover2 = 0,
                        mu_param = 1/7,
                        omega_param_trajectory,
                        ksi_param = 1/18,
                        nu_param = 1/4,
                        gene_rho_param = 0.01,
                        sd_param = 0.69,
                        frac_infec_param= 0.64,
                        case_rho = 9E-5, 
                        case_kappa = 5,
                        obs_dist = "normal",
                        ODE = TRUE){
  # set.seed(seed)
  # r0_trajectory = r0
  # num_sim=1
  # seed = 225
  # pop_size = 300000
  # init_incid = 10
  # init_latent = 0
  # init_recover1 = 0
  # init_recover2 = 0
  # mu_param = 1/10
  # omega_param = 0
  # ksi_param = 1/24
  # nu_param = 1/4
  # rhon1_param = 0.003
  # rhon2_param = 0.0006
  # frac_infec_param= 0.64
  # kappa_param = 1.31
  # case_rho = 9E-5 
  # case_kappa = 5
  
  pop_size <- pop_size
  
  R0 <- r0_trajectory
  strata <- NULL
  # in R1, we still have emission, in R2, there is no emission
  compartments <- c("S", "E",  "I", "R1", "R2")
  rates <- list(rate("beta_t * I", "S", "E", incidence = T),
                rate("nu", "E", "I", incidence = T),
                rate("mu", "I", "R1", incidence = T),
                rate("ksi", "R1", "R2", incidence = T),
                rate("omega_t", "R2", "S", incidence = T))
  
  state_initializer <- list(stem_initializer(c(S = pop_size - init_latent - init_incid  - init_recover1 - init_recover2, 
                                               E = init_latent, 
                                               I = init_incid, 
                                               R1 = init_recover1,
                                               R2 = init_recover2),
                                             fixed = T, 
                                             prior = c(pop_size - init_latent - init_incid - init_recover1 - init_recover2,
                                                       init_latent, 
                                                       init_incid,
                                                       init_recover1,
                                                       init_recover2)))
  adjacency <- NULL
  
  
  # setting up time varying R0, translate through Beta which is the actual parameter in the model
  
  parameters = c(mu = mu_param, 
                 nu = nu_param, 
                 gene_rho = gene_rho_param,
                 sd = sd_param,
                 ksi = ksi_param,
                 pop_size = pop_size,
                 frac_infec = frac_infec_param,
                 case_rho = case_rho, 
                 case_kappa = case_kappa)
  
  time <- 0:(length(r0_trajectory))
  tcovar <- cbind(time = time,
                  beta_t = R0 * parameters[["mu"]] / pop_size,
                  omega_t = omega_param_trajectory,
                  tests = tests)
  
  constants <- c(t0 = 0)
  
  t0 <- 0; 
  tmax <- length(r0_trajectory);
  
  dynamics <-
    stem_dynamics(
      rates = rates,
      tmax = tmax,
      parameters = parameters,
      state_initializer = state_initializer,
      compartments = compartments,
      constants = constants,
      strata = strata,
      adjacency = adjacency,
      tcovar = tcovar,
      messages = T,
      compile_ode = T,
      compile_rates = T,
      compile_lna = F,
      rtol = 1e-6,
      atol = 1e-6,
      step_size = 1e-6
    )
  
  if (obs_dist == "normal"){
    emissions <-
      list(emission(meas_var = "log_gene_copies1", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies2", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies3", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies4", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies5", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies6", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies7", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies8", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies9", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "log_gene_copies10", 
                    distribution = "gaussian", 
                    emission_params = c("log((I*frac_infec + R1*(1-frac_infec)) * gene_rho)", "sd"),
                    incidence = FALSE,
                    obstimes = seq(1, tmax, by = 1)),
           emission(meas_var = "cases", 
                    distribution = "negbinomial", 
                    emission_params = c("case_kappa", "E2I * case_rho * tests"),
                    incidence = TRUE,
                    obstimes = seq(1, tmax, by = 1)))
    
  }
  
  if (obs_dist == "student") {
    emissions <-
      list(emission(meas_var = "cases", 
                    distribution = "negbinomial", 
                    emission_params = c("case_kappa", "E2I * case_rho * tests"),
                    incidence = TRUE,
                    obstimes = seq(1, tmax, by = 1)))
    
  }
  measurement_process <- stem_measure(emissions = emissions, 
                                      dynamics = dynamics, 
                                      messages = T)
  
  stem_object <- make_stem(dynamics = dynamics, 
                           measurement_process = measurement_process)
  
  if (ODE == TRUE) {
    stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                          method       = "ode",
                                          paths        = TRUE,
                                          observations = T,
                                          nsim         = num_sim,
                                          census_times = unique(c(0:tmax)))
    
  } else {
    stem_data_stochastic <- simulate_stem(stem_object  = stem_object,
                                          method       = "gillespie",
                                          paths        = TRUE,
                                          observations = T,
                                          nsim         = num_sim,
                                          census_times = unique(c(0:tmax)))
    
  }
  
  outs <- list(stem_data_stochastic, tcovar)
  
}



# create stemr data -------------------------------------------------------
# sim = ww_data_stoch
# r0_trajectory = r0
# tests = tests
# pop_size = pop_size
# gene_rho_param = gene_rho
# frac_infec_param = lambda
# index = 11
# lump = 2
# t_df = 2.99
# t_sd = 0.5
# obs_dist = "student"
# ODE = TRUE
create_ww_data <- function(sim, 
                              r0_trajectory, 
                              tests,
                              pop_size,
                              gene_rho_param,
                              frac_infec_param,
                              index = 1,
                              lump = 2,
                              t_df = 1,
                              t_sd = 1,
                              obs_dist = "normal",
                              ODE = TRUE) {
  # 
  # sim <- ww_data_stoch
  # index = 1
  # lump = 2
  if (ODE == TRUE) {
    sim_paths <- as.data.frame(sim[[1]][["paths"]][1])
    natural_paths <- as.data.frame(sim[[1]][["natural_paths"]][1])
    
  } else {
    sim_paths <- as.data.frame(sim[[1]][["paths"]][index])
    
  }
  #view(sim_paths)
  if (ODE == TRUE) {
    true_incidence <- sim_paths %>% 
                      left_join(natural_paths, by = "time") %>% 
                      dplyr::select(time, S2E, E2I, S, E, I, R1, R2) %>% 
                      rename("incidence" = "S2E")
    
  } else {
    true_incidence <- sim_paths %>% 
      dplyr::select(time, S2E, E2I, S, E, I, R1, R2) %>% 
      rename("incidence" = "S2E")
    
  }
  

  #data
  sim_data_sets <- as.data.frame(sim[[1]][["datasets"]][index])
  sim_data_sets <- rbind(c(0,0), sim_data_sets)
  
  
  time <- 0:(length(r0))
  tests_frame <- data.frame(time, c(0,tests))
  
  
  true_r0 <- data.frame(time = time, r0 = c(r0[1],r0))
  
  if (obs_dist == "normal") {
    data <- sim_data_sets %>%
      mutate(keep = time %% lump,
             week = floor((time -1) / 7),
             lumptime = floor(time/lump)) %>%
      left_join(true_r0, by = "time") %>% 
      left_join(true_incidence, by = "time") %>%
      left_join(tests_frame, by = "time") %>%
      mutate(percent_s = S/pop_size,
             true_rt = r0 * percent_s,
             log_gene_mean = log(I * frac_infec_param + (1-frac_infec_param)*R1) + log(gene_rho_param),
             log_mean_copiesten = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3) +
                                         exp(log_gene_copies4) + exp(log_gene_copies5) + exp(log_gene_copies6) +
                                         exp(log_gene_copies7) + exp(log_gene_copies8) + exp(log_gene_copies9) + 
                                         exp(log_gene_copies10))/10),
             log_mean_copiesthree = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3))/3)) %>%
      group_by(week) %>%
      mutate(order = row_number())
    
  }

  if (obs_dist == "student") {
    data <- sim_data_sets %>%
      mutate(keep = time %% lump,
             week = floor((time -1) / 7),
             lumptime = floor(time/lump)) %>%
      left_join(true_r0, by = "time") %>% 
      left_join(true_incidence, by = "time") %>%
      left_join(tests_frame, by = "time") %>%
      mutate(percent_s = S/pop_size,
             lag_percent_s = lag(percent_s),
             true_rt = r0 * lag_percent_s, # we are defining R_t_1 = R_0_1 * S(t-1)/N
             log_gene_mean = log(I * frac_infec_param + (1-frac_infec_param)*R1) + log(gene_rho_param)) %>%
      rowwise() %>%
      mutate(
             log_gene_copies1 = log_gene_mean + (t_sd * rt(1,t_df)),
             log_gene_copies2 = log_gene_mean + (t_sd * rt(1,t_df)),
             log_gene_copies3 = log_gene_mean + (t_sd * rt(1,t_df)),
             log_gene_copies4 = log_gene_mean + (t_sd * rt(1,t_df)),
             log_gene_copies5 = log_gene_mean + t_sd * rt(1,t_df),
             log_gene_copies6 = log_gene_mean + t_sd * rt(1,t_df),
             log_gene_copies7 = log_gene_mean + t_sd * rt(1,t_df),
             log_gene_copies8 = log_gene_mean + t_sd * rt(1,t_df),
             log_gene_copies9 = log_gene_mean + t_sd * rt(1,t_df),
             log_gene_copies10 = log_gene_mean + t_sd * rt(1,t_df)
      )
    
    
    
    
    data <- data %>%
             mutate(log_mean_copiesten = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3) +
                                         exp(log_gene_copies4) + exp(log_gene_copies5) + exp(log_gene_copies6) +
                                         exp(log_gene_copies7) + exp(log_gene_copies8) + exp(log_gene_copies9) +
                                         exp(log_gene_copies10))/10),
             log_mean_copiesthree = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3))/3)) %>%
      group_by(week) %>%
      mutate(order = row_number())
    
  }

  
  full_lump_data <- data %>%
    ungroup() %>% 
    filter(keep == 1)  %>%
    ungroup() 

  res <- list("lump_data" =  full_lump_data, "full_data"=data)
  return(res)
  
}

# create ww data for SEIIRR -------------------------------------------------------
# sim = ww_data_stoch
# r0_trajectory = r0
# tests = tests
# pop_size = pop_size
# gene_rho_param = gene_rho
# frac_infec_param = lambda
# index = 11
# lump = 2
# t_df = 2.99
# t_sd = 0.5
# obs_dist = "student"
# ODE = TRUE
create_ww_data_seiirr <- function(sim, 
                           r0_trajectory, 
                           tests,
                           pop_size,
                           gene_rho_param,
                           frac_total_infec_param, #what fraction of shedding is in infectiousness total
                           frac_infec1_param, # what fraction of the shedding in infectiousness is in I1
                           index = 1,
                           lump = 2,
                           t_df = 1,
                           t_sd = 1,
                           ODE = TRUE) {
  # 
  # sim <- ww_data_stoch
  # index = 1
  # lump = 2
  if (ODE == TRUE) {
    sim_paths <- as.data.frame(sim[[1]][["paths"]][1])
    natural_paths <- as.data.frame(sim[[1]][["natural_paths"]][1])
    
  } else {
    sim_paths <- as.data.frame(sim[[1]][["paths"]][index])
    
  }
  #view(sim_paths)
  if (ODE == TRUE) {
    true_incidence <- sim_paths %>% 
      left_join(natural_paths, by = "time") %>% 
      dplyr::select(time, S2E, E2I1, S, E, I1, I2, R1, R2) %>% 
      rename("incidence" = "S2E")
    
  } else {
    true_incidence <- sim_paths %>% 
      dplyr::select(time, S2E, E2I1, S, E, I1, I2, R1, R2) %>% 
      rename("incidence" = "S2E")
    
  }
  
  
  #data
  sim_data_sets <- as.data.frame(sim[[1]][["datasets"]][index])
  sim_data_sets <- rbind(c(0,0), sim_data_sets)
  
  
  time <- 0:(length(r0))
  tests_frame <- data.frame(time, c(0,tests))
  
  
  true_r0 <- data.frame(time = time, r0 = c(r0[1],r0))
  
    data <- sim_data_sets %>%
      mutate(keep = time %% lump,
             week = floor((time -1) / 7),
             lumptime = floor(time/lump)) %>%
      left_join(true_r0, by = "time") %>% 
      left_join(true_incidence, by = "time") %>%
      left_join(tests_frame, by = "time") %>%
      mutate(percent_s = S/pop_size,
             lag_percent_s = lag(percent_s),
             true_rt = r0 * lag_percent_s, # we are defining R_t_1 = R_0_1 * S(t-1)/N
             log_gene_mean = log((I1 * frac_infec1_param * frac_total_infec_param) + 
                                   (I2 * (1-frac_infec1_param) * frac_total_infec_param) + 
                                   (1-frac_total_infec_param)*R1) + 
                             log(gene_rho_param)) %>%
      rowwise() %>%
      mutate(
        log_gene_copies1 = log_gene_mean + (t_sd * rt(1,t_df)),
        log_gene_copies2 = log_gene_mean + (t_sd * rt(1,t_df)),
        log_gene_copies3 = log_gene_mean + (t_sd * rt(1,t_df)),
        log_gene_copies4 = log_gene_mean + (t_sd * rt(1,t_df)),
        log_gene_copies5 = log_gene_mean + t_sd * rt(1,t_df),
        log_gene_copies6 = log_gene_mean + t_sd * rt(1,t_df),
        log_gene_copies7 = log_gene_mean + t_sd * rt(1,t_df),
        log_gene_copies8 = log_gene_mean + t_sd * rt(1,t_df),
        log_gene_copies9 = log_gene_mean + t_sd * rt(1,t_df),
        log_gene_copies10 = log_gene_mean + t_sd * rt(1,t_df)
      )
    
    
    
    
    data <- data %>%
      mutate(log_mean_copiesten = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3) +
                                         exp(log_gene_copies4) + exp(log_gene_copies5) + exp(log_gene_copies6) +
                                         exp(log_gene_copies7) + exp(log_gene_copies8) + exp(log_gene_copies9) +
                                         exp(log_gene_copies10))/10),
             log_mean_copiesthree = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3))/3)) %>%
      group_by(week) %>%
      mutate(order = row_number())
    

  
  
  full_lump_data <- data %>%
    ungroup() %>% 
    filter(keep == 1)  %>%
    ungroup() 
  
  res <- list("lump_data" =  full_lump_data, "full_data"=data)
  return(res)
  
}

# create_case_data --------------------------------------------------------
create_case_data <- function(sim, 
                           r0_trajectory, 
                           tests,
                           pop_size,
                           index = 1,
                           ODE = TRUE) {
  # 
  # sim <- ww_data_stoch
  # index = 1
  # lump = 2
  if (ODE == TRUE) {
    sim_paths <- as.data.frame(sim[[1]][["paths"]][1])
    natural_paths <- as.data.frame(sim[[1]][["natural_paths"]][1])
    
  } else {
    sim_paths <- as.data.frame(sim[[1]][["paths"]][index])
    
  }
  #view(sim_paths)
  if (ODE == TRUE) {
    true_incidence <- sim_paths %>% 
      left_join(natural_paths, by = "time") %>% 
      dplyr::select(time, S2E, E2I) %>% 
      rename("incidence" = "S2E")
    
  } else {
    true_incidence <- sim_paths %>% 
      dplyr::select(time, S2E, E2I) %>% 
      rename("incidence" = "S2E")
    
  }
  
  
  #data
  sim_data_sets <- as.data.frame(sim[[1]][["datasets"]][index])
  
  
  time <- sim_data_sets$time
  tests_frame <- data.frame(time, tests)
  
  
  weekly_data <- sim_data_sets %>% 
                 left_join(tests_frame, by = "time") %>%
                 left_join(true_incidence, by = "time") %>%
                 mutate(week = floor((time - 1)/7)) %>%
                 group_by(week) %>%
                 summarise(total_tests = sum(tests),
                           total_cases = sum(cases),
                           total_E2I = sum(E2I)) %>%
                 mutate(new_week = week + 1)
                 
  return(weekly_data)
  
}


# create new lump data from full data -------------------------------------

new_lump_data <- function(full_data, lump) {
  data <- full_data %>%
    mutate(keep = time %% lump,
           week = floor(time / 7),
           lumptime = floor(time/lump))
  
  lump_case_data <- data %>% 
    group_by(lumptime) %>%
    summarise(total_cases = sum(cases),
              total_tests = sum(tests))
  
  
  full_lump_data <- data %>%
    ungroup() %>% 
    filter(keep == 1)  %>%
    ungroup() %>%
    left_join(lump_case_data, by = "lumptime")
  
  return(full_lump_data)
  
}


# make prior samples ------------------------------------------------------

make_prior_samples <- function(prior_gq) {
  prior_gq_samples <-
    prior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]", negate = T))

  return(prior_gq_samples)
  
}


# make quantiles for time varying parameters ------------------------------
make_timevarying_prior_quantiles <- function(prior_gq) {
  prior_gq_samples <-
    prior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]")) %>%
    mutate(time = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(.,tibble(time = 0:max(.$time)))
  
  return(prior_gq_samples)
  
}


# make prior quantiles ----------------------------------------------------
make_prior_quantiles <- function(prior_samples) {
  prior_quantiles <- prior_samples %>%
    group_by(name) %>%
    summarise(quantiles = quantile(value, c(0.025, 0.5, 0.975)),
              q = c(0.025, 0.5, 0.975)) %>%
    ungroup() %>%
    pivot_wider(id_cols = name, names_from = q, values_from = quantiles) 
  
  return(prior_quantiles)
}


# plot simulated data -----------------------------------------------------

plot_sim_data <- function(full_sim_data) {
  
  maxI <- max(full_sim_data$I) 
  maxRe <- max(full_sim_data$R1)
  maxgene <- max(c(full_sim_data$log_gene_copies1, 
                   full_sim_data$log_gene_copies2, 
                   full_sim_data$log_gene_copies3))
  
  maxcase <- max(full_sim_data$cases)
  maxtest <- max(full_sim_data$tests)
  
  sim_prev <- full_sim_data %>%
    ggplot(aes(x = time, y = I)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxI + 10)) +
    ggtitle("Prev I") +
    ylab("Prev I") +
    xlab("") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  sim_R1 <- full_sim_data %>%
    ggplot(aes(x = time, y = R1)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxRe)) +
    scale_x_continuous(expand = c(0,0)) +
    ggtitle("Prev R1") +
    ylab("Prev R1") +
    xlab("") +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  
  sim_rt <- full_sim_data %>%
    ggplot(aes(x = time, y = true_rt)) + 
    geom_line() +
    ggtitle("Rt") +
    ylab("Rt") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  sim_data_long <- full_sim_data %>%
    ungroup() %>%
    dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
    pivot_longer(cols = -time, values_to = "log_genes") 
  
  
  full_sim_data <- full_sim_data %>% 
      mutate(mean_log = (log_gene_copies1 + log_gene_copies2 + log_gene_copies3)/3)
  loggene_plot <- sim_data_long %>%
    ggplot(aes(x = time, y = log_genes, color = name)) + 
    geom_point() + 
    geom_line() + 
    geom_line(data = full_sim_data, aes(x = time, y = mean_log), color = "black", size = 1.2) +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
    scale_x_continuous(expand = c(0,0)) +
    
    ggtitle("Log Gene Counts") + 
    ylab("Log RNA") + 
    xlab("Time") +
    theme_bw()+
    theme(text = element_text(size = 14),
          legend.position = "none")
  
  sim_cases <- full_sim_data %>%
    ggplot(aes(x = time, y = cases)) +
    geom_point() +
    geom_line() +
    ggtitle("Observed Daily Cases") +
    ylab("Cases") +
    xlab("Time") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxcase)) +
    
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  sim_tests <- full_sim_data %>%
    ggplot(aes(x = time, y = tests)) +
    geom_point() +
    geom_line() +
    ggtitle("Daily Tests") +
    ylab("Tests") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxtest + 20)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  full_plot <- (sim_rt + sim_prev) /
    (sim_R1 + sim_tests) /
    (loggene_plot + sim_cases) +
    plot_layout(guides = 'collect')
  
}

plot_sim_gene <- function(full_sim_data,
                          mean_line = TRUE) {
  
  maxgene <- max(c(full_sim_data$log_gene_copies1, 
                   full_sim_data$log_gene_copies2, 
                   full_sim_data$log_gene_copies3))
  


  sim_data_long <- full_sim_data %>%
    ungroup() %>%
    dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
    pivot_longer(cols = -time, values_to = "log_genes") 
  full_sim_data <- full_sim_data %>% 
    mutate(mean_log = (log_gene_copies1 + log_gene_copies2 + log_gene_copies3)/3)
  
  if (mean_line == TRUE) {
    loggene_plot <- sim_data_long %>%
      ggplot(aes(x = time, y = log_genes, color = name)) + 
      geom_point() + 
      geom_line() + 
      geom_line(data = full_sim_data, aes(x = time, y = mean_log), color = "black", size = 1.2) +
      scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
      scale_x_continuous(expand = c(0,0)) +
      
      ggtitle("Log Gene Counts") + 
      ylab("Log RNA") + 
      xlab("Time") +
      theme_bw()+
      theme(text = element_text(size = 14),
            legend.position = "none")
    
  }
  
  if (mean_line == FALSE) {
    loggene_plot <- sim_data_long %>%
      ggplot(aes(x = time, y = log_genes, color = name)) + 
      geom_point() + 
      geom_line() + 
      scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
      scale_x_continuous(expand = c(0,0)) +
      
      ggtitle("Log Gene Counts") + 
      ylab("Log RNA") + 
      xlab("Time") +
      theme_bw()+
      theme(text = element_text(size = 14),
            legend.position = "none")
    
  }
  
  loggene_plot
  

}

plot_sim_data_short <- function(full_sim_data) {
  
  maxI <- max(full_sim_data$I) 
  maxRe <- max(full_sim_data$R1)
  maxgene <- max(c(full_sim_data$log_gene_copies1, 
                   full_sim_data$log_gene_copies2, 
                   full_sim_data$log_gene_copies3))
  
  maxcase <- max(full_sim_data$cases)
  maxtest <- max(full_sim_data$tests)
  

  
  sim_rt <- full_sim_data %>%
    ggplot(aes(x = time, y = true_rt)) + 
    geom_line() +
    ggtitle("Rt") +
    ylab("Rt") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  sim_data_long <- full_sim_data %>%
    ungroup() %>%
    dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
    pivot_longer(cols = -time, values_to = "log_genes") 
  
  
  full_sim_data <- full_sim_data %>% 
    mutate(mean_log = (log_gene_copies1 + log_gene_copies2 + log_gene_copies3)/3)
  loggene_plot <- sim_data_long %>%
    ggplot(aes(x = time, y = log_genes, color = name)) + 
    geom_point() + 
    geom_line() + 
    geom_line(data = full_sim_data, aes(x = time, y = mean_log), color = "black", size = 1.2) +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
    scale_x_continuous(expand = c(0,0)) +
    
    ggtitle("Log Gene Counts") + 
    ylab("Log RNA") + 
    xlab("Time") +
    theme_bw()+
    theme(text = element_text(size = 14),
          legend.position = "none")
  
  sim_cases <- full_sim_data %>%
    ggplot(aes(x = time, y = cases)) +
    geom_point() +
    geom_line() +
    ggtitle("Observed Daily Cases") +
    ylab("Cases") +
    xlab("Time") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxcase)) +
    
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  

  full_plot <- (sim_rt + loggene_plot + sim_cases) +
    plot_layout(guides = 'collect')
  
}


# make fixed posterior samples --------------------------------------------

make_fixed_posterior_samples <- function(posterior_gq) {
  posterior_gq_samples <- posterior_gq %>%
  filter(str_detect(name, "\\[\\d+\\]", negate = T))
  
  return(posterior_gq_samples)
}


# make time varying posterior quantiles -----------------------------------
make_timevarying_posterior_quantiles <- function(posterior_gq) {
    timevarying_posterior_quantiles <-
    posterior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]")) %>%
    mutate(time = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(.,tibble(time = 0:max(.$time)))
    
    return(timevarying_posterior_quantiles)
  
}


# make rt plot ------------------------------------------------------------
make_rt_plot <- function(timevarying_quantiles, 
                         simdata, 
                         initial_conds, 
                         varnames, 
                         varnum = 8,
                         SEIRR = TRUE){
  
  if (SEIRR == TRUE) {
    rt_plot <- timevarying_quantiles %>%
      filter(name == var_names[varnum]) %>%
      mutate(new_time = time) %>%
      ggplot() +
      geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper,fill = fct_rev(ordered(.width)))) +
      geom_point(data = simdata, mapping = aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept=1, linetype="dashed", color = "black") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw()+
      theme(legend.position = "none") + 
      ylab("Rt") +
      xlab("Time") +
      ggtitle("Rt Wastewater") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
      scale_x_continuous(expand = c(0,0))
    
  }
  
  if (SEIRR == FALSE) {
    timevarying_quantiles <- timevarying_quantiles %>% 
      filter(name == var_names[varnum]) %>%
      mutate(new_time = time * 7)
                             
    last_quantiles <- timevarying_quantiles %>% 
                      filter(new_time == max(new_time)) %>%
                      mutate(new_time = max(simdata$new_time))
    
    rt_quantiles <- bind_rows(timevarying_quantiles, last_quantiles)
    rt_plot <- rt_quantiles %>%
      ggplot() +
      geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper,fill = fct_rev(ordered(.width)))) +
      geom_point(data = simdata, mapping = aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept=1, linetype="dashed", color = "black") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw()+
      theme(legend.position = "none") + 
      ylab("Rt") +
      xlab("Time") +
      ggtitle("Rt Wastewater") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
      scale_x_continuous(expand = c(0,0))
    
    
  }
  
  return(rt_plot)
  
}


# make trajectory plot ----------------------------------------------------

make_comp_plot <- function(posterior_timevarying_quantiles, 
                           initial_conds, 
                           traj_names,
                           true_trajectory, 
                           cases = FALSE){
  if (cases == FALSE) {
    true_init_trajectory <- initial_conds %>% dplyr::select(S, E, I , R1, log_genes_mean, new_time) %>%
      pivot_longer(cols = -(new_time)) %>%
      rename("true_value" = "value")
  } else {
    true_init_trajectory <- initial_conds %>% dplyr::select(S, E, I, new_time) %>%
      pivot_longer(cols = -(new_time)) %>%
      rename("true_value" = "value")
    
  }

  if (cases == FALSE) {
    posterior_trajectory <- posterior_timevarying_quantiles %>%
      filter(name %in% traj_names) %>%
      left_join(true_trajectory, by = c("time" = "new_time", 
                                        "name"= "name"))
  }

  else {
    posterior_trajectory <- posterior_timevarying_quantiles %>%
      filter(name %in% traj_names) %>%
      left_join(true_trajectory, by = c("time" = "new_week", 
                                        "name"= "name"))
    }
  
  traj_plot <- posterior_trajectory %>%
    ggplot() +
    geom_lineribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper)) +
    geom_point(data = posterior_trajectory, mapping = aes(x = time, y = true_value), color = "coral1") +
    geom_point(data = true_init_trajectory, mapping = aes(x = new_time, y = true_value), shape= 0, color = "goldenrod2") +
    scale_fill_brewer(name = "Credible Interval Width") +
    theme_bw() +
    facet_wrap(.~name, scales = "free_y") +
    ggtitle("Posteriors of ODE compartments")

  return(traj_plot)
}


# make plot of prior and posteriors for fixed params ----------------------
make_fixed_param_plot <- function(posterior_gq_samples, prior_gq_samples){
  priors_and_posteriors <- rbind(posterior_gq_samples %>% dplyr::select(name, value, type, true_value), prior_gq_samples)
  
  
  param_plot <- priors_and_posteriors %>%
    ggplot(aes(value, type, fill = type)) +
    stat_halfeye(normalize = "xy")  +
    geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
    facet_wrap(. ~ name, scales = "free_x") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    ggtitle("Fixed parameter prior and posteriors")
  
  return(param_plot)
  
}


# make posterior predictive intervals -------------------------------------
# posterior_predictive = eirr_post_pred
# sim_data = simdata
# cases = FALSE
# ten_sim = ten_sim_val
# three_mean = three_mean_val
make_post_pred_intervals <- function(posterior_predictive, sim_data, cases = FALSE, ten_sim = FALSE, three_mean = FALSE){
  if (cases == FALSE & ten_sim == FALSE & three_mean == FALSE) {
    obs_time <- sim_data %>%
    dplyr::select(log_gene_copies1, log_gene_copies2, log_gene_copies3, new_time) %>%
    pivot_longer(-(new_time)) %>%
    filter(value > 0) %>%
    arrange(name) %>%
    mutate(obs_index = row_number()) %>%
    dplyr::select(obs_index, new_time)
  
  posterior_predictive_samples <- posterior_predictive %>%
    pivot_longer(-c(iteration, chain)) %>%
    mutate(obs_index = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                summarize(value = sum(value),
                          .groups = "drop"))  %>%
    left_join(obs_time, by = "obs_index")
  
  posterior_predictive_intervals <- posterior_predictive_samples %>%
    select(new_time, name, value) %>%
    group_by(new_time, name) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    select(new_time, name, value, starts_with("."))
  } else if (cases == FALSE & ten_sim == TRUE & three_mean == FALSE) {
    obs_time <- sim_data %>%
      dplyr::select(log_gene_copies1, 
                    log_gene_copies2, 
                    log_gene_copies3, 
                    log_gene_copies4, 
                    log_gene_copies5, 
                    log_gene_copies6, 
                    log_gene_copies7,
                    log_gene_copies8,
                    log_gene_copies9,
                    log_gene_copies10,
                    new_time) %>%
      pivot_longer(-(new_time)) %>%
      filter(value > 0) %>%
      arrange(name) %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_time)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
    posterior_predictive_intervals <- posterior_predictive_samples %>%
      select(new_time, name, value) %>%
      group_by(new_time, name) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)) %>%
      select(new_time, name, value, starts_with("."))
    
  } else if (cases == FALSE & ten_sim == FALSE & three_mean == TRUE) {
    obs_time <- sim_data %>%
      dplyr::select(log_mean_copiesthree,
                    new_time) %>%
      pivot_longer(-(new_time)) %>%
      filter(value > 0) %>%
      arrange(name) %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_time)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
    posterior_predictive_intervals <- posterior_predictive_samples %>%
      select(new_time, name, value) %>%
      group_by(new_time, name) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)) %>%
      select(new_time, name, value, starts_with("."))
    
  }
  else {
    obs_time <- sim_data %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_week)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
    posterior_predictive_intervals <- posterior_predictive_samples %>%
      select(new_week, name, value) %>%
      group_by(new_week, name) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)) %>%
      select(new_week, name, value, starts_with("."))
    
  }
  
  return(posterior_predictive_intervals)
}


# make posterior predictive plot ------------------------------------------
make_post_pred_plot <- function(posterior_predictive_intervals, 
                                sim_data, 
                                cases = FALSE,
                                ten_sim = FALSE,
                                three_mean = FALSE) {

  if (cases == FALSE & ten_sim == FALSE & three_mean == FALSE) {
    true_data <- sim_data %>%
    dplyr::select(new_time, 
                  log_gene_copies1, 
                  log_gene_copies2, 
                  log_gene_copies3) %>%
    rename("log_copies1" = "log_gene_copies1",
           "log_copies2" = "log_gene_copies2",
           "log_copies3" = "log_gene_copies3") %>%
    pivot_longer(cols = - new_time) %>%
    rename("true_value" = "value") %>%
    filter(true_value > 0)
  
  posterior_predictive_intervals <- posterior_predictive_intervals %>%
    left_join(true_data, by = c("new_time"))
  
  posterior_predictive_plot <- posterior_predictive_intervals %>%
    ggplot() +
    geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
    geom_line(aes(x = new_time, y = value)) + 
    geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
    scale_fill_brewer(name = "Credible Interval Width") +
    # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
    theme_bw() + 
    ggtitle("Posterior Predictive (ODE)")
  } else if (cases == FALSE & ten_sim == TRUE & three_mean == FALSE) {
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_gene_copies1, 
                    log_gene_copies2, 
                    log_gene_copies3,
                    log_gene_copies4,
                    log_gene_copies5,
                    log_gene_copies6,
                    log_gene_copies7,
                    log_gene_copies8,
                    log_gene_copies9,
                    log_gene_copies10,
      ) %>%
      rename("log_copies1" = "log_gene_copies1",
             "log_copies2" = "log_gene_copies2",
             "log_copies3" = "log_gene_copies3",
             "log_copies4" = "log_gene_copies4",
             "log_copies5" = "log_gene_copies5",
             "log_copies6" = "log_gene_copies6",
             "log_copies7" = "log_gene_copies7",
             "log_copies8" = "log_gene_copies8",
             "log_copies9" = "log_gene_copies9",
             "log_copies10" = "log_gene_copies10") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  } else if (cases == FALSE & ten_sim == FALSE & three_mean == TRUE) {
    # not sure what is going to happen here, wait until we have a posterior to work with before finishing
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_mean_copiesthree
      ) %>%
      rename("log_mean_copies" = "log_mean_copiesthree") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  }
  else {
    true_data <- sim_data %>%
      dplyr::select(new_week, 
                    total_cases) %>%
      rename("true_value" = "total_cases")
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_week"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_week, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_week, y = value)) + 
      geom_point(mapping = aes(x = new_week, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      # scale_fill_manual(values=c("skyblue1", "skyblue2", "skyblue3"), name="fill") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  }
  
  return(posterior_predictive_plot)
  
}


# calculate huisman rt ----------------------------------------------------

calculate_huisman_rt <- function(simdata, sim = TRUE) {
  
  if (sim == TRUE){
    date_length <- length(simdata$time)
    fake_date_one <- ymd("20220110")
    fake_date <- fake_date_one
    for (i in 2:date_length) {
      fake_date <- append(fake_date, fake_date[i-1] + ddays(2))
    }
    
    mean_data <- simdata %>%
      mutate(expone = exp(log_gene_copies1),
             exptwo = exp(log_gene_copies2),
             expthree = exp(log_gene_copies3)) %>%
      group_by(new_time) %>%
      summarise(mean_gene = (expone + exptwo + expthree)/3)
    
    
    
    
    data <- mean_data
    dates <- fake_date
    data <- cbind(data, date = dates)
    
  }
  
  if (sim == FALSE) {
    data <- simdata %>%
            mutate(mean_gene = (gene_copy_1 + gene_copy_2 + gene_copy_3)/3)
  }
  
    data <- data %>%
      mutate(orig_data = TRUE) %>%
      complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
      mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
      mutate(region = 'LA') %>%
      distinct() 
    
    ###########################################################
    ## Normalisation of WW data ####
    
    # Min observed
    norm_min <- min(data$mean_gene)
    
    ### ALL Normalised WW DATA ####
    ww_data = data %>%
      mutate(norm_mean_gene = mean_gene/norm_min)
    
    ###########################################################
    ##### Deconvolve and Estimate WW Re #####
    
    config_df = expand.grid("region" = c('LA'),  
                            'incidence_var' = c('norm_mean_gene'),
                            'FirstGamma' = 'sim_latent',
                            'SecondGamma' = 'sim_sld' )
    
    
    deconv_ww_data <- data.frame()
    Re_ww <- data.frame()
    mean_si <- 11.00176
    std_si <- 8.070758
    
    for(row_i in 1:nrow(config_df)){
      row_i = 1
      new_deconv_data = deconvolveIncidence(ww_data %>% filter(region == config_df[row_i, 'region']), 
                                            incidence_var = config_df[row_i, 'incidence_var'],
                                            getCountParams(as.character(config_df[row_i, 'FirstGamma'])), 
                                            getCountParams(as.character(config_df[row_i, 'SecondGamma'])),
                                            smooth_param = TRUE, n_boot = 50)
      
      new_deconv_data <- new_deconv_data %>%
        mutate(incidence_var = config_df[row_i, 'incidence_var'])
      
      ##### Get Re #####
      new_Re_ww = getReBootstrap(new_deconv_data, mean_si, std_si)
      new_Re_ww <- new_Re_ww %>%
        mutate(variable = config_df[row_i, 'incidence_var'],
               region = config_df[row_i, 'region'])
      
      deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
      Re_ww = bind_rows(Re_ww, new_Re_ww)
      
    }
  
  
  return(Re_ww)
  
}
plot_huisman_rt <- function(Re_ww, simdata, full_simdata, sim = TRUE) {
  if (sim == TRUE){
    
    
    date_ranges <- Re_ww %>%
      group_by(region) %>%
      summarise(min_date = min(date),
                max_date = max(date)) 
    
    date_length <- length(simdata$time)
    fake_date_one <- ymd("20220110")
    fake_date <- fake_date_one
    for (i in 2:date_length) {
      fake_date <- append(fake_date, fake_date[i-1] + ddays(2))
    }
    
    true_rt <- simdata %>% dplyr::select(time, new_time, true_rt) %>% cbind(date = fake_date)
    initial_conds <- full_simdata %>% mutate(new_time = time - start_time + 7) %>% filter(new_time == 0)
    
    ## Plot Re ####
    
    Re_ww <- Re_ww %>%
      left_join(true_rt, by = "date")
    huisman_Re_plot <- Re_ww %>%
      filter(new_time >= 15) %>%
      mutate(.width = "0.95") %>%
      ggplot() + 
      geom_lineribbon(aes(x = new_time, 
                          y= median_R_mean, 
                          ymin = median_R_lowHPD, 
                          ymax = median_R_highHPD)) + 
      geom_point(data = true_rt, aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept = 1, linetype = "dotted") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw() + 
      theme(legend.position = c(0.6, 0.8), legend.background = element_rect("transparent")) +
      ylab("Rt") +
      xlab("time") +
      ggtitle("Rt wastewater (state-of-art)") +     
      scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
      scale_x_continuous(expand = c(0,0))
    
  }
  
  return(huisman_Re_plot)
}
# generation time discretization functions --------------------------------

zero_epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, (y+1))
  pmf[1] <- pgamma(0.5, alpha, rate = beta)
  for (i in 2:(y+1)) {
    pmf[i] <- pgamma(i -1 +.5, alpha, rate = beta) - pgamma(i -1 -.5, alpha, rate = beta)
  }
  
  pmf
}


epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, y)
  pmf[1] <- phypoexp(1.5, rates)
  for (i in 2:y) {
    pmf[i] <- phypoexp(i+.5, rates) - phypoexp(i-.5, rates)
  }
  
  pmf
}


epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i+.5, meanlog = params[1], sdlog = params[2]) - 
      plnorm(i-.5, meanlog = params[1], sdlog = params[2])
  }
  
  pmf
}

epidemia_weibull <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- pweibull(1.5, shape = params[1], scale = params[2])
  for (i in 2:y) {
    pmf[i] <- pweibull(i+.5, shape = params[1], scale = params[2]) - 
      pweibull(i-.5, shape = params[1], scale = params[2])
  }
  
  pmf
}



# create weekly CA data ---------------------------------------------------
# create data sets from CA data
#values are weekly cases, weekly tests, and the first date of the week
create_weekly_data <- function(ca_data, 
                               county, 
                               start_sunday = "2020-08-02",
                               end_sunday = "2022-01-09"){
  weekly_data <- ca_data %>%
    filter(area == county) %>%
    mutate(week = epiweek(date), 
           year = epiyear(date)) %>%
    group_by(year, week ) %>%
    summarise(total_cases = sum(cases),
              total_tests = sum(total_tests),
              min_date = min(date)) %>%
    ungroup() %>%
    filter(min_date >= start_sunday & min_date <= end_sunday) %>%
    mutate(time = row_number(),
           epidemia_time = time + 1)
  
}


# Create automated process for choosing kappa parameter from a spl --------

run_nb_spline <- function(data, 
                          seed = 17,
                          iter = 4000,
                          warmup = 1000,
                          thin = 10, 
                          refresh = 0,
                          adapt_delta = 0.99) {
  spline_model <- brm(bf(total_cases ~ s(time)),
                      data = data, family = negbinomial(), cores = 4, seed = seed,
                      iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                      control = list(adapt_delta = adapt_delta))
  
  return(spline_model)
}


compare_kappa_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- qlnorm(c(0.025, 0.975),
                                meanlog = candidate_params[1], 
                                sdlog = candidate_params[2])
  
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 + (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}

choose_kappa_params <- function(spline_posterior) {
  posterior_pars <- summary(spline_posterior)
  start_mean <- log(posterior_pars[["spec_pars"]][[1]])
  start_sd <- 0.3
  true_lb <- posterior_pars[["spec_pars"]][[3]]
  true_ub <- posterior_pars[["spec_pars"]][[4]]
  
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  
  optim_params <- optim(par = start_params,
                        fn = compare_kappa_quantiles, 
                        true_quantiles = true_quantiles)
}

# function for using epiestim to choose starting points for log rt --------
get_logrtstart <- function(data,
                           window = 1, 
                           GI_mean = 11.5/7
) {
  
  window = window
  GI_mean = GI_mean
  GI_var = 2*(GI_mean/2)^2
  
  ts <- data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  estimate_R(
    incid = data$total_cases,
    method = "uncertain_si",
    config = make_config(
      list(
        mean_si = GI_mean,
        min_mean_si = 1,
        max_mean_si = GI_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(GI_var),
        min_std_si = sqrt(GI_var)*.8,
        max_std_si = sqrt(GI_var)*1.2,
        n1 = 50,
        n2 = 100, 
        t_start=ts,
        t_end=te
      )
    )
  ) -> ee_outs
  
  ee_quantile <- ee_outs[["R"]] %>%
    dplyr::select(t_start, 
                  rt_mean = `Mean(R)`, 
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start) 
  
  
  log_ee_median <- log(ee_quantile %>% pull(rt_median))
  
  first_one <- log_ee_median[1]
  rt_start <- c(first_one, log_ee_median)
  
  return(rt_start)
}
# fit exp_seed model ------------------------------------------------------
# this is the current version of the model
# assuming generation time is a hypo-expo distribution
# assuming delay time is a gamma distribution, admitting the zero case
fit_estimgamma_model <- function(data,
                                 gen_params,
                                 delay_params,
                                 prev_vals,
                                 log_nu_mean = -2,
                                 log_nu_sd = 0.7,
                                 log_sigma_mean = -0.6,
                                 log_sigma_sd = 0.6,
                                 log_rho_mean,
                                 log_rho_sd,
                                 log_r0_mean = log(1),
                                 log_r0_sd = 0.75,
                                 kappa_mean,
                                 kappa_sd,
                                 seed_exp_rate = 0.3,
                                 init_func,
                                 iterations = 2000,
                                 thin = 2,
                                 adapt_delta = 0.99,
                                 treedepth = 12,
                                 seed = 45,
                                 chain = 4,
                                 gen_dist = "hypo-exp",
                                 delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  
  if (gen_dist == "hypo-exp") {
    gen_weights <- epidemia_hypoexp(data_length, gen_params)
    
  }
  
  
  
  if (gen_dist == "log-normal") {
    gen_weights <- epidemia_lognormal(data_length, gen_params)
  }
  
  if (gen_dist == "weibull") {
    gen_weights <- epidemia_weibull(data_length, gen_params)
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w = gen_weights,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd,
                       seed_exp_rate = seed_exp_rate)
  
  
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = "src/rt_model_gamma_expseed.stan",
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    thin = thin,
                    chain = chain,
                    init = init_func,
                    control = control_list)
  
  return(model_fit)
}

# create rt posterior for real data, no truth involved --------------------
summarise_realdata_rt_estimgamma <- function(stan_posterior,
                                             weekly_data,
                                             start_date,
                                             include_chains){
  
  
  time_week <- weekly_data %>%
    dplyr::select( min_date, max_date, time)
  
  rt_posterior <- stan_posterior %>%
    spread_draws(log_rt[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + 0 + start_date -1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(date, rt) %>%
    group_by(date) %>% 
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    dplyr::select(date, rt_median = rt, .lower, .upper, .width) %>%
    left_join(time_week, by = c("date" = "time"))
  
  return(rt_posterior)
}

# plot_estimgamma_rt ------------------------------------------------------
# first function, creates graphable summary of posteriors
# true_data = data 
# start_date = 1
# good_chains = c(1,4)
summarise_estimgamma_rtposterior <- function(true_data, 
                                posterior, 
                                start_date,
                                good_chains = c(1,2,3,4)) {
  
  truth <- true_data %>%
    dplyr::select(new_time, true_rt)
  
  draws <- posterior %>%
    spread_draws(log_rt[i]) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    filter(.chain %in% good_chains) %>%
    mutate(time = i*7 - 1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(time, rt) %>%
    group_by(time) %>% 
    median_qi(.width = c(0.5, 0.8, 0.95))%>%
    left_join(truth, by = c("time" = "new_time"))
  
  draws
  
}

graph_rt_posterior <- function(true_data,
                               initial_conds,
                               posterior, 
                               start_date) {
  # draws = test
  # true_data = scenario2_weekly_simdata
  # initial_conds = scenario2_initial_conds
  draws <- summarise_estimgamma_rtposterior(true_data, posterior, start_date)

    new_plot <- draws %>% 
      ggplot() +
      geom_lineribbon(aes(x = time, y = rt, ymin = .lower, ymax = .upper,fill = fct_rev(ordered(.width)))) +
      geom_point(data = true_data, mapping = aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept=1, linetype="dashed", color = "black") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw()+
      theme(legend.position = "none") + 
      ylab("Rt") +
      xlab("Time") +
      ggtitle("Rt estimgamma") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
      scale_x_continuous(expand = c(0,0))
    
}

# run negative binomial spline on cases --------

run_nb_spline <- function(data, 
                          seed = 17,
                          iter = 4000,
                          warmup = 1000,
                          thin = 10, 
                          refresh = 0,
                          adapt_delta = 0.99) {
  spline_model <- brm(bf(total_cases ~ s(time)),
                      data = data, family = negbinomial(), cores = 4, seed = seed,
                      iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                      control = list(adapt_delta = adapt_delta))
  
  return(spline_model)
}


# create weekly CA data ---------------------------------------------------
# create data sets from CA data
#values are weekly cases, weekly tests, and the first date of the week
create_weekly_data <- function(ca_data, 
                               county, 
                               start_sunday = "2020-08-02",
                               end_sunday = "2022-01-09"){
  weekly_data <- ca_data %>%
    filter(area == county) %>%
    mutate(week = epiweek(date), 
           year = epiyear(date)) %>%
    group_by(year, week ) %>%
    summarise(total_cases = sum(cases),
              total_tests = sum(total_tests),
              min_date = min(date)) %>%
    ungroup() %>%
    filter(min_date >= start_sunday & min_date <= end_sunday) %>%
    mutate(time = row_number(),
           epidemia_time = time + 1)
  
}




# rt_metrics --------------------------------------------------------------
# used for calculating frequentist characteristics of inference for rt
rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - true_rt),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag(true_rt),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs(true_rt - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}


# freqentist metrics ---------------------------------------------------------------

# used for calculating frequentist characteristics of inference for rt
freq_metrics<- function(data, value, true_value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - {{ true_value }}),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = {{ true_value }} >= {{ lower }} & {{ true_value }} <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag({{ true_value }}),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs({{ true_value }} - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}






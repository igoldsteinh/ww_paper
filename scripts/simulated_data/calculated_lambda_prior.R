# calculate priors for normalized lambda 
# by simulating outbreak under different infectious/shedding but recovered period distributions
# and then fitting glmnet to each outbreak
# use the priors of our model to draw the different means of the distributions

library(tidyverse)
library(greybox)
source("src/simulate_stochastic_seirr.R")

library(glmnet)

num_sims = 1000
counts = vector(mode='list', length=num_sims)
nozeros = vector(mode = 'list', length = num_sims)
x = vector(mode = 'list', length = num_sims)
y = vector(mode = 'list', length = num_sims)
loglogx = vector(mode = 'list', length = num_sims)
loglogy = vector(mode = 'list', length = num_sims)

set.seed(1234)
gamma_vec = rlnorm(num_sims, log(1/4), 0.2)
nu_vec = rlnorm(num_sims,log(1/7), 0.2)
eta_vec = rlnorm(num_sims,log(1/18), 0.2)
for (i in 1:1000) {
  gamma = gamma_vec[i]
  nu = nu_vec[i]
  eta = eta_vec[i]
  
  counts[[i]] = sim_and_calc_genes(N = 1000, 
                                   I_init = 5, 
                                   beta = 2 * nu, 
                                   gamma = gamma, 
                                   nu = nu, 
                                   eta = eta, 
                                   num_sims = 1)[["counts"]]
}

for (i in 1:1000) {
  nozeros[[i]] = counts[[i]] %>% filter(total_count > 0) %>% mutate(logcounts = log(total_count),
                                                                    logI = log(num_I), 
                                                                    logR1 = log(num_R1))
  x[[i]] = as.matrix(nozeros[[i]] %>% dplyr::select(num_I, num_R1))
  y[[i]] = nozeros[[i]] %>% pull(total_count)
}

for (i in 1:1000) {
  loglogx[[i]] = as.matrix(nozeros[[i]] %>% filter(num_I > 0 & num_R1 > 0) %>% dplyr::select(logI, logR1))
  loglogy[[i]] = as.matrix(nozeros[[i]] %>% filter(num_I > 0 & num_R1 >0) %>% pull(logcounts))
  
}
          
write_rds(nozeros, here::here("data", "sim_data", "lambda_sims.rds"))


# create average concentrations -------------------------------------------
popsize = 1000
num_sims = 1000
nozeros <- read_rds( here::here("data", "sim_data", "lambda_sims.rds"))
x = vector(mode = 'list', length = num_sims)
y = vector(mode = 'list', length = num_sims)

for (i in 1:1000) {
  x[[i]] = as.matrix(nozeros[[i]] %>% dplyr::select(num_I, num_R1))
  y[[i]] = nozeros[[i]] %>% pull(total_count)/popsize
}


# calculate coef for linear relationship ----------------------------------
linear_coefs = map2(x,y, ~coef(glmnet(.x, .y, lambda = 0, lower.limits = 0, intercept = FALSE))[2:3,1])



linear_combined_coefs = bind_rows(linear_coefs) %>%
  mutate(total_num = num_I + num_R1,
         normalized_I = num_I/total_num)

linear_true_quantiles = quantile(linear_combined_coefs$normalized_I, probs = c(0.025, 0.975))

# calculate coef for logY = I + R1 ----------------------------------
log_coefs = map(nozeros, ~summary(lm(logcounts ~ num_I + num_R1 - 1, data =.x))$coefficients[1:2, 1])

log_combined_coefs = bind_rows(log_coefs) %>%
  mutate(total_num = num_I + num_R1,
         normalized_I = num_I/total_num)

log_true_quantiles = quantile(log_combined_coefs$normalized_I, probs = c(0.025, 0.975))

# calculate coef for logY = logI + logR1 ----------------------------------
loglog_coefs = map2(loglogx, loglogy, ~coef(glmnet(.x, .y, lambda = 0, lower.limits = 0, intercept = FALSE))[2:3,1])

loglog_combined_coefs = bind_rows(loglog_coefs) %>%
  mutate(total_num = logI + logR1,
         normalized_I = logI/total_num)

loglog_true_quantiles = quantile(loglog_combined_coefs$normalized_I, probs = c(0.025, 0.975))

# functions for calculating priors ----------------------------------------

compare_lambda_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- qlogitnorm(c(0.025, 0.975),
                                mu = candidate_params[1], 
                                sigma = candidate_params[2])
  
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 + (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}

choose_lambda_params <- function(true_quantiles, start_mean, start_sd) {
  true_lb <- true_quantiles[1]
  true_ub <- true_quantiles[2]
  
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  
  optim_params <- optim(par = start_params,
                        fn = compare_lambda_quantiles, 
                        true_quantiles = true_quantiles)
}


# find linear lambda prior -------------------------------------------------------
start_mean <- qlogis(0.98)
start_sd <- 0.25

linear_lambda_params <- choose_lambda_params(linear_true_quantiles, start_mean, start_sd)
# mu = 5.685528
# sigma = 2.178852

linear_lambda_prior <- plogis(rnorm(10000, 5.685528, 2.178852))
quantile(linear_lambda_prior, c(0.025, 0.5, 0.975))

quantile(linear_combined_coefs$normalized_I, probs = c(0.025, 0.5, 0.975))

# find log lambda prior ---------------------------------------------------
start_mean <- qlogis(0.4)
start_sd <- 0.25

log_lambda_params <- choose_lambda_params(log_true_quantiles, start_mean, start_sd)
# mu = -0.507017
# sigma = 0.4278839

log_lambda_prior <- plogis(rnorm(10000, -0.507017, 0.4278839))
quantile(log_lambda_prior, c(0.025, 0.5, 0.975))

# find loglog lambda prior ------------------------------------------------

start_mean <- qlogis(0.4)
start_sd <- 0.25

loglog_lambda_params <- choose_lambda_params(loglog_true_quantiles, start_mean, start_sd)
# mu = -0.7956065
# sigma = 0.6253220

loglog_lambda_prior <- plogis(rnorm(10000, -0.7956065 , 0.6253220))
quantile(log_lambda_prior, c(0.025, 0.5, 0.975))


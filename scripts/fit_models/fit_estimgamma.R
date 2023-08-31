# fit rt estim gamma to real data
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
library(lubridate)
source("src/wastewater_functions.R")
source("src/rt_estim_gamma_functions.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


data <- read_csv(here::here("data", "LA_EIR_data.csv"))
data <- data %>%
        rename("time" = "new_week")


# run spline for kappa priors
county_spline <- run_nb_spline(data)

# calculate kappa priors
county_kappa <- choose_kappa_params(county_spline)

# calculate quantile for tests
test_quantile <- quantile(data$total_tests)

# fit model
data_length = dim(data)[1]

# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(data)

#next choose incidence starting points
incid_start <- 1/0.2 * data$total_cases

init_func <- function() list(log_incid_rate_raw = 0,
                             log_rt0_raw = 0,
                             rho = 0.2/test_quantile[2],
                             kappa = county_kappa$par[1],
                             seed_incid_one_raw =1,
                             incid = incid_start,
                             log_rt = logrt_start)

county_posterior <- fit_estimgamma_model(data,
                                         gen_params = c(log(7.872346) + log(1/7), 
                                                        0.642713),
                                         delay_params = c(4.05, 7*0.74),
                                         prev_vals = 4,
                                         log_rho_mean = log(0.066/test_quantile[2]),
                                         log_rho_sd = 0.3,
                                         kappa_mean = county_kappa$par[1],
                                         kappa_sd = county_kappa$par[2],
                                         iterations = 6000,
                                         init_func = init_func,
                                         gen_dist = "log-normal",
                                         seed = 225,
                                         thin = 3)
file_name <- "LA_estimgamma_posterior.rds"
write_rds(county_posterior, here::here("results", "estimgamma", file_name))


# alternate seed distr ----------------------------------------------------
la_cases <-read_csv("data/covid19cases_test.csv") %>% 
  mutate(case_date = as.Date("2021-07-04") - days(28)) %>%
  filter(area == "Los Angeles") %>%
  filter(date < as.Date("2021-07-04") & date >= case_date) %>% pull(cases)


plot_la_cases <-read_csv("data/covid19cases_test.csv") %>% 
  mutate(case_date = as.Date("2021-07-04") - days(28)) %>%
  filter(area == "Los Angeles") %>%
  filter(date < as.Date("2021-09-01") & date >= case_date)

looksee <- plot_la_cases %>% 
         ggplot(aes(x = date, y = cases)) +
         geom_point() + 
        theme_bw()


rough_guess = sum(la_cases) * (1/0.2)/4
rough_param = 1/rough_guess
alt_county_posterior <- fit_estimgamma_model(data,
                                         gen_params = c(log(7.872346) + log(1/7), 
                                                        0.642713),
                                         delay_params = c(4.05, 7*0.74),
                                         prev_vals = 4,
                                         log_rho_mean = log(0.2/test_quantile[2]),
                                         log_rho_sd = 0.3,
                                         log_r0_mean = log(2),
                                         log_r0_sd = 0.1,
                                         kappa_mean = county_kappa$par[1],
                                         kappa_sd = county_kappa$par[2],
                                         seed_exp_rate = rough_param,
                                         iterations = 6000,
                                         init_func = init_func,
                                         gen_dist = "log-normal",
                                         seed = 225,
                                         thin = 3)
file_name <- "LA_estimgamma_posterior_alt.rds"
write_rds(alt_county_posterior, here::here("results", "estimgamma", file_name))

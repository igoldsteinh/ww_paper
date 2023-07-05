# fit epidemia to real LA data
library(tidyverse)
library(tidybayes)
library(patchwork)
library(epidemia)
library(rstanarm)
library(lubridate)
library(rstan)
library(GGally)
library(scales)
library(gridExtra)
library(cowplot)
library(posterior)
source("src/epidemia_functions.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# read in data ------------------------------------------------------------
dat <- read_csv(here::here("data", "LA_EIR_data.csv"))

dat$start_date <- as.Date(dat$min_date, format = "%m/%d/%y")



# let's try it with epidemia ----------------------------------------------
data_length <- dim(dat)[1]

gen_params = c(log(7.872346) + log(1/7), 
               0.642713)
delay_params = c(4.05, 7*0.74)

gen_weights <- epidemia_lognormal(10, gen_params)

sum <- sum(gen_weights)

gen_weights <- gen_weights/ sum

delay_weights <- epidemia_gamma(10, delay_params[1], delay_params[2])

delay_sum <- sum(delay_weights)

delay_weights <- delay_weights/delay_sum


date <- seq(from = ymd('2020-07-04'), length.out = length(dat$total_cases) + 1, by = "days")

data <- data.frame(
  city = "LA",
  cases = c(NA, dat$total_cases),
  date = date, 
  day = weekdays(date)
)


rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.1),
  prior_intercept = normal(log(2), 0.1),
  link = 'log'
)

obs <-  epiobs(
  formula = cases ~ 1,
  prior_intercept = rstanarm::normal(location=0.066, scale=0.05),
  link = "identity",
  i2o = delay_weights
)


args <- list(
  rt = rt,
  inf = epiinf(gen = gen_weights),
  obs = obs,
  data = data,
  iter = 5e3,
  seed = 225
)


#to redo with incidence as parameter
args$inf <- epiinf(gen = gen_weights, 
                   latent=TRUE, 
                   prior_aux = normal(10,2))
fm2 <- do.call(epim, args)

posteriors <- as.data.frame(fm2)

subset_samples <- as_draws_df(fm2[["stanfit"]])

mcmc_summary <- summarise_draws(subset_samples)


write_rds(fm2, here::here("results", "estimnormal", "estimnormal_LA_fit.rds"))

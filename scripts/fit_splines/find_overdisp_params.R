# find overdispersion priors
# find parameters for overdispersion priors
# based on posteriors from NB spline
library(tidyverse)
library(rstan)
library(brms)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("src/spline_functions.R")
sim = 1

if (sim == 2) {
  data <- read_csv("data/sim_data/scenario2_overdispdata.csv") 
  
}

if (sim == 9) {
  data <- read_csv("data/sim_data/scenario9_lump7data.csv") 
  
}

if (sim == 1) {
  data <- read_csv("data/sim_data/scenario1_overdispdata.csv") 
  
}

if (sim == 3) {
  data <- read_csv("data/sim_data/scenario3_lump7data.csv") 
  
}

if (sim == 111) {
  data <- read_csv("data/sim_data/scenario111_overdispdata.csv") 
  
}

if (sim == "weekly") {
  data <- read_csv(here::here("data", "sim_data", "test_weekly_data.csv")) 
}

if (sim == "real") {
  data <- read_csv(here::here("data", "LA_EIR_data.csv"))
}
# find case priors ---------------------------------------------
set.seed(225)
case_spline <- run_nb_spline(data = data,
                             response = "cases")

case_overdisp <- choose_kappa_params(case_spline)

case_params <- case_overdisp$par
# save the results --------------------------------------------------------

labels <- c("case")

priors <- t(case_params) %>%
  cbind(labels)

priors <- data.frame(priors)

colnames(priors) <- c("mean", "sd", "labels")
rownames(priors) <- NULL

file_name <- paste("overdisp_priors_sim", sim, ".csv", sep = "")
write_csv(priors, paste("data/sim_data/", file_name, sep = ""))

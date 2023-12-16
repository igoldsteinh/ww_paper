### Process Results eirrc_closed
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
source("src/wastewater_functions.R")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  sim = 102
  seed = 1
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}

dir_create(path("results", "eirrc_closed", "mcmc_summaries"))

priors_only = sim == 0

# a hack for dealing with re-using scenario 1 data
repeat_scenario1s = c(3,4,5,6,7,8,9,10, 12, "frw", 102)
repeat_scenario41s = c(411, 4111)

if (sim %in% repeat_scenario1s) {
  scenario_sim = 1
} else if (sim %in% repeat_scenario41s) {
  scenario_sim = 41
} else {
  scenario_sim = sim
}

# priors only -------------------------------------------------------------
if(priors_only == TRUE) {

  priorname <- paste0("prior_generated_quantities_scenario", sim, "_seed", seed, ".csv")
  prior_gq_samples_all <- read_csv(here::here("results",
                                              "eirrc_closed",
                                              priorname)) %>%
    pivot_longer(-c(iteration, chain)) %>%
    select( name, value)
  
  
  priors <- make_prior_samples(prior_gq_samples_all)
  
  prior_quantiles <- make_prior_quantiles(priors)
  
  prior_timevarying_quantiles <- make_timevarying_prior_quantiles(prior_gq_samples_all)
  
  prior_samp_name <- paste0("prior_samples_scenario", sim , "_seed", seed, ".csv")
  prior_timevarying_name <- paste0("prior_timevaryingquantiles_scenario", sim, "_seed", seed, ".csv")
  prior_quant_name <- paste0("prior_quantiles_scenario", sim, "_seed", seed, ".csv")
  write_csv(priors, here::here("results", "eirrc_closed", prior_samp_name))
  write_csv(prior_quantiles, here::here("results", "eirrc_closed", prior_quant_name))
  write_csv(prior_timevarying_quantiles, here::here("results", "eirrc_closed", prior_timevarying_name))
  
  quit()
}

# posterior ---------------------------------------------------------------

# calculate MCMC diagnostics after burnin
gq_address <- paste0("results/eirrc_closed/generated_quantities/generated_quantities_scenario", 
                     sim, 
                     "_seed", 
                     seed,
                     ".csv")

posterior_samples <- read_csv(gq_address) %>%
  rename(.iteration = iteration,
         .chain = chain) %>%
  as_draws()

max_iteration = max(posterior_samples$.iteration)
min_iteration = round(max_iteration/2)


subset_samples <- subset_draws(posterior_samples)

mcmc_summary <- summarise_draws(subset_samples)

mcmc_summary_address <- paste0("results/eirrc_closed/mcmc_summaries/mcmc_summary_scenario", 
                               sim, 
                               "_seed",
                               seed,
                               ".csv")
write_csv(mcmc_summary, mcmc_summary_address)

# lp trace plot -----------------------------------------------------------
lp_df <-read_csv(here::here("results", "eirrc_closed", "generated_quantities", paste0("posterior_df_scenario",
                                                                                   sim, 
                                                                                   "_seed",
                                                                                   seed,
                                                                                   ".csv"))) 


trace_plot <- lp_df %>%
  ggplot(aes(x = iteration, y = lp, color = as.factor(chain))) + 
  geom_line() +
  theme_bw() + 
  ggtitle("EIRR LP Trace")

ggsave(here::here("results", "eirrc_closed", "mcmc_summaries", paste0("trace_scenario", sim, "_seed", seed, ".png" )), trace_plot, width = 5, height = 5)

# create long format fixed samples and time-varying quantiles -----------------
posterior_gq_samples_all <- subset_samples  %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  select( name, value)


posterior_fixed_samples <- make_fixed_posterior_samples(posterior_gq_samples_all)

fixed_samples_address <- paste0("results/eirrc_closed/generated_quantities/posterior_fixed_samples_scenario",
                                sim, 
                                "_seed",
                                seed,
                                ".csv")

write_csv(posterior_fixed_samples, fixed_samples_address)

rm(posterior_fixed_samples)

posterior_timevarying_quantiles <- make_timevarying_posterior_quantiles(posterior_gq_samples_all)

timevarying_quantiles_address <- paste0("results/eirrc_closed/generated_quantities/posterior_timevarying_quantiles_scenario",
                                        sim,
                                        "_seed",
                                        seed,
                                        ".csv")

write_csv(posterior_timevarying_quantiles, timevarying_quantiles_address)

rm(posterior_timevarying_quantiles)

rm(posterior_gq_samples_all)

# create posterior predictive quantiles -----------------------------------
post_pred_address <- paste0("results/eirrc_closed/posterior_predictive/posterior_predictive_scenario",
                            sim,
                            "_seed",
                            seed,
                            ".csv")
eirr_post_pred <- read_csv(post_pred_address)

if (sim != "real" & sim != "uci" & sim!= "uci_region1" & sim != "uci_region2" & sim != "uci_region3" & sim != 101) {
  if (sim != 12 & scenario_sim != 41) {
    simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_fitted_genecount_obsdata.csv")
    
    simdata <- read_csv(simdata_address) %>%
      dplyr::filter(seed == seed)
    
  } else if (sim == 12) {
    simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_lump2data_seiirr_100sims.csv")
    
    simdata <- read_csv(simdata_address) %>%
      dplyr::filter(seed == seed)
    
  } else if (scenario_sim == 41) {
    simdata <- read_csv(here::here("data", "sim_data", "scenario41_fitted_obsdata.csv")) 
  }

} else if (sim == "real") {
  simdata_address <- "data/LA_daily_data_feb2022.csv"

  simdata <- read_csv(simdata_address)

} else if (sim == "uci") {

  simdata <- read_csv(here::here("data", "uci_data", "uci_fitting_data.csv"))
  
} else if (sim == "uci_region1") {
  
  simdata <- read_csv(here::here("data", "uci_data", "uci_fitting_data_region1.csv"))
  
} else if (sim == "uci_region2") {
  
  simdata <- read_csv(here::here("data", "uci_data", "uci_fitting_data_region2.csv"))
  
} else if (sim == "uci_region3") {
  
  simdata <- read_csv(here::here("data", "uci_data", "uci_fitting_data_region3.csv"))
  
} else if (sim == 101) {
  simdata <- read_csv(here::here("data", "sim_data", paste0("scenario101_seed", seed, "_fitted_genecount_obsdata.csv")))
}

if (sim == 3) {
  ten_sim_val = TRUE
} else {
  ten_sim_val = FALSE
}

if (sim == 4) {
  three_mean_val = TRUE
} else {
  three_mean_val = FALSE
}
eirr_post_pred_intervals <- make_post_pred_intervals(eirr_post_pred, simdata, ten_sim = ten_sim_val, three_mean = three_mean_val)

post_pred_interval_address <- paste0("results/eirrc_closed/posterior_predictive/posterior_predictive_intervals_scenario",
                              sim,
                              "_seed",
                              seed,
                              ".csv")

write_csv(eirr_post_pred_intervals, post_pred_interval_address)


### Process Results seir cases
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
source("src/wastewater_functions.R")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  sim = 1
  seed = 2
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}

dir_create(path("results", "seir_cases", "mcmc_summaries"))

repeat_scenario1s = c(11,13, 14, 1111)
if (sim %in% repeat_scenario1s) {
  scenario_sim = 1
} else {
  scenario_sim = sim
}

priors_only = sim == 0
# priors only -------------------------------------------------------------
if(priors_only == TRUE) {
  # sim = 1
  priorname <- paste0("prior_generated_quantities_scenario", sim, ".csv")
  prior_gq_samples_all <- read_csv(here::here("results",
                                              "seir_cases",
                                              priorname)) %>%
    pivot_longer(-c(iteration, chain)) %>%
    select( name, value)
  
  
  priors <- make_prior_samples(prior_gq_samples_all)
  
  prior_quantiles <- make_prior_quantiles(priors)
  
  prior_timevarying_quantiles <- make_timevarying_prior_quantiles(prior_gq_samples_all)
  
  prior_samp_name <- paste0("prior_samples_scenario", sim , ".csv")
  prior_timevarying_name <- paste0("prior_timevaryingquantiles_scenario", sim, ".csv")
  prior_quant_name <- paste0("prior_quantiles_scenario", sim , ".csv")
  write_csv(priors, here::here("results", "seir_cases", prior_samp_name))
  write_csv(prior_quantiles, here::here("results", "seir_cases", prior_quant_name))
  write_csv(prior_timevarying_quantiles, here::here("results", "seir_cases", prior_timevarying_name))
  
  quit()
}

# calculate MCMC diagnostics after burnin
gq_address <- paste0("results/seir_cases/generated_quantities/generated_quantities_scenario", 
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

mcmc_summary_address <- paste0("results/seir_cases/mcmc_summaries/mcmc_summary_scenario", 
                               sim,
                               "_seed", 
                               seed,
                               ".csv")
write_csv(mcmc_summary, mcmc_summary_address)


# lp trace plot -----------------------------------------------------------
lp_df <-read_csv(here::here("results", "seir_cases", "generated_quantities", paste0("posterior_df_scenario",
                                                                                      sim, 
                                                                                      "_seed",
                                                                                      seed,
                                                                                      ".csv"))) 


trace_plot <- lp_df %>%
  ggplot(aes(x = iteration, y = lp, color = as.factor(chain))) + 
  geom_line() +
  theme_bw() + 
  ggtitle("SEIR Cases LP Trace")

ggsave(here::here("results", "seir_cases", "mcmc_summaries", paste0("trace_scenario", sim, "_seed", seed, ".png" )), trace_plot, width = 5, height = 5)


# create long format fixed samples and time-varying quantiles -----------------
posterior_gq_samples_all <- subset_samples %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  select(.chain, .iteration, name, value)

  
posterior_fixed_samples <- make_fixed_posterior_samples(posterior_gq_samples_all)

fixed_samples_address <- paste0("results/seir_cases/generated_quantities/posterior_fixed_samples_scenario",
                                sim,
                                "_seed", 
                                seed,
                                ".csv")

write_csv(posterior_fixed_samples, fixed_samples_address)

rm(posterior_fixed_samples)

posterior_timevarying_quantiles <- make_timevarying_posterior_quantiles(posterior_gq_samples_all)


timevarying_quantiles_address <- paste0("results/seir_cases/generated_quantities/posterior_timevarying_quantiles_scenario",
                                        sim,
                                        "_seed", 
                                        seed,
                                        ".csv")

write_csv(posterior_timevarying_quantiles, timevarying_quantiles_address)

rm(posterior_timevarying_quantiles)

rm(posterior_gq_samples_all)


# create posterior predictive quantiles -----------------------------------
post_pred_address <- paste0("results/seir_cases/posterior_predictive/posterior_predictive_scenario",
                            sim,
                            "_seed", 
                            seed,
                            ".csv")
seir_post_pred <- read_csv(post_pred_address)

simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_fitted_cases_obsdata.csv")

simdata <- read_csv(simdata_address) %>%
  dplyr::filter(seed == seed)

seir_post_pred_intervals <- make_post_pred_intervals(seir_post_pred, simdata, cases = TRUE)

post_pred_interval_address <- paste0("results/seir_cases/posterior_predictive/posterior_predictive_intervals_scenario",
                                     sim,
                                     "_seed", 
                                     seed,
                                     ".csv")

write_csv(seir_post_pred_intervals, post_pred_interval_address)
 

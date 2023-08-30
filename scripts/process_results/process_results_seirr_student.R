### Process Results SEIRR_student
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
source("src/wastewater_functions.R")

args <- commandArgs(trailingOnly=TRUE)
sim <- as.integer(args[1])

dir_create(path("results", "seirr_student", "mcmc_summaries"))

if (length(args) == 0) {
  sim = "own"
  seed = 2
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}

priors_only = sim == 0

# a hack for dealing with re-using scenario 1 data
repeat_scenario1s = c(3, "alt_prior", "fixed", "frw")
if (sim %in% repeat_scenario1s) {
  scenario_sim = 1
} else {
  scenario_sim = sim
}

# priors only -------------------------------------------------------------
if(priors_only == TRUE) {

  priorname <- paste0("prior_generated_quantities_scenario", sim, ".csv")
  prior_gq_samples_all <- read_csv(here::here("results",
                                                        "seirr_student",
                                                        priorname)) %>%
    pivot_longer(-c(iteration, chain)) %>%
    select( name, value)


  priors <- make_prior_samples(prior_gq_samples_all)
  
  prior_quantiles <- make_prior_quantiles(priors)
  
  prior_timevarying_quantiles <- make_timevarying_prior_quantiles(prior_gq_samples_all)
  
  prior_samp_name <- paste0("prior_samples_scenario", sim , ".csv")
  prior_timevarying_name <- paste0("prior_timevaryingquantiles_scenario", sim, ".csv")
  prior_quant_name <- paste0("prior_quantiles_scenario", sim , ".csv")
  write_csv(priors, here::here("results", "seirr_student", prior_samp_name))
  write_csv(prior_quantiles, here::here("results", "seirr_student", prior_quant_name))
  write_csv(prior_timevarying_quantiles, here::here("results", "seirr_student", prior_timevarying_name))
  
  quit()
}

# posterior ---------------------------------------------------------------
# calculate MCMC diagnostics after burnin
gq_address <- paste0("results/seirr_student/generated_quantities/generated_quantities_scenario", 
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

mcmc_summary_address <- paste0("results/seirr_student/mcmc_summaries/mcmc_summary_scenario", 
                               sim, 
                               "_seed",
                               seed,
                               ".csv")
write_csv(mcmc_summary, mcmc_summary_address)
# lp trace plot -----------------------------------------------------------
lp_df <-read_csv(here::here("results", "seirr_student", "generated_quantities", paste0("posterior_df_scenario",
                                                                                     sim, 
                                                                                     "_seed",
                                                                                     seed,
                                                                                     ".csv"))) 


trace_plot <- lp_df %>%
  ggplot(aes(x = iteration, y = lp, color = as.factor(chain))) + 
  geom_line() +
  theme_bw() + 
  ggtitle("SEIRR LP Trace")

ggsave(here::here("results", "seirr_student", "mcmc_summaries", paste0("trace_scenario", sim, "_seed", seed, ".png" )), trace_plot, width = 5, height = 5)

# create long format fixed samples and time-varying quantiles -----------------
posterior_gq_samples_all <- subset_samples %>%
  pivot_longer(-c(.iteration, .chain, .draw)) %>%
  select( name, value)

  
posterior_fixed_samples <- make_fixed_posterior_samples(posterior_gq_samples_all)

fixed_samples_address <- paste0("results/seirr_student/generated_quantities/posterior_fixed_samples_scenario",
                                sim, 
                                "_seed",
                                seed,
                                ".csv")

write_csv(posterior_fixed_samples, fixed_samples_address)

rm(posterior_fixed_samples)

posterior_timevarying_quantiles <- make_timevarying_posterior_quantiles(posterior_gq_samples_all)


timevarying_quantiles_address <- paste0("results/seirr_student/generated_quantities/posterior_timevarying_quantiles_scenario",
                                        sim, 
                                        "_seed",
                                        seed,
                                        ".csv")

write_csv(posterior_timevarying_quantiles, timevarying_quantiles_address)

rm(posterior_timevarying_quantiles)

rm(posterior_gq_samples_all)

# create posterior predictive quantiles -----------------------------------
post_pred_address <- paste0("results/seirr_student/posterior_predictive/posterior_predictive_scenario",
                            sim, 
                            "_seed",
                            seed,
                            ".csv")
seirr_post_pred <- read_csv(post_pred_address)

simdata_address <-  paste0("data/sim_data/scenario", scenario_sim, "_fitted_genecount_obsdata.csv")

simdata <- read_csv(simdata_address) %>%
  dplyr::filter(seed == seed)

if (sim == 3) {
  ten_sim_val = TRUE
} else {
  ten_sim_val = FALSE
}

seirr_post_pred_intervals <- make_post_pred_intervals(seirr_post_pred, simdata, ten_sim = ten_sim_val)

post_pred_interval_address <- paste0("results/seirr_student/posterior_predictive/posterior_predictive_intervals_scenario",
                              sim, 
                              "_seed",
                              seed,
                              ".csv")

write_csv(seirr_post_pred_intervals, post_pred_interval_address)


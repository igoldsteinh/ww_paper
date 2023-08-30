### Summarise processed results seir_cases
library(tidyverse)
library(tidybayes)
library(posterior)
library(fs)
library(gridExtra)
library(ggplot2)
library(scales)
library(cowplot)
source("src/wastewater_functions.R")

args <- commandArgs(trailingOnly=TRUE)


if (length(args) == 0) {
  sim = 1
} else {
  sim <- as.integer(args[1])
}

repeat_scenario1s = c(11,13, 14, 1111)
if (sim %in% repeat_scenario1s) {
  scenario_sim = 1
} else {
  scenario_sim = sim
}

# read in data and results ---------------------------------------------------------
num_sims = 100
seed = 1:num_sims
full_stan_diag <- map(seed, ~read_csv(paste0("results/seir_cases/mcmc_summaries/mcmc_summary_scenario", 
                                                             sim, 
                                                             "_seed",
                                                             .x,
                                                             ".csv"))) %>%
  bind_rows(.id = "seed") %>%
  group_by(seed) %>% 
  filter(variable != "R[0]") %>%
  summarise(min_rhat = min(rhat),
            max_rhat = max(rhat),
            min_ess_bulk = min(ess_bulk),
            max_ess_bulk = max(ess_bulk),
            min_ess_tail = min(ess_tail),
            max_ess_tail = max(ess_tail))

write_csv(full_stan_diag, here::here("results", "seir_cases", paste0("seir_cases_scenario",sim, "_allseeds_stan_diag.csv")))
# create final rt frame ---------------------------------------------------
full_simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_full_genecount_obsdata.csv")

full_simdata <- read_csv(full_simdata_address) %>% filter(seed == 1) %>% rename("true_rt" = "Rt",
                                                                                "r0" = "R0")

# the data set is simulated for a certain period of time
# but then based on how we choose to space apart observations (every two days, every seven etc)
# there is a max observed time in the data set, we should not judge the model beyond the fitted data (for now)
# this time should be the same across models (it will be the same for the case models even though they're slightly different)
gene_fitted_simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_fitted_genecount_obsdata.csv")

gene_fitted_simdata <- read_csv(gene_fitted_simdata_address)
max_time <- max(gene_fitted_simdata$time)

timevarying_quantiles <- map(seed, ~read_csv(paste0("results/seir_cases/generated_quantities/posterior_timevarying_quantiles_scenario",
                                                                    sim,
                                                                    "_seed",
                                                                    .x,
                                                                    ".csv"))) 

rt_quantiles <- timevarying_quantiles %>%
                map(~.x %>% filter(name == "rt_t_values") %>%
                      rename(week = time) %>%
                      right_join(full_simdata, by = "week") %>%
                dplyr::select(week, time, true_rt, value, .lower, .upper, .width,.point, .interval) %>%
                filter(time <= max_time,
                       week >= 0)
                ) %>%
  bind_rows(.id = "seed")

r0_quantiles <- timevarying_quantiles %>%
  map(~.x %>% filter(name == "r0_t_values") %>%
        rename(week = time) %>%
        right_join(full_simdata, by = "week") %>%
        dplyr::select(week, time, r0, value, .lower, .upper, .width,.point, .interval) %>%
        filter(time <= max_time,
               week >= 0)
  ) %>%
  bind_rows(.id = "seed")

I_quantiles <- timevarying_quantiles %>%
  map(~.x %>% filter(name == "I") %>%
        rename(week = time) %>%
        right_join(full_simdata, by = "week") %>%
        dplyr::select(week, time, num_I, value, .lower, .upper, .width,.point, .interval) %>%
        filter(time <= max_time,
               week >= 0)
  ) %>%
  bind_rows(.id = "seed")

write_csv(rt_quantiles, here::here("results", "seir_cases", paste0("seir_cases_scenario", sim,  "_allseeds_rt_quantiles.csv")))
write_csv(rt_quantiles, here::here("results", "seir_cases", paste0("seir_cases_scenario", sim,  "_allseeds_r0_quantiles.csv")))
write_csv(I_quantiles, here::here("results", "seir_cases", paste0("seir_cases", sim,  "_allseeds_prevI_quantiles.csv")))

# visualize results
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

make_rt_plot <- function(seed_val) {
  rt_quantiles %>%
    filter(seed == seed_val) %>%
    ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    geom_point(aes(time, true_rt), color = "coral1") + 
    scale_y_continuous("Rt", label = comma) +
    scale_x_continuous(name = "Time") +
    ggtitle(str_c("SEIR Cases Posterior Rt Scenario ", sim, " Seed ", seed_val)) +
    my_theme
}

make_r0_plot <- function(seed_val) {
  r0_quantiles %>%
    filter(seed == seed_val) %>%
    ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    geom_point(aes(time, r0), color = "coral1") + 
    scale_y_continuous("R0", label = comma) +
    scale_x_continuous(name = "Time") +
    ggtitle(str_c("SEIR Cases Posterior Rt Scenario ", sim, " Seed ", seed_val)) +
    my_theme
}

ggsave2(filename = here::here("results", "seir_cases", paste0("seir_cases_rt_plots_scenario", sim, ".pdf")),
        plot = rt_quantiles %>%
          distinct(seed) %>%
          arrange(seed) %>%
          pull(seed) %>%
          map(make_rt_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)

ggsave2(filename = here::here("results", "seir_cases", paste0("seir_cases_r0_plots_scenario", sim, ".pdf")),
        plot = r0_quantiles %>%
          distinct(seed) %>%
          arrange(seed) %>%
          pull(seed) %>%
          map(make_r0_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)



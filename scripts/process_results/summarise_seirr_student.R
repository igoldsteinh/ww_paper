### Summarise processed results seirr_student
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

scenario_sim = sim
# read in data and results ---------------------------------------------------------
num_sims = 100
seed = 1:num_sims
full_stan_diag <- map(seed, ~read_csv(paste0("results/seirr_student/mcmc_summaries/mcmc_summary_scenario", 
                                                             sim, 
                                                             "_seed",
                                                             .x,
                                                             ".csv"))) %>%
  bind_rows(.id = "seed") %>%
  group_by(seed) %>% 
  filter(variable != "R2[0]") %>%
  summarise(min_rhat = min(rhat),
            max_rhat = max(rhat),
            min_ess_bulk = min(ess_bulk),
            max_ess_bulk = max(ess_bulk),
            min_ess_tail = min(ess_tail),
            max_ess_tail = max(ess_tail))

write_csv(full_stan_diag, here::here("results", "seirr_student", paste0("seirr_scenario", sim, "_allseeds_stan_diag.csv")))
# create final rt frame ---------------------------------------------------
full_simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_full_genecount_obsdata.csv")

full_simdata <- read_csv(full_simdata_address) %>% filter(seed == 1) %>% rename("true_rt" = "Rt")

print(head(full_simdata))

# the data set is simulated for a certain period of time
# but then based on how we choose to space apart observations (every two days, every seven etc)
# there is a max observed time in the data set, we should not judge the model beyond the fitted data (for now)
# this time should be the same across models (it will be the same for the case models even though they're slightly different)
fitted_simdata_address <-  paste0("data/sim_data/scenario", scenario_sim, "_fitted_genecount_obsdata.csv")

fitted_simdata <- read_csv(fitted_simdata_address)
max_time <- max(fitted_simdata$time)

timevarying_quantiles <- map(seed, ~read_csv(paste0("results/seirr_student/generated_quantities/posterior_timevarying_quantiles_scenario",
                                                                    sim,
                                                                    "_seed",
                                                                    .x,
                                                                    ".csv"))) 
var_names <- unique(timevarying_quantiles[[1]]$name)

rt_quantiles <- timevarying_quantiles %>%
                map(~.x %>% filter(name == "rt_t_values") %>%
                      mutate(time = time + 1) %>%
                      right_join(full_simdata, by = c("time" = "new_time")) %>%
                dplyr::select(time, week, true_rt, value, .lower, .upper, .width,.point, .interval) %>%
                filter(time <= max_time,
                       week >= 0)
                ) %>%
  bind_rows(.id = "seed")

I_quantiles <- timevarying_quantiles %>%
  map(~.x %>% filter(name == "I") %>%
        right_join(full_simdata, by = c("time" = "new_time")) %>%
        dplyr::select(time, week, num_I, value, .lower, .upper, .width,.point, .interval) %>%
        filter(time <= max_time,
               week >= 0)
  ) %>%
  bind_rows(.id = "seed")

write_csv(rt_quantiles, here::here("results", "seirr_student", paste0("seirr_scenario", sim,  "_allseeds_rt_quantiles.csv")))
write_csv(I_quantiles, here::here("results", "seirr_student", paste0("seirr_scenario", sim,  "_allseeds_prevI_quantiles.csv")))

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
    ggtitle(str_c("SEIRR Posterior Rt Scenario ", sim, " Seed ", seed_val)) +
    my_theme
}

ggsave2(filename = here::here("results", "seirr_student", paste0("seirr_rt_plots_scenario", sim, ".pdf")),
        plot = rt_quantiles %>%
          distinct(seed) %>%
          arrange(seed) %>%
          pull(seed) %>%
          map(make_rt_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)



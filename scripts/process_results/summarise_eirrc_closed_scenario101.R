# summarise eirrc results for stochastic Rt simulation 
# scenario 101
# summarise results eirrc_closed 
### Summarise processed results eirrc_closed
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
  sim = 57
} else {
  sim <- as.integer(args[1])
}

scenario_sim = sim


print(sim)
print(scenario_sim)
# read in data and results ---------------------------------------------------------
num_sims = 100
seed = 1:num_sims


timevarying_suffix <- paste0("posterior_timevarying_quantiles_scenario",
                             sim)
timevarying_list <- list.files(path = path("results", "eirrc_closed", "generated_quantities"),  pattern = timevarying_suffix)

mcmc_suffix <- paste0("mcmc_summary_scenario", 
                      sim)
mcmc_list <- list.files(path = path("results", "eirrc_closed", "mcmc_summaries"),  pattern = mcmc_suffix)

full_stan_diag <- map(mcmc_list, ~read_csv(here::here("results", "eirrc_closed",  "mcmc_summaries", .x)) %>% 
                        mutate(address = .x)) %>%
  bind_rows() %>%
  mutate(seed = as.numeric(stringr::str_extract(address, stringr::regex("(\\d+)(?!.*\\d)")))) %>%  
  group_by(seed) %>% 
  filter(variable != "R2[0]" & variable != "C[0]") %>%
  summarise(min_rhat = min(rhat),
            max_rhat = max(rhat),
            min_ess_bulk = min(ess_bulk),
            max_ess_bulk = max(ess_bulk),
            min_ess_tail = min(ess_tail),
            max_ess_tail = max(ess_tail))

write_csv(full_stan_diag, here::here("results", "eirrc_closed", paste0("eirrc_closed_scenario", sim,  "_allseeds_stan_diag.csv")))
# create final rt frame ---------------------------------------------------
full_simdata <- vector(mode='list', length=100)
for (i in 1:100) {
  full_simdata[[i]] <- read_csv(here::here("data", "sim_data", paste0("scenario101_seed", i, "_full_genecount_obsdata.csv")))
}


full_simdata <- full_simdata %>% bind_rows()%>% rename("true_rt" = "Rt") %>% 
  group_by(seed) %>% 
  fill(true_rt, .direction = "down")  # if Rt is missing, it is because no conditions changed in the day that it was missing, in such a case, the previous value is the true Rt

# the data set is simulated for a certain period of time
# but then based on how we choose to space apart observations (every two days, every seven etc)
# there is a max observed time in the data set, we should not judge the model beyond the fitted data (for now)
# this time should be the same across models (it will be the same for the case models even though they're slightly different)
fitted_simdata <- vector(mode='list', length=100)
for (i in 1:100) {
  fitted_simdata[[i]] <- read_csv(here::here("data", "sim_data", paste0("scenario101_seed", i, "_fitted_genecount_obsdata.csv")))
}



fitted_simdata <- fitted_simdata %>% bind_rows()
max_time <- fitted_simdata %>% group_by(seed) %>% summarise(max_time = max(time))

timevarying_quantiles <- map(timevarying_list, ~read_csv(here::here("results", "eirrc_closed", "generated_quantities", .x)) %>% 
                               mutate(address = .x,
                                      seed = as.numeric(stringr::str_extract(address, stringr::regex("(\\d+)(?!.*\\d)"))),
                                      scenario = stringr::str_match(address, "scenario(\\w+)_seed")[,2]))

rt_quantiles <- timevarying_quantiles %>%
  map(~.x %>% filter(name == "rt_t_values") %>%
        rename(week = time) %>% 
        dplyr::select(seed, week, scenario, value, .lower, .upper, .width,.point, .interval)) %>%
  bind_rows() %>% 
  filter(scenario == sim) %>%
  right_join(full_simdata, by = c("week", "seed")) %>%
  left_join(max_time, by = "seed") %>%
  filter(time <= max_time,
         week >= 0)

# I_quantiles <- timevarying_quantiles %>% 
#   map(~.x %>% filter(name == "I") %>%
#         right_join(full_simdata, by = c("time" = "new_time")) %>%
#         dplyr::select(week, time, num_I, value, .lower, .upper, .width,.point, .interval) %>%
#         filter(time <= max_time,
#                week >= 0)
#   ) %>%
#   bind_rows(.id = "seed")

write_csv(rt_quantiles, here::here("results", "eirrc_closed", paste0("eirrc_scenario", sim,  "_allseeds_rt_quantiles.csv")))

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
    ggtitle(str_c("EIRR-WW Posterior Rt Scenario ", sim, " Seed ", seed_val)) +
    my_theme
}
# 
# testing =   rt_quantiles %>%
#   filter(seed == seed_val) %>% 
#   dplyr::select(seed, week, time,value, .lower, .upper, true_rt)
ggsave2(filename = here::here("results", "eirrc_closed", paste0("eirrc_closed_rt_plots_scenario", sim, ".pdf")),
        plot = rt_quantiles %>%
          distinct(seed) %>%
          arrange(seed) %>%
          pull(seed) %>%
          map(make_rt_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)



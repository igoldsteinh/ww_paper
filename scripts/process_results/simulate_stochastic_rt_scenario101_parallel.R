# finish simulating scenario 101
# in parallel

library(tidyverse)
library(glmnet)
library(patchwork)
library(GGally)
library(scales)
library(cowplot)
library(gridExtra)


args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  seed = 1
} else {
  seed <- as.integer(args[1])
}

source("src/simulate_stochastic_seirr.R")
# make a changing vector for r0
r01 = rep(1.75, 7 * 5)  
r02 = rep(0.7, 7 * 4)
r022 = exp(seq(log(0.7), log(0.9), length.out = (7 * 2)))
r03 = exp(seq(log(0.9), log(2.5), length.out = (7 * 5))) # then a five week jump up to 2.5
r0_init = r01[1]
r0_vec = c(r01, r02, r022, r03)
r01_val = unique(r01)
r02_change = unique(r02)
r022_change = unique(r022)
unique_r0 = c(r02_change, r022_change, r03)
change1 = min(which(r0_vec == r02_change))
# we start at 0.95 at the start of r03, so that is the first time point
next_changes = (max(which(r0_vec == r02_change))):(length(r0_vec))
change_points = c(change1, next_changes)

N = 100000
I_init = 200
E_init = 200
gamma = 1/4
nu = 1/7
eta = 1/18
beta_init = r0_init * nu
beta_vec = unique_r0 * nu 

for (seed in 1:100) {
set.seed(seed)
# we have checked and the first 100 all have enough data, so we will just use these
individ_data <- read_csv(paste0("data/sim_data/scenario101_seed",
                                           seed,
                                           "individ_data.csv")) %>% 
                      mutate(seed = seed,
                             max_time = max(t),
                             ok = max_time >= (7*30))


state_data <- read_csv(paste0("data/sim_data/scenario101_seed",
                                          seed,
                                         "truecurve.csv")) %>% 
                    mutate(seed = seed)


# plot the data -----------------------------------------------------------
# creating true gene count data ------------------------------------------------------
# write_csv(individ_data, here::here("data", "sim_data", "scenario101_individ_data.csv"))
# individ_data = read_csv(here::here("data", "sim_data", "scenario101_individ_data.csv"))

wide_format = individ_data %>% group_by(labels) %>% 
                    arrange(labels, type) %>%
                    dplyr::select(t, type, labels) %>%
                    filter(type != 5) %>%
                    pivot_wider(id_cols = labels, names_from = type, values_from = t ) %>%
                    rename("infection_time" = `1`, "infectious_time" = `2`, "recover_time" = `3`, "stopshed_time" = `4`) %>%
                    dplyr::select(labels, infection_time, infectious_time, recover_time, stopshed_time) %>%
                    mutate(infectious_period = recover_time - infectious_time,
                           r1_period = stopshed_time - recover_time) 

obs_start = max(which(r0_vec == 0.9)) 
obs_end = obs_start + (19 * 7) - 1


times = obs_start:obs_end


# create true RNA conc data ----------------------------------------------

  true_conc_data = map(times, ~calc_total_gene_counts(wide_format, time = .x, N = N)) %>%
    bind_rows()
  



# simulate data sets from the true values (Scenario 1)------------------------
obs_true_data = true_conc_data %>% filter(time >= obs_start & time <= obs_end)

obs_data = simulate_gene_data(obs_true_data,
                              seed = seed,
                              rho = exp(-16),
                              t_sd = 0.5,
                              t_df = 2.99)
# take every other day 
final_obs_data = obs_data %>%
                       mutate(new_time = time - obs_start + 1,
                              everyother = new_time %% 2 == 1)  %>% 
                       mutate(week = floor((new_time - 1)/7)) %>%
                       filter(everyother == TRUE)


full_obs_data <- obs_data %>%
                        left_join(state_data, by = c("time" = "integer_day")) %>%
                        mutate(new_time = time - obs_start + 1) %>% 
                        mutate(week = floor((new_time - 1)/7))

obs_true_data <- obs_true_data %>% 
                        left_join(state_data, by = c("time" = "integer_day"))


full_obs_data <- full_obs_data %>% 
  mutate(seed = seed)

final_obs_data <- final_obs_data %>% 
  mutate(seed = seed)

write_csv(obs_true_data, here::here("data", "sim_data", paste0("scenario101_seed", seed, "truegenecounts.csv")))
write_csv(full_obs_data, here::here("data", "sim_data", paste0("scenario101_seed", seed, "_full_genecount_obsdata.csv")))
write_csv(final_obs_data, here::here("data", "sim_data", paste0("scenario101_seed", seed, "_fitted_genecount_obsdata.csv")))


# initial conditions ------------------------------------------------------
initial_states = state_data %>% 
                       mutate(diff = (obs_start-1) - integer_day) %>%
                       filter(diff > 0) %>%
                       filter(diff == min(diff)) %>% 
                       mutate(seed = seed)

write_csv(initial_states, here::here("data", "sim_data", paste0("scenario101_seed", seed, "_initstates.csv")))
}

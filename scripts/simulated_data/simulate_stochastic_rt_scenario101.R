# simulate observed gene concentrations from individual seirr
# where each simulation comes from a different Rt trajectory
# simulate observed gene counts from our individual seirr engine
# labeled scenario 101
library(tidyverse)
library(glmnet)
library(patchwork)
library(GGally)
library(scales)
library(cowplot)
library(gridExtra)

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

set.seed(1)
# we have checked and the first 100 all have enough data, so we will just use these
seed = 1:100
individ_data <- map(seed, ~read_csv(paste0("data/sim_data/scenario101_seed",
                                                    .x,
                                                    "individ_data.csv")) %>% 
                     mutate(seed = .x,
                            max_time = max(t),
                            ok = max_time >= (7*30)))


state_data <- map(seed, ~read_csv(paste0("data/sim_data/scenario101_seed",
                                           .x,
                                           "truecurve.csv")) %>% 
                      mutate(seed = .x))

write_rds(individ_data, here::here("data", "sim_data", "scenario101_individ_data.rds"))
write_rds(state_data, here::here("data", "scenario101_truecurve.rds"))

# plot the data -----------------------------------------------------------
# creating true gene count data ------------------------------------------------------
# write_csv(individ_data, here::here("data", "sim_data", "scenario101_individ_data.csv"))
# individ_data = read_csv(here::here("data", "sim_data", "scenario101_individ_data.csv"))

wide_format = map(individ_data, ~.x %>% group_by(labels) %>% 
                    arrange(labels, type) %>%
                    dplyr::select(t, type, labels) %>%
                    filter(type != 5) %>%
                    pivot_wider(id_cols = labels, names_from = type, values_from = t ) %>%
                    rename("infection_time" = `1`, "infectious_time" = `2`, "recover_time" = `3`, "stopshed_time" = `4`) %>%
                    dplyr::select(labels, infection_time, infectious_time, recover_time, stopshed_time) %>%
                    mutate(infectious_period = recover_time - infectious_time,
                           r1_period = stopshed_time - recover_time)) 

times = map(individ_data, ~1:floor(max(.x$t)))


# create true RNA conc data ----------------------------------------------

true_conc_data <- vector(mode='list', length=100)
for (i in (1:100)) {
  current_times = times[[i]]
  
  true_conc_data[[i]] = map(current_times, ~calc_total_gene_counts(wide_format[[i]], time = .x, N = N)) %>%
    bind_rows()
  
}



# simulate data sets from the true values (Scenario 1)------------------------
obs_start = max(which(r0_vec == 0.9)) 
obs_end = obs_start + (19 * 7) - 1
seeds = 1:100

obs_true_data = map(true_conc_data, ~.x %>% filter(time >= obs_start & time <= obs_end))

obs_data = map2(seeds, obs_true_data, ~simulate_gene_data(.y, 
                                                          seed = .x, 
                                                          rho = exp(-16), 
                                                          t_sd = 0.5,
                                                          t_df = 2.99))
# take every other day 
final_obs_data = map(obs_data, ~.x %>%
                       mutate(new_time = time - obs_start + 1,
                              everyother = new_time %% 2 == 1)  %>% 
                       mutate(week = floor((new_time - 1)/7)) %>%
                       filter(everyother == TRUE))


full_obs_data <- map2(obs_data, state_data, ~.x %>%
                        left_join(.y, by = c("time" = "integer_day")) %>%
                        mutate(new_time = time - obs_start + 1) %>% 
                        mutate(week = floor((new_time - 1)/7)))

obs_true_data <- map2(obs_true_data, state_data, ~.x %>% 
                        left_join(.y, by = c("time" = "integer_day")))


full_obs_data <- full_obs_data %>% 
  bind_rows(.id = "seed")

final_obs_data <- final_obs_data %>% 
  bind_rows(.id = "seed")

write_rds(obs_true_data, here::here("data", "sim_data", "scenario101_truegenecounts.rds"))
write_csv(full_obs_data, here::here("data", "sim_data", "scenario101_full_genecount_obsdata.csv"))
write_csv(final_obs_data, here::here("data", "sim_data", "scenario101_fitted_genecount_obsdata.csv"))


# initial conditions ------------------------------------------------------
initial_states = map(state_data, ~.x %>% mutate(diff = (obs_start-1) - integer_day) %>%
                                         filter(diff > 0) %>%
                                         filter(diff == min(diff))) %>% bind_rows(.id = "seed")

testing = state_data %>% bind_rows(.id = "seed") %>% filter(seed == 3)
write_csv(initial_states, here::here("data", "sim_data", "scenario101_initstates.csv"))


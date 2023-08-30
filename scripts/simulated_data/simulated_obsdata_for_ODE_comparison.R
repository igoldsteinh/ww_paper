# resimulating data at half-daily resolution in order to compare ODE solver 
# vs closed form solution model 
library(tidyverse)
library(glmnet)
library(patchwork)
library(GGally)
library(scales)

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

individ_data = read_csv(here::here("data", "sim_data", "scenario1_individ_data.csv"))
state_data = read_csv(here::here("data", "sim_data", "scenario1_truecurve.csv"))


wide_format = individ_data %>% group_by(labels) %>% 
  arrange(labels, type) %>%
  dplyr::select(t, type, labels) %>%
  filter(type != 5) %>%
  pivot_wider(id_cols = labels, names_from = type, values_from = t ) %>%
  rename("infection_time" = `1`, "infectious_time" = `2`, "recover_time" = `3`, "stopshed_time" = `4`) %>%
  dplyr::select(labels, infection_time, infectious_time, recover_time, stopshed_time) %>%
  mutate(infectious_period = recover_time - infectious_time,
         r1_period = stopshed_time - recover_time) 

times = seq(from = 2, to = floor(max(individ_data$t)), by = 0.5)

# create true RNA count data ----------------------------------------------
true_count_data = map(times, ~calc_total_gene_counts(wide_format, time = .x, N = N)) %>%
  bind_rows()

# simulate multiple data sets from the true values (Scenario 1)------------------------
obs_start = max(which(r0_vec == 0.9)) 
obs_end = obs_start + (19 * 7) - 1
seeds = 1:100

obs_true_data = true_count_data %>% filter(time >= obs_start & time <= obs_end)

obs_data = map(seeds, ~simulate_gene_data(obs_true_data, 
                                          seed = .x, 
                                          rho = exp(-16), # we should choose rho so that we're in the ballpark of the real data I think 
                                          t_sd = 0.5,
                                          t_df = 2.99))
final_obs_data = map(obs_data, ~.x %>%
                       mutate(new_time = time - obs_start + 1,
                              everyother = new_time %% 2 == 1)  %>% 
                       mutate(week = floor((new_time - 1)/7)))



full_obs_data <- map(obs_data, ~.x %>%
                       left_join(state_data, by = c("time" = "integer_day")) %>%
                       mutate(new_time = time - obs_start + 1) %>% 
                       mutate(week = floor((new_time - 1)/7)))

obs_true_data <- obs_true_data %>% 
  left_join(state_data, by = c("time" = "integer_day"))


full_obs_data <- full_obs_data %>% 
  bind_rows(.id = "seed")

final_obs_data <- final_obs_data %>% 
  bind_rows(.id = "seed")

write_csv(obs_true_data, here::here("data", "sim_data", "ODE_comp_truegenecounts.csv"))
write_csv(full_obs_data, here::here("data", "sim_data", "ODE_comp_full_genecount_obsdata.csv"))
write_csv(final_obs_data, here::here("data", "sim_data", "ODE_comp_fitted_genecount_obsdata.csv"))





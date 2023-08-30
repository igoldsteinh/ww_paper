# simulate observed gene counts from our individual seirr engine
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

N = 100000
I_init = 200
E_init = 200
gamma = 1/4
nu = 1/7
eta = 1/18
beta_init = r0_init * nu
beta_vec = unique_r0 * nu 
set.seed(3)
testing = sim_SEIRR_nonconst(N = N, 
                             I_init = I_init, 
                             E_init,
                             beta_init = beta_init,
                             beta_vec = beta_vec,
                             change_points = change_points,
                             gamma = gamma,
                             nu = nu, 
                             eta = eta)

individ_data = testing[[1]]
state_data = create_daily_data(testing[[2]]) %>% 
  ungroup() %>% 
  mutate(Rt = R0 * lag(S)/N) 

write_csv(individ_data, here::here("data", "sim_data", "scenario1_individ_data.csv"))
write_csv(state_data, here::here("data", "sim_data", "scenario1_truecurve.csv"))

# plot the data -----------------------------------------------------------

plot_rt = state_data %>% 
          ggplot(aes(x = time, y = Rt)) +
          geom_point() + 
          theme_bw()


plot_states <- state_data %>% 
               dplyr::select(time,  I) %>% 
               pivot_longer(-time) %>%
               ggplot(aes(x = time, y = value, color = name)) + 
               geom_point()
# creating true gene count data ------------------------------------------------------
# write_csv(individ_data, here::here("data", "sim_data", "scenario1_individ_data.csv"))
# individ_data = read_csv(here::here("data", "sim_data", "scenario1_individ_data.csv"))

wide_format = individ_data %>% group_by(labels) %>% 
  arrange(labels, type) %>%
  dplyr::select(t, type, labels) %>%
  filter(type != 5) %>%
  pivot_wider(id_cols = labels, names_from = type, values_from = t ) %>%
  rename("infection_time" = `1`, "infectious_time" = `2`, "recover_time" = `3`, "stopshed_time" = `4`) %>%
  dplyr::select(labels, infection_time, infectious_time, recover_time, stopshed_time) %>%
  mutate(infectious_period = recover_time - infectious_time,
         r1_period = stopshed_time - recover_time) 

times = 2:floor(max(individ_data$t))

# creating true case data -------------------------------------------------
true_E2I_data = map(times, ~calc_E2I_counts(wide_format, time = .x)) %>% 
  bind_rows()

plot_E2I <- true_E2I_data %>% 
  filter(time >= 92 & time <= 120) %>%
  ggplot(aes(x = time, y = E2I_transitions)) + 
  geom_point()

plot_E2I

# create simulated case data ----------------------------------------------
obs_start = max(which(r0_vec == 0.9)) 
obs_end = obs_start + (19 * 7) - 1
seeds = 1:100

true_E2I_data <- true_E2I_data %>% filter(time >= obs_start & time <= obs_end)



case_data = map(seeds, ~simulate_case_data(true_E2I_data, 
                                           seed = .x, 
                                           rho = 0.2, 
                                           phi = 57.55))
# aggregate to the week
final_case_data = map(case_data, ~ .x %>%
                        mutate(new_time = time - obs_start + 1) %>% 
                        mutate(week = floor((new_time - 1)/7)) %>%
                        group_by(week) %>%
                        summarise(total_cases = sum(cases),
                                  total_E2I = sum(E2I_transitions),
                                  num_days = n()) %>%
                        mutate(new_week = week + 1))


case_data <- case_data %>% 
  bind_rows(.id = "seed")

final_case_data <- final_case_data %>% 
  bind_rows(.id = "seed")

write_csv(true_E2I_data, here::here("data", "sim_data", "scenario1_truecasecounts.csv"))
write_csv(case_data, here::here("data", "sim_data", "scenario1_full_cases_obsdata.csv"))
write_csv(final_case_data, here::here("data", "sim_data", "scenario1_fitted_cases_obsdata.csv"))

# create true RNA count data ----------------------------------------------
true_count_data = map(times, ~calc_total_gene_counts(wide_format, time = .x, N = N)) %>%
  bind_rows()


# simulate multiple data sets from the true values (Scenario 1)------------------------
obs_true_data = true_count_data %>% filter(time >= obs_start & time <= obs_end)

obs_data = map(seeds, ~simulate_gene_data(obs_true_data, 
                              seed = .x, 
                              rho = exp(-16), # we should choose rho so that we're in the ballpark of the real data I think 
                              t_sd = 0.5,
                              t_df = 2.99))
# take every other day 
final_obs_data = map(obs_data, ~.x %>%
  mutate(new_time = time - obs_start + 1,
         everyother = new_time %% 2 == 1)  %>% 
    mutate(week = floor((new_time - 1)/7)) %>%
  filter(everyother == TRUE))


# make a copy of state data 

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

write_csv(obs_true_data, here::here("data", "sim_data", "scenario1_truegenecounts.csv"))
write_csv(full_obs_data, here::here("data", "sim_data", "scenario1_full_genecount_obsdata.csv"))
write_csv(final_obs_data, here::here("data", "sim_data", "scenario1_fitted_genecount_obsdata.csv"))

# create overdisp data ----------------------------------------------------
true_E2I_data <- true_E2I_data %>% filter(time >= obs_start & time <= obs_end)
case_data = simulate_case_data(true_E2I_data, 
                                           seed = 111, 
                                           rho = 0.2, 
                                           phi = 57.55)
# aggregate to the week
final_case_data = case_data %>%
                        mutate(new_time = time - obs_start + 1) %>% 
                        mutate(week = floor((new_time - 1)/7)) %>%
                        group_by(week) %>%
                        summarise(total_cases = sum(cases),
                                  total_E2I = sum(E2I_transitions),
                                  num_days = n()) %>%
                        mutate(new_week = week + 1)



write_csv(final_case_data, here::here("data", "sim_data", "scenario1_overdisp.csv"))

# initial conditions ------------------------------------------------------
# state_data = read_csv(here::here("data", "sim_data", "scenario1_truecurve.csv"))

initial_states = state_data %>% filter(integer_day == obs_start - 1)

# WNAR Data --------------------------------
obs_start = min(next_changes)
obs_end = 225

obs_true_data = true_count_data %>% filter(time >= obs_start & time <= obs_end)

obs_data = simulate_gene_data(obs_true_data, 
                              seed = 1, 
                              rho = exp(-16),
                              t_sd = 0.5,
                              t_df = 2.99)
# take every other day 
final_obs_data = obs_data %>%
                 mutate(new_time = time - obs_start + 1,
                        everyother = new_time %% 2 == 1) %>%
                 filter(everyother == TRUE)

state_data <- state_data %>% dplyr::select(-time)
full_obs_data <- obs_data %>% 
                 left_join(state_data, by = c("time" = "integer_day"))
                
obs_true_data <- obs_true_data %>% 
                 left_join(state_data, by = c("time" = "integer_day"))
write_csv(obs_true_data, here::here("data", "sim_data", "scenario41_truegenecounts.csv"))
write_csv(full_obs_data, here::here("data", "sim_data", "scenario41_full_obsdata.csv"))
write_csv(final_obs_data, here::here("data", "sim_data", "scenario41_fitted_obsdata.csv"))



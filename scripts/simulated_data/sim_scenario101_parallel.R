# simulating takes a while with 100K population
# lets get the simulation running in parralel
# check if its long enough after the fact
# simulate observed gene concentrations from individual seirr
# where each simulation comes from a different Rt trajectory
# simulate observed gene counts from our individual seirr engine
# labeled scenario 101
library(tidyverse)

source("src/simulate_stochastic_seirr.R")
# make a changing vector for r0

args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  seed = 1
} else {
  seed <- as.integer(args[1])
}

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

# simulate the epidemic
set.seed(seed)
potential_epidemic = sim_SEIRR_nonconst(N = N, 
                                          I_init = I_init, 
                                          E_init,
                                          beta_init = beta_init,
                                          beta_vec = beta_vec,
                                          change_points = change_points,
                                          gamma = gamma,
                                          nu = nu, 
                                          eta = eta)

# save the individual times, and also the aggregated state counts
individ_data = potential_epidemic[[1]]
state_data = create_daily_data(potential_epidemic[[2]]) %>% 
  ungroup() %>% 
  mutate(Rt = R0 * lag(S)/N) 

write_csv(individ_data, here::here("data", "sim_data", paste0("scenario101_seed", seed, "individ_data.csv")))
write_csv(state_data, here::here("data", "sim_data", paste0("scenario101_seed", seed, "truecurve.csv")))

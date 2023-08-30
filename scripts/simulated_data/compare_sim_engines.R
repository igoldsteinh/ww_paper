# sanity checking our seirr sim vs agent based
library(tidyverse)
library(patchwork)
library(tidybayes)
library(GGally)
library(gridExtra)
library(cowplot)
library(scales)
library(latex2exp)

source("src/simulate_stochastic_seirr.R")

# our sim 
pop_size = 100
r0 = 1.5
nu = 1/7
beta =  (r0 * nu)
beta_indiv = (r0 * nu)/pop_size
gamma = 1/4
eta = 1/18
I_init = 5
N = pop_size

set.seed(1)
our_sim <- sim_SEIRR_nonconst(N = N, 
                              I_init = I_init, 
                              E_init = 0,
                              beta_init = beta,
                              beta_vec = -3,
                              change_points = Inf,
                              gamma = gamma,
                              nu = nu, 
                              eta = eta)
our_daily <- create_daily_data(our_sim[[2]])


agent_sim <- sim_agent_SEIRR(pop_size, I_init = I_init, beta = beta_indiv, gamma = gamma, nu = nu, eta = eta)

agent_daily <- create_daily_from_agent(agent_sim)

# visualize ---------------------------------------------------------------
truecurve <- our_daily %>% 
  dplyr::select(integer_day, S, E, I, R1) %>% 
  pivot_longer(-integer_day) %>% 
  filter(integer_day < 210) %>%
  rename(Compartment = name)

level_list <- c("S", "E", "I", "R1")

truecurve$Compartment <- factor(truecurve$Compartment, levels=level_list)

epi_curve <- truecurve %>%
  ggplot(aes(x = integer_day, y = value, color = Compartment)) + 
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("My Engine") +
  geom_point() + 
  theme_bw() + 
  theme(legend.position = c(0.63, 0.78),
        legend.background = element_rect("transparent"),
        text = element_text(size = 18)) 

#redo for agent based
agentcurve <- agent_daily %>% 
  dplyr::select(integer_day, S, E, I, R1) %>% 
  pivot_longer(-integer_day) %>% 
  filter(integer_day < 210) %>%
  rename(Compartment = name)

level_list <- c("S", "E", "I", "R1")

agentcurve$Compartment <- factor(agentcurve$Compartment, levels=level_list)

agent_epi_curve <- agentcurve %>%
  ggplot(aes(x = integer_day, y = value, color = Compartment)) + 
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("Agent Engine") +
  geom_point() + 
  theme_bw() + 
  theme(legend.position = c(0.63, 0.78),
        legend.background = element_rect("transparent"),
        text = element_text(size = 18)) 

compare <- agent_epi_curve + epi_curve

compare

# lets do it for 100 ------------------------------------------------------
num_sims = 10000
mysims = vector(mode='list', length=num_sims)
mydaily = vector(mode='list', length=num_sims)
agentsims = vector(mode='list', length=num_sims)
agentdaily = vector(mode='list', length=num_sims)

# params
pop_size = 100
r0 = 1.5
nu = 1/7
beta =  (r0 * nu)
beta_indiv = (r0 * nu)/pop_size
gamma = 1/4
eta = 1/18
I_init = 5
N = pop_size

for (i in 1:num_sims) {
  set.seed(i)
  mysims[[i]] = sim_SEIRR_nonconst(N = N, 
                                   I_init = I_init, 
                                   E_init = 0,
                                   beta_init = beta,
                                   beta_vec = -3,
                                   change_points = Inf,
                                   gamma = gamma,
                                   nu = nu, 
                                   eta = eta)
  mydaily[[i]] = create_daily_data(mysims[[i]][[2]])
}

for (i in 1:num_sims) {
  agentsims[[i]] = sim_agent_SEIRR(pop_size, I_init = I_init, beta = beta_indiv, gamma = gamma, nu = nu, eta = eta)
  
  
  
  agentdaily[[i]] = create_daily_from_agent(agentsims[[i]])
  
}

write_rds(agentsims, here::here("data", "sim_data", "agentsims.rds"))
allmy = bind_rows(mydaily, .id = "id") 
allagent = bind_rows(agentdaily, .id = "id")

avgmy = allmy %>% 
  dplyr::select(integer_day,  S, E, I, R1) %>% 
  pivot_longer(-integer_day) %>% 
  group_by(name, integer_day) %>%
  mean_qi(.width = c(0.5, 0.8, 0.95))   %>%
  rename(Compartment = name) 

  

avgagent <- allagent %>% 
  dplyr::select(integer_day,  S, E, I, R1) %>% 
  pivot_longer(-integer_day) %>% 
  group_by(name, integer_day) %>%
  mean_qi(.width = c(0.5, 0.8, 0.95))   %>%
  rename(Compartment = name) 

  

level_list <- c("S", "E", "I", "R1")

avgmy$Compartment <- factor(avgmy$Compartment, levels=level_list)

my_theme <- list(
  scale_fill_brewer(name = "Quantiles",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme())


avgmy_curve <- avgmy %>%
  filter(integer_day <= 100) %>% 
  ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() +
  facet_wrap(vars(Compartment), scales = "free") +
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("My Engine") +
  my_theme + 
  theme(text = element_text(size = 18)) 

#redo for agent based

avgagent$Compartment <- factor(avgagent$Compartment, levels=level_list)

avgagent_curve <- avgagent %>%
  filter(integer_day <= 100) %>% 
  ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() +
  facet_wrap(vars(Compartment), scales = "free") +
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("Agent Engine") +
  my_theme + 
  theme(text = element_text(size = 18)) 

avg_compare <- avgmy_curve + avgagent_curve + plot_layout(guides = "collect")

avg_compare

# calculate regular gillespie states --------------------------------------
num_sims = 10000
gillespiesims = vector(mode='list', length=num_sims)
gillespiedaily = vector(mode='list', length=num_sims)

# params
pop_size = 100
r0 = 1.5
nu = 1/7
beta =  (r0 * nu)
beta_indiv = (r0 * nu)/pop_size
gamma = 1/4
eta = 1/18
I_init = 5
N = pop_size

for (i in 1:num_sims) {
  set.seed(i)
  gillespiesims[[i]] = gillespie_seirr(pop_size = N, 
                                   I_init = I_init, 
                                   beta = beta,
                                   gamma = gamma,
                                   nu = nu, 
                                   eta = eta)
  
  gillespiedaily[[i]] = create_gillespie_daily_data(gillespiesims[[i]])
}


allgillespie = bind_rows(gillespiedaily, .id = "id") 

avggillespie = allgillespie %>% 
  dplyr::select(integer_day,  S, E, I, R1) %>% 
  pivot_longer(-integer_day) %>% 
  group_by(name, integer_day) %>%
  mean_qi(.width = c(0.5, 0.8, 0.95))   %>%
  rename(Compartment = name) 

level_list <- c("S", "E", "I", "R1")

avggillespie$Compartment <- factor(avggillespie$Compartment, levels=level_list)

avggillespie_curve <- avggillespie %>%
  filter(integer_day <= 100) %>%
  ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() +
  facet_wrap(vars(Compartment), scales = "free") +
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("Gillespie Engine") +
  my_theme + 
  theme(text = element_text(size = 18)) 

full_compare <- avg_compare + avggillespie_curve + plot_layout(guides = "collect")

full_compare

# compare by compartment --------------------------------------------------
allmy <- allmy %>% dplyr::select(id, integer_day, S, E, I, R1) %>% mutate(engine = "Mine")
allgillespie <- allgillespie %>% dplyr::select(id, integer_day, S, E, I, R1) %>% mutate(engine = "Gillespie")
allagent <- allagent %>% dplyr::select(id, integer_day, S, E, I, R1) %>% mutate(engine = "Agent")

write_rds(allmy, here::here("data", "sim_data", "allmy_sims.rds"))
write_rds(allgillespie, here::here("data", "sim_data", "allgillespie_sims.rds"))
write_rds(allagent, here::here("data", "sim_data", "allagent.rds"))
combined_sims <- bind_rows(allmy, allgillespie, allagent) %>%
                 dplyr::select(-id) %>% 
                 pivot_longer(-c(engine, integer_day)) %>%
                 group_by(engine, name, integer_day) %>% 
                 mean_qi(.width = c(0.5, 0.8, .95))

mean_sd_calcs <- bind_rows(allmy, allgillespie, allagent) %>%
  dplyr::select(-id) %>% 
  pivot_longer(-c(engine, integer_day)) %>%
  group_by(engine, name, integer_day) %>% 
  summarise(mean = mean(value),
            sd = sd(value),
            num = n(),
            se = sd/sqrt(num),
            lb = mean - 2 * se, 
            ub = mean + 2 * se) %>% 
  pivot_wider(id_cols = c(name, integer_day), names_from = engine, values_from = c(mean, sd, num, se, lb, ub)) %>% 
  group_by(name, integer_day) %>% 
  mutate(diff_agent_my = abs(mean_Agent - mean_Mine),
         diff_agent_gill = abs(mean_Agent - mean_Gillespie),
         diff_my_gill = abs(mean_Mine - mean_Gillespie),
         cover_agent_my = (ub_Mine >= lb_Agent & ub_Mine <= ub_Agent) | (ub_Agent >= lb_Mine & ub_Agent <= ub_Mine),
         cover_agent_gill = (ub_Agent >= lb_Gillespie & ub_Agent <= ub_Gillespie) | (ub_Gillespie >= lb_Agent & ub_Gillespie <= ub_Agent), 
         cover_my_gill = (ub_Mine >= lb_Gillespie & ub_Mine <= ub_Gillespie) | (ub_Gillespie >= lb_Mine & ub_Gillespie <= ub_Mine),
         all_cover = (cover_agent_my == TRUE) & (cover_agent_gill == TRUE) & (cover_my_gill == TRUE),
         min_num = min(num_Mine, num_Gillespie, num_Agent))


summarise_means <- mean_sd_calcs %>% 
                   group_by(name) %>%
                   filter(integer_day <= 100) %>% 
                   summarise(max_diff = max(diff_agent_my, diff_agent_gill, diff_my_gill),
                             prop_cover = sum(all_cover)/n(),
                             max_mc_error = max(se_Mine, se_Agent), 
                             min_mc_error = min(se_Mine, se_Agent))

looksee = mean_sd_calcs %>% filter(max(diff_agent_my, diff_agent_gill, diff_my_gill) >= 1.21)
S_alone <- mean_sd_calcs %>% filter(name == "S")

example_day <- S_alone %>% filter(integer_day == 50) 
combined_sims$engine[combined_sims$engine == "Mine"] <- "Variation"
combined_sims$engine[combined_sims$engine == "Gillespie"] <- "Compartment"

S_plot <- combined_sims %>% 
          filter(name == "S") %>% 
          filter(integer_day <= 100) %>%
          ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
          geom_lineribbon() +
          facet_wrap(vars(engine)) +
          xlab("Time") + 
          ylab("Compartment Counts") + 
          ggtitle("S Compartment") +
          my_theme + 
          theme(text = element_text(size = 18)) 

E_plot <- combined_sims %>% 
  filter(name == "E") %>% 
  filter(integer_day <= 100) %>%
  ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() +
  facet_wrap(vars(engine)) +
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("E Compartment") +
  my_theme + 
  theme(text = element_text(size = 18)) 

I_plot <- combined_sims %>% 
  filter(name == "I") %>% 
  filter(integer_day <= 100) %>%
  ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() +
  facet_wrap(vars(engine)) +
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("I Compartment") +
  my_theme + 
  theme(text = element_text(size = 18)) 


R1_plot <- combined_sims %>% 
  filter(name == "R1") %>% 
  filter(integer_day <= 100) %>%
  ggplot(aes(x = integer_day, y = value, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() +
  facet_wrap(vars(engine)) +
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("R1 Compartment") +
  my_theme + 
  theme(text = element_text(size = 18)) 

ggsave(here::here("figures", "S_mc_compplot.pdf"), S_plot, width = 10, height = 6)
ggsave(here::here("figures", "E_mc_compplot.pdf"), E_plot, width = 10, height = 6)
ggsave(here::here("figures", "I_mc_compplot.pdf"), I_plot, width = 10, height = 6)
ggsave(here::here("figures", "R1_mc_compplot.pdf"), R1_plot, width = 10, height = 6)

# calculate pathogen shedding for both engines ----------------------------
num_sims = 10000
N = 100
my_conc = vector(mode='list', length=num_sims)

for (i in 1:num_sims) {
  # print(i)
  wide_format = mysims[[i]][[1]] %>% group_by(labels) %>% 
    arrange(labels, type) %>%
    dplyr::select(t, type, labels) %>%
    filter(type != 5) 
  
  wide_format$type <- as.factor(wide_format$type)
  
  wide_format <- wide_format %>%
    pivot_wider(id_cols = labels, names_from = type, values_from = t, names_expand = TRUE) %>%
    rename("infectious_time" = `2`, "recover_time" = `3`, "stopshed_time" = `4`) %>%
    dplyr::select(labels, infectious_time, recover_time, stopshed_time) %>%
    mutate(infectious_period = recover_time - infectious_time,
           r1_period = stopshed_time - recover_time) 
  
  times = 2:floor(max(mysims[[i]][[1]]$t))
  
  true_count_data = map(times, ~calc_total_gene_counts(wide_format, time = .x, N = N)) %>%
    bind_rows() %>% 
    mutate(total_count = total_count/N)
  
  obs_data = simulate_gene_data(true_count_data, 
                                seed = i, 
                                rho = exp(-16) * 100000,
                                t_sd = 0.5,
                                t_df = 2.99)
  
  
  my_conc[[i]] = obs_data
  
}

# agent engine ------------------------------------------------------------
num_sims = 10000
agent_conc = vector(mode='list', length=num_sims)

for (i in 1:num_sims) {
  # print(i)
  dataframe = as.data.frame(agentsims[[i]])
  dataframe$V1 <- as.numeric(dataframe$V1)
  testing = dataframe %>% 
    rename("time" = "V1") %>% 
    pivot_longer(-time)
  
  testing$time <- as.numeric(testing$time)
  testing = testing %>% 
    group_by(name) %>%
    arrange(name, time) %>% 
    mutate(lag_state = lag(value),
           change = lag_state != value) %>% 
    filter(change == TRUE)
  
  testing$value <- as.factor(testing$value)
  
  if(!("I" %in% unique(testing$value))) {
    new_rows <- dataframe %>% ungroup() %>% filter(V1 == 0) %>% 
      rename("time" = "V1") %>% 
      pivot_longer(-time) %>% 
      filter(name %in% testing$name)  %>% 
      mutate(lag_state = NA,
             change = NA)
    
    testing = bind_rows(testing, new_rows)
    
  }
  
  testing = testing %>% 
    pivot_wider(id_cols = name, names_from = value, values_from = time, names_expand = TRUE) %>% 
    rename("infectious_time" = "I",
           "recover_time" = "R1",
           "stopshed_time" = "R2") 
  
  testing$infectious_time[is.na(testing$infectious_time)] <- 0
  
  testing = testing %>% 
    mutate(infectious_period = recover_time - infectious_time,
           r1_period = stopshed_time - recover_time) 
  
  times = 2:floor(max(dataframe$V1))
  N = 100
  true_count_data = map(times, ~calc_total_gene_counts(testing, time = .x, N = N)) %>%
    bind_rows() %>%
    mutate(total_count = total_count/N)
  
  obs_data = simulate_gene_data(true_count_data, 
                                seed = i, 
                                rho = exp(-16) * 100000,
                                t_sd = 0.5,
                                t_df = 2.99)
  
  
  agent_conc[[i]] = obs_data
  
}

write_rds(agent_conc, here::here("data", "sim_data", "agent_conc.rds"))
write_rds(my_conc, here::here("data", "sim_data", "my_conc.rds"))

allagent_conc <- agent_conc %>% bind_rows(.id = "id")
allmy_conc <- my_conc %>% bind_rows(.id = "id")

allagent_conc <- allagent_conc %>% mutate(engine = "Agent")
allmy_conc <- allmy_conc %>% mutate(engine = "Mine")

all_conc <- bind_rows(allagent_conc, allmy_conc)


all_conc <- all_conc %>% 
            dplyr::select(id, engine, time, log_gene_copies1) %>% 
            group_by(engine, id) %>% 
            pivot_longer(-c(time, engine, id))

all_conc$value[all_conc$value == -Inf] <- 0
mean_concs <- all_conc %>% 
  group_by(engine, time) %>% 
  summarise(mean = mean(value),
            sd = sd(value),
            num = n(),
            se = sd/sqrt(num),
            lb = mean - 2 * se, 
            ub = mean + 2 * se) %>% 
  pivot_wider(id_cols = c(time), names_from = engine, values_from = c(mean, sd, num, se, lb, ub)) %>% 
  group_by(time) %>% 
  mutate(diff_agent_my = abs(mean_Agent - mean_Mine),
         cover_agent_my = (ub_Mine >= lb_Agent & ub_Mine <= ub_Agent) | (ub_Agent >= lb_Mine & ub_Agent <= ub_Mine),
         min_num = min(num_Mine, num_Agent))

summarise_mean_concs <- mean_concs %>% 
  ungroup() %>% 
  filter(time <= 100) %>% 
  summarise(max_diff = max(diff_agent_my),
            prop_cover = sum(cover_agent_my)/n(),
            max_mc_error = max(se_Mine, se_Agent), 
            min_mc_error = min(se_Mine, se_Agent))

not_covered <- mean_concs %>% filter(cover_agent_my == FALSE & min_num >= 1000)

all_conc$engine[all_conc$engine == "Mine"] <- "Variation"

conc_plot_data <- all_conc %>% ungroup() %>% dplyr::select(engine, time, value)
conc_plot <- all_conc %>% 
             ungroup() %>% 
             dplyr::select(engine, time, value) %>%
             group_by(engine, time) %>% 
             filter(time <= 100) %>% 
             mean_qi(.width = c(0.5, 0.8, 0.95)) %>% 
             ggplot(aes(x = time, y = value, ymin = .lower, ymax = .upper)) + 
             geom_lineribbon() +
             facet_wrap(vars(engine)) +
             xlab("Time") + 
             ylab("Log Genetic Conc.") + 
             ggtitle("Log Genetic Concentrations by Engine") +
             my_theme + 
             theme(text = element_text(size = 18)) 
             
ggsave(here::here("figures", "conc_mc_compplot.pdf"), conc_plot, width = 10, height = 6)

# visualize individual models 
# scenario 1
library(tidyverse)
library(GGally)
library(gridExtra)
library(cowplot)
library(scales)
library(latex2exp)
source(here::here("src", "wastewater_functions.R"))

# scenario 1 --------------------------------------------------------------
# read in data, lets just look at 80% CI
seirr_rt_scenario1 <- read_csv(here::here("results", "seirr_student", "seirr_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "SEIRR") 
eirr_rt_scenario1 <- read_csv(here::here("results", "eirr_closed", "eirr_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR")

seir_rt_scenario1 <- read_csv(here::here("results", "seir_cases", "seir_cases_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "SEIR")

eir_rt_scenario1 <- read_csv(here::here("results", "eir_cases", "eir_cases_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate("EIR")

my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme())

seed = 1
seirr_scenario1_rt_plot <- seirr_rt_scenario1 %>%
  filter(seed == seed) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(time, true_rt), color = "coral1") + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ylim(c(0,3.5)) +
  ggtitle(str_c("SEIRR-ww")) +
  my_theme +
  ylab(TeX('$R_{t}$')) +
  theme(legend.position = "none",
        text = element_text(size = 18))

# swap for eirrc fit ------------------------------------------------------
scenario_sim = 1
sim = 1
full_simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_full_genecount_obsdata.csv")

full_simdata <- read_csv(full_simdata_address) %>% filter(seed == 1) %>% rename("true_rt" = "Rt",
                                                                                "r0" = "R0")
fitted_simdata_address <-  paste0("data/sim_data/scenario", scenario_sim, "_fitted_genecount_obsdata.csv")

fitted_simdata <- read_csv(fitted_simdata_address)
max_time <- max(fitted_simdata$time)

timevarying_quantiles <- read_csv(paste0("results/eirrc_closed/generated_quantities/posterior_timevarying_quantiles_scenario",
                                         sim,
                                         "_seed",
                                         seed,
                                         ".csv")) 

eirr_rt_scenario1 <- timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(week = time) %>%
  right_join(full_simdata, by = "week") %>%
  dplyr::select(week,time, new_time, true_rt, value, .lower, .upper, .width,.point, .interval) %>%
  filter(time <= max_time,
         week >= 0)

eirr_scenario1_rt_plot <- eirr_rt_scenario1 %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(time, true_rt), color = "coral1") + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ylim(c(0,3.5)) +
  ggtitle(str_c("EIRR-ww")) +
  my_theme + 
  theme(legend.position = c(0.6, 0.8),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size = 18)) +
  ylab(TeX('$R_{t}$')) 

# back to normal things ---------------------------------------------------
seir_scenario1_rt_plot <- seir_rt_scenario1 %>%
  filter(seed == seed) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(time, true_rt), color = "coral1") + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ylim(c(0,3.5)) +
  ggtitle(str_c("SEIR-cases")) +
  my_theme + 
  theme(legend.position = "none",
        text = element_text(size = 18)) +
  ylab(TeX('$R_{t}$')) 
  
eir_scenario1_rt_plot <- eir_rt_scenario1 %>%
  filter(seed == seed) %>%
  mutate(Truth = "") %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(time, true_rt, color = Truth)) + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ylim(c(0,3.5)) +
  ggtitle(str_c("EIR-cases")) +
  my_theme + 
  theme(legend.position = "none",
        text = element_text(size = 18)) +
  ylab(TeX('$R_{t}$')) 
  
point_plot <- eir_rt_scenario1 %>% 
                mutate(Truth = " ") %>%
  ggplot() +
  geom_point(aes(time, true_rt, color = Truth)) + 
  scale_color_manual(values = c( # actual line types used in the plot
    " " = "coral1"),
    drop = F) + 
  theme_bw() +
  theme(legend.background = element_rect(fill="transparent"),
    text = element_text(size = 18))

point_legend <- get_legend(point_plot)

line_plot <- eir_rt_scenario1 %>% 
             mutate(Median = " ") %>%
             ggplot() + 
             geom_line(aes(x = time, y = value, linetype = Median), size = 2) + 
            scale_linetype_manual(values = c( # actual line types used in the plot
    " " = "solid"),
    drop = F) + 
  
            theme_bw()+
  theme(legend.background = element_rect(fill="transparent"),
    text = element_text(size = 18))

pointline_plot <- eir_rt_scenario1 %>% 
  mutate(Median = " ") %>%
  mutate(Truth = " ") %>%
  ggplot() + 
  geom_line(aes(x = time, y = value, linetype = Median)) + 
  geom_point(aes(x = time, y = value, color = Truth)) +
  scale_linetype_manual(values = c( # actual line types used in the plot
    " " = "solid"),
    drop = F) + 
    scale_color_manual(values = c( # actual line types used in the plot
      " " = "coral1"),
      drop = F) + 
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent",
                                          colour = NA_character_), # necessary to avoid drawing panel outline
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          plot.background = element_rect(fill = "transparent",
                                         colour = NA_character_), # necessary to avoid drawing plot outline
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent"),
        text = element_text(size = 18))

pointline_legend <- get_legend(pointline_plot)

line_legend <- get_legend(line_plot)

fill_plot <- seirr_rt_scenario1 %>%
  filter(seed == 1) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(time, true_rt), color = "coral1") + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ylim(c(0,3.5)) +
  ggtitle(str_c("SEIRR-ww")) +
  my_theme +
  ylab("R[t]") +
  theme(
        text = element_text(size = 18))

pointline_legend <- plot_grid(
                              line_legend, 
                              point_legend,
                              nrow = 2
) 

combined_plot <- (seirr_scenario1_rt_plot + eirr_scenario1_rt_plot) / 
                (seir_scenario1_rt_plot + eir_scenario1_rt_plot + inset_element(
                  pointline_legend, 
                  left = 0.5, 
                  bottom = 0.5, 
                  right = unit(1, 'npc') - unit(1, 'cm'), 
                  top = unit(1, 'npc') - unit(1, 'cm')
                )) 


ggsave(here::here("figures", "scenario1_example_rt_plot.pdf"), combined_plot, width = 14, height =10)
                 
# visualize scenario data -------------------------------------------------

seed = 1
scenario1_simdata <- read_csv("data/sim_data/scenario1_fitted_genecount_obsdata.csv") %>%
  filter(seed == seed)
scenario1_truecurve <- read_csv(here::here("data", "sim_data", "scenario1_truecurve.csv"))

scenario1_casedata <- read_csv("data/sim_data/scenario1_fitted_cases_obsdata.csv") %>%
  filter(seed == seed)

start_time <- min(scenario1_simdata$time)
end_time <- max(scenario1_simdata$time)

truecurve <- scenario1_truecurve %>% 
             dplyr::select(integer_day, S, E, I, R1) %>% 
             pivot_longer(-integer_day) %>% 
             mutate(time = integer_day - start_time) %>% 
             filter(integer_day < 210) %>%
             rename(Compartment = name)

level_list <- c("S", "E", "I", "R1")

truecurve$Compartment <- factor(truecurve$Compartment, levels=level_list)

epi_curve <- truecurve %>%
  ggplot(aes(x = time, y = value, color = Compartment)) + 
  xlab("Time") + 
  ylab("Compartment Counts") + 
  ggtitle("Simulated Epidemic") +
  geom_point() + 
  theme_bw() +
  theme(legend.position = c(0.83, 0.78),
        legend.background = element_blank(),
        text = element_text(size = 18)) 


gene_plot <- scenario1_simdata %>% 
  dplyr::select(new_time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>% 
  pivot_longer(-new_time) %>% 
  ggplot(aes(x = new_time, y = value)) + 
  geom_point() + 
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Time") + 
  ylab("Log Genome Conc.") +
  ggtitle("Wastewater Data")

case_plot <- scenario1_casedata %>% 
  mutate(time = new_week * 7) %>% 
  ggplot(aes(x = time, y = log(total_cases))) + 
  geom_point() + 
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Time") + 
  ylab("Log Weekly Cases") +
  ggtitle("Case Data")

rt_plot <- scenario1_truecurve %>%
  mutate(new_time = time - start_time) %>% 
  filter(time <= end_time & time >= start_time) %>% 
  ggplot(aes(x = new_time, y = Rt)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Time") +
  ylab("Rt") +
  ggtitle("True Rt")


scenario1_sim_plot_paper <- (epi_curve + rt_plot) / (gene_plot + case_plot)
scenario1_sim_plot_paper

ggsave(here::here("figures", "scenario1_sim_plot_paper.pdf"), scenario1_sim_plot_paper, width = 11, height = 11)

# presentation plot -------------------------------------------------------
eirr_scenario1_rt_plot <- eirr_rt_scenario1 %>%
  mutate(new_time = time - start_time) %>% 
  ggplot(aes(new_time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(new_time, true_rt), color = "coral1") + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ylim(c(0,3.5)) +
  ggtitle(str_c("EIRR-ww")) +
  my_theme + 
  theme(legend.position = c(0.6, 0.8),
        legend.background = element_rect(fill = "transparent"),
        text = element_text(size = 18)) +
  ylab(TeX('$R_{t}$')) 

presentation_plot <- epi_curve / gene_plot / eirr_scenario1_rt_plot
ggsave(here::here("figures", "scenario1_presentationplot.pdf"), presentation_plot, width = 8, height = 10)

# compartment plot --------------------------------------------------------
true_E2I_data <- read_csv(here::here("data", "sim_data", "scenario1_truecasecounts.csv")) %>% 
                 mutate(new_time = time - min(time) + 1) %>% 
                 rename(
                        "true_value" = "E2I_transitions") %>% 
                 mutate(name = "Incidence") %>% 
                 dplyr::select(new_time, name, true_value)

scenario1_true_trajectory <- full_simdata %>% dplyr::select(E, I , R1, new_time) %>%
  pivot_longer(cols = -(new_time)) %>%
  rename("true_value" = "value") %>%
  bind_rows(true_E2I_data)

traj_names <- c("E", "I", "R1", "Incidence")

timevarying_quantiles$name[timevarying_quantiles$name == "incid"] = "Incidence"
posterior_trajectory <- timevarying_quantiles %>%
  filter(name %in% traj_names) %>%
  left_join(scenario1_true_trajectory, by = c("time" = "new_time", 
                                              "name"= "name"))



traj_plot <- posterior_trajectory %>%
  ggplot() +
  geom_lineribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper)) +
  geom_point(data = posterior_trajectory, mapping = aes(x = time, y = true_value), color = "coral1") +
  scale_fill_brewer(name = "Credible Interval Width") +
  theme_bw() +
  facet_wrap(.~name, scales = "free_y") +
  xlab("Time") + 
  ggtitle("EIRR Scenario 1 Latent Trajectories") + 
  ylab("Counts")

traj_plot

ggsave(here::here("figures", "scenario1_trajectoryplot.pdf"), traj_plot, width = 10, height = 8)

# fixed parameters --------------------------------------------------------

fixed_names <- c("gamma",
                 "nu",
                 "eta",
                 "rho_gene",
                 "rt_sigma",
                 "rt_init",
                 "tau",
                 "lambda",
                 "df",
                 "I_init",
                 "E_init",   
                 "R1_init")

true_vals <- c(1/4,
               1/7,
               1/18,
               NA,
               NA,
               0.88,
               0.5,
               NA,
               2.99,
               489,
               225,
               2075)

true_frame <- data.frame(name = fixed_names, true_value = true_vals)

eirr_scenario1_fixed_posterior_samples <- read_csv(here::here("results",
                                                              "eirrc_closed",
                                                              "generated_quantities",
                                                              paste0("posterior_fixed_samples_scenario", sim, "_seed", seed, ".csv"))) %>%
  left_join(true_frame, by = "name") %>%
  mutate(type = "posterior") %>%
  filter(name != ".draw")

eirr_prior_samples <- read_csv(here::here("results", "eirrc_closed", paste0("prior_samples_scenario", sim, "_seed", seed, ".csv"))) %>%
  mutate(type = "prior") %>%
  left_join(true_frame, by = "name")


eirr_scenario1_fixed_param_plot <- make_fixed_param_plot(eirr_scenario1_fixed_posterior_samples,
                                                         eirr_prior_samples) +
                                   ggtitle("EIRR Scenario 1 Fixed Parameters") +
                                   ylab("Type") + 
                                   xlab("Value")


eirr_scenario1_fixed_param_plot

ggsave(here::here("figures", "scenario1_fixedparamplot.pdf"), eirr_scenario1_fixed_param_plot, width = 10, height = 10)

# posterior predictive ----------------------------------------------------

eirr_scenario1_post_pred_intervals <- read_csv(here::here("results",
                                                          "eirrc_closed",
                                                          "posterior_predictive",
                                                          paste0("posterior_predictive_intervals_scenario", sim, "_seed", seed, ".csv")))

true_data <- scenario1_simdata %>%
  dplyr::select(new_time, 
                log_gene_copies1, 
                log_gene_copies2, 
                log_gene_copies3) 

eirr_scenario1_post_pred_intervals <- eirr_scenario1_post_pred_intervals %>%
  left_join(true_data, by = c("new_time"))

posterior_predictive_plot <- eirr_scenario1_post_pred_intervals %>%
  ggplot() +
  geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
  geom_line(aes(x = new_time, y = value)) +
  geom_point(mapping = aes(x = new_time, y = log_gene_copies1), color = "coral1") +
  geom_point(mapping = aes(x = new_time, y = log_gene_copies2), color = "coral1") +
  geom_point(mapping = aes(x = new_time, y = log_gene_copies3), color = "coral1") +
  scale_fill_brewer(name = "Credible Interval Width") +
  theme_bw() + 
  ggtitle("EIRR Scenario 1 Posterior Predictive") +
  xlab("Time") +
  ylab("Log Genetic Conc.")


posterior_predictive_plot
ggsave(here::here("figures", "scenario1_posteriorpredictive.pdf"), posterior_predictive_plot, width = 8, height = 5)

# huisman visualization ---------------------------------------------------

# scenario 1 --------------------------------------------------------------
# read in data, lets just look at 80% CI
huisman_rt_scenario1 <- read_csv(here::here("results", "huisman", "huisman_scenario1_allseeds_rt_quantiles.csv"))

my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme())


seed = 1
huisman_scenario1_rt_plot <- huisman_rt_scenario1 %>%
  filter(seed == seed) %>%
  mutate(time = time - min(time) + 1) %>%
  ggplot(aes(time, median_R_mean, ymin = median_R_lowHPD, ymax = median_R_highHPD)) +
  geom_lineribbon() +
  geom_point(aes(time, true_rt), color = "coral1") + 
  scale_y_continuous("Rt", label = comma) +
  scale_x_continuous(name = "Time") +
  ggtitle(str_c("Huisman (WW)")) +
  my_theme +
  ylab(TeX('$R_{t}$')) +
  theme(legend.position = "none",
        text = element_text(size = 18))
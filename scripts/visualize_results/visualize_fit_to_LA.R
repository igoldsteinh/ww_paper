# visualize various model fits to LA Data
library(tidyverse)
library(GGally)
library(gridExtra)
library(cowplot)
library(scales)
library(epidemia)
source(here::here("src", "wastewater_functions.R"))
source(here::here("src", "rt_estim_gamma_functions.R"))
# read in LA data ---------------------------------------------------------

seed = 1
sim = "real"

real_data <- read_csv("data/LA_daily_data_feb2022.csv")

sim = "real"
seed = 1
date_week_crosswalk <- real_data %>% 
  dplyr::select(date, epi_week, new_time) %>%
  mutate(time = epi_week - 27)

# visualized data ----------------------------------------------------------
ww_data_plot <- real_data %>% 
  dplyr::select(date, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>% 
  pivot_longer(-date) %>% 
  filter(value > 0) %>%
  ggplot(aes(x = date, y = value)) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 18)) +
  scale_y_continuous("Log Genome Conc.", breaks = 0:20) +
  scale_x_date("Date", breaks = "months") +
  ggtitle("LA Wastewater Jul 21 - Feb 22")

case_data <- read_csv(here::here("data", "LA_weekly_case_data.csv")) 

case_plot <- case_data %>% 
  filter(max_date <= max(real_data$date) & max_date >= min(real_data$date)) %>%
  ggplot(aes(x = max_date, y = log(total_cases))) + 
  geom_point() + 
  geom_line() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 18)) +
  scale_y_continuous("Log Weekly Cases", labels = comma) +
  scale_x_date("Date", breaks = "months") +
  ggtitle("LA Cases Jul 21 - Feb 22")

combined_data_plot <- ww_data_plot + case_plot

combined_data_plot

ggsave(here::here("figures", "LA_data_jul21_feb22.pdf"), combined_data_plot, width = 10, height = 4)

real_data$epi_week[real_data$year == 2022] <- real_data$epi_week[real_data$year == 2022] + 52
date_week_crosswalk <- real_data %>% 
  dplyr::select(date, epi_week, new_time) %>%
  mutate(time = epi_week - 27)

# EIR ---------------------------------------------------------------------
timevarying_quantiles_eir <- read_csv(paste0("results/eir_cases/generated_quantities/posterior_timevarying_quantiles_scenario",
                                         sim,
                                         "_seed",
                                         seed,
                                         ".csv")) 


rt_quantiles_eir <- timevarying_quantiles_eir %>%
  filter(name == "rt_t_values") %>%
  left_join(date_week_crosswalk, by = "time") %>%
  dplyr::select(time, date, epi_week, value, .lower, .upper, .width,.point, .interval) 

fill_date = as.Date("2021-10-11")
fill_week = 41

rt_quantiles_eir$date[is.na(rt_quantiles_eir$date)] <- fill_date
rt_quantiles_eir$epi_week[is.na(rt_quantiles_eir$epi_week)] <- fill_week
# visualize results
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw()+
  theme(legend.position = "bottom"))

eir_realdata_rt_plot_seed1 <- rt_quantiles_eir %>%
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  scale_y_continuous("Rt", label = comma, breaks = c(0:5), limits = c(0,5)) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIR-cases Posterior Rt") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        text = element_text(size = 18)) 

eir_realdata_rt_plot_seed1

# Huisman -----------------------------------------------------------------
la_huisman_rt <- read_csv(here::here("results", "huisman", "huisman_la_rt_quantiles.csv"))

missing_dates <- la_huisman_rt %>% slice(1) %>% mutate(date = as.Date("2022-02-23"),
                                                       median_R_mean = NA,
                                                       median_R_lowHPD = NA,
                                                       median_R_highHPD = NA)
huisman_realdata_rt_plot <- la_huisman_rt %>%
  mutate(.width = 0.95) %>%
  ggplot(aes(date, median_R_mean, ymin = median_R_lowHPD, ymax = median_R_highHPD)) +
  geom_lineribbon() +
  scale_y_continuous("Rt", label = comma) +
  scale_x_date(name = "Date", date_breaks = "month", limits = as.Date(c("2021-07-04", "2022-02-23"))) +
  ggtitle("Huisman Posterior Rt") +
  my_theme + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        ,
        text = element_text(size = 18)) 

# EIRRC --------------------------------------------------------------------
seed = 1
sim = "real"

timevarying_quantiles <- read_csv(paste0("results/eirrc_closed/generated_quantities/posterior_timevarying_quantiles_scenario",
                                         sim,
                                         "_seed",
                                         seed,
                                         ".csv")) 

rt_quantiles_eirr <- timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  left_join(date_week_crosswalk, by = "time") %>%
  dplyr::select(time, date, epi_week, value, .lower, .upper, .width,.point, .interval) 

fill_date = as.Date("2021-10-11")
fill_week = 41

rt_quantiles_eirr$date[is.na(rt_quantiles_eirr$date)] <- fill_date
rt_quantiles_eirr$epi_week[is.na(rt_quantiles_eirr$epi_week)] <- fill_week
# visualize results
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme(legend.position = "bottom"))

eirrc_realdata_rt_plot_seed1 <- rt_quantiles_eirr %>%
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  scale_y_continuous("Rt", label = comma, breaks = c(0:5), limits = c(0,5)) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIRR-ww Posterior Rt") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.4, 0.75),
        text = element_text(size = 18),
        legend.background = element_blank())

eirrc_realdata_rt_plot_seed1

# estimnormal (epidemia) --------------------------------------------------
dat <- read_csv(here::here("data", "LA_EIR_data.csv"))

dat$start_date <- as.Date(dat$min_date, format = "%m/%d/%y")
fm2 <- read_rds(here::here("results", "estimnormal", "estimnormal_LA_fit.rds"))

dat <- dat %>% 
  mutate(time = row_number())
start_date <- min(dat$time)
max_date <- max(dat$time)

estimnormal_posterior_rt <-
  posterior_rt(fm2)[["draws"]] %>%
  data.frame() %>%
  `colnames<-`(start_date:(max_date + 1)) %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws,
               names_to = "epidemia_time",
               values_to = "value",
               names_transform = list(epidemia_time = as.integer)
  ) %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(
    method = "estim_normal",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

plot_data <- dat %>% 
  left_join(estimnormal_posterior_rt, by = "time")

rt_posterior = plot_data %>%
  mutate(year = year(max_date),
         epi_week = ifelse(year == 2021, epiweek(max_date), epiweek(max_date) + 52)) %>%
  left_join(date_week_crosswalk %>% dplyr::select(-time), by = "epi_week") %>%
  dplyr::select(time, date, epi_week, value, .lower, .upper, .width) 

fill_date = as.Date("2021-10-11")
fill_week = 41

rt_posterior$date[is.na(rt_posterior$date)] <- fill_date
rt_posterior$epi_week[is.na(rt_posterior$epi_week)] <- fill_week

estimnormal_ladata_plot <- rt_posterior %>%
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous("Rt", label = comma, breaks = 0:5, limits = c(0, 5)) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("Epidemia Posterior Rt ") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        text = element_text(size = 18)) 

estimnormal_ladata_plot

# rt estim gamma fit to real data -----------------------------------------
file_name <- "LA_estimgamma_posterior_alt.rds"
posterior <- read_rds(here::here("results", "estimgamma", file_name))

data <- read_csv(here::here("data", "LA_EIR_data.csv"))
data <- data %>%
  rename("time" = "new_week")

trace <- rstan::traceplot(posterior, pars = "lp__") +
  ggtitle("lp traceplot")


rt_posterior = summarise_realdata_rt_estimgamma(posterior,
                                                data,
                                                1,
                                                include_chains= c(1,2,3))

# visualize results
# all credit to Damon Bayer for plot functions 

real_data <- read_csv("data/LA_daily_data_feb2022.csv")

rt_posterior = rt_posterior %>%
  dplyr::select(-date) %>%
  mutate(year = year(max_date),
         epi_week = ifelse(year == 2021, epiweek(max_date), epiweek(max_date) + 52)) %>%
  left_join(date_week_crosswalk, by = "epi_week") %>%
  dplyr::select(time, date, epi_week, rt_median, .lower, .upper, .width) 

fill_date = as.Date("2021-10-11")
fill_week = 41

rt_posterior$date[is.na(rt_posterior$date)] <- fill_date
rt_posterior$epi_week[is.na(rt_posterior$epi_week)] <- fill_week
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_bw(),
  theme(legend.position = "bottom"))

estimgamma_ladata_plot <- rt_posterior %>%
  ggplot(aes(date, rt_median, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous("Rt", label = comma, breaks = c(0:5), limits = c(0, 5)) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("Rt-estim-gamma Posterior Rt ") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        text = element_text(size = 18)) 

estimgamma_ladata_plot

# paper main plot  --------------------------------------------------------
main_plot <- huisman_realdata_rt_plot + eir_realdata_rt_plot_seed1 + eirrc_realdata_rt_plot_seed1
ggsave(here::here("figures", "main_real_plot.pdf"), main_plot, width = 14, height = 6)

# paper appendix plot -----------------------------------------------------
app_eir <- eir_realdata_rt_plot_seed1 + theme(
  text = element_text(size = 18),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) +
  xlab("")

app_eirr <-  eirrc_realdata_rt_plot_seed1 + theme(
  text = element_text(size = 18),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) +
  xlab("")
appendix_plot <- (app_eir + app_eirr) / (estimgamma_ladata_plot + estimnormal_ladata_plot)

ggsave(here::here("figures", "app_real_plot.pdf"), appendix_plot, width = 10, height = 10)

# WNAR presentation plot --------------------------------------------------
WNAR_estimgamma <- estimgamma_ladata_plot + ggtitle("LA Posterior Rt (Cases)") 
WNAR_eirr <- eirrc_realdata_rt_plot_seed1 + ggtitle("LA Posterior Rt (Wastewater)")
WNAR_plot <- WNAR_estimgamma + WNAR_eirr

WNAR_eir <- eir_realdata_rt_plot_seed1 + ggtitle("LA Posterior Rt (Cases no Tests)")
WNAR_case_plot <- WNAR_estimgamma + WNAR_eir
ggsave(here::here("figures", "WNAR_real_plot.pdf"), WNAR_plot, width = 10, height = 6)
ggsave(here::here("figures", "WNAR_case_plot.pdf"), WNAR_case_plot, width = 10.5, height = 6)

# calculate case detection rate from incidence ----------------------------

daily_case_data <- read_csv(here::here("data", "LA_daily_case_data.csv"))

# relative detection ------------------------------------------------------

gq_address <- paste0("results/eirrc_closed/generated_quantities/generated_quantities_scenario", 
                     sim, 
                     "_seed", 
                     seed,
                     ".csv")

posterior_samples <- read_csv(gq_address) 

max_iteration = max(posterior_samples$.iteration)
min_iteration = round(max_iteration/2)

posterior_gq_samples_all <- posterior_samples  %>%
  pivot_longer(-c(iteration, chain))

timevarying_posterior <-
  posterior_gq_samples_all %>%
  filter(str_detect(name, "\\[\\d+\\]")) %>%
  mutate(time = name %>%
           str_extract("(?<=\\[)\\d+(?=\\])") %>%
           as.numeric(),
         name = name %>%
           str_extract("^.+(?=\\[)") %>%
           str_remove("data_")) %>%
  group_by(name, time)

C_vals <- timevarying_posterior %>% filter(name == "C") %>%
  group_by(chain, iteration) %>% 
  left_join(daily_case_data, by = c("time" = "new_time")) %>%
  mutate(detect = (cases * 0.48)/value,
         prev_detect = lag(detect),
         rel_change = detect/prev_detect,
         daily_incidence = value - lag(value)
  )

# lets do this at a weekly scale tho
weekly_rel_change <- C_vals %>%
  group_by(chain, iteration, year, epi_week) %>% 
  summarise(total_val = sum(value),
            total_incid= sum(daily_incidence),
            total_cases = sum(cases) * 0.48,
            total_tests = sum(tests) * 0.48,
            total_detect = total_cases/total_incid,
            rho = total_cases/(total_tests * total_incid)) 

first_detects <- weekly_rel_change %>% 
                 group_by(chain, iteration) %>% 
                 filter(epi_week == 27) %>% 
                 dplyr::select(chain, iteration, first_detect = total_detect, first_rho = rho)
  

weekly_rel_change <- weekly_rel_change %>% 
  group_by(chain, iteration) %>% 
  left_join(first_detects, by = c("chain", "iteration")) %>% 
  mutate(rel_change = total_detect/first_detect,
         rel_change_rho = rho/first_rho) 

weekly_rel_change_quantiles <- weekly_rel_change %>%
  group_by(year, epi_week) %>% 
  median_qi(rel_change, .width = c(0.5, 0.8, 0.95)) %>%
  mutate(epi_week = ifelse(year == 2021, epi_week, epi_week + 52)) %>%
  right_join(date_week_crosswalk, by = "epi_week")

weekly_detect_quantiles <- weekly_rel_change %>%
  group_by(year, epi_week) %>% 
  median_qi(total_detect, .width = c(0.5, 0.8, 0.95))  %>%
  mutate(epi_week = ifelse(year == 2021, epi_week, epi_week + 52)) %>%
  right_join(date_week_crosswalk, by = "epi_week")

weekly_rho_quantiles <- weekly_rel_change %>%
  group_by(year, epi_week) %>% 
  median_qi(rho, .width = c(0.5, 0.8, 0.95)) %>%
  mutate(epi_week = ifelse(year == 2021, epi_week, epi_week + 52)) %>%
  right_join(date_week_crosswalk, by = "epi_week")


weekly_relrho_quantiles <- weekly_rel_change %>%
  group_by(year, epi_week) %>% 
  median_qi(rel_change_rho, .width = c(0.5, 0.8, 0.95)) %>%
  mutate(epi_week = ifelse(year == 2021, epi_week, epi_week + 52)) %>%
  right_join(date_week_crosswalk, by = "epi_week")

weekly_incid_quantiles <- weekly_rel_change %>%
  group_by(year, epi_week) %>% 
  median_qi(total_incid, .width = c(0.5, 0.8, 0.95)) %>%
  mutate(epi_week = ifelse(year == 2021, epi_week, epi_week + 52)) %>%
  right_join(date_week_crosswalk, by = "epi_week")

plot_weekly_change <- weekly_rel_change_quantiles %>% 
  drop_na() %>%
  ggplot(aes(x = date, y = rel_change, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() + 
  scale_y_continuous("Case Detection Ratio") +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIRRC Posterior Case Detection Ratio") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.25, 0.75)) 

plot_weekly_detect <- weekly_detect_quantiles %>% 
  drop_na() %>%
  ggplot(aes(x = date, y = total_detect, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() + 
  scale_y_continuous("Case Detection Rate") +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIRR Posterior Detection Rate") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") 

plot_rho_quantiles <- weekly_rho_quantiles %>% 
  drop_na() %>%
  ggplot(aes(x = date, y = rho, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() + 
  scale_y_continuous("Noramlized CDR", labels = comma) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIRR Posterior Normalized CDR") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.42, 0.75),
        legend.background = element_blank()) 

plot_relrho_quantiles <- weekly_relrho_quantiles %>% 
  drop_na() %>%
  ggplot(aes(x = date, y = rel_change_rho, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() + 
  scale_y_continuous("Rho Ratio") +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIRR Posterior Rho Ratio") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.25, 0.75),
        legend.background = element_blank()) 

plot_incid_quantiles <- weekly_incid_quantiles %>% 
  drop_na() %>%
  ggplot(aes(x = date, y = total_incid, ymin = .lower, ymax = .upper)) + 
  geom_lineribbon() + 
  scale_y_continuous("Weekly Incidence") +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("EIRRC Posterior Weekly Incidence") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(0.25, 0.75),
        legend.background = element_blank()) 

case_plus_rho <- plot_weekly_detect + plot_rho_quantiles  
ggsave(here::here("figures", "la_case_detect.pdf"), case_plus_rho, width = 12, height = 4)

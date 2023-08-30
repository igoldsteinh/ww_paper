# Calculate realistic noise parameters for simulation purposes by fitting 
# student-t splines to real data
library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
source("src/wastewater_functions.R")
LA_data <- read_csv(here::here("data", "LA_daily_data_feb2022.csv"))


model_data <- LA_data %>% 
            dplyr::select(date, gene_copy_1, gene_copy_2, gene_copy_3) %>%
            pivot_longer(cols = all_of(c("gene_copy_1", 
                             "gene_copy_2",
                             "gene_copy_3")), 
             names_to = "copy") %>%
            filter(value != -1) %>%
            mutate(value_scaled = value/1000,
                   log_values = log(value_scaled),
                   days_since = date - min(LA_data$date))

model_data$days_since <- as.numeric(model_data$days_since)

# plot data ---------------------------------------------------------------
data_plot <- model_data %>%
  ggplot(aes(x = date, y = log_values, color = copy)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  ylab("Log concentration (copies/ml)") + 
  xlab("Time")



# lets just fit a spline to it and see what happens -----------------------



seed = 456
iter = 4000
warmup = 1000
thin = 10
refresh = 0
adapt_delta = 0.99

norm_spline <- brm(bf(log_values ~ s(days_since)),
                      data = model_data, family = gaussian(), cores = 4, seed = seed,
                      iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                      control = list(adapt_delta = adapt_delta))

write_rds(norm_spline, here::here("results", 
                                  "splines", 
                                  "la_normspline_2021_2022.rds"))

normspline_output <- read_rds(here::here("results", 
                                         "splines",
                                         "la_normspline_2021_2022.rds"))


normpp_draws <- add_predicted_draws(model_data, normspline_output, cores = 4) %>%
  ungroup() %>%
  dplyr::select(days_since, .prediction) %>%
  group_by(days_since) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

normpp_graph_data <- model_data %>%
  left_join(normpp_draws, by = "days_since")

normspline_case_plot <- normpp_graph_data %>%
  ggplot(aes(x = date, y = .prediction, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(aes(x = date,fill = fct_rev(ordered(.width)))) +
  theme_bw() +
  ylab("Log Concentration") +
  xlab("Date") +
  geom_point(aes(x = date, y = log_values),  col = "#F8766D") +
  geom_line(aes(x = date, y = log_values, col = "#F8766D")) +
  ggtitle("N1 Norm Spline PP") +
  theme(text = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_date(date_labels = "%b %Y") +
  xlab("") +
  scale_fill_brewer(name = "Credible Interval Width") +
  scale_linetype(name = NULL) +
  theme(legend.background = element_rect(fill="transparent")) +
  xlab("Date")



norm_e_draws <- add_linpred_draws(model_data, normspline_output, cores = 4) %>%
  ungroup() %>%
  dplyr::select(days_since, .linpred) %>%
  group_by(days_since) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

norm_fitted_graph_data <- model_data %>%
  left_join(norm_e_draws, by = "days_since")

norm_fitted_plot <- norm_fitted_graph_data %>%
  ggplot(aes(x=date, y=.linpred, ymax = .upper, ymin = .lower)) +
  geom_lineribbon(aes(x = date,fill = fct_rev(ordered(.width)))) +
  theme_bw()+
  ylab("Log Counts")+
  xlab("Date")+
  geom_point(aes(x = date, y = log_values), color = "Red") +
  scale_fill_brewer(name = "Credible Interval Width") +
  theme(legend.position = "None")+
  ggtitle("Norm Spline (not posterior predictive)") 

# student t spline --------------------------------------------------------


seed = 456
iter = 4000
warmup = 1000
thin = 10
refresh = 0
adapt_delta = 0.99

student_spline <- brm(bf(log_values ~ s(days_since)),
                   data = model_data, family = student(), cores = 4, seed = seed,
                   iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                   control = list(adapt_delta = adapt_delta))

write_rds(student_spline, here::here("results", 
                                  "splines", 
                                  "la_studentspline_2021_2022.rds"))

studentspline_output <- read_rds(here::here("results", 
                                         "splines",
                                         "la_studentspline_2021_2022.rds"))


studentpp_draws <- add_predicted_draws(model_data, studentspline_output, cores = 4) %>%
  ungroup() %>%
  dplyr::select(days_since, .prediction) %>%
  group_by(days_since) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

studentpp_graph_data <- model_data %>%
  left_join(studentpp_draws, by = "days_since")

studentpline_case_plot <- studentpp_graph_data %>%
  ggplot(aes(x = date, y = .prediction, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(aes(x = date,fill = fct_rev(ordered(.width)))) +
  theme_bw() +
  ylab("Log Concentration") +
  xlab("Date") +
  geom_point(aes(x = date, y = log_values),  col = "#F8766D") +
  geom_line(aes(x = date, y = log_values, col = "#F8766D")) +
  ggtitle("N1 Student Spline PP") +
  theme(text = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_date(date_labels = "%b %Y") +
  xlab("") +
  scale_fill_brewer(name = "Credible Interval Width") +
  scale_linetype(name = NULL) +
  theme(legend.background = element_rect(fill="transparent")) +
  xlab("Date")



student_e_draws <- add_linpred_draws(model_data, studentspline_output, cores = 4) %>%
  ungroup() %>%
  dplyr::select(days_since, .linpred) %>%
  group_by(days_since) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

student_fitted_graph_data <- model_data %>%
  left_join(student_e_draws, by = "days_since")

student_fitted_plot <- student_fitted_graph_data %>%
  ggplot(aes(x=date, y=.linpred, ymax = .upper, ymin = .lower)) +
  geom_lineribbon(aes(x = date,fill = fct_rev(ordered(.width)))) +
  theme_bw()+
  ylab("Log Counts")+
  xlab("Date")+
  geom_point(aes(x = date, y = log_values), color = "Red") +
  scale_fill_brewer(name = "Credible Interval Width") +
  theme(legend.position = "None")+
  ggtitle("Student Spline (not posterior predictive)") 

compare_pp <- studentpline_case_plot + normspline_case_plot

compare_posterior_spline <- student_fitted_plot + norm_fitted_plot

# case spline
ca_data <- read_csv(here::here("data", "covid19cases_test.csv"))

# assume la has 10 million
# then the sim has 100K so we need to know the kappa at 0.01 of cases
la_case_data <- create_weekly_data(ca_data, 
                              county = "Los Angeles",
                              start_sunday = "2021-07-04",
                              end_sunday = "2022-02-27" ) %>%
                mutate(total_cases = round(0.01 * total_cases))

case_spline <- run_nb_spline(la_case_data)

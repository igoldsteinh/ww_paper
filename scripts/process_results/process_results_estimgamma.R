# process results of rt-estim-gamma
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
library(lubridate)
library(scales)
library(GGally)
library(cowplot)
library(gridExtra)

source("src/wastewater_functions.R")
source("src/rt_estim_gamma_functins.R")
file_name <- "LA_estimgamma_posterior.rds"
posterior <- read_rds(here::here("results", "estimgamma", file_name))

data <- read_csv(here::here("data", "LA_EIR_data.csv"))
data <- data %>%
  rename("time" = "new_week")

trace <- rstan::traceplot(posterior, pars = "lp__") +
  ggtitle("lp traceplot")

ggsave(here::here("results", "estimgamma", str_c("scenario", sim, "_trace", ".pdf", sep = "")), 
       plot = trace, 
       width = 5, 
       height = 5)

rt_posterior = summarise_realdata_rt_estimgamma(posterior,
                                 data,
                                 1,
                                 include_chains= c(1,2,3))

# visualize results
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

estimgamma_ladata_plot <- rt_posterior %>%
  ggplot(aes(max_date, rt_median, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous("Rt", label = comma) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("LA Jul 21 - Feb 22 Posterior Rt Seed 1") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90)) 

estimgamma_ladata_plot
write_csv(rt_posterior, here::here("results", "estimgamma", paste0("rt_posterior_sim", sim, ".csv")))

# alt prior ---------------------------------------------------------------
file_name <- "LA_estimgamma_posterior_alt.rds"
alt_county_posterior <- read_rds(here::here("results", "estimgamma", file_name))

data <- read_csv(here::here("data", "LA_EIR_data.csv"))
data <- data %>%
  rename("time" = "new_week")


trace <- rstan::traceplot(alt_county_posterior, pars = "lp__") +
  ggtitle("lp traceplot")

ggsave(here::here("results", "estimgamma", str_c("scenario", sim, "_trace", ".pdf", sep = "")), 
       plot = trace, 
       width = 5, 
       height = 5)

rt_posterior = summarise_realdata_rt_estimgamma(alt_county_posterior,
                                                data,
                                                1,
                                                include_chains= c(1,2,3))

# visualize results
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

estimgamma_ladata_plot <- rt_posterior %>%
  ggplot(aes(max_date, rt_median, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous("Rt", label = comma) +
  scale_x_date(name = "Date", date_breaks = "month") +
  ggtitle("LA Jul 21 - Feb 22 Posterior Rt Seed 1") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90)) 

estimgamma_ladata_plot

write_csv(rt_posterior, here::here("results", "estimgamma", paste0("rt_posterior_sim", sim, ".csv")))

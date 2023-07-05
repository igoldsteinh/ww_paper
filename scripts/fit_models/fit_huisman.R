# fit huisman model 
library(tidyverse)
library(lubridate)
source(here::here("src", "huisman_functions.R"))
source(here::here("src", "wastewater_functions.R"))


source(here::here("src", "2_utils_getInfectionIncidence.R"))
source(here::here("src", "3_utils_doReEstimates.R"))

args <- commandArgs(trailingOnly=TRUE)


if (length(args) == 0) {
  snum = 1
  seed_val = 1
  scenario_snum = snum
} else {
  snum <- as.integer(args[1])
  seed_val <- as.integer(args[2])
  
}

all_data <- read_csv("data/sim_data/scenario1_fitted_genecount_obsdata.csv")

data <- all_data %>% filter(seed == seed_val)

start_time <- min(data$time)


# set up true data --------------------------------------------------------
fake_date_zero <- ymd("20220109")

fake_dates <- seq.Date(fake_date_zero, length.out = 133, by = "day")

full_simdata_address <- paste0("data/sim_data/scenario", scenario_snum, "_full_genecount_obsdata.csv")

full_simdata <- read_csv(full_simdata_address) %>% filter(seed == 1) %>% rename("true_rt" = "Rt",
                                                                                "r0" = "R0")

full_simdata <- cbind(full_simdata, date = fake_dates)




rt_quantiles <- NULL
# attach results to data ---------------------------------------------
for (i in 1:100){
  seed_val = i
  
  data <- all_data %>% filter(seed == seed_val)
  
  scenario1_huisman_rt <- calculate_huisman_rt(data, sim = TRUE)
  
  sim_rt_quantiles <- scenario1_huisman_rt %>% 
    dplyr::select(date, median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
    right_join(full_simdata, by = "date") %>%
    dplyr::select(date,time, true_rt, median_R_mean, median_R_highHPD, median_R_lowHPD) %>% 
    drop_na() %>% 
    mutate(seed = seed_val)
  
  rt_quantiles <- bind_rows(rt_quantiles, sim_rt_quantiles)
  
}

write_csv(rt_quantiles, here::here("results", "huisman", "huisman_scenario1_allseeds_rt_quantiles.csv"))


# fit huisman to real data ------------------------------------------------

la_data <- read_csv("data/LA_daily_data_feb2022.csv")
set.seed(1234)
la_huisman_rt <- calculate_huisman_rt(la_data, sim = FALSE)

write_csv(la_huisman_rt, here::here("results", "huisman", "huisman_la_rt_quantiles.csv"))
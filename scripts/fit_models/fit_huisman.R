# fit huisman model 
library(tidyverse)
library(lubridate)
source(here::here("src", "huisman_functions.R"))
source(here::here("src", "wastewater_functions.R"))


source(here::here("src", "2_utils_getInfectionIncidence.R"))
source(here::here("src", "3_utils_doReEstimates.R"))

args <- commandArgs(trailingOnly=TRUE)


if (length(args) == 0) {
  sim = 1
  seed = 1
  scenario_sim = sim
} else {
  sim <- as.integer(args[1])
  seed <- as.integer(args[2])
  
}

all_data <- read_csv("data/sim_data/scenario1_fitted_genecount_obsdata.csv")

data <- all_data %>% filter(seed == seed)

start_time <- min(data$time)


# set up true data --------------------------------------------------------
fake_date_zero <- ymd("20220109")

fake_dates <- seq.Date(fake_date_zero, length.out = 133, by = "day")

full_simdata_address <- paste0("data/sim_data/scenario", scenario_sim, "_full_genecount_obsdata.csv")

full_simdata <- read_csv(full_simdata_address) %>% filter(seed == 1) %>% rename("true_rt" = "Rt",
                                                                                "r0" = "R0")

full_simdata <- cbind(full_simdata, date = fake_dates)




rt_quantiles <- NULL
# attach results to data ---------------------------------------------
for (i in 1:100){
  seed = i
  data <- all_data %>% filter(seed == seed)
  
  scenario1_huisman_rt <- calculate_huisman_rt(data, sim = TRUE, mean_si = 11, std_si = sqrt(65))
  
  sim_rt_quantiles <- scenario1_huisman_rt %>% 
    dplyr::select(date, median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
    right_join(full_simdata, by = "date") %>%
    dplyr::select(date,time, true_rt, median_R_mean, median_R_highHPD, median_R_lowHPD) %>% 
    drop_na() %>% 
    mutate(seed = seed)
  
  rt_quantiles <- bind_rows(rt_quantiles, sim_rt_quantiles)
  
}

write_csv(rt_quantiles, here::here("results", "huisman", "huisman_scenario1_allseeds_rt_quantiles.csv"))


# fit huisman to real data ------------------------------------------------

la_data <- read_csv("data/LA_daily_data_feb2022.csv")
set.seed(1234)
la_huisman_rt <- calculate_huisman_rt(la_data, sim = FALSE,  mean_si = 9.7, std_si = sqrt(2*(9.7/2)^2))

write_csv(la_huisman_rt, here::here("results", "huisman", "huisman_la_rt_quantiles.csv"))



# sanity checking huisman on their own data -------------------------------
## Wastewater data - SCAN Pilot Project ####


raw_data <- read.csv(here::here("data", "CA_ww_data.csv")) %>%
  mutate(n_orig = n_gene,
         s_orig = s_gene,
         orf1a_orig = ORF1a)

## Normalisation ####
norm_min <- raw_data %>%
  filter(n_gene != 0) %>%
  pull(n_gene) %>% min()

## Final Cleaning ####
#  gaps in the data prior to 15 Nov.
# stop on 18.3.2021
# 03 January, 18 February do not meet quality control
raw_data$date <- as.Date(raw_data$date)
clean_data <- raw_data %>%
  mutate(orig_data = TRUE) %>%
  filter(date >= as_date('2020-11-15'),
         date <= as_date('2021-03-18'),
         date != as_date('2021-01-03'),
         date != as_date('2021-02-18') ) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
  mutate(across(c(n_gene, s_gene, ORF1a), ~ ./norm_min ) ) %>%
  mutate(county = 'Santa Clara',
         region = county)

###### Deconvolve #####

deconv_result <- data.frame()
result <- data.frame()
set.seed(50)
mean_si = 3
std_si = 2.4
for (incidence_var_i in c('n_gene', 's_gene', 'ORF1a')){
  new_deconv_result = deconvolveIncidence(clean_data, 
                                          incidence_var = incidence_var_i,
                                          getCountParams('incubation'), 
                                          getCountParams('benefield'),
                                          smooth_param = TRUE, n_boot = 50) %>% #n_boot = 1000 in paper
    mutate(data_type = incidence_var_i)
  
  new_result = getReBootstrap(new_deconv_result, mean_si, std_si)
  new_result = new_result %>%
    mutate(data_type = incidence_var_i)
  new_result['variable'] = incidence_var_i
  
  deconv_result = bind_rows(deconv_result, new_deconv_result)
  result = bind_rows(result, new_result)
  
}



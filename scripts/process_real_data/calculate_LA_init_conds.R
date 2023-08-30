# file for showing work for choosing LA initial conditions
# and also for estimating rho_gene
library(tidyverse)
library(lubridate)
# Initial conditons for LA ------------------------------------------------


# rough rule of thumb for initial conditions
# S = (popsize - # vaccinated) / 2
#        E = Number of observed cases in the previous 2 weeks * 4 * 1 / 3
#      I = Number of observed cases in the previous 2 weeks * 4 * 2 / 3


popsize = 4.8E6
# http://publichealth.lacounty.gov/media/Coronavirus/vaccine/vaccine-dashboard.htm
# About 60% of pop was back to being susceptible
S = (popsize - 0.4 * popsize)

# http://dashboard.publichealth.lacounty.gov/covid19_surveillance_dashboard/
la_cases <-read_csv("data/covid19cases_test.csv") %>% 
  mutate(case_date = as.Date("2021-07-04") - days(11),
         new_case_date = as.Date("2021-07-04") - days(9),
         recover_date = case_date - days(18),
         new_recover_date = case_date - days(36),
         newnew_recover_date = case_date - days(27),
         recover_date_17 = new_case_date - days(27),
         recover_date_19 = new_case_date - days(18),
         recover_date_21 = case_date - days(10)) %>%
  filter(date >= new_recover_date & date <= "2021-10-30") %>%
  filter(area == "Los Angeles") %>%
  mutate(year = year(date),
         epi_week = epiweek(date)) 
la_cases <-  la_cases %>%
  mutate(case_date = as.Date("2021-07-04") - days(11),
         new_case_date = as.Date("2021-07-04") - days(9),
         recover_date = case_date - days(18),
         new_recover_date = case_date - days(36),
         newnew_recover_date = case_date - days(27),
         recover_date_19 = new_case_date - days(18),
         recover_date_21 = case_date - days(10),
         current_cases = date >= case_date & date < "2021-07-04",
         new_current_cases = date >= new_case_date & date < "2021-07-04",
         recovered_cases = date >= recover_date & date < case_date,
         new_recovered_cases = date >= new_recover_date & date < case_date,
         newnew_recovered_cases = date >= newnew_recover_date & date < case_date,
         newnewnew_recovered_cases = date >= recover_date_17 & date < new_case_date,
         recovered_cases19 = date >= recover_date_19 & date < new_case_date,
         recovered_cases21 = date >= recover_date_21 & date < case_date) 

cases_lasttwo <- la_cases %>% 
  filter(current_cases == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)
cases_extra24 <- la_cases %>% 
  filter(recovered_cases == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)

cases_extra36 <- la_cases %>% 
  filter(new_recovered_cases == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)

cases_extra27 <- la_cases %>% 
  filter(newnew_recovered_cases == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)


cases_seed17 <- la_cases %>% 
  filter(new_current_cases == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)

r1_seed17 <- la_cases %>% 
  filter(newnewnew_recovered_cases == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)

r1_seed19 <- la_cases %>% 
  filter(recovered_cases19 == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)

r1_seed21 <- la_cases %>% 
  filter(recovered_cases21 == TRUE) %>%
  ungroup() %>%
  summarise(total_cases = sum(cases)) %>%
  pull(total_cases)
# detecting 1 in 5 seems plausible, though it is probably low balling 
E = sum(cases_lasttwo) * 5 * 1/3 * 0.48
I = sum(cases_lasttwo) * 5 * 2/3 * 0.48
R1 = sum(cases_extra24) * 5 * 0.48
R2 = popsize - S - E - I - R1

N = S + E + I + R1 + R2

new_R1 = cases_extra36 * 5 * 0.48
newnew_R1 = cases_extra27 * 5 * 0.48


# seed17 ------------------------------------------------------------------
E = sum(cases_seed17) * 5 * 2/9 * 0.48
I = sum(cases_seed17) * 5 * 7/9 * 0.48
R1 = sum(r1_seed17) * 5 * 0.48


# seed19 ------------------------------------------------------------------

E = sum(cases_seed17) * 5 * 2/9 * 0.48
I = sum(cases_seed17) * 5 * 7/9 * 0.48
R1 = sum(r1_seed19) * 5 * 0.48



# seed 20 -----------------------------------------------------------------
# detecting 1 in 5 seems plausible, though it is probably low balling 
E = sum(cases_lasttwo) * 5 * 4/11 * 0.48
I = sum(cases_lasttwo) * 5 * 7/11 * 0.48
R1 = sum(cases_extra24) * 5 * 0.48
R2 = popsize - S - E - I - R1

shedding_sum = 0.4 * I + 0.6 * R1

rho = mean_copies/shedding_sum


# seed21 ------------------------------------------------------------------
E = sum(cases_lasttwo) * 5 * 1/3 * 0.48
I = sum(cases_lasttwo) * 5 * 2/3 * 0.48
R1 = sum(r1_seed21) * 5 * 0.48


# now calculate rho given these values ------------------------------------
# current LA data
LA_ww_data <- read_csv(here::here("data", "LA_daily_data_feb2022.csv"))

mean_copies <- LA_ww_data %>% 
               filter(date == "2021-07-04") %>% 
               summarise(mean_copies = (gene_copy_1 + gene_copy_2 + gene_copy_3)/3000) %>% pull()

shedding_sum = 0.2 * I + 0.8 * R1

rho = mean_copies/shedding_sum


# lets redo assuming more like 1 in 50 rates of detection -----------------
E = sum(cases_lasttwo) * 200 * 4/11 * 0.48
I = sum(cases_lasttwo) * 200 * 7/11 * 0.48
R1 = sum(cases_extra24) * 200 * 0.48
R2 = popsize - S - E - I - R1


# now calculate rho given these values ------------------------------------
# current LA data
LA_ww_data <- read_csv(here::here("data", "LA_daily_data_feb2022.csv"))

mean_copies <- LA_ww_data %>% 
  filter(date == "2021-07-04") %>% 
  summarise(mean_copies = (gene_copy_1 + gene_copy_2 + gene_copy_3)/3000) %>% pull()

shedding_sum = 0.2 * I + 0.8 * R1

rho = mean_copies/shedding_sum



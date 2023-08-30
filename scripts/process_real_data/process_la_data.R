# this file is for processing la data for use in julia ode model
library(tidyverse)
library(lubridate)

# here the counts are in copies/L
raw_dat <- read_csv(here::here("data", "ladata_2021_2022.csv"))

raw_dat$date <- as.Date(raw_dat$date, format = "%m/%d/%Y")
dat <- raw_dat %>%
       dplyr::select(1:4) %>%
       filter(date >= "2021-07-04" & date <= "2021-10-30") %>%
       mutate(
              year = year(date),
              epi_week = epiweek(date),
              year_day = yday(date))


# filter out the days with no values
dat[is.na(dat)] <- -1
LA_dat <- dat %>% 
          mutate(sum_genes = gene_copy_1 + gene_copy_2 + gene_copy_3) %>% 
          filter(sum_genes >= 0) %>%
          mutate(new_time = year_day - min(year_day) + 7,
                 log_gene_copies1 = log(gene_copy_1/1000),
                 log_gene_copies2 = log(gene_copy_2/1000),
                 log_gene_copies3 = log(gene_copy_3/1000)) 

LA_dat[is.na(LA_dat)] <- -1
          
# filter out 8/23 which is a dumb day
LA_dat <- LA_dat %>% filter(date != "2021-08-23")

write_csv(LA_dat, here::here("data", "LA_daily_data.csv"))


# make a longer version of LA data all the way to the end for alpha model
raw_dat$date <- as.Date(raw_dat$date, format = "%m/%d/%Y")
dat <- raw_dat %>%
  dplyr::select(1:4) %>%
  filter(date >= "2021-07-04") %>%
  mutate(
    year = year(date),
    epi_week = epiweek(date),
    year_day = yday(date))

dat$year_day[dat$year == 2022] <- dat$year_day[dat$year == 2022] + 365
# filter out the days with no values
dat[is.na(dat)] <- -1
LA_dat <- dat %>% 
  mutate(sum_genes = gene_copy_1 + gene_copy_2 + gene_copy_3) %>% 
  filter(sum_genes >= 0) %>%
  mutate(new_time = year_day - min(year_day) + 1,
         log_gene_copies1 = log(gene_copy_1/1000),
         log_gene_copies2 = log(gene_copy_2/1000),
         log_gene_copies3 = log(gene_copy_3/1000)) 

LA_dat[is.na(LA_dat)] <- -1

# filter out 8/23 which is a weird day
# and filter our 7/7 which also seems weird
# filter out 7/11 which has a note saying its below LOD
LA_dat <- LA_dat %>% filter(date != "2021-08-23") %>% filter(date != "2021-07-07") %>% filter(date != "2021-07-011")

write_csv(LA_dat,here::here("data", "LA_daily_data_feb2022.csv"))

# creating week data (defunct, no longer doing this) ----------------------


week_data <- dat %>%
             filter(gene_copy_1 != -1) %>%
             group_by(year, epi_week) %>%
             sample_n(1) %>%
             select(-order,  -date_diff)

week_41mean1 <- week_data %>% filter(epi_week == 40 | epi_week == 42) %>% pull(gene_copy_1) %>% mean()
week_41mean2 <- week_data %>% filter(epi_week == 40 | epi_week == 42) %>% pull(gene_copy_2) %>% mean()
week_41mean3 <- week_data %>% filter(epi_week == 40 | epi_week == 42) %>% pull(gene_copy_3) %>% mean()

week_41 <- data.frame(date = as.Date("2021-10-10"), 
                      gene_copy_1 = week_41mean1,
                      gene_copy_2 = week_41mean2,
                      gene_copy_3 = week_41mean3,
                      year = 2021, 
                      epi_week = 41)
week_data <- rbind(week_data, week_41) %>%
             arrange(year, epi_week)


week_data <- week_data %>%
             ungroup() %>%
             mutate(log_gene_copies1 = log(gene_copy_1),
                    log_gene_copies2 = log(gene_copy_2), 
                    log_gene_copies3 = log(gene_copy_3),
                    new_time = row_number()) 

avg <- week_data %>% mutate(mean_copies = (exp(log_gene_copies1) + 
                                          exp(log_gene_copies2) + 
                                          exp(log_gene_copies3))/3)

write_csv(week_data, "LA_week_data.csv")


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
  mutate(case_date = as.Date("2021-07-04") - days(14),
         recover_date = case_date - days(24)) %>%
  filter(date >= recover_date & date <= "2021-10-30") %>%
  filter(area == "Los Angeles") %>%
  mutate(year = year(date),
         epi_week = epiweek(date)) 
la_cases <-  la_cases %>%
                  mutate(case_date = as.Date("2021-07-04") - days(14),
                         recover_date = case_date - days(24),
                         current_cases = date >= case_date & date < "2021-07-04",
                         recovered_cases = date >= recover_date & date <= case_date) 

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

  
  
E = sum(cases_lasttwo) * 4 * 1/3 * 0.48
I = sum(cases_lasttwo) * 4 * 2/3 * 0.48
R1 = sum(cases_extra24) * 5 * 0.48
R2 = popsize - S - E - I - R1

N = S + E + I + R1 + R2

# let's futz with the S compartment so that S/N give us closer to 2
S = 2 * N / 5
E = sum(cases_lasttwo) * 15 * 1/3 * 0.48
I = sum(cases_lasttwo) * 15 * 2/3 * 0.48
R1 = sum(cases_extra24) * 15 * 0.48
R2 = popsize - S - E - I - R1

# rough calculations for rho in LA on Aug 23rd
la_data <- read_csv("data/covid19cases_test.csv") %>% 
           filter(area == "Los Angeles") %>%
           mutate(case_date = as.Date("2021-08-23") - days(14),
                  recover_date = case_date - days(24)) %>%
           filter(date >= recover_date & date <= "2021-08-23") %>%
           mutate(current_cases = date >= case_date,
                  recovered_cases = date >= recover_date & date <= case_date)

current_cases <- la_data %>% 
                 filter(current_cases == TRUE) %>%
                 summarise(sum_cases = sum(cases) * 4 * (2/3) * 0.48) %>%
                 pull(sum_cases)
                

current_recovered <- la_data %>% 
  filter(recovered_cases == TRUE) %>%
  summarise(sum_recov = sum(cases) * 4 * 0.48) %>%
  pull(sum_recov)

total_emit <- current_cases + current_recovered

# minimal average value is 8/23 
# estimated rho is then minimal value divided by estimated I + R1
estim_rho <- 51530.79/total_emit

# estimating tau across time
raw_dat <- read_csv("ladata_2021_2022.csv")

raw_dat$date <- as.Date(raw_dat$date, format = "%m/%d/%Y")
dat <- raw_dat %>%
  dplyr::select(1:4) %>%
  filter(date >= "2021-07-04") %>%
  mutate(date_diff = difftime(date, lag(date)),
         year = year(date),
         epi_week = epiweek(date)) %>%
  group_by(year, epi_week) %>%
  mutate(order = row_number())

dat$gene_copy_1[dat$date == "2021-09-27"] <- (dat$gene_copy_2[dat$date == "2021-09-27"] + 
                                                dat$gene_copy_3[dat$date == "2021-09-27"])/2
dat$gene_copy_1[dat$date == "2021-10-03"] <- (dat$gene_copy_2[dat$date == "2021-10-03"] + 
                                                dat$gene_copy_3[dat$date == "2021-10-03"])/2
dat$gene_copy_2[dat$date == "2021-10-04"] <- (dat$gene_copy_1[dat$date == "2021-10-04"] + 
                                                dat$gene_copy_3[dat$date == "2021-10-04"])/2
month_year_data <- dat %>%
  filter(gene_copy_1 != -1) %>%
  mutate(month_year = floor_date(date, "month"),
         log_gene_copies1 = log(gene_copy_1),
                log_gene_copies2 = log(gene_copy_2), 
                log_gene_copies3 = log(gene_copy_3),
                new_time = row_number()) %>%
  ungroup() %>%
  dplyr::select(year, month_year, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
  pivot_longer(cols = c(-year, -month_year), names_to = "name", values_to = "log_copies") %>%
  group_by(year, month_year) %>%
  summarise(mom_tau = sd(log_copies),
            num_data = n())

plot_tau <- month_year_data %>% 
            ggplot(aes(x = month_year, y = mom_tau)) +
            geom_line() + 
            geom_point() + 
            theme_bw() + 
            ylab("MoM Tau") + 
            ggtitle("Tau estimates by Month")

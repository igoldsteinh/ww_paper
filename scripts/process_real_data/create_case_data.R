# file for creating LA case data
library(tidyverse)
library(ckanr)
library(lubridate)
library(fs)



quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


ckanr_setup(url="https://data.ca.gov")
ckan <- quiet(ckanr::src_ckan("https://data.ca.gov"))

# get resources
resources <- rbind(resource_search("name:covid-19", as = "table")$results,
                   resource_search("name:hospitals by county", as = "table")$results)


cases_deaths_url <- resources %>% filter(name == "Statewide COVID-19 Cases Deaths Tests") %>% pull(url)

cases <-
  read_csv(cases_deaths_url) %>%
  mutate(date = lubridate::ymd(date),
         deaths = as.integer(deaths),
         reported_cases = as.integer(reported_cases),
         cases = as.integer(cases),
         positive_tests = as.integer(positive_tests),
         total_tests = as.integer(total_tests)) %>%
  select(date,
         cases = cases,
         tests = total_tests,
         deaths,
         county = area) %>%
  arrange(date, county) %>%  
  mutate(year = year(date),
         epi_week = epiweek(date),
         year_day = yday(date)) %>% 
  filter(county == "Los Angeles",
         date >= "2021-06-01")

weekly_data <- cases %>% 
               group_by(year, epi_week) %>%
               summarise(total_cases = sum(cases),
                         total_tests = sum(tests),
                         min_date = min(date),
                         max_date = max(date))

write_csv(weekly_data, here::here("data", "LA_weekly_case_data.csv"))
            

# create data to use with EIR ---------------------------------------------
weekly_data <- read_csv(here::here("data", "LA_weekly_case_data.csv"))

EIR_data <- weekly_data %>% filter(max_date >= "2021-07-10" & max_date <= "2022-02-26") %>%
            mutate(new_week = row_number())

write_csv(EIR_data, here::here("data", "LA_EIR_data.csv"))


# create daily case data to estimate case detection -----------------------
daily_case_data <- cases %>% 
                   filter(date >= "2021-07-04" & date <= "2022-02-23") %>%
                   mutate(new_time = row_number())

write_csv(daily_case_data, here::here("data", "LA_daily_case_data.csv"))


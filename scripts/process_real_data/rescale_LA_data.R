library(tidyverse)


la_data <- read_csv(here::here("data", "LA_daily_data.csv"))

model_data <- la_data %>% 
  dplyr::select(date, new_time, gene_copy_1, gene_copy_2, gene_copy_3) %>%
  pivot_longer(cols = all_of(c("gene_copy_1", 
                               "gene_copy_2",
                               "gene_copy_3")), 
               names_to = "copy") %>%
  filter(value != -1) %>%
  mutate(value_scaled = value/1000,
         log_values = log(value_scaled),
         days_since = date - min(la_data$date)) %>%
  filter(date != "2021-07-07")


write_csv(model_data, here::here("data", "formatted_LA_daily_data.csv"))

la_data <- la_data %>% 
           mutate(scaled_1 = gene_copy_1/1000,
                  scaled_2 = gene_copy_2/1000,
                  scaled_3 = gene_copy_3/1000,
                  log_scaled_1 = log(scaled_1),
                  log_scaled_2 = log(scaled_2),
                  log_scaled_3 = log(scaled_3))

la_data[is.na(la_data)] <- -1

la_data <- la_data %>% 
           filter(date != "201-07-07")

write_csv(la_data, here::here("data", "LA_daily_data_scaled.csv"))
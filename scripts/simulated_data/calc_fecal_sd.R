# calculate mean standard deviation of concentrations
# this number is used when simulation data
library(tidyverse)
library(jsonlite)
library(httr)
#studies are from https://github.com/tillahoffmann/shedding 
# paper is Hoffman and Alsing Faecal shedding models for SARS-CoV-2 RNA among hos- pitalised patients and implications for wastewater-based epidemiology
han = read_json(here::here("raw_data", "Han2020.json"))
lui = read_json(here::here("raw_data", "Lui2020.json"))
woelfel = read_json(here::here("raw_data", "Woelfel2020.json"))

woelfel_frame = woelfel %>% pluck("loads") %>% enframe() %>% unnest_wider(value) %>% mutate(study = "woelfel")
lui_frame = lui %>% pluck("loads") %>% enframe() %>% unnest_wider(value) %>% mutate(study = "lui")
han_frame = han %>% pluck("loads") %>% enframe() %>% unnest_wider(value) %>% mutate(study = "han")

full_data = bind_rows(woelfel_frame, lui_frame, han_frame)

dot_plot = full_data %>% 
           ggplot(aes(x = day, y = value)) + 
           geom_point() + 
           theme_bw()

sd_by_day = full_data %>% 
            group_by(day) %>%
            drop_na() %>%
            summarise(mean_val = mean(value),
                      sd_val = sd(value),
                      num_points = n())
mean(sd_by_day$sd_val, na.rm = TRUE)
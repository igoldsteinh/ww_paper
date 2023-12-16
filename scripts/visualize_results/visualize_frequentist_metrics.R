# visualize frequentist performance of various models on various scenarios
library(tidyverse)
source(here::here("src", "wastewater_functions.R"))

# scenario 1 --------------------------------------------------------------
# read in data, lets just look at 80% CI
seirr_rt_scenario1 <- read_csv(here::here("results", "seirr_student", "seirr_scenario1_allseeds_rt_quantiles.csv")) %>% 
                      mutate(model = "SEIRR") %>% 
                      filter(.width == 0.8)

eirr_rt_scenario1 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario1_allseeds_rt_quantiles.csv")) %>% 
                     mutate(model = "EIRR-ww")%>% 
  filter(.width == 0.8)

seir_rt_scenario1 <- read_csv(here::here("results", "seir_cases", "seir_cases_scenario1_allseeds_rt_quantiles.csv")) %>% 
                     mutate(model = "SEIR-cases")%>% 
  filter(.width == 0.8)

eir_rt_scenario1 <- read_csv(here::here("results", "eir_cases", "eir_cases_scenario1_allseeds_rt_quantiles.csv")) %>% 
                    mutate("EIR-cases")%>% 
  filter(.width == 0.8)

checking_eirr <- eirr_rt_scenario1 %>% 
                 group_by(seed) %>% 
                 summarise(num = n()) %>% 
                 pull(num) %>% 
                 unique()

checking_seir <- seir_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()
checking_eir <- eir_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking_seirr <- seirr_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

# create metrics

eirr_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "EIRR-ww")
  eirr_rt_metrics = bind_rows(eirr_rt_metrics, metrics)
}

seirr_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = seirr_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "SEIRR-ww")
  seirr_rt_metrics = bind_rows(seirr_rt_metrics, metrics)
}

eir_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = eir_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "EIR-cases")
  eir_rt_metrics = bind_rows(eir_rt_metrics, metrics)
}

seir_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = seir_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "SEIR-cases")
  seir_rt_metrics = bind_rows(seir_rt_metrics, metrics)
}

all_metrics <- bind_rows(eirr_rt_metrics, seirr_rt_metrics, eir_rt_metrics, seir_rt_metrics) %>% 
               pivot_longer(cols = -c(seed, model), names_to = "metric")


level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")

label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")
all_metrics$metric <- factor(all_metrics$metric, levels=level_list, labels=label_list)

model_level_list <-c ("EIRR-ww", "SEIRR-ww", "EIR-cases", "SEIR-cases")
all_metrics$model <- factor(all_metrics$model, levels = model_level_list)

all_metrics$useful_value <- 0
all_metrics$useful_value[all_metrics$metric == "Envelope"] <- 0.8
all_metrics$useful_value[all_metrics$metric == "MCIW"] <- NA
all_metrics$useful_value[all_metrics$metric == "MASV"] <- 0.02361307
all_metric_plot <- all_metrics %>% 
                   filter(metric != "True MASV") %>% 
                   ggplot(aes(x = model, y = value)) + 
                   geom_boxplot() + 
                   geom_hline(aes(yintercept = useful_value)) + 
                   theme_bw() + 
                   theme(text = element_text(size = 18)) +
                   facet_wrap(vars(metric),
                              scales = "free") +
                   ggtitle("Frequentist Metrics Across Models") + 
                   ylab("") + 
                   xlab("Model")

all_metric_dev_plot <- all_metrics %>% 
            filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x=element_blank(),
        ) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")

all_metric_env_plot <- all_metrics %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")


all_metric_mciw_plot <- all_metrics %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")

all_metric_masv_plot <- all_metrics %>% 
  filter(metric == "MASV") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")

new_all_metric_plot <- all_metric_dev_plot + all_metric_env_plot + all_metric_mciw_plot + all_metric_masv_plot + plot_annotation(
  title = 'Frequentist Metrics Across Models'
) &
theme(text = element_text(size = 18))

ggsave(here::here("figures", "scenario1_frequentist_metrics.pdf"), new_all_metric_plot, width = 11, height =11)

# calculating MCIW --------------------------------------------------------

mciw_calc <- all_metrics %>% filter(metric == "MCIW") %>%
             group_by(model) %>%
             summarise(median_MCIW = median(value))

masv_calc <- all_metrics %>% filter(metric == "MASV") %>%
  group_by(model) %>%
  summarise(median_MCIW = median(value))
# presentation version ----------------------------------------------------

presentation_all_metric_plot <- all_metrics %>% 
  filter(metric != "True MASV" & metric != "MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ylab("") +
  xlab("Model") +
  facet_wrap(vars(metric),
             scales = "free") +
  ggtitle("Frequentist Metrics")

ggsave(here::here("figures", "scenario1_presentation_frequentist_metrics.png"), 
       presentation_all_metric_plot, 
       width = 10, 
       height =4)

# Scenario 1 vs 3 vs 4 vs 5 vs 6 (data type plot) -----------------------------------
# read in data, sanity check we have the same number of comparison points
# ten replicates
eirr_rt_scenario3 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario3_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Ten)")%>% 
  filter(.width == 0.8)
# three mean
eirr_rt_scenario4 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario4_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Mean)")%>% 
  filter(.width == 0.8)
# 10 mean
eirr_rt_scenario5 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario5_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Mean)")%>% 
  filter(.width == 0.8)
# 1 replicate
eirr_rt_scenario6 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario6_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Mean)")%>% 
  filter(.width == 0.8)

checking3 <- eirr_rt_scenario3 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking4 <- eirr_rt_scenario3 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking5 <- eirr_rt_scenario5 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking6 <- eirr_rt_scenario6 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

# create metrics

eirr_rt_metrics3 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario3 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "10-rep")
  eirr_rt_metrics3 = bind_rows(eirr_rt_metrics3, metrics)
}

eirr_rt_metrics4 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario4 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "3-mean")
  eirr_rt_metrics4 = bind_rows(eirr_rt_metrics4, metrics)
}

eirr_rt_metrics5 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario5 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "10-mean")
  eirr_rt_metrics5 = bind_rows(eirr_rt_metrics5, metrics)
}

eirr_rt_metrics6 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario6 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "1-rep")
  eirr_rt_metrics6 = bind_rows(eirr_rt_metrics6, metrics)
}

# make graph

all_metrics2 <- bind_rows(eirr_rt_metrics, 
                          eirr_rt_metrics3, 
                          eirr_rt_metrics4,
                          eirr_rt_metrics5,
                          eirr_rt_metrics6) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")

level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")

label_list <- c("EIRR-ww Deviation", "EIRR-ww Envelope", "EIRR-ww MCIW", "EIRR-ww MASV", "True MASV")

all_metrics2$metric <- factor(all_metrics2$metric, levels=level_list, labels=label_list)

all_metrics2$model[all_metrics2$model == "EIRR-ww"] <- "3-rep"
model_level_list <-c ("3-rep", "1-rep", "10-rep", "3-mean", "10-mean")
all_metrics2$model <- factor(all_metrics2$model, levels = model_level_list)

all_metrics2$useful_value <- 0
all_metrics2$useful_value[all_metrics2$metric == "EIRR-ww Envelope"] <- 0.8
all_metrics2$useful_value[all_metrics2$metric == "EIRR-ww MCIW"] <- NA
all_metrics2$useful_value[all_metrics2$metric == "EIRR-ww MASV"] <- 0.02361307

all_metric_plot2 <- all_metrics2 %>% 
  filter(metric != "True MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  facet_wrap(vars(metric), scales = "free") +
  ggtitle("Frequentist Metrics Across Data Sources") + 
  ylab("") +
  xlab("Fitted Data")


all_metric2_dev_plot <- all_metrics2 %>% 
  filter(metric == "EIRR-ww Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x=element_blank(),
  ) +
  ggtitle("EIRR-ww Deviation") + 
  ylab("Deviation") + 
  xlab("")

all_metric2_env_plot <- all_metrics2 %>% 
  filter(metric == "EIRR-ww Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  ggtitle("EIRR-ww Envelope") + 
  ylab("Envelope") + 
  xlab("")


all_metric2_mciw_plot <- all_metrics2 %>% 
  filter(metric == "EIRR-ww MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("EIRR-ww MCIW") + 
  ylab("MCIW") + 
  xlab("Fitted Data")

all_metric2_masv_plot <- all_metrics2 %>% 
  filter(metric == "EIRR-ww MASV") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("EIRR-ww MASV") + 
  ylab("MASV") + 
  xlab("Fitted Data")

new_all_metric2_plot <- all_metric2_dev_plot + all_metric2_env_plot + all_metric2_mciw_plot + all_metric2_masv_plot + plot_annotation(
  title = 'Frequentist Metrics Across Data Sources'
) &
  theme(text = element_text(size = 18))

ggsave(here::here("figures", "otherscenarios_frequentist_metrics.pdf"), new_all_metric2_plot, width = 11, height =11)

# presentation version ----------------------------------------------------
presentation_all_metric2_plot <- all_metrics2 %>% 
  filter(metric != "True MASV" & metric != "MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 20)) +
  ylab("") +
  xlab("Model") +
  facet_wrap(vars(metric), scales = "free") +
  ggtitle("Frequentist Metrics (Changing WW Data)")

ggsave(here::here("figures", "otherscenarios_presentation_frequentist_metrics.png"), 
       presentation_all_metric2_plot, 
       width = 12, 
       height =6)

# scenario 1 vs 7-10 (sensitivity analysis) -------------------------------
# lambda centered at 0.8
eirr_rt_scenario8 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario8_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Mean)")%>% 
  filter(.width == 0.8)
# E,I are at 75%
eirr_rt_scenario9 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario9_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Mean)")%>% 
  filter(.width == 0.8)

# E,I are at 133%
eirr_rt_scenario10 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario10_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR (Mean)")%>% 
  filter(.width == 0.8)

# sanity check

checking8 <- eirr_rt_scenario8 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking9 <- eirr_rt_scenario9 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking10 <- eirr_rt_scenario10 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

# create metrics

eirr_rt_metrics8 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario8 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "Low Prop")
  eirr_rt_metrics8 = bind_rows(eirr_rt_metrics8, metrics)
}

eirr_rt_metrics9 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario9 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "Low Init")
  eirr_rt_metrics9 = bind_rows(eirr_rt_metrics9, metrics)
}

eirr_rt_metrics10 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario10 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "High Init")
  eirr_rt_metrics10 = bind_rows(eirr_rt_metrics10, metrics)
}

# make graph

all_metrics3 <- bind_rows(eirr_rt_metrics, 
                          eirr_rt_metrics8,
                          eirr_rt_metrics9,
                          eirr_rt_metrics10) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")


level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")

label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")


all_metrics3$metric <- factor(all_metrics3$metric, levels=level_list, labels=label_list)

model_level_list <-c ("EIRR-ww", "Low Prop", "High Prop", "Low Init", "High Init")
all_metrics3$model <- factor(all_metrics3$model, levels = model_level_list)

all_metrics3$useful_value <- 0
all_metrics3$useful_value[all_metrics3$metric == "Envelope"] <- 0.8
all_metrics3$useful_value[all_metrics3$metric == "MCIW"] <- NA
all_metrics3$useful_value[all_metrics3$metric == "MASV"] <- 0.02361307

all_metric_plot3 <- all_metrics3 %>% 
  filter(metric != "True MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  facet_wrap(vars(metric),
             scales = "free") +
  ggtitle("Frequentist Metrics (Sensitivity Analysis)") + 
  ylab("") + 
  xlab("Sensitivity Setting")

all_metric3_dev_plot <- all_metrics3 %>% 
  filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x=element_blank(),
  ) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")

all_metric3_env_plot <- all_metrics3 %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")

all_metric3_mciw_plot <- all_metrics3 %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Sensitivity Setting")

all_metric3_masv_plot <- all_metrics3 %>% 
  filter(metric == "MASV") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Sensitivity Setting")

new_all_metric3_plot <- all_metric3_dev_plot + all_metric3_env_plot + all_metric3_mciw_plot + all_metric3_masv_plot + plot_annotation(
  title = 'Frequentist Metrics (Sensitivity Analysis)'
) &
  theme(text = element_text(size = 18))

ggsave(here::here("figures", "sensitivity_frequentist_metrics.pdf"), new_all_metric3_plot, width = 11, height =11)

# presentation version ----------------------------------------------------
presentation_all_metric3_plot <- all_metrics3 %>% 
  filter(metric != "True MASV" & metric != "MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 90)) +
  xlab("Setting") +
  ylab("") +
  facet_wrap(vars(metric),
             scales = "free") +
  ggtitle("Frequentist Metrics (Sensitivity Analysis)")

ggsave(here::here("figures", "sensitivity_presentation_frequentist_metrics.png"), 
       presentation_all_metric3_plot, 
       width = 10, 
       height =6)

# comparing against huisman -----------------------------------------------

huisman_rt_scenario1 <- read_csv(here::here("results", "huisman", "huisman_scenario1_allseeds_rt_quantiles.csv"))

huisman_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = huisman_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = median_R_mean, upper = median_R_highHPD, lower = median_R_lowHPD) %>% 
            mutate(seed = i,
                   model = "Huisman")
  huisman_rt_metrics = bind_rows(huisman_rt_metrics, metrics)
}

eirr_rt_scenario1_0.95 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR-ww")%>% 
  filter(.width == 0.95)


abbreviated_eirr = NULL

for (i in 1:100) {
  sub_frame = eirr_rt_scenario1_0.95 %>% filter(seed == i) %>% filter(time <= max(huisman_rt_scenario1$time) & time >= min(huisman_rt_scenario1$time))
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = i,
           model = "EIRR-ww")
  abbreviated_eirr = bind_rows(abbreviated_eirr, metrics)
}

# make graph

all_metrics4 <- bind_rows(huisman_rt_metrics, 
                          abbreviated_eirr) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")


level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")

label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")


all_metrics4$metric <- factor(all_metrics4$metric, levels=level_list, labels=label_list)

model_level_list <-c ("EIRR-ww", "Huisman")
all_metrics4$model <- factor(all_metrics4$model, levels = model_level_list)

all_metrics4$useful_value <- 0
all_metrics4$useful_value[all_metrics4$metric == "Envelope"] <- 0.95
all_metrics4$useful_value[all_metrics4$metric == "MCIW"] <- NA
all_metrics4$useful_value[all_metrics4$metric == "MASV"] <- 0.02361307

all_metric_plot4 <- all_metrics4 %>% 
  filter(metric != "True MASV") %>% 
  filter(value <= 600) %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  facet_wrap(vars(metric), scales = "free") +
  ggtitle("Frequentist Metrics (State of the Art)") + 
  ylab("") +
  xlab("Model")

all_metric4_dev_plot <- all_metrics4 %>% 
  filter(metric == "Deviation") %>%
  filter(value <= 600) %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x=element_blank(),
  ) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")

all_metric4_env_plot <- all_metrics4 %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")

all_metric4_mciw_plot <- all_metrics4 %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")

all_metric4_masv_plot <- all_metrics4 %>% 
  filter(metric == "MASV") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")

new_all_metric4_plot <- all_metric4_dev_plot + all_metric4_env_plot + all_metric4_mciw_plot + all_metric4_masv_plot + plot_annotation(
  title = 'Frequentist Metrics (State of the Art)'
) &
  theme(text = element_text(size = 18))

ggsave(here::here("figures", "stateoftheart_frequentist_metrics.pdf"), new_all_metric4_plot, width = 10, height =10)

# presentation version ----------------------------------------------------
presentation_all_metric4_plot <- all_metrics4 %>% 
  filter(metric != "True MASV" & metric != "MASV") %>% 
  filter(value <= 2) %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  xlab("Model") +
  ylab("") +
  facet_wrap(vars(metric), scales = "free") +
  ggtitle("Frequentist Metrics (State of the Art)")

ggsave(here::here("figures", "stateoftheart_presentation_frequentist_metrics.png"), 
       presentation_all_metric4_plot, 
       width = 10, 
       height =4)

# repeating model comparison figure with 95% CI ---------------------------
seirr_rt_scenario1 <- read_csv(here::here("results", "seirr_student", "seirr_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "SEIRR-ww") %>% 
  filter(.width == 0.95)

eirr_rt_scenario1 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR-ww")%>% 
  filter(.width == 0.95)

seir_rt_scenario1 <- read_csv(here::here("results", "seir_cases", "seir_cases_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "SEIR-cases")%>% 
  filter(.width == 0.95)

eir_rt_scenario1 <- read_csv(here::here("results", "eir_cases", "eir_cases_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate("EIR-cases")%>% 
  filter(.width == 0.95)

checking_eirr <- eirr_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking_seir <- seir_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()
checking_eir <- eir_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

checking_seirr <- seirr_rt_scenario1 %>% 
  group_by(seed) %>% 
  summarise(num = n()) %>% 
  pull(num) %>% 
  unique()

# create metrics

eirr_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "EIRR-ww")
  eirr_rt_metrics = bind_rows(eirr_rt_metrics, metrics)
}

seirr_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = seirr_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "SEIRR-ww")
  seirr_rt_metrics = bind_rows(seirr_rt_metrics, metrics)
}

eir_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = eir_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "EIR-cases")
  eir_rt_metrics = bind_rows(eir_rt_metrics, metrics)
}

seir_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = seir_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "SEIR-cases")
  seir_rt_metrics = bind_rows(seir_rt_metrics, metrics)
}

all_metrics <- bind_rows(eirr_rt_metrics, seirr_rt_metrics, eir_rt_metrics, seir_rt_metrics) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")


level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")

label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")
all_metrics$metric <- factor(all_metrics$metric, levels=level_list, labels=label_list)


model_level_list <-c ("EIRR-ww", "SEIRR-ww", "EIR-cases", "SEIR-cases")
all_metrics$model <- factor(all_metrics$model, levels = model_level_list)

all_metrics$useful_value <- 0
all_metrics$useful_value[all_metrics$metric == "Envelope"] <- 0.95
all_metrics$useful_value[all_metrics$metric == "MCIW"] <- NA
all_metrics$useful_value[all_metrics$metric == "MASV"] <- 0.02361307
all_metric_plot_95 <- all_metrics %>% 
  filter(metric != "True MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  facet_wrap(vars(metric),
             scales = "free") +
  ggtitle("Frequentist Metrics Across Models (95% CI)") + 
  ylab("") + 
  xlab("Model")

all_metric_dev_plot_95 <- all_metrics %>% 
  filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x=element_blank(),
  ) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")

all_metric_env_plot_95 <- all_metrics %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")


all_metric_mciw_plot_95 <- all_metrics %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")

all_metric_masv_plot_95 <- all_metrics %>% 
  filter(metric == "MASV") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")

new_all_metric_plot_95 <- all_metric_dev_plot_95 + all_metric_env_plot_95 + all_metric_mciw_plot_95 + all_metric_masv_plot_95 + plot_annotation(
  title = 'Frequentist Metrics Across Models (95% CI)'
) &
  theme(text = element_text(size = 18))

ggsave(here::here("figures", "scenario1_frequentist_metrics_95CI.pdf"), new_all_metric_plot_95, width = 11, height =11)

# scenario 101 vs 1 --------------------------------------------------------------
# read in data, lets just look at 80% CI

eirr_rt_scenario1 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario1_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR-ww")%>% 
  filter(.width == 0.8)

eirr_rt_scenario101 <- read_csv(here::here("results", "eirrc_closed", "eirrc_scenario101_allseeds_rt_quantiles.csv")) %>% 
  mutate(model = "EIRR-ww (Stoch Rt)")%>% 
  filter(.width == 0.8)


# create metrics

eirr_rt_metrics = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario1 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "EIRR-ww")
  eirr_rt_metrics = bind_rows(eirr_rt_metrics, metrics)
}

eirr_rt_metrics_101 = NULL
for (i in 1:100) {
  sub_frame = eirr_rt_scenario101 %>% filter(seed == i) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% mutate(seed = i,
                                                                                            model = "EIRR-ww (Stoch Rt)")
  eirr_rt_metrics_101 = bind_rows(eirr_rt_metrics_101, metrics)
}


all_metrics <- bind_rows(eirr_rt_metrics, eirr_rt_metrics_101) %>% 
  pivot_longer(cols = -c(seed, model), names_to = "metric")


level_list <- c("mean_dev", "mean_env", "MCIW", "MASV", "true_MASV")

label_list <- c("Deviation", "Envelope", "MCIW", "MASV", "True MASV")
all_metrics$metric <- factor(all_metrics$metric, levels=level_list, labels=label_list)

model_level_list <-c ("EIRR-ww", "EIRR-ww (Stoch Rt)")
all_metrics$model <- factor(all_metrics$model, levels = model_level_list)

all_metrics$useful_value <- 0
all_metrics$useful_value[all_metrics$metric == "Envelope"] <- 0.8
all_metrics$useful_value[all_metrics$metric == "MCIW"] <- NA
all_metrics$useful_value[all_metrics$metric == "MASV"] <- 0.02361307
all_metric_plot <- all_metrics %>% 
  filter(metric != "True MASV") %>% 
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  facet_wrap(vars(metric),
             scales = "free") +
  ggtitle("Frequentist Metrics Across Models") + 
  ylab("") + 
  xlab("Model")

all_metric_dev_plot <- all_metrics %>% 
  filter(metric == "Deviation") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x=element_blank(),
  ) +
  ggtitle("Deviation") + 
  ylab("Deviation") + 
  xlab("")

all_metric_env_plot <- all_metrics %>% 
  filter(metric == "Envelope") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18),
        axis.text.x = element_blank()) +
  ggtitle("Envelope") + 
  ylab("Envelope") + 
  xlab("")


all_metric_mciw_plot <- all_metrics %>% 
  filter(metric == "MCIW") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MCIW") + 
  ylab("MCIW") + 
  xlab("Model")

all_metric_masv_plot <- all_metrics %>% 
  filter(metric == "MASV") %>%
  ggplot(aes(x = model, y = value)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = useful_value)) + 
  theme_bw() + 
  theme(text = element_text(size = 18)) +
  ggtitle("MASV") + 
  ylab("MASV") + 
  xlab("Model")

all_metric_plot_stochRt <- all_metric_dev_plot + all_metric_env_plot + all_metric_mciw_plot + all_metric_masv_plot + plot_annotation(
  title = 'Frequentist Metrics (Stochastic vs Fixed Rt)'
) &
  theme(text = element_text(size = 18))

ggsave(here::here("figures", "stochvsfixed_frequentist_metrics.pdf"), all_metric_plot_stochRt, width = 11, height =11)

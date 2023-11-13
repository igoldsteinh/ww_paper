# General utility functions
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)
library(EpiEstim)
library(zoo)
library(tidybayes)
set.seed(1234)

# make prior samples ------------------------------------------------------

make_prior_samples <- function(prior_gq) {
  prior_gq_samples <-
    prior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]", negate = T))

  return(prior_gq_samples)
  
}


# make quantiles for time varying parameters ------------------------------
make_timevarying_prior_quantiles <- function(prior_gq) {
  prior_gq_samples <-
    prior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]")) %>%
    mutate(time = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(.,tibble(time = 0:max(.$time)))
  
  return(prior_gq_samples)
  
}


# make prior quantiles ----------------------------------------------------
make_prior_quantiles <- function(prior_samples) {
  prior_quantiles <- prior_samples %>%
    group_by(name) %>%
    summarise(quantiles = quantile(value, c(0.025, 0.5, 0.975)),
              q = c(0.025, 0.5, 0.975)) %>%
    ungroup() %>%
    pivot_wider(id_cols = name, names_from = q, values_from = quantiles) 
  
  return(prior_quantiles)
}


# plot simulated data -----------------------------------------------------

plot_sim_data <- function(full_sim_data) {
  
  maxI <- max(full_sim_data$I) 
  maxRe <- max(full_sim_data$R1)
  maxgene <- max(c(full_sim_data$log_gene_copies1, 
                   full_sim_data$log_gene_copies2, 
                   full_sim_data$log_gene_copies3))
  
  maxcase <- max(full_sim_data$cases)
  maxtest <- max(full_sim_data$tests)
  
  sim_prev <- full_sim_data %>%
    ggplot(aes(x = time, y = I)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxI + 10)) +
    ggtitle("Prev I") +
    ylab("Prev I") +
    xlab("") +
    theme_bw() +
    theme(text = element_text(size = 14))
  
  sim_R1 <- full_sim_data %>%
    ggplot(aes(x = time, y = R1)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxRe)) +
    scale_x_continuous(expand = c(0,0)) +
    ggtitle("Prev R1") +
    ylab("Prev R1") +
    xlab("") +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  
  sim_rt <- full_sim_data %>%
    ggplot(aes(x = time, y = true_rt)) + 
    geom_line() +
    ggtitle("Rt") +
    ylab("Rt") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  sim_data_long <- full_sim_data %>%
    ungroup() %>%
    dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
    pivot_longer(cols = -time, values_to = "log_genes") 
  
  
  full_sim_data <- full_sim_data %>% 
      mutate(mean_log = (log_gene_copies1 + log_gene_copies2 + log_gene_copies3)/3)
  loggene_plot <- sim_data_long %>%
    ggplot(aes(x = time, y = log_genes, color = name)) + 
    geom_point() + 
    geom_line() + 
    geom_line(data = full_sim_data, aes(x = time, y = mean_log), color = "black", size = 1.2) +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
    scale_x_continuous(expand = c(0,0)) +
    
    ggtitle("Log Gene Counts") + 
    ylab("Log RNA") + 
    xlab("Time") +
    theme_bw()+
    theme(text = element_text(size = 14),
          legend.position = "none")
  
  sim_cases <- full_sim_data %>%
    ggplot(aes(x = time, y = cases)) +
    geom_point() +
    geom_line() +
    ggtitle("Observed Daily Cases") +
    ylab("Cases") +
    xlab("Time") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxcase)) +
    
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  sim_tests <- full_sim_data %>%
    ggplot(aes(x = time, y = tests)) +
    geom_point() +
    geom_line() +
    ggtitle("Daily Tests") +
    ylab("Tests") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxtest + 20)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  full_plot <- (sim_rt + sim_prev) /
    (sim_R1 + sim_tests) /
    (loggene_plot + sim_cases) +
    plot_layout(guides = 'collect')
  
}

plot_sim_gene <- function(full_sim_data,
                          mean_line = TRUE) {
  
  maxgene <- max(c(full_sim_data$log_gene_copies1, 
                   full_sim_data$log_gene_copies2, 
                   full_sim_data$log_gene_copies3))
  


  sim_data_long <- full_sim_data %>%
    ungroup() %>%
    dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
    pivot_longer(cols = -time, values_to = "log_genes") 
  full_sim_data <- full_sim_data %>% 
    mutate(mean_log = (log_gene_copies1 + log_gene_copies2 + log_gene_copies3)/3)
  
  if (mean_line == TRUE) {
    loggene_plot <- sim_data_long %>%
      ggplot(aes(x = time, y = log_genes, color = name)) + 
      geom_point() + 
      geom_line() + 
      geom_line(data = full_sim_data, aes(x = time, y = mean_log), color = "black", size = 1.2) +
      scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
      scale_x_continuous(expand = c(0,0)) +
      
      ggtitle("Log Gene Counts") + 
      ylab("Log RNA") + 
      xlab("Time") +
      theme_bw()+
      theme(text = element_text(size = 14),
            legend.position = "none")
    
  }
  
  if (mean_line == FALSE) {
    loggene_plot <- sim_data_long %>%
      ggplot(aes(x = time, y = log_genes, color = name)) + 
      geom_point() + 
      geom_line() + 
      scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
      scale_x_continuous(expand = c(0,0)) +
      
      ggtitle("Log Gene Counts") + 
      ylab("Log RNA") + 
      xlab("Time") +
      theme_bw()+
      theme(text = element_text(size = 14),
            legend.position = "none")
    
  }
  
  loggene_plot
  

}

plot_sim_data_short <- function(full_sim_data) {
  
  maxI <- max(full_sim_data$I) 
  maxRe <- max(full_sim_data$R1)
  maxgene <- max(c(full_sim_data$log_gene_copies1, 
                   full_sim_data$log_gene_copies2, 
                   full_sim_data$log_gene_copies3))
  
  maxcase <- max(full_sim_data$cases)
  maxtest <- max(full_sim_data$tests)
  

  
  sim_rt <- full_sim_data %>%
    ggplot(aes(x = time, y = true_rt)) + 
    geom_line() +
    ggtitle("Rt") +
    ylab("Rt") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  
  
  sim_data_long <- full_sim_data %>%
    ungroup() %>%
    dplyr::select(time, log_gene_copies1, log_gene_copies2, log_gene_copies3) %>%
    pivot_longer(cols = -time, values_to = "log_genes") 
  
  
  full_sim_data <- full_sim_data %>% 
    mutate(mean_log = (log_gene_copies1 + log_gene_copies2 + log_gene_copies3)/3)
  loggene_plot <- sim_data_long %>%
    ggplot(aes(x = time, y = log_genes, color = name)) + 
    geom_point() + 
    geom_line() + 
    geom_line(data = full_sim_data, aes(x = time, y = mean_log), color = "black", size = 1.2) +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxgene)) +
    scale_x_continuous(expand = c(0,0)) +
    
    ggtitle("Log Gene Counts") + 
    ylab("Log RNA") + 
    xlab("Time") +
    theme_bw()+
    theme(text = element_text(size = 14),
          legend.position = "none")
  
  sim_cases <- full_sim_data %>%
    ggplot(aes(x = time, y = cases)) +
    geom_point() +
    geom_line() +
    ggtitle("Observed Daily Cases") +
    ylab("Cases") +
    xlab("Time") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxcase)) +
    
    scale_x_continuous(expand = c(0,0)) +
    theme(text = element_text(size = 14)) +
    theme_bw()+
    theme(text = element_text(size = 14))
  

  full_plot <- (sim_rt + loggene_plot + sim_cases) +
    plot_layout(guides = 'collect')
  
}


# make fixed posterior samples --------------------------------------------

make_fixed_posterior_samples <- function(posterior_gq) {
  posterior_gq_samples <- posterior_gq %>%
  filter(str_detect(name, "\\[\\d+\\]", negate = T))
  
  return(posterior_gq_samples)
}


# make time varying posterior quantiles -----------------------------------
make_timevarying_posterior_quantiles <- function(posterior_gq) {
    timevarying_posterior_quantiles <-
    posterior_gq %>%
    filter(str_detect(name, "\\[\\d+\\]")) %>%
    mutate(time = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    left_join(.,tibble(time = 0:max(.$time)))
    
    return(timevarying_posterior_quantiles)
  
}


# make rt plot ------------------------------------------------------------
make_rt_plot <- function(timevarying_quantiles, 
                         simdata, 
                         initial_conds, 
                         varnames, 
                         varnum = 8,
                         SEIRR = TRUE){
  
  if (SEIRR == TRUE) {
    rt_plot <- timevarying_quantiles %>%
      filter(name == var_names[varnum]) %>%
      mutate(new_time = time) %>%
      ggplot() +
      geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper,fill = fct_rev(ordered(.width)))) +
      geom_point(data = simdata, mapping = aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept=1, linetype="dashed", color = "black") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw()+
      theme(legend.position = "none") + 
      ylab("Rt") +
      xlab("Time") +
      ggtitle("Rt Wastewater") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
      scale_x_continuous(expand = c(0,0))
    
  }
  
  if (SEIRR == FALSE) {
    timevarying_quantiles <- timevarying_quantiles %>% 
      filter(name == var_names[varnum]) %>%
      mutate(new_time = time * 7)
                             
    last_quantiles <- timevarying_quantiles %>% 
                      filter(new_time == max(new_time)) %>%
                      mutate(new_time = max(simdata$new_time))
    
    rt_quantiles <- bind_rows(timevarying_quantiles, last_quantiles)
    rt_plot <- rt_quantiles %>%
      ggplot() +
      geom_lineribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper,fill = fct_rev(ordered(.width)))) +
      geom_point(data = simdata, mapping = aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept=1, linetype="dashed", color = "black") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw()+
      theme(legend.position = "none") + 
      ylab("Rt") +
      xlab("Time") +
      ggtitle("Rt Wastewater") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
      scale_x_continuous(expand = c(0,0))
    
    
  }
  
  return(rt_plot)
  
}


# make trajectory plot ----------------------------------------------------

make_comp_plot <- function(posterior_timevarying_quantiles, 
                           initial_conds, 
                           traj_names,
                           true_trajectory, 
                           cases = FALSE){
  if (cases == FALSE) {
    true_init_trajectory <- initial_conds %>% dplyr::select(S, E, I , R1, log_genes_mean, new_time) %>%
      pivot_longer(cols = -(new_time)) %>%
      rename("true_value" = "value")
  } else {
    true_init_trajectory <- initial_conds %>% dplyr::select(S, E, I, new_time) %>%
      pivot_longer(cols = -(new_time)) %>%
      rename("true_value" = "value")
    
  }

  if (cases == FALSE) {
    posterior_trajectory <- posterior_timevarying_quantiles %>%
      filter(name %in% traj_names) %>%
      left_join(true_trajectory, by = c("time" = "new_time", 
                                        "name"= "name"))
  }

  else {
    posterior_trajectory <- posterior_timevarying_quantiles %>%
      filter(name %in% traj_names) %>%
      left_join(true_trajectory, by = c("time" = "new_week", 
                                        "name"= "name"))
    }
  
  traj_plot <- posterior_trajectory %>%
    ggplot() +
    geom_lineribbon(aes(x = time, y = value, ymin = .lower, ymax = .upper)) +
    geom_point(data = posterior_trajectory, mapping = aes(x = time, y = true_value), color = "coral1") +
    geom_point(data = true_init_trajectory, mapping = aes(x = new_time, y = true_value), shape= 0, color = "goldenrod2") +
    scale_fill_brewer(name = "Credible Interval Width") +
    theme_bw() +
    facet_wrap(.~name, scales = "free_y") +
    ggtitle("Posteriors of ODE compartments")

  return(traj_plot)
}


# make plot of prior and posteriors for fixed params ----------------------
make_fixed_param_plot <- function(posterior_gq_samples, prior_gq_samples){
  priors_and_posteriors <- rbind(posterior_gq_samples %>% dplyr::select(name, value, type, true_value), prior_gq_samples)
  
  
  param_plot <- priors_and_posteriors %>%
    ggplot(aes(value, type, fill = type)) +
    stat_halfeye(normalize = "xy")  +
    geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
    facet_wrap(. ~ name, scales = "free_x") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    ggtitle("Fixed parameter prior and posteriors")
  
  return(param_plot)
  
}


# make posterior predictive intervals -------------------------------------
make_post_pred_intervals <- function(posterior_predictive, sim_data, cases = FALSE, ten_sim = FALSE, three_mean = FALSE){
  if (cases == FALSE & ten_sim == FALSE & three_mean == FALSE) {
    obs_time <- sim_data %>%
    dplyr::select(log_gene_copies1, log_gene_copies2, log_gene_copies3, new_time) %>%
    pivot_longer(-(new_time)) %>%
    filter(value > 0) %>%
    arrange(name) %>%
    mutate(obs_index = row_number()) %>%
    dplyr::select(obs_index, new_time)
  
  posterior_predictive_samples <- posterior_predictive %>%
    pivot_longer(-c(iteration, chain)) %>%
    mutate(obs_index = name %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.numeric(),
           name = name %>%
             str_extract("^.+(?=\\[)") %>%
             str_remove("data_")) %>%
    bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                summarize(value = sum(value),
                          .groups = "drop"))  %>%
    left_join(obs_time, by = "obs_index")
  
  posterior_predictive_intervals <- posterior_predictive_samples %>%
    select(new_time, name, value) %>%
    group_by(new_time, name) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    select(new_time, name, value, starts_with("."))
  } else if (cases == FALSE & ten_sim == TRUE & three_mean == FALSE) {
    obs_time <- sim_data %>%
      dplyr::select(log_gene_copies1, 
                    log_gene_copies2, 
                    log_gene_copies3, 
                    log_gene_copies4, 
                    log_gene_copies5, 
                    log_gene_copies6, 
                    log_gene_copies7,
                    log_gene_copies8,
                    log_gene_copies9,
                    log_gene_copies10,
                    new_time) %>%
      pivot_longer(-(new_time)) %>%
      filter(value > 0) %>%
      arrange(name) %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_time)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
    posterior_predictive_intervals <- posterior_predictive_samples %>%
      select(new_time, name, value) %>%
      group_by(new_time, name) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)) %>%
      select(new_time, name, value, starts_with("."))
    
  } else if (cases == FALSE & ten_sim == FALSE & three_mean == TRUE) {
    obs_time <- sim_data %>%
      dplyr::select(log_mean_copiesthree,
                    new_time) %>%
      pivot_longer(-(new_time)) %>%
      filter(value > 0) %>%
      arrange(name) %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_time)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
    posterior_predictive_intervals <- posterior_predictive_samples %>%
      select(new_time, name, value) %>%
      group_by(new_time, name) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)) %>%
      select(new_time, name, value, starts_with("."))
    
  }
  else {
    obs_time <- sim_data %>%
      mutate(obs_index = row_number()) %>%
      dplyr::select(obs_index, new_week)
    
    posterior_predictive_samples <- posterior_predictive %>%
      pivot_longer(-c(iteration, chain)) %>%
      mutate(obs_index = name %>%
               str_extract("(?<=\\[)\\d+(?=\\])") %>%
               as.numeric(),
             name = name %>%
               str_extract("^.+(?=\\[)") %>%
               str_remove("data_")) %>%
      bind_rows(., group_by(., chain, iteration, obs_index, name) %>%
                  summarize(value = sum(value),
                            .groups = "drop"))  %>%
      left_join(obs_time, by = "obs_index")
    
    posterior_predictive_intervals <- posterior_predictive_samples %>%
      select(new_week, name, value) %>%
      group_by(new_week, name) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)) %>%
      select(new_week, name, value, starts_with("."))
    
  }
  
  return(posterior_predictive_intervals)
}


# make posterior predictive plot ------------------------------------------
make_post_pred_plot <- function(posterior_predictive_intervals, 
                                sim_data, 
                                cases = FALSE,
                                ten_sim = FALSE,
                                three_mean = FALSE) {

  if (cases == FALSE & ten_sim == FALSE & three_mean == FALSE) {
    true_data <- sim_data %>%
    dplyr::select(new_time, 
                  log_gene_copies1, 
                  log_gene_copies2, 
                  log_gene_copies3) %>%
    rename("log_copies1" = "log_gene_copies1",
           "log_copies2" = "log_gene_copies2",
           "log_copies3" = "log_gene_copies3") %>%
    pivot_longer(cols = - new_time) %>%
    rename("true_value" = "value") %>%
    filter(true_value > 0)
  
  posterior_predictive_intervals <- posterior_predictive_intervals %>%
    left_join(true_data, by = c("new_time"))
  
  posterior_predictive_plot <- posterior_predictive_intervals %>%
    ggplot() +
    geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
    geom_line(aes(x = new_time, y = value)) + 
    geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
    scale_fill_brewer(name = "Credible Interval Width") +
    theme_bw() + 
    ggtitle("Posterior Predictive (ODE)")
  } else if (cases == FALSE & ten_sim == TRUE & three_mean == FALSE) {
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_gene_copies1, 
                    log_gene_copies2, 
                    log_gene_copies3,
                    log_gene_copies4,
                    log_gene_copies5,
                    log_gene_copies6,
                    log_gene_copies7,
                    log_gene_copies8,
                    log_gene_copies9,
                    log_gene_copies10,
      ) %>%
      rename("log_copies1" = "log_gene_copies1",
             "log_copies2" = "log_gene_copies2",
             "log_copies3" = "log_gene_copies3",
             "log_copies4" = "log_gene_copies4",
             "log_copies5" = "log_gene_copies5",
             "log_copies6" = "log_gene_copies6",
             "log_copies7" = "log_gene_copies7",
             "log_copies8" = "log_gene_copies8",
             "log_copies9" = "log_gene_copies9",
             "log_copies10" = "log_gene_copies10") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  } else if (cases == FALSE & ten_sim == FALSE & three_mean == TRUE) {
    true_data <- sim_data %>%
      dplyr::select(new_time, 
                    log_mean_copiesthree
      ) %>%
      rename("log_mean_copies" = "log_mean_copiesthree") %>%
      pivot_longer(cols = - new_time) %>%
      rename("true_value" = "value") %>%
      filter(true_value > 0)
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_time"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_time, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_time, y = value)) + 
      geom_point(mapping = aes(x = new_time, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  }
  else {
    true_data <- sim_data %>%
      dplyr::select(new_week, 
                    total_cases) %>%
      rename("true_value" = "total_cases")
    
    posterior_predictive_intervals <- posterior_predictive_intervals %>%
      left_join(true_data, by = c("new_week"))
    
    posterior_predictive_plot <- posterior_predictive_intervals %>%
      ggplot() +
      geom_ribbon(aes(x = new_week, y = value, ymin = .lower, ymax = .upper, fill = fct_rev(ordered(.width)))) +
      geom_line(aes(x = new_week, y = value)) + 
      geom_point(mapping = aes(x = new_week, y = true_value), color = "coral1") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw() + 
      ggtitle("Posterior Predictive (ODE)")
    
  }
  
  return(posterior_predictive_plot)
  
}


# calculate huisman rt ----------------------------------------------------
calculate_huisman_rt <- function(simdata, sim = TRUE, mean_si, std_si) {
  
  if (sim == TRUE){
    date_length <- length(simdata$time)
    fake_date_one <- ymd("20220110")
    fake_date <- fake_date_one
    for (i in 2:date_length) {
      fake_date <- append(fake_date, fake_date[i-1] + ddays(2))
    }
    
    mean_data <- simdata %>%
      mutate(expone = exp(log_gene_copies1),
             exptwo = exp(log_gene_copies2),
             expthree = exp(log_gene_copies3)) %>%
      group_by(new_time) %>%
      summarise(mean_gene = (expone + exptwo + expthree)/3)
    
    
    
    
    data <- mean_data
    dates <- fake_date
    data <- cbind(data, date = dates)
    
  }
  
  if (sim == FALSE) {
    data <- simdata %>%
            mutate(mean_gene = (gene_copy_1 + gene_copy_2 + gene_copy_3)/3)
  }
  
    data <- data %>%
      mutate(orig_data = TRUE) %>%
      complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
      mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
      mutate(region = 'LA') %>%
      distinct() 
    
    ###########################################################
    ## Normalisation of WW data ####
    
    # Min observed
    norm_min <- min(data$mean_gene)
    
    ### ALL Normalised WW DATA ####
    ww_data = data %>%
      mutate(norm_mean_gene = mean_gene/norm_min)
    
    ###########################################################
    ##### Deconvolve and Estimate WW Re #####
    
    config_df = expand.grid("region" = c('LA'),  
                            'incidence_var' = c('norm_mean_gene'),
                            'FirstGamma' = 'sim_latent',
                            'SecondGamma' = 'sim_sld' )
    
    
    deconv_ww_data <- data.frame()
    Re_ww <- data.frame()

    for(row_i in 1:nrow(config_df)){
      row_i = 1
      new_deconv_data = deconvolveIncidence(ww_data %>% filter(region == config_df[row_i, 'region']), 
                                            incidence_var = config_df[row_i, 'incidence_var'],
                                            getCountParams(as.character(config_df[row_i, 'FirstGamma'])), 
                                            getCountParams(as.character(config_df[row_i, 'SecondGamma'])),
                                            smooth_param = TRUE, n_boot = 1000)
      
      new_deconv_data <- new_deconv_data %>%
        mutate(incidence_var = config_df[row_i, 'incidence_var'])
      
      ##### Get Re #####
      new_Re_ww = getReBootstrap(new_deconv_data, mean_si = mean_si, std_si = std_si)
      new_Re_ww <- new_Re_ww %>%
        mutate(variable = config_df[row_i, 'incidence_var'],
               region = config_df[row_i, 'region'])
      
      deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
      Re_ww = bind_rows(Re_ww, new_Re_ww)
      
    }
  
  
  return(Re_ww)
  
}
plot_huisman_rt <- function(Re_ww, simdata, full_simdata, sim = TRUE) {
  if (sim == TRUE){
    
    
    date_ranges <- Re_ww %>%
      group_by(region) %>%
      summarise(min_date = min(date),
                max_date = max(date)) 
    
    date_length <- length(simdata$time)
    fake_date_one <- ymd("20220110")
    fake_date <- fake_date_one
    for (i in 2:date_length) {
      fake_date <- append(fake_date, fake_date[i-1] + ddays(2))
    }
    
    true_rt <- simdata %>% dplyr::select(time, new_time, true_rt) %>% cbind(date = fake_date)
    initial_conds <- full_simdata %>% mutate(new_time = time - start_time + 7) %>% filter(new_time == 0)
    
    ## Plot Re ####
    
    Re_ww <- Re_ww %>%
      left_join(true_rt, by = "date")
    huisman_Re_plot <- Re_ww %>%
      filter(new_time >= 15) %>%
      mutate(.width = "0.95") %>%
      ggplot() + 
      geom_lineribbon(aes(x = new_time, 
                          y= median_R_mean, 
                          ymin = median_R_lowHPD, 
                          ymax = median_R_highHPD)) + 
      geom_point(data = true_rt, aes(x = new_time, y = true_rt), color = "coral1") +
      geom_point(data = initial_conds, mapping = aes(x = new_time, y = true_rt), shape = 0, color = "goldenrod2") + 
      geom_hline(yintercept = 1, linetype = "dotted") +
      scale_fill_brewer(name = "Credible Interval Width") +
      theme_bw() + 
      theme(legend.position = c(0.6, 0.8), legend.background = element_rect("transparent")) +
      ylab("Rt") +
      xlab("time") +
      ggtitle("Rt wastewater (state-of-art)") +     
      scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) +
      scale_x_continuous(expand = c(0,0))
    
  }
  
  return(huisman_Re_plot)
}


# freqentist metrics ---------------------------------------------------------------

# used for calculating frequentist characteristics of inference for rt
freq_metrics<- function(data, value, true_value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - {{ true_value }}),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = {{ true_value }} >= {{ lower }} & {{ true_value }} <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag({{ true_value }}),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs({{ true_value }} - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}

# rt_metrics --------------------------------------------------------------
# used for calculating frequentist characteristics of inference for rt
rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - true_rt),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag(true_rt),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs(true_rt - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}






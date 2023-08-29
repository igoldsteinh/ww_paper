# Functions for simulating data 
library(fields)
library(tidyverse)
# stochastic seirr --------------------------------------------------------
# start with the case of only infectious_individuals
# based on code by Phil O'Neil and Theodore Kypraios
sim_SEIRR <- function(N, I_init, beta, gamma, nu, eta) {
  # initial number of infectives and susceptibles;
  R1 <- 0
  E <- 0
  I <- I_init
  S <- N-I_init;
  
  # recording time;
  # create initial infectious, recover and stop shedding times
  init_infectious <- rep(rexp(1, rate = gamma), I)
  init_recover <- init_infectious + rexp(I, rate = nu)
  init_stopshed <- init_recover + rexp(I, rate = eta) 
  
  t <- unique(init_infectious);
  times <- c(rep(t,I_init));
  
  

  recover_times <- c(init_recover)
  stopshed_times <- c(init_stopshed)
  infectious_times <- c(rep(NA, I))
  # a vector which records the type of event (1=infection, 2=infectious, 3 = recovered, 4 = stopshedding)
  type <- c(rep(2, length(times)));
  
  # a counter for labelling the individuals
  lambda <- I;
  
  # a vector to store the labels
  labels <- c(1:lambda);
  
  while (I > 0 | R1 > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0 & I > 0) {
      infec_time  <- rexp(1, (beta/N)*I*S)
    }
    else {
      infec_time <- Inf;
    }
    
    # next become infectious 
    next_infectious <- min(infectious_times, na.rm = TRUE)
    
    # time to next removal
    next_recover <- min(recover_times, na.rm=TRUE)
    
    # next stop shed
    next_stopshed <- min(stopshed_times, na.rm = TRUE)

    if (t + infec_time < min(next_recover, next_stopshed, next_infectious, na.rm = TRUE)) {
      # infection occurs
      E <- E + 1
      S <- S-1;
      
      # simulate all of the next event times for this individual
      new_infectious_period <- rexp(1, gamma)
      new_recover_period <- rexp(1,nu)
      new_stopshed_period <- rexp(1, eta)
      infectious_times <- append(infectious_times, t + infec_time + new_infectious_period)
      recover_times <- append(recover_times, t + infec_time + new_infectious_period + new_recover_period)
      stopshed_times <- append(stopshed_times, t + infec_time + new_infectious_period + new_recover_period + new_stopshed_period)

      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      type <- append(type, 1);
      times <- append(times, t + infec_time);
      t <- t + infec_time
    }
    else if (next_infectious < min(next_recover, next_stopshed, na.rm = TRUE)){
      #transition to I occurs
      I <- I + 1
      E <- E - 1
      type <- append(type, 2);
      index_infectious_time <- which(infectious_times==next_infectious)
      infectious_times[index_infectious_time] <- NA
      labels <- append(labels, index_infectious_time)
      times <- append(times, next_infectious);
      t <- next_infectious     
    }
    
    else if (next_recover < next_stopshed) {
      # transition from I to R1 occurs
      I <- I - 1
      R1 <- R1 + 1
      
      type <- append(type, 3)
      index_recover_time <- which(recover_times == next_recover)
      recover_times[index_recover_time] <- NA
      labels <- append(labels, index_recover_time)
      times <- append(times, next_recover)
      t <- next_recover
    }
    
    else {
      # transition out of R1 occurs
      R1 <- R1 - 1
      
      type <- append(type, 4)
      index_stopshed_time <- which(stopshed_times == next_stopshed)
      stopshed_times[index_stopshed_time] <- NA
      labels <- append(labels, index_stopshed_time)
      times <- append(times, next_stopshed)
      t <- next_stopshed
    }
  }
  
  
  res <- data.frame("t"=times, "type"=type, "labels"=labels);
  res
}

# stochastic seirr with non-constant beta--------------------------------------------------------
# start with the case of only infectious_individuals
# based on code by Phil O'Neil and Theodore Kypraios
sim_SEIRR_nonconst <- function(N, I_init, E_init, beta_init, beta_vec, change_points, gamma, nu, eta) {
  # initial number of infectives and susceptibles;
  R1 <- 0
  E <- E_init
  I <- I_init
  S <- N-I_init;
  
  # ok for everyone who is already in a compartment, their times will be the same
  # but then their subsequent times will be different 
  #  the E compartment
  E_infection <- rep(0, E) # infected at time 0
  E_infectious <- E_infection + rexp(E, rate = gamma)
  E_recover <- E_infectious + rexp(E, rate = nu)
  E_stopshed <- E_recover + rexp(E, rate = eta)
  
  E_types = rep(1, E)
  # do the same for I compartment
  I_infectious <- rep(0, I)
  I_recover <- I_infectious + rexp(I, rate = nu)
  I_stopshed <- I_recover + rexp(I, rate = eta)
  
  I_types = rep(2, I)
  
  # create initial infectious, recover and stop shedding times
  beta_t = beta_init
  t <- 0;
  times <- c(rep(0,E + I + R1));
  
  
  
  recover_times <- c(E_recover, I_recover)
  stopshed_times <- c(E_stopshed, I_stopshed)
  infectious_times <- c(E_infectious, I_infectious)
  
  # set the zeros to NA so that they are not valid times to use
  infectious_times[which(infectious_times == 0)] <- NA
  # a vector which records the type of event (1=infection, 2=infectious, 3 = recovered, 4 = stopshedding, 5 = change point)
  type <- c(E_types, I_types);
  
  # a counter for labelling the individuals
  lambda <- E + I + R1; # lambda is the last individual
  
  # a vector to store the labels, i.e., the ids of the individuals corresponding to the time and type of action which has taken place
  labels <- c(1:lambda); # this is each individual from 1 to the last
  
  # a dataframe to store the states as well as beta and r0
  states <- as.matrix(t(c(t,S,E,I,R1, beta_t, beta_t/nu)))  
    while (I > 0 | R1 > 0) {
    
    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0 & I > 0) {
      infec_time  <- rexp(1, (beta_t/N)*I*S)
    }
    else {
      infec_time <- Inf;
    }
    
    # change points
    next_cp <- min(change_points, na.rm = TRUE)
    
    # next become infectious 
    next_infectious <- min(infectious_times, na.rm = TRUE)
    
    # time to next removal
    next_recover <- min(recover_times, na.rm=TRUE)
    
    # next stop shed
    next_stopshed <- min(stopshed_times, na.rm = TRUE)

    if (next_cp < min(next_recover, next_stopshed, next_infectious, t+ infec_time, na.rm = TRUE)) {
      
      index = which(change_points == next_cp)
      beta_t = beta_vec[index]
      change_points[index] = NA
      # a label of 2 times the population size (ie impossible for an individual) is code for this was a change point
      labels <- append(labels, NA)
      type = append(type, 5)
      times <- append(times, next_cp)
      t <- next_cp
      
    } else if (t + infec_time < min(next_recover, next_stopshed, next_infectious, na.rm = TRUE)) {
      # infection occurs
      E <- E + 1
      S <- S-1;
      
      # simulate all of the next event times for this individual
      new_infectious_period <- rexp(1, gamma)
      new_recover_period <- rexp(1,nu)
      new_stopshed_period <- rexp(1, eta)
      infectious_times <- append(infectious_times, t + infec_time + new_infectious_period)
      recover_times <- append(recover_times, t + infec_time + new_infectious_period + new_recover_period)
      stopshed_times <- append(stopshed_times, t + infec_time + new_infectious_period + new_recover_period + new_stopshed_period)
      
      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      type <- append(type, 1);
      times <- append(times, t + infec_time);
      t <- t + infec_time
      new_states = t(c(t, S, E, I, R1, beta_t, beta_t/nu))
      states = rbind(states, new_states)
    }
    else if (next_infectious < min(next_recover, next_stopshed, na.rm = TRUE)){
      #transition to I occurs
      I <- I + 1
      E <- E - 1
      type <- append(type, 2);
      index_infectious_time <- which(infectious_times==next_infectious)
      infectious_times[index_infectious_time] <- NA
      labels <- append(labels, index_infectious_time)
      times <- append(times, next_infectious);
      t <- next_infectious   
      new_states = t(c(t, S, E, I, R1, beta_t, beta_t/nu))
      states = rbind(states, new_states)
      
    }
    
    else if (next_recover < next_stopshed) {
      # transition from I to R1 occurs
      I <- I - 1
      R1 <- R1 + 1
      
      type <- append(type, 3)
      index_recover_time <- which(recover_times == next_recover)
      recover_times[index_recover_time] <- NA
      labels <- append(labels, index_recover_time)
      times <- append(times, next_recover)
      t <- next_recover
      new_states = t(c(t, S, E, I, R1, beta_t, beta_t/nu))
      states = rbind(states, new_states)
      
    }
    
    else {
      # transition out of R1 occurs
      R1 <- R1 - 1
      
      type <- append(type, 4)
      index_stopshed_time <- which(stopshed_times == next_stopshed)
      stopshed_times[index_stopshed_time] <- NA
      labels <- append(labels, index_stopshed_time)
      times <- append(times, next_stopshed)
      t <- next_stopshed
      new_states = t(c(t, S, E, I, R1, beta_t, beta_t/nu))
      states = rbind(states, new_states)
      
    }
  }
  
  
  res <- list(data.frame("t"=times, "type"=type, "labels"=labels), states)
  res
}


# function for creating day level state data ------------------------------
create_daily_data <- function(sim_states) {
  day_data <- data.frame(sim_states) %>%
    rename("time" = "X1",
           "S" = "X2",
           "E" = "X3", 
           "I" = "X4", 
           "R1" = "X5", 
           "beta_t" = "X6",
           "R0" = "X7") %>%
    mutate(integer_day = ceiling(time),
           time_diff = integer_day - time) %>%
    group_by(integer_day) %>% 
    filter(time_diff == min(time_diff)) 
  
  return(day_data)
}

# function for calculating individual gene counts at time t ---------------
x <- c(-3, -1, 1, 3, 5, 7, 9, 13, 17, 21, 25, 29)
y <- c(5, 6.9, 6.7, 6.5, 6.2, 5.9, 5.5, 4.75, 3.9, 3.3, 2.1, 1.3)
x_adj <- x + 3

#translate from log base 10 scale to real scale
exp_y <- 10^(y)
#tps cant do this for some reason, so instead do the spline on log base 10
data <- data.frame(x_adj, exp_y, y)
# tp_spline <- Tps(data$x_adj, data$y)
set.seed(1234)
tp_spline <- Tps(data$x_adj, data$y)

calc_individ_counts <- function(t) {
  mean_log10_prediction<- predict(tp_spline, t)
  
  # hard cutoff at zero
  if (mean_log10_prediction < 0){
    prediction <- 0
  } else {
    prediction <- 10^(rnorm(1, mean_log10_prediction, sd = 1.09)) # should be 1.09
  }
  
  return(prediction)
}
# function for calculating gene counts at time t --------------------------
calc_total_gene_counts <- function(epi_curve, time, N) {
  # first collect all active individuals
  # epi_curve = wide_format
  # time = 5
  active_individuals <- epi_curve %>% 
                      ungroup() %>% 
                        mutate(infectious_individuals = infectious_time <= time & recover_time > time,
                               r1_individuals = recover_time <= time & stopshed_time > time) %>%
                        filter(infectious_individuals == TRUE | r1_individuals == TRUE)
  
  if (dim(active_individuals)[1] > 0) {
      active_individuals <- active_individuals %>%
                            mutate(time_since_infectious = time - infectious_time) %>% 
                            rowwise() %>% 
                            mutate(gene_counts = calc_individ_counts(time_since_infectious))
      
  # report the total gene counts, and the total number in each compartment
  res <- data.frame("time" = time,
                    "total_count" = sum(active_individuals$gene_counts), 
                    "total_count_I" = sum(active_individuals$gene_counts[active_individuals$infectious_individuals == TRUE]),
                    "total_count_R1" = sum(active_individuals$gene_counts[active_individuals$r1_individuals == TRUE]),
                    "num_I" = sum(active_individuals$infectious_individuals),
                    "num_R1" = sum(active_individuals$r1_individuals))
  
  return(res)
  } else {
    res <- data.frame("time" = time,
                      "total_count" = 0, 
                      "total_count_I" = 0,
                      "total_count_R1" = 0,
                      "num_I" = 0,
                      "num_R1" = 0
                    )
    
    return(res)
    
  }
}


# function for counting transitions from E to I in time (t-1,t] -------------
calc_E2I_counts <- function(epi_curve, time) {
  E2I <- epi_curve %>% 
         filter(infectious_time > (time -1) & infectious_time <= time) %>%
         ungroup() %>% 
         summarise(total = n()) %>% pull()
  res <- data.frame("time" = time,
                    "E2I_transitions" = E2I)
  
  return(res)
}

# function for creating data set from sim results -------------------------
create_individ_data <- function(sims) {
  wide_forms = sims %>% group_by(labels) %>% 
    arrange(labels, type) %>%
    dplyr::select(t, type, labels) %>%
    pivot_wider(id_cols = labels, names_from = type, values_from = t ) %>%
    rename("infection_time" = `1`, "infectious_time" = `2`, "recover_time" = `3`, "stopshed_time" = `4`) %>%
    dplyr::select(labels, infection_time, infectious_time, recover_time, stopshed_time) %>%
    mutate(infectious_period = recover_time - infectious_time,
           r1_period = stopshed_time - recover_time) 
  
  return(wide_forms)
}


# function for finding the gene counts across many simulations ----------------
sim_and_calc_genes <- function(N, I_init, beta, gamma, nu, eta, num_sims) {
  seed = sample.int(num_sims)
  # set up output lists
  sims <- vector(mode='list', length=length(seed))
  wide_forms <- vector(mode='list', length=length(seed))
  counts <- vector(mode='list', length=length(seed))

  # generate simulations, record counts across simulations
  for (i in 1:length(seed)){
    set.seed(seed[i])
    # simulate epidemic
    sims[[i]] = sim_SEIRR(N = N, 
                          I_init = I_init, 
                          beta = beta, 
                          gamma = gamma, 
                          nu = nu, 
                          eta = eta)
    
    wide_forms[[i]] = create_individ_data(sims[[i]])

    # format as needed
    # find valid integer times for the simulation
    min_int = ceiling(min(sims[[i]]$t))
    max_int = floor(max(sims[[i]]$t))
    
    # create times to check the counts, bespoke to each simulation 
    times = min_int:max_int
    
    # calculate the counts
    counts[[i]] = map(times, ~calc_total_gene_counts(wide_forms[[i]], time = .x)) %>%
      bind_rows()
    
  
  }


  # combine counts, plot findings
  combined_counts <- bind_rows(counts, .id = "seed") 
  res <- list("sims" = sims, "wide_data" = wide_forms, "counts" = combined_counts, "seeds_used" = seed)
  return(res)
}


# simulate gene data from a realization of SEIRR --------------------------
simulate_gene_data <- function(true_gene_counts, seed, rho, t_sd, t_df){
  set.seed(seed)
  sim_data <- true_gene_counts %>% 
              mutate(log_genes_mean = log(total_count) + log(rho)) %>%
              rowwise() %>%
              mutate(
                      log_gene_copies1 = log_genes_mean  + (t_sd * rt(1,t_df)),
                      log_gene_copies2 = log_genes_mean  + (t_sd * rt(1,t_df)),
                      log_gene_copies3 = log_genes_mean  + (t_sd * rt(1,t_df)),
                      log_gene_copies4 = log_genes_mean  + (t_sd * rt(1,t_df)),
                      log_gene_copies5 = log_genes_mean  + t_sd * rt(1,t_df),
                      log_gene_copies6 = log_genes_mean  + t_sd * rt(1,t_df),
                      log_gene_copies7 = log_genes_mean  + t_sd * rt(1,t_df),
                      log_gene_copies8 = log_genes_mean  + t_sd * rt(1,t_df),
                      log_gene_copies9 = log_genes_mean  + t_sd * rt(1,t_df),
                      log_gene_copies10 = log_genes_mean + t_sd * rt(1,t_df)) %>%
              mutate(log_mean_copiesten = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3) +
                                       exp(log_gene_copies4) + exp(log_gene_copies5) + exp(log_gene_copies6) +
                                       exp(log_gene_copies7) + exp(log_gene_copies8) + exp(log_gene_copies9) +
                                       exp(log_gene_copies10))/10),
                     log_mean_copiesthree = log((exp(log_gene_copies1) + exp(log_gene_copies2) + exp(log_gene_copies3))/3))
  
  return(sim_data)
      
}

# simulate case data ------------------------------------------------------
simulate_case_data <- function(true_E2I_counts, seed, rho, phi) {
  set.seed(seed)
  size = phi
  sim_data <- true_E2I_counts %>%
              rowwise() %>% 
              mutate(cases = rnbinom(1, size, prob = 1/(1 + ((E2I_transitions * rho)/phi))))
}


# simulate agent-based seirr ----------------------------------------------

sim_agent_SEIRR <- function(pop_size, I_init, beta, gamma, nu, eta) {
  # setting up initial states and frames
  pop_vec <- rep("S", length.out = pop_size)
  pop_vec[1:I_init] <- "I"
  
  rate_vec <- rep(0, length.out = pop_size)
  id_vec <- 1:pop_size
  
  # copying what I did in the time-varying model
  t <- 0
  
  rate_frame <- data.frame(id = id_vec, 
                           state = pop_vec,
                           rate = rate_vec)
  
  # setting up the output frame
  init_state <- t(c(t, rate_frame$state))
  state_frame <- NULL
  state_frame <- rbind(state_frame, init_state)
  
  # have to count num infectious and num R1 to copy other simulation
  num_infectious <- rate_frame %>% 
    mutate(isI = state == "I") %>% 
    summarise(total = sum(isI)) %>%
    pull(total)
  
  num_r1 <- rate_frame %>% 
    mutate(isR1 = state == "R1") %>% 
    summarise(total = sum(isR1)) %>%
    pull(total)
  
  # frame of individual rates for each individual
  rate_frame <- rate_frame %>% 
    mutate(rate = ifelse(state == "S", beta * num_infectious, 
                         ifelse(state == "E", gamma, 
                                ifelse(state == "I", nu, 
                                       ifelse (state == "R1", eta, 0)))))
  
  while (num_infectious > 0 | num_r1 > 0) {
    # time to next event is the min of all possible events
    next_event = rexp(1, rate = sum(rate_frame$rate))
    
    # choose which event happens proportional to the rates
    which_id = sample(rate_frame$id, 1, prob = rate_frame$rate)
    
    # update the state of the individual who changes
    rate_frame$state[which(rate_frame$id == which_id)] <- ifelse(rate_frame$state[which(rate_frame$id == which_id)] == "S", "E", 
                                                                   ifelse(rate_frame$state[which(rate_frame$id == which_id)] == "E", "I", 
                                                                            ifelse(rate_frame$state[which(rate_frame$id == which_id)] == "I", 
                                                                                   "R1", "R2")))
    # update total infections and r1 counts
    num_infectious <- rate_frame %>% 
      mutate(isI = state == "I") %>% 
      summarise(total = sum(isI)) %>%
      pull(total)
    
    num_r1 <- rate_frame %>% 
      mutate(isR1 = state == "R1") %>% 
      summarise(total = sum(isR1)) %>%
      pull(total)
    
    #re-calculate rate frame (changes with num_infectious, as well state changes)
    rate_frame <- rate_frame %>% 
      mutate(rate = ifelse(state == "S", beta * num_infectious, 
                           ifelse(state == "E", gamma, 
                                  ifelse(state == "I", nu, 
                                         ifelse (state == "R1", eta, 0)))))
    
    # update time 
    t <- t + next_event 
    current_state <- c(t, rate_frame$state)
    state_frame <- rbind(state_frame, current_state)
  }
  return(state_frame)
  
}


# create daily data from agent model --------------------------------------
create_daily_from_agent <- function(res) {
  day_data <- data.frame(res) %>%
    rename("time" = "X1") %>%
    pivot_longer(-time, names_to = "id", values_to = "state") 
  
  day_data$time <- as.numeric(day_data$time) 
  
  day_data <- day_data %>% 
    mutate(integer_day = ceiling(time),
           time_diff = integer_day - time) %>%
    group_by(integer_day) %>% 
    filter(time_diff == min(time_diff)) %>%
    group_by(integer_day) %>% 
    summarise(S = sum(state == "S"),
              E = sum(state == "E"),
              I = sum(state == "I"), 
              R1 = sum(state == "R1"),
              R2 = sum(state == "R2"))
  
  return(day_data)
}


# simulate classic gillespie seirr ----------------------------------------
gillespie_seirr <- function(pop_size, I_init, beta, gamma, nu, eta) {
  # initial number of infectives and susceptibles;
  R1 <- 0
  E <- 0
  I <- I_init
  S <- pop_size-I_init
  
  t <- 0;

  # a dataframe to store the states as well as beta and r0
  states <- as.matrix(t(c(t,S,E,I,R1)))  
  
  while (I >0 | R1 >0) {
    total_rate <- (beta/pop_size * I * S) + gamma*E + nu*I + eta*R1
    
    next_event <- rexp(1, total_rate)
    
    type_event <- sample(c("infection", "infectious", "recover", "stop_shed"), size = 1, prob = c((beta/pop_size * I * S),
                                                                                         gamma*E,
                                                                                         nu*I, 
                                                                                         eta*R1))
    
    t <- t + next_event
    if (type_event == "infection") {
      S <- S - 1
      E <- E + 1
    } else if (type_event == "infectious") {
      E <- E - 1
      I <- I + 1
    } else if (type_event == "recover") {
      I <- I - 1
      R1 <- R1 + 1
    } else {
      R1 <- R1 - 1
    }
    
    current_state <- t(c(t,S,E,I,R1))
    states <- rbind(states, current_state)
  }
  return(states)
}

# function for creating day level state data ------------------------------
create_gillespie_daily_data <- function(sim_states) {
  day_data <- data.frame(sim_states) %>%
    rename("time" = "X1",
           "S" = "X2",
           "E" = "X3", 
           "I" = "X4", 
           "R1" = "X5") %>%
    mutate(integer_day = ceiling(time),
           time_diff = integer_day - time) %>%
    group_by(integer_day) %>% 
    filter(time_diff == min(time_diff)) 
  
  return(day_data)
}


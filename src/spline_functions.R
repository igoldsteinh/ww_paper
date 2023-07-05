# Create automated process for choosing overdispersion parameter from a spline

run_nb_spline <- function(data,
                          response = c("hosp", "icu", "cases"),
                          seed = 17,
                          iter = 4000,
                          warmup = 1000,
                          thin = 10,
                          refresh = 0,
                          adapt_delta = 0.99) {

  if (response == "hosp") {
    spline_model <- brm(bf(hospitalizations ~ s(time)),
                        data = data, family = negbinomial(), cores = 4, seed = seed,
                        iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                        control = list(adapt_delta = adapt_delta))

  }

  if (response == "icu") {
    spline_model <- brm(bf(icu~ s(time)),
                        data = data, family = negbinomial(), cores = 4, seed = seed,
                        iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                        control = list(adapt_delta = adapt_delta))

  }
  
  if (response == "cases") {
    spline_model <- brm(bf(total_cases~ s(new_week)),
                        data = data, family = negbinomial(), cores = 4, seed = seed,
                        iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                        control = list(adapt_delta = adapt_delta))
    
  }

  return(spline_model)
}

compare_kappa_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- qlnorm(c(0.025, 0.975),
                                meanlog = candidate_params[1], 
                                sdlog = candidate_params[2])
  
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 + (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}


choose_kappa_params <- function(spline_posterior) {
  posterior_pars <- summary(spline_posterior)
  start_mean <- log(posterior_pars[["spec_pars"]][[1]])
  start_sd <- 0.3
  true_lb <- posterior_pars[["spec_pars"]][[3]]
  true_ub <- posterior_pars[["spec_pars"]][[4]]
  
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  
  optim_params <- optim(par = start_params,
                        fn = compare_kappa_quantiles, 
                        true_quantiles = true_quantiles)
}

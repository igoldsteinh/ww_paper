# epidemia functions
epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, y)
  pmf[1] <- pgamma(1.5, alpha, rate = beta)
  for (i in 2:y) {
    pmf[i] <- pgamma(i+.5, alpha, rate = beta) - pgamma(i-.5, alpha, rate = beta)
  }
  
  pmf
}

epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i+.5, meanlog = params[1], sdlog = params[2]) - 
      plnorm(i-.5, meanlog = params[1], sdlog = params[2])
  }
  
  pmf
}

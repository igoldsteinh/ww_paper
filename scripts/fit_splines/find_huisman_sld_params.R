# calculate Huisman SLD gamma parameters
library(fields)
library(tidyverse)
library(DescTools)

# function for calculating individual gene counts at time t based on Nourbakhsh
x <- c(-3, -1, 1, 3, 5, 7, 9, 13, 17, 21, 25, 29)
y <- c(5, 6.9, 6.7, 6.5, 6.2, 5.9, 5.5, 4.75, 3.9, 3.3, 2.1, 1.3)
x_adj <- x + 3.01


#translate from log base 10 scale to real scale, add 0,0 to anchor it
exp_y <- 10^(y)
x_adj <- c(0, x_adj)
exp_y <- c(0, exp_y)
data <- data.frame(x_adj, exp_y)
# tp_spline <- Tps(data$x_adj, data$y)
set.seed(1234)
tp_spline <- Tps(data$x_adj, data$exp_y)


# create values of the sld to match to a gamma -------------------
x_val = seq(0, 25, by = 0.1)
expy_val = rep(-1, length.out = length(x_val))
for (i in 1:length(x_val)) {
  expy_val[i] <- predict(tp_spline, x_val[i])
  
}

norm_y_val = (expy_val / AUC(x_val, expy_val, method = "spline")) 


# find a gamma using optim ------------------------------------------------

compare_y_vals <- function(candidate_params, x_vals, true_vals) {
  candidate_vals <- dgamma(x_vals, 
                           shape = candidate_params[1],
                           scale = candidate_params[2])
  
  loss <- sum((candidate_vals - true_vals)^2)
  return(loss)
}

choose_gamma_params <- function(true_vals, x_vals, start_shape, start_scale) {
  start_params <- c(start_shape, start_scale)

  optim_params <- optim(par = start_params,
                        fn = compare_y_vals, 
                        x_vals = x_vals,
                        true_vals = true_vals)
}



testing <- dgamma(x_val, shape = 9, scale = 0.5)
plot(x_val, testing)

gamma_params <- choose_gamma_params(norm_y_val, x_val, start_shape = 9, start_scale = 0.5)$par



test <- dgamma(x_val, shape = gamma_params[1], scale = gamma_params[2])

my_data <- data.frame(x = x_val, true_val = norm_y_val, test_val = test) %>% pivot_longer(-x)
my_plot <- my_data %>% 
           ggplot(aes(x = x, y = value, color = name)) + 
           geom_point() +
           theme_bw()

my_plot
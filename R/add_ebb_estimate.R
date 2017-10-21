library(tidyverse)
library(broom)
source(file = "R/deconvolveR.R")

add_ebb_estimate <- function(dat, total, x){
  delta <- 0.01
  df <- as.data.frame(list(n=dat$total,x=dat$x))
  tau <- seq(from = 0.01, to = 0.99, by = delta)
  result <- deconv(tau = tau, X = df, family = "Binomial", c0 = 1,pDegree = 5)
  theta <- result$stats[, 'theta']
  gTheta <- result$stats[, 'g']
  f_alpha <- function(n_k, x_k) {
    ## .01 is the delta_theta in the Riemann sum
    sum(dbinom(x = x_k, size = n_k, prob = theta) * gTheta) * .01
  }
  
  g_theta_hat <- function(n_k, x_k) {
    gTheta * dbinom(x = x_k, size = n_k, prob = theta) / f_alpha(n_k, x_k)
  }
  col1 <- enquo(total)
  col2 <- enquo(x)
  dat %>%
    rowwise() %>%
    mutate(.fitted=sum(theta*g_theta_hat(x_k = !!col1, n_k = !!col2))*delta)
}
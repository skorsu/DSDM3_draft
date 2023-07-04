rm(list = ls())
library(tidyverse)

### Data Simulation
data_sim <- function(N, K, J, clus_prop, pi_g, pi_w, s2, lb, ub){
  
  ### cluster assignment
  ci <- sample(1:K, N, replace = TRUE, prob = clus_prop)
  
  ### gamma for each observation based on its cluster (at-risk indicator)
  gm_mat <- t(apply(matrix(pi_g[ci]), 1, rbinom, n = J, size = 1))
  
  ### w (important variable)
  w <- rbinom(J, 1, pi_w)
  
  ### Generate xi = exp(beta) where beta ~ N(0, s2).
  beta <- matrix(rnorm(J * K, 0, sqrt(s2)), ncol = J)
  xi <- exp(beta)
  
  ### alpha matrix
  gm_w <- sweep(gm_mat, MARGIN = 2, w, `*`)
  w_xi <- sweep(xi[ci, ], MARGIN = 2, w, `*`)
  alpha_mat <- matrix(rgamma(N * J, gm_w * w_xi, 1), ncol = J)
  
  ### normalize the alpha matrix
  alpha_norm <- t(apply(alpha_mat, 1, function(x){x/sum(x)}))
  
  ### simulate zi
  sum_zi <- round(runif(N, lb, ub))
  zi <- matrix(NA, ncol = J, nrow = N)
  for(i in 1:N){
    zi[i, ] <- rmultinom(1, sum_zi[i], alpha_norm[i, ])
  }
  
  ### Return the simulated data
  list(ci = ci, gamma = gm_mat, w = w, beta = beta, z = zi)
}

set.seed(12)
test_dat <- data_sim(10, 2, 20, c(0.4, 0.6), c(0.4, 0.6), 0.5, 1, 50, 100)

new_gamma <- update_gamma(z = test_dat$z, clus_assign = test_dat$ci, w = test_dat$w, 
                          old_gamma = test_dat$gamma, xi = exp(test_dat$beta), 
                          b0g = 1, b1g = 1)

cbind(test_dat$z[5, ], test_dat$gamma[5, ], new_gamma[5, ])

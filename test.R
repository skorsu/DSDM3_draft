rm(list = ls())
library(tidyverse)

### Data Simulation
data_sim <- function(N, K, J, clus_prop, pi_g, pi_w, s2, lb_a, ub_a, lb_ia, ub_ia){
  
  ### cluster assignment
  ci <- sample(1:K, N, replace = TRUE, prob = clus_prop)
  
  ### gamma for each observation based on its cluster (at-risk indicator)
  gm_mat <- t(apply(matrix(pi_g[ci]), 1, rbinom, n = J, size = 1))
  
  ### w (important variable)
  w <- sort(rbinom(J, 1, pi_w), decreasing = TRUE)
  n_active <- sum(w)
  
  ### Generate xi = exp(beta) where beta ~ N(0, s2).
  beta <- matrix(rnorm(J * K, 0, sqrt(s2)), ncol = J)
  xi <- exp(beta)
  
  ### Generate the data separately for active and inactive cluster
  alpha_mat <- matrix(rgamma(K * J, xi, 1), ncol = J)
  
  ### active variable
  active_alpha <- alpha_mat[ci, (1:n_active)]
  gm_alpha_active <- gm_mat[, (1:n_active)] * active_alpha
  alpha_norm_active <- t(apply(gm_alpha_active, 1, function(x){x/sum(x)}))
  sum_zi_active <- round(runif(N, lb_a, ub_a))
  zi_active <- matrix(NA, ncol = n_active, nrow = N)
  for(i in 1:N){
    zi_active[i, ] <- rmultinom(1, sum_zi_active[i], alpha_norm_active[i, ])
  }
  
  ### inactive variable
  inactive_alpha <- alpha_mat[ci, (n_active+1):J]
  gm_alpha_inactive <- gm_mat[, (n_active+1):J] * inactive_alpha
  alpha_norm_inactive <- t(apply(gm_alpha_inactive, 1, function(x){x/sum(x)}))
  sum_zi_inactive <- round(runif(N, lb_ia, ub_ia))
  zi_inactive <- matrix(NA, ncol = J - n_active, nrow = N)
  for(i in 1:N){
    zi_inactive[i, ] <- rmultinom(1, sum_zi_inactive[i], alpha_norm_inactive[i, ])
  }
  
  zi <- cbind(zi_active, zi_inactive)
  
  ### Rearrange the variable order
  variable_order <- sample(1:J, J)
  gm_mat <- gm_mat[, variable_order]
  w <- w[variable_order]
  beta <- beta[, variable_order]
  zi <- zi[, variable_order]
  
  ### Return the simulated data
  list(ci = ci, gamma = gm_mat, w = w, beta = beta, z = zi)
}

set.seed(12)
test_dat <- data_sim(N = 100, K = 3, J = 50, clus_prop = c(0.3, 0.4, 0.3), 
                     pi_g = c(0.75, 0.9, 0.5), pi_w = 0.15, s2 = 1, lb_a = 20, ub_a = 30, lb_ia = 25, ub_ia = 35)


test_bet <- log_beta_k(beta_k = test_dat$beta[1, ], ci = 2 - 1, z = test_dat$z, 
                       gamma_mat = test_dat$gamma, beta = test_dat$beta, 
                       w = test_dat$w, clus_assign = test_dat$ci - 1, s2 = 1)


identical(round(test_bet, 1), round(log(dnorm(test_dat$beta[1, test_dat$w == 1], 0, sqrt(10))), 1))


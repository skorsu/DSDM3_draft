rm(list = ls())
library(tidyverse)

### Function: Simulate the data
### Note: the cluster index must start from 0
### Note: set the number of active variable = 2*K
data_sim <- function(K, J, r, pi_g, beta_base, lb, ub){

  ### Actual Assignment
  ci_actual <- rep(1:K, r)
  N <- length(ci_actual)
  
  ### Active Variables
  conc_mat <- 1.5 * diag(K)
  conc_mat[conc_mat == 0] <- 0.1
  conc_mat <- apply(conc_mat, 1, rep, r)
  conc_mat <- cbind(conc_mat, conc_mat)
  
  xi_mat <- matrix(NA, ncol = 2*K, nrow = K)
  z <- matrix(NA, ncol = 2*K, nrow = N)
  
  for(k in 1:K){
    xi_mat[k, ] <- exp(rnorm(2*K, beta_base[k], 0.1))
  }
  
  for(i in 1:N){
    sum_zi <- round(runif(1, lb, ub))
    z[i, ] <- rmultinom(1, sum_zi, (conc_mat[i, ] * xi_mat[ci_actual[i], ])/sum(conc_mat[i, ] * xi_mat[ci_actual[i], ]))
  }
  
  ### Non-active variable
  non_active <- J - (2*K)
  z_non_active <- matrix(NA, ncol = non_active, nrow = N)
  gm_mat <- matrix(NA, ncol = non_active, nrow = N)
  for(i in 1:N){
    gm_vec <- rbinom(non_active, 1, pi_g)
    gm_mat[i, ] <- gm_vec
    alp <- exp(rnorm(non_active, 0, 1)) * gm_vec
    sum_zi <- round(runif(1, lb, ub))
    z_non_active[i, ] <- rmultinom(1, sum_zi, alp/sum(alp))
  }
  
  list(conc_mat = conc_mat, ci = ci_actual - 1, xi = xi_mat, 
       z = cbind(z, z_non_active), 
       gamma = cbind(matrix(1, ncol = 2*K, nrow = N), gm_mat))
  
}

test_dat <- data_sim(K = 5, J = 50, r = 10, pi_g = 0.25,
                     beta_base = c(2, 1, 0, -1, -2), lb = 50, ub = 60)
table(test_dat$gamma[2, ], test_dat$z[2, ])


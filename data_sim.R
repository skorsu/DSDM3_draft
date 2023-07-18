rm(list = ls())
library(tidyverse)

### Function: Simulate the data
### Note: the cluster index must start from 0
### Note: we have one constraint: K <= J_imp.
data_sim <- function(N, K, J, J_imp, z_lb, z_ub, z_lb_ova, z_ub_ova, 
                     pi_gk, pi_g_ova){
  
  ### Cluster Assignment
  ci_actual <- sample(1:K, N, replace = TRUE)
  
  ### Important Variables (wj = 1 for the first J_imp variables)
  ### Note: K <= J_imp
  gm_base <- cbind(diag(J_imp)[1:K, ])
  xi_base <- gm_base
  xi_base[xi_base == 1] <- 10
  xi_dum <- runif(sum(xi_base == 0))
  xi_base[xi_base == 0] <- xi_dum
  
  ### Unimportant variable (wj = 0 for the last J - J_imp variables)
  xi_ova <- matrix(runif((J - J_imp) * K), nrow = K, ncol = J - J_imp)
  
  zi_active <- matrix(0, nrow = N, ncol = J_imp)
  gm_active <- matrix(0, nrow = N, ncol = J_imp)
  zi_inactive <- matrix(0, nrow = N, ncol = J - J_imp)
  gm_inactive <- matrix(0, nrow = N, ncol = J - J_imp)
  
  for(i in 1:N){
    ### Active variables
    gm_i <- gm_base[ci_actual[i], ]
    gm_i[gm_i == 0] <- rbinom(sum(gm_i == 0), 1, pi_gk[ci_actual[i]])
    gm_active[i, ] <- gm_i
    gwx_i <- gm_i * xi_base[ci_actual[i], ]
    zi_active[i, ] <- rmultinom(1, round(runif(1, z_lb, z_ub)), gwx_i/sum(gwx_i))
    
    ### Inactive variables
    gm_inc <- rbinom(J - J_imp, 1, pi_g_ova)
    gm_inactive[i, ] <- gm_inc
    gwx_inc <- gm_inc * xi_ova[ci_actual[i], ]
    zi_inactive[i, ] <- rmultinom(1, round(runif(1, z_lb_ova, z_ub_ova)), gwx_inc/sum(gwx_inc))
  }
  
  list(gm = cbind(gm_active, gm_inactive), xi = cbind(xi_base, xi_ova), 
       z = cbind(zi_active, zi_inactive), ci = ci_actual - 1)  
  
}
library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)
library(latex2exp)
library(salso)

### Import the external function
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(1291)
sim_list <- data_sim(n = 100, K = 5, J_imp = 5, 
                     pi_gm_mat = matrix(c(1), ncol = 10, nrow = 5),
                     xi_scale = 10, sum_zi = 500)
sim_tau <- c(rgamma(5, 1, 1))

### Function for adjusting tau and beta
adjust_tau_beta(beta_mat = matrix(1:10, ncol = 5, nrow = 10), 
                tau_vec = rep(1, 10), clus_assign = 0:9)

adjust_tau_beta(beta_mat = matrix(1:10, ncol = 5, nrow = 10), 
                tau_vec = rep(1, 10), clus_assign = c(1))

adjust_tau_beta(beta_mat = matrix(1:10, ncol = 5, nrow = 10), 
                tau_vec = rep(1, 10), clus_assign = c(1, 4, 2, 8))

### Function: Reallocation
tt <- realloc(z = sim_list$z, clus_assign = rep(0:1, 50), 
              gamma_mat = sim_list$gamma, w = c(rep(1, 4), rep(0, 6)), 
              beta_mat = matrix(1, nrow = 5, ncol = 10), 
              tau_vec = sim_tau, theta_vec = rep(1, 5))

### Function: Split-Merge
tt1 <- sm(K_max = 5, z = sim_list$z, clus_assign = sim_list$ci,
         gamma_mat = sim_list$gamma, w = c(rep(1, 5), rep(0, 5)), 
         beta_mat = sim_list$beta, tau_vec = sim_tau, theta_vec = rep(1, 5), 
         launch_iter = 10, s2 = 1, r0c = 1, r1c = 1)
table(sim_list$ci, tt1$clus_assign)
tt1$split_ind
tt1$accept_MH
tt1$logA

set.seed(1)
bet_test <- t(matrix(c(rnorm(10), rep(0, 40)), nrow = 10))
tau_test <- c(rgamma(1, 1, 1), rep(0, 4))
tt2 <- sm(K_max = 5, z = sim_list$z, clus_assign = rep(0, 100),
          gamma_mat = sim_list$gamma, w = c(rep(1, 5), rep(0, 5)), 
          beta_mat = bet_test, tau_vec = tau_test, theta_vec = rep(1, 5), 
          launch_iter = 10, s2 = 1, r0c = 1, r1c = 1)
table(sim_list$ci, tt2$clus_assign)
c(tt2$split_ind, tt2$accept_MH, tt2$logA)

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

tau_vec <- rgamma(5, 1, 1)
U <- rgamma(1, 1, 1)
ttt <- update_tau_u(clus_assign = sim_list$ci, tau_vec = tau_vec, 
                    theta_vec = rep(1, 5), old_U = U)
cbind(tau_vec, ttt$new_tau)
cbind(U, ttt$new_U)

### Test: Full Function
start_time <- Sys.time()
test_result <- full_func(iter = 10000, K_max = 5, z = sim_list$z,
                         gamma_mat = sim_list$gamma, w = c(rep(1, 5), rep(0, 5)),
                         theta = rep(1, 5), launch_iter = 10,
                         MH_var = 0.01, s2 = 1, r0c = 1, r1c = 1)
Sys.time() - start_time
table(test_result$ci)
test_result$beta[, , 2000]
table(test_result$sm_result)
table(test_result$accept_result)
table(test_result$sm_result, test_result$accept_result)

table(salso(t(test_result$ci)[5001:10000, ]), sim_list$ci)

### Test with the same simulated dataset
rm(list = ls())
load("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/result_c1.Rdata")

KK <- 5
run_time <- Sys.time()
result <- full_func(iter = 10000, K_max = KK, z = result_list$data$z,
                    gamma_mat = result_list$data$gamma, w = c(rep(1, 5), rep(0, 2)),
                    theta = rep(1, KK), launch_iter = 10, MH_var = 0.01, 
                    s2 = 1, r0c = 1, r1c = 1)
Sys.time() - run_time
ci_result <- salso(t(result$ci)[-c(1:5000), ], maxNClusters = KK)
table(actual = result_list$data$ci, result = ci_result)

table(result$accept_result)

plot(apply(t(result$ci), 1, function(x){length(unique(x))}), type = "l")

KK <- 30
run_time <- Sys.time()
result <- full_func_1(iter = 10000, K_max = KK, z = result_list$data$z,
                    gamma_mat = result_list$data$gamma, w = c(rep(1, 5), rep(0, 2)),
                    theta = rep(1, KK), launch_iter = 10, MH_var = 0.01, 
                    s2 = 1, r0c = 1, r1c = 1)
Sys.time() - run_time
ci_result <- salso(t(result$ci)[-c(1:5000), ], maxNClusters = KK)
table(actual = result_list$data$ci, result = ci_result)

plot(apply(t(result$ci), 1, function(x){length(unique(x))}), type = "l")

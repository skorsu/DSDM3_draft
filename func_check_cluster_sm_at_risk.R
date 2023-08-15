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
set.seed(21356)
sim_list <- data_sim(n = 100, K = 5, J_imp = 10, 
                     pi_gm_mat = matrix(c(0.25), ncol = 10, nrow = 5),
                     xi_scale = 10, sum_zi = 1500)
sim_list$z
sim_list$gamma

sum(sim_list$z == 0)

### Test: Full Function
KK <- 15
start_time <- Sys.time()
test_result <- full_func_at_risk(iter = 50000, K_max = KK, z = sim_list$z,
                                 w = rep(1, 10), theta = rep(1, KK), launch_iter = 10,
                                 MH_var = 0.01, s2 = 1, r0g = 3, r1g = 1, 
                                 r0c = 1, r1c = 1)
Sys.time() - start_time
table(test_result$sm_result)
table(test_result$accept_result)
table(test_result$sm_result, test_result$accept_result)

plot(apply(test_result$ci, 2, function(x){length(unique(x))}))

table(salso(t(test_result$ci)[30001:50000, ]), sim_list$ci)

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

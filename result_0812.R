rm(list = ls())

### Required Libraries
library(salso)
library(foreach)
library(doParallel)
library(mclustcomp)

### Import the function for simulated data
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")

### Scenario 1
# sim_list <- data_sim(n = 500, K = 3, J_imp = 7, 
#                      pi_gm_mat = matrix(1, nrow = 3, ncol = 15),
#                      xi_scale = 10, sum_zi = 1000)

### Scenario 2
# sim_list <- data_sim(n = 500, K = 5, J_imp = 5, 
#                      pi_gm_mat = matrix(1, nrow = 5, ncol = 10),
#                      xi_scale = 10, sum_zi = 500)

### Hyperparameters
set.seed(124, kind = "L'Ecuyer-CMRG")
KK <- 10

registerDoParallel(5)

run_ova <- Sys.time()
result_para <- foreach(t = 1:10) %dopar% {
  
  ### Simulate the data
  sim_list <- data_sim(n = 500, K = 3, J_imp = 7, 
                       pi_gm_mat = matrix(1, nrow = 3, ncol = 15),
                       xi_scale = 10, sum_zi = 1000)
  
  ### Run the model
  runtime_1 <- Sys.time()
  result_1 <- full_func(iter = 10000, K_max = KK, z = sim_list$z,
                        gamma_mat = sim_list$gamma, w = c(rep(1, 5), rep(0, 5)),
                        theta = rep(1, 10), launch_iter = 10,
                        MH_var = 0.01, s2 = 1, r0c = 1, r1c = 1)
  runtime_1 <- as.numeric(Sys.time() - runtime_1)
  
  runtime_all <- Sys.time()
  result_all <- full_func_1(iter = 10000, K_max = KK, z = sim_list$z,
                            gamma_mat = sim_list$gamma, w = c(rep(1, 5), rep(0, 5)),
                            theta = rep(1, 10), launch_iter = 10,
                            MH_var = 0.01, s2 = 1, r0c = 1, r1c = 1)
  runtime_all <- as.numeric(Sys.time() - runtime_all)
  
  ### Result
  ci_mcmc_1 <- as.numeric(salso(t(result_1$ci)[-(1:5000), ], maxNClusters = KK))
  ci_mcmc_all <- as.numeric(salso(t(result_all$ci)[-(1:5000), ], maxNClusters = KK))
  result_list <- list(n_init_all = result_all$init_clus, actual_ci = sim_list$ci,
                      ci_mcmc_1 = ci_mcmc_1, ci_mcmc_all = ci_mcmc_all,
                      n_clus_1 = apply(result_1$ci, 2, function(x){length(unique(x))}),
                      n_clus_all = apply(result_all$ci, 2, function(x){length(unique(x))}),
                      comp_1 = mclustcomp(sim_list$ci, ci_mcmc_1)[c(1, 5, 22), 2],
                      comp_all = mclustcomp(sim_list$ci, ci_mcmc_all)[c(1, 5, 22), 2],
                      runtime_1 = runtime_1, runtime_all = runtime_all)

  return(result_list)
  
}
Sys.time() - run_ova

stopImplicitCluster()

### Report the result
mcomp <- matrix(NA, ncol = 6, nrow = 10)
rtime <- matrix(NA, ncol = 2, nrow = 10)
for(i in 1:10){
  mcomp[i, ] <- c(result_para[[i]]$comp_1, result_para[[i]]$comp_all)
  rtime[i, ] <- c(result_para[[i]]$runtime_1, result_para[[i]]$runtime_all)
}

apply(mcomp, 2, mean)
apply(mcomp, 2, sd)
apply(rtime, 2, mean)
apply(rtime, 2, sd)


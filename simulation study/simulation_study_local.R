###: ---------------------------------------------------------------------------
###: Change the settings here
n_core <- 5
dat_rep <- 30

base_path <- "/Users/kevin-imac/Desktop/"
case_name <- "weak_sig_with_overlap_zi"

#### Data Hyperparameters
xi_c <- 100
pig <- 0.65

#### Model Hyperparameters
MH_v <- 1
burn_in <- 5000

###: ---------------------------------------------------------------------------

### Required Libraries
library(foreach)
library(doParallel)
library(sparseMbClust)

### Import the function for simulating the data
source(paste0(base_path, "Github - Repo/ClusterZI/data/data_sim_DM.R"))

### Simulate the data
registerDoParallel(n_core)
data_sim <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  pat_mat <- matrix(0, nrow = 3, ncol = 20)
  pat_mat[1, 1:10] <- 1
  pat_mat[2, 8:15] <- 1
  pat_mat[3, 15:20] <- 1
  sim_dat <- simDM(n = 200, pattern = pat_mat, xi_conc = xi_c, pi_gm = c(pig), 
                   pi_c = c(1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                   theta = 0.01)
  sim_dat_sum <- simDM_sum(sim_dat)$summary_zero
  return(list(data = sim_dat, zero = sim_dat_sum))
}
stopImplicitCluster()


### Run the models: ------------------------------------------------------------

start_overall <- Sys.time()

### ZIDM-ZIDM
registerDoParallel(n_core)
zz_result <- foreach(t = 1:dat_rep) %dopar% {
  
  set.seed((10 * t) + 21)
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  zz_mod <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$data$z, 
                      theta_vec = rep(1, 10), launch_iter = 5, MH_var = MH_v, 
                      mu = 0, s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                      print_iter = 2500)
  zz_time <- difftime(Sys.time(), start_time, units = "secs")
  
  return(list(zz_mod = zz_mod, zz_time = zz_time))
  
}
stopImplicitCluster()

### DM-ZIDM
registerDoParallel(n_core)
dz_result <- foreach(t = 1:dat_rep) %dopar% {
  
  set.seed((10 * t) + 21)
  
  ### DM-ZIDM
  start_time <- Sys.time()
  dz_mod <- DM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$data$z, 
                    theta_vec = rep(1, 10), launch_iter = 5, MH_var = MH_v, mu = 0, 
                    s2 = 1, r0c = 1, r1c = 1, print_iter = 2500)
  dz_time <- difftime(Sys.time(), start_time, units = "secs")
  

  return(list(dz_mod = dz_mod, dz_time = dz_time))
  
}
stopImplicitCluster()

### Shi Model
registerDoParallel(n_core)
shi_result <- foreach(t = 1:dat_rep) %dopar% {
  
  start_time <- Sys.time()
  shi_mod <- PerformClustering(t(data_sim[[t]]$data$z), "DP", w = 1, 
                               beta1 = 1, beta2 = 1, totaliter = 10000, 
                               burnin = 5000, thin = 1)
  shi_time <- difftime(Sys.time(), start_time, units = "secs")

  return(list(shi_mod = shi_mod, shi_time = shi_time))
  
}
stopImplicitCluster()

### DM-DM
registerDoParallel(n_core)
dd_result <- foreach(t = 1:dat_rep) %:%
  foreach(k = 2:10) %dopar% {
    set.seed((k * t) + 12)
    start_time <- Sys.time()
    dd_mod <- DM_DM(iter = 10000, K_max = k, z = data_sim[[t]]$data$z,
                    theta_vec = rep(1, k), MH_var = MH_v, mu = 0, s2 = 1,
                    print_iter = 2500)
    dd_time <- difftime(Sys.time(), start_time, units = "secs")
    return(list(t = t, k = k, dd_mod = dd_mod, dd_time = dd_time))
  }
stopImplicitCluster()

print(Sys.time() - start_overall)

### Save Result: ---------------------------------------------------------------
#### Data
saveRDS(list(dat = data_sim), paste0(base_path, "simu_result/data_", case_name, ".RData"))
### readRDS(paste0(base_path, "simu_result/data_", case_name, ".RData"))
#### Result
saveRDS(list(zz_result = zz_result, dz_result = dz_result,
             dd_result = dd_result, shi_result = shi_result), 
        paste0(base_path, "simu_result/result_", case_name, ".RData"))

### Fix Shi's Model: -----------------------------------------------------------
#### First, Import the data
base_path <- "/Users/kevin-imac/Desktop/"
case_name <- "strong_sig_with_overlap_zi"
n_core = 5
dat_rep = 30
dat <- readRDS(paste0(base_path, "simu_result/data_", case_name, ".RData"))

start_shi_fix <- Sys.time()
registerDoParallel(n_core)
shi_result <- foreach(t = 1:dat_rep) %dopar% {
  
  start_time <- Sys.time()
  shi_mod <- PerformClustering(t(dat$dat[[t]]$data$z), "DP", w = 1, 
                               beta1 = 1, beta2 = 1, totaliter = 10000, 
                               burnin = 5000, thin = 1)
  shi_time <- difftime(Sys.time(), start_time, units = "secs")
  
  return(list(shi_mod = shi_mod, shi_time = shi_time))
  
}
stopImplicitCluster()
Sys.time() - start_shi_fix

saveRDS(list(fix_shi_result = shi_result), 
        paste0(base_path, "simu_result/fix_shi_result_", case_name, ".RData"))

rm(list = ls())



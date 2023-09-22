###: Change the settings here
n_core <- 5
dat_rep <- 2

base_path <- "/Users/kevin-imac/Desktop/"
case_name <- "case3a_zi"

#### Data Hyperparameters
xi_c <- 100
pig <- 0.85

#### Model Hyperparameters
MH_v <- 1

#### Check
##### () Normal: pi_gm = 0.95 and x_conc = 10
##### () ZI: pi_gm = 0.85 and xi_conc = 100

###: ---------------------------------------------------------------------------

### Required Libraries
library(foreach)
library(doParallel)
library(sparseMbClust)

### Import the function for simulating the data
source(paste0(base_path, "Github - Repo/ClusterZI/data/data_sim_DM.R"))

start_overall <- Sys.time()

### Simulate the data
registerDoParallel(n_core)
data_sim <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  pat_mat <- matrix(0, nrow = 3, ncol = 30)
  pat_mat[1, c(1, 3:5, 7:10)] <- 1
  pat_mat[2, c(8:10, 12, 14:15)] <- 1
  pat_mat[3, c(15:16, 18, 20)] <- 1
  sim_dat <- simDM(n = 200, pattern = pat_mat, xi_conc = xi_c, pi_gm = c(pig), 
                   pi_c = c(1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                   theta = 0.01)
  return(sim_dat)
}
stopImplicitCluster()

### ZIDM-ZIDM
zz_overall <- Sys.time()
registerDoParallel(n_core)
zz_model <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  result <- ZIDM_ZIDM(iter = 20000, K_max = 10, z = data_sim[[t]]$z, 
                      theta_vec = rep(1, 10), launch_iter = 5, MH_var = MH_v, 
                      mu = 0, s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                      print_iter = 2500)
  time <- difftime(Sys.time(), start_time, "sec")
  return(list(zz_assign = result$assign, zz_time = time))
}
stopImplicitCluster()
print(paste0("ZIDM-ZIDM: ", difftime(Sys.time(), zz_overall, "sec")))

### DM-ZIDM
dz_overall <- Sys.time()
registerDoParallel(n_core)
dz_model <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  result <- DM_ZIDM(iter = 20000, K_max = 10, z = data_sim[[t]]$z, 
                    theta_vec = rep(1, 10), launch_iter = 5, MH_var = MH_v, 
                    mu = 0, s2 = 1, r0c = 1, r1c = 1, print_iter = 2500)
  time <- difftime(Sys.time(), start_time, "sec")
  return(list(dz_assign = result$assign, dz_time = time))
}
stopImplicitCluster()
print(paste0("DM-ZIDM: ", difftime(Sys.time(), dz_overall, "sec")))

### DM-DM
dd_overall <- Sys.time()
registerDoParallel(n_core)
dd_model <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  result <- DM_DM(iter = 20000, K_max = 5, z = data_sim[[t]]$z,
                  theta_vec = rep(1, 5), MH_var = MH_v, mu = 0, s2 = 1,
                  print_iter = 2500)
  time <- difftime(Sys.time(), start_time, "sec")
  return(list(dd_assign = result, dd_time = time))
}
stopImplicitCluster()
print(paste0("DM-DM: ", difftime(Sys.time(), dd_overall, "sec")))

### DM-sDM
dsd_overall <- Sys.time()
registerDoParallel(n_core)
dsd_model <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  result <- DM_DM(iter = 20000, K_max = 10, z = data_sim[[t]]$z,
                  theta_vec = rep(0.001, 5), MH_var = MH_v, mu = 0, s2 = 1,
                  print_iter = 2500)
  time <- difftime(Sys.time(), start_time, "sec")
  return(list(dsd_assign = result, dsd_time = time))
}
stopImplicitCluster()
print(paste0("DM-sDM: ", difftime(Sys.time(), dsd_overall, "sec")))

### Shi's (DP)
shi_dp_overall <- Sys.time()
registerDoParallel(n_core)
shi_dp_model <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  result <- PerformClustering(t(data_sim[[t]]$z), "DP", w = 1, 
                              beta1 = 1, beta2 = 1, totaliter = 20000, 
                              burnin = 10000, thin = 1)
  time <- difftime(Sys.time(), start_time, "sec")
  return(list(shi_dp_assign = result$crec, shi_dp_time = time))
}
stopImplicitCluster()
print(paste0("Shi (DP): ", difftime(Sys.time(), shi_dp_overall, "sec")))

### Shi's (MFM)
shi_mfm_overall <- Sys.time()
registerDoParallel(n_core)
shi_mfm_model <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  result <- PerformClustering(t(data_sim[[t]]$z), "MFM", w = 1, 
                              totaliter = 20000, burnin = 10000, thin = 1)
  time <- difftime(Sys.time(), start_time, "sec")
  return(list(shi_mfm_assign = result$crec, shi_mfm_time = time))
}
stopImplicitCluster()
print(paste0("Shi (MFM): ", difftime(Sys.time(), shi_mfm_overall, "sec")))

mod <- list(zz_model = zz_model, dz_model = dz_model, dd_model = dd_model, 
            dsd_model = dsd_model, shi_dp_model = shi_dp_model, 
            shi_mfm_model = shi_mfm_model)

### Save the data and result
saveRDS(data_sim, paste0(base_path, "data_", case_name, ".RData"))
saveRDS(mod, paste0(base_path, "mod_", case_name, ".RData"))




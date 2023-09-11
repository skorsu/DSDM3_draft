###: ---------------------------------------------------------------------------
###: Change the settings here
n_core <- 5
dat_rep <- 2

base_path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
case_name <- "weak_sig_no_overlap_zi"

#### Data Hyperparameters
xi_c <- 100
pig <- 0.85

#### Model Hyperparameters
MH_v <- 10
burn_in <- 5000

#### Check
##### () Normal: pi_gm = 1 and x_conc = 10
##### () ZI: pi_gm = 0.65 and xi_conc = 100

###: ---------------------------------------------------------------------------

### Required Libraries
library(ClusterZI)
library(salso)
library(foreach)
library(doParallel)

### Import the function for simulating the data
source(paste0(base_path, "data/data_sim_DM.R"))

start_overall <- Sys.time()

### Simulate the data
registerDoParallel(n_core)
data_sim <- foreach(t = 1:dat_rep) %dopar% {
  set.seed(t)
  pat_mat <- matrix(0, nrow = 3, ncol = 20)
  pat_mat[1, 1:7] <- 1
  pat_mat[2, 8:14] <- 1
  pat_mat[3, 15:20] <- 1
  sim_dat <- simDM(n = 200, pattern = pat_mat, xi_conc = xi_c, pi_gm = c(pig), 
                   pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                   theta = 0.01)
  sim_dat_sum <- simDM_sum(sim_dat)$summary_zero
  return(list(data = sim_dat, zero = sim_dat_sum))
}
stopImplicitCluster()

print("Perform x-ZIDM")

### ZIDM-ZIDM, DM-ZIDM
ss <- Sys.time()
registerDoParallel(n_core)
shi_zidm <- foreach(t = 1:dat_rep) %dopar% {
  
  set.seed((10 * t) + 21)
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  zz_mod <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$data$z, 
                      theta_vec = rep(1, 10), launch_iter = 5, MH_var = MH_v, 
                      mu = 0, s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                      print_iter = 2500)
  zz_time <- difftime(Sys.time(), start_time, units = "secs")
  
  ### DM-ZIDM
  start_time <- Sys.time()
  dz_mod <- DM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$data$z, 
                    theta_vec = rep(1, 10), launch_iter = 5, MH_var = MH_v, mu = 0, 
                    s2 = 1, r0c = 1, r1c = 1, print_iter = 2500)
  dz_time <- difftime(Sys.time(), start_time, units = "secs")
  
  ### Shi (Error due to cannot install in goose)
  ### start_time <- Sys.time()
  ### shi_mod <- PerformClustering(t(data_sim[[t]]$data$z), "DP", w = 1, 
  ###                              beta1 = 1000, beta2 = 1, totaliter = 10000, 
  ###                              burnin = 5000, thin = 1)
  ### shi_time <- difftime(Sys.time(), start_time, units = "secs")
  ### shi_assign <- as.numeric(salso(shi_mod$crec, maxNClusters = 10))

  return(list(zz_mod = zz_mod, dz_mod = dz_mod, 
              zz_time = zz_time, dz_time = dz_time))
  
}
stopImplicitCluster()

### DM-DM
registerDoParallel(n_core)
dm_dm_result <- foreach(t = 1:dat_rep) %:%
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

### Analyze the result by using salso

#### Actual Cluster Assignment and Zero Count
registerDoParallel(n_core)
actual_ci <- foreach(t = 1:dat_rep) %dopar% {
  list(actual_ci = data_sim[[t]]$data$ci, 
       zero = data_sim[[t]]$zero)
}
stopImplicitCluster()

#### Assign: ZIDM-ZIDM
registerDoParallel(n_core)
zz_salso <- foreach(t = 1:dat_rep, .combine = "cbind") %dopar% {
  as.numeric(salso(shi_zidm[[t]]$zz_mod$assign[-(1:burn_in), ]))
}
stopImplicitCluster()

#### Assign: DM-ZIDM
registerDoParallel(n_core)
dz_salso <- foreach(t = 1:dat_rep, .combine = "cbind") %dopar% {
  as.numeric(salso(shi_zidm[[t]]$dz_mod$assign[-(1:burn_in), ]))
}
stopImplicitCluster()

#### Assign: DM-DM
registerDoParallel(n_core)
dd_salso <- foreach(t = 1:dat_rep) %:%
  foreach(k = 1:9, .combine = "cbind") %dopar% {
    salso(dm_dm_result[[t]][[k]]$dd_mod[-(1:burn_in), ], maxNClusters = (k + 1))
  }
stopImplicitCluster()

#### Time
registerDoParallel(n_core)
time_mat <- foreach(t = 1:dat_rep, .combine = "rbind") %dopar% {
  c(shi_zidm[[t]]$zz_time, shi_zidm[[t]]$dz_time)
}
stopImplicitCluster()

colnames(time_mat) <- c("ZIDM-ZIDM", "DM-ZIDM")

registerDoParallel(n_core)
time_dd <- foreach(t = 1:dat_rep, .combine = "rbind") %:%
  foreach(k = 1:9, .combine = "cbind") %dopar% {
    dm_dm_result[[t]][[k]]$dd_time
  }
stopImplicitCluster()

colnames(time_dd) <- paste0("DM-DM", 2:10)
time_mat <- cbind(time_mat, time_dd)

### Save the data
saveRDS(data_sim, file = paste0(base_path, "data/", case_name, ".RData"))
### readRDS(paste0(base_path, "data/", case_name, ".RData"))

list_model <- list(shi_zidm = shi_zidm, dm_dm_result = dm_dm_result)
saveRDS(list_model, paste0(base_path, "model/", case_name, ".RData"))

list_result <- list(actual_ci = actual_ci, zz_salso = zz_salso, 
                    dz_salso = dz_salso, dd_salso = dd_salso, time_mat = time_mat)
saveRDS(list_result, paste0(base_path, "result/", case_name, ".RData"))
### readRDS(paste0(base_path, "result/", case_name, ".RData"))

print(Sys.time() - start_overall)

print("END")

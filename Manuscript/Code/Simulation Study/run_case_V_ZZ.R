library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_V.rds"))

### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
foreach(t = 1:20) %dopar% {
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 25, 
                      z = as.matrix(dat[[t]]$dat), atrisk_init = matrix(1, nrow = 50, ncol = 250), 
                      beta_init = matrix(0, nrow = 10, ncol = 250), 
                      ci_init = rep(0, 50), theta = 1, mu = 0, s2 = 0.1, 
                      s2_MH = 1e-3, t_thres = 1000, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 9, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_V_chain_", t, "_ZZ.rds"))
  
}
stopImplicitCluster()

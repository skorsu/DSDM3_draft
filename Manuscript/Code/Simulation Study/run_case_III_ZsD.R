library(foreach)
library(doParallel)
library(DTMM)
library(ape)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_III.rds"))

### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(3)
foreach(t = 1:20) %dopar% {
  
  ### ZIDM-sDM
  start_time <- Sys.time()
  mod <- ZIDM_DM(iter = 5000, Kmax = 10, nbeta_split = 20, 
                 z = as.matrix(dat[[t]]$dat), atrisk_init = matrix(1, nrow = 50, ncol = 100), 
                 beta_init = matrix(0, nrow = 10, ncol = 100), 
                 ci_init = rep(0, 50), theta = 1e-2, mu = 0, s2 = 0.1, 
                 s2_MH = 1e-3, t_thres = 1000, r0g = 1, r1g = 1, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_III_chain_", t, "_ZsD.rds"))
  
}
stopImplicitCluster()

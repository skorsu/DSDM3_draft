library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_III.rds"))

### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(2)
foreach(t = 1:20) %dopar% {
  
  ### DM-DM
  start_time <- Sys.time()
  mod <- DM_DM(iter = 5000, Kmax = 10, z = as.matrix(dat[[t]]$dat),
               beta_init = matrix(0, nrow = 10, ncol = 100), 
               ci_init = rep(0, 50), theta = 1e-5, mu = 0, 
               s2 = 0.1, s2_MH = 1e-3, t_thre = 1000, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_III_chain_", t, "_DsD.rds"))
  
}
stopImplicitCluster()
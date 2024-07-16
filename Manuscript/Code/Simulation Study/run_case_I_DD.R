library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_I.rds"))

result <- foreach(t = 1:20) %:%
  foreach(k = 2:9) %dopar% {
    
    ### DM-DM
    start_time <- Sys.time()
    mod <- DMDM(iter = 2500, Kmax = k, z = as.matrix(dat[[t]]$dat),
                beta_mat = matrix(0, nrow = k, ncol = 250), 
                ci_init = rep(0, 50), theta = 1, thin = 1)
    tot_time <- difftime(Sys.time(), start_time, units = "secs")
    list(mod = mod, time = tot_time) 
            
}

saveRDS(result, 
        file = paste0(path, "Result/Simulation Study/simu_data_case_I_DD.rds"))

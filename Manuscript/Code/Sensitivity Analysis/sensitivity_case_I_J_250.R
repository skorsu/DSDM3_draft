library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_I.rds"))

### Set of the hyperparameter: -------------------------------------------------
hyperParam <- list(c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 20, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 0.1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 0.1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 10, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-5, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 5, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 20, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 0.1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 10, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 0.1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 10, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 4, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 9, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 4, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 9, t_thres = 5000),
                   c(Kmax = 5, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 20, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 2500),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1e-3, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 7500))

### Run all models and save the result -----------------------------------------
set.seed(1, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(3)
result <- foreach(h = 1:length(hyperParam)) %:% ### Each set of the hyperparameter
  foreach(r = 1:20) %dopar% { ### The number of replicated data
    start_time <- Sys.time()
    mod <- mod_adaptive(iter = 10000, Kmax = hyperParam[[h]]["Kmax"], 
                        nbeta_split = hyperParam[[h]]["nbeta_split"], 
                                z = dat[[r]]$dat, 
                                atrisk_init = matrix(1, ncol = 250, nrow = 50), 
                                beta_init = matrix(0, ncol = 250, nrow = hyperParam[[h]]["Kmax"]), 
                                ci_init = rep(0, 50), theta = hyperParam[[h]]["theta"], 
                                mu = 0, s2 = hyperParam[[h]]["s2"], 
                                s2_MH = hyperParam[[h]]["s2_MH"], 
                                t_thres = hyperParam[[h]]["t_thres"],
                                launch_iter = hyperParam[[h]]["launch_iter"], 
                                r0g = hyperParam[[h]]["r0g"], 
                                r1g = hyperParam[[h]]["r1g"], 
                                r0c = hyperParam[[h]]["r0c"], 
                                r1c = hyperParam[[h]]["r1c"], 
                                thin = 1)
    tot_time <- difftime(Sys.time(), start_time, units = "secs")
    saveRDS(list(mod = mod, time = tot_time), 
            file = paste0(path, "Result/Sensitivity Analysis/simu_data_case_I_chain_", r, "_hyperparamset_", h, ".rds"))
    
    list(time = tot_time, result = clus_result$ci_result)
  }
stopImplicitCluster()
difftime(Sys.time(), start_ova)
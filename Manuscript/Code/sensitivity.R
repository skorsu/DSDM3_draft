library(Rcpp)
library(salso)
library(foreach)
library(doParallel)

### Directory path: ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

sourceCpp(paste0(path, "src/clusterZI.cpp"))
caseName <- "diffindex_3_K_5_simDat"
dat <- readRDS(paste0(path, "Manuscript/Data/", caseName, ".RData"))

### Set of the hyperparameter: -------------------------------------------------
hyperParam <- list(c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 20, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 0.1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 0.1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 10, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 0.1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 5, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 20, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 0.1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 10, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 0.1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 10, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 4, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 9, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 4, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 9, t_thres = 5000),
                   c(Kmax = 5, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 20, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 5000),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 2500),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, t_thres = 7500))

### Run all models and save the result -----------------------------------------
set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
result <- foreach(h = 1:length(hyperParam)) %:% ### Each set of the hyperparameter
  foreach(r = 1:length(dat)) %dopar% { ### The number of replicated data
    start_time <- Sys.time()
    clus_result <- mod_adaptive(iter = 10000, Kmax = hyperParam[[h]]["Kmax"], 
                                nbeta_split = hyperParam[[h]]["nbeta_split"], 
                                z = dat[[r]]$dat, 
                                atrisk_init = matrix(1, ncol = 50, nrow = 100), 
                                beta_init = matrix(0, ncol = 50, nrow = hyperParam[[h]]["Kmax"]), 
                                ci_init = rep(0, 100), theta = hyperParam[[h]]["theta"], 
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
    list(time = tot_time, result = clus_result$ci_result)
  }
stopImplicitCluster()
difftime(Sys.time(), start_ova)

### Save the result: -----------------------------------------------------------
for(i in 1:nHyperSet){
  
  file_case <- paste(paste0(names(hyperParam[[i]]), "_", hyperParam[[i]]), 
                     collapse = "_")
  save_result <- list(hyper = hyperParam[[i]], 
                      result = result[[i]])
  saveRDS(save_result, 
          paste0(path, caseName, "_", file_case, "_MB_Adaptive.RData"))
  
}

### ----------------------------------------------------------------------------

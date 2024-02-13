### Load libraries
library(Rcpp)
library(foreach)
library(doParallel)

### Change the settings --------------------------------------------------------
case_name <- "diffindex_3_K_2"
nData <- 20

### Import the data ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

dat <- readRDS(paste0(path, "ClusterZI/Manuscript/Data/", case_name, "_simDat.RData"))
dat_J <- ncol(dat[[1]]$dat) 
dat_n <- nrow(dat[[1]]$dat) 

save_path <- paste0(path, "ClusterZI/Manuscript/Result/Simulation Study/")
sourceCpp(paste0(path, "ClusterZI/src/clusterZI.cpp"))

### Run the models -------------------------------------------------------------
#### ZIDM-ZIDM
set.seed(3124, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
resultZZ <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- mod(iter = 100000, Kmax = 10, nbeta_split = 5, z = dat[[t]]$dat, 
                     atrisk_init = matrix(1, ncol = dat_J, nrow = dat_n), 
                     beta_init = matrix(0, ncol = dat_J, nrow = 10), 
                     ci_init = rep(0, dat_n), theta = 1, mu = 0, s2 = 1, s2_MH = 1, 
                     launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                     thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(resultZZ, paste0(save_path, case_name, "_ZZ.RData"))

#### DM-ZIDM
start_ova <- Sys.time()
registerDoParallel(5)
resultDZ <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- DMZIDM(iter = 100000, Kmax = 10, z = dat[[t]]$dat, 
                        beta_mat = matrix(0, ncol = dat_J, nrow = 10), 
                        ci_init = rep(0, dat_n), 
                        theta = 1, launch_iter = 10, r0c = 1, r1c = 1, 
                        thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDZ, paste0(save_path, case_name, "_DZ.RData"))

#### DM-DM
start_ova <- Sys.time()
registerDoParallel(5)
resultDD <- foreach(t = 1:nData) %:%
  foreach(k = 2:10) %dopar% {
    start_time <- Sys.time()
    clus_result <- DMDM(iter = 100000, Kmax = k, z = dat[[t]]$dat, 
                        beta_mat = matrix(0, ncol = dat_J, nrow = k), 
                        ci_init = rep(0, dat_n), 
                        theta = 1, thin = 100)
    tot_time <- difftime(Sys.time(), start_time, units = "secs")
    list(time = tot_time, result = clus_result)
  }
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDD, paste0(save_path, case_name, "_DD.RData"))

#### DM-sDM
start_ova <- Sys.time()
registerDoParallel(5)
resultDsD <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- DMDM(iter = 100000, Kmax = 10, z = dat[[t]]$dat, 
                      beta_mat = matrix(0, ncol = dat_J, nrow = 10), 
                      ci_init = rep(0, dat_n), 
                      theta = 1e-10, thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDsD, paste0(save_path, case_name, "_DsD.RData"))

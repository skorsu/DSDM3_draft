### Setting
nCores <- 5
path_data <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/result_1224/"
path_result <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/result_1224/"
case_name <- "diffindex_3"

### Required Library
library(Rcpp)
library(sparseMbClust)
library(foreach)
library(doParallel)

### Analysis
### Import the data
dat <- readRDS(paste0(path_data, case_name, "_simDat.RData"))
n <- length(dat$dat)

### Run the model with DP
start_ova <- Sys.time()
registerDoParallel(nCores)
resultDP <- foreach(t = 1:n) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- PerformClustering(otutable = t(dat$dat[[t]]), ClusteringMethod = "DP",
                                   alpha1 = 1, alpha2 = 1, nu = 1,
                                   w = 1, beta1 = 1, beta2 = 1, a = 1, b = 1, pargamma = 1, lambda = 1,
                                   totaliter = 100000, burnin = 0, thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result$crec)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDP, paste0(path_result, case_name, "_DP_goose.RData"))

### Run the model with MFM
start_ova <- Sys.time()
registerDoParallel(nCores)
resultMFM <- foreach(t = 1:n) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- PerformClustering(otutable = t(dat$dat[[t]]), ClusteringMethod = "MFM",
                                   alpha1 = 1, alpha2 = 1, nu = 1,
                                   w = 1, beta1 = 1, beta2 = 1, a = 1, b = 1, pargamma = 1, lambda = 1,
                                   totaliter = 100000, burnin = 0, thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result$crec)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultMFM, paste0(path_result, case_name, "_MFM_goose.RData"))
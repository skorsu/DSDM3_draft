library(DTMM)
library(ape)
library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_II.rds"))

### Adjusted the DTMM function 
adjDTMM <- function(Y, tree, tau_vec = 10 ^ seq(-1, 4, 0.5), nu_vec = 1, 
                    theta_vec = seq(0.01, 0.99, 0.08), init_c, 
                    beta = "default", mcmc_iter = 2500){
  
  if(beta == "default"){
    temp = runif(1)
    alpha = temp/(1 - temp)
  }
  
  tree <- reorder(tree, order = "postorder")
  edge <- apply(tree$edge, 2, rev)
  p = dim(Y)[2]
  gamma_sample = rbinom(p - 1, 1, 0.5)
  
  HDTMcpp(Y, edge, tau_vec, nu_vec, theta_vec,
          init_c, gamma_sample, alpha, mcmc_iter, select = FALSE)
}

#### DTMM: with no structures
set.seed(1, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(3)
resultDTMM <- foreach(t = 1:20) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- adjDTMM(dat[[t]]$dat, 
                         tree = read.tree(text = paste("(", paste(1:250, collapse = ", "), ");")), 
                         init_c = rep(1, 100), mcmc_iter = 10000)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = t(clus_result$post_c))
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_II_chain_", t, "_DMDP.rds"))
  
}
stopImplicitCluster()



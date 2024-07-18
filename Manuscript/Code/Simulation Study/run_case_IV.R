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
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_IV.rds"))

### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(3)
foreach(t = 1:20) %dopar% {
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 30, 
                      z = as.matrix(dat[[t]]$dat), atrisk_init = matrix(1, nrow = 150, ncol = 150), 
                      beta_init = matrix(0, nrow = 10, ncol = 150), 
                      ci_init = rep(0, 150), theta = 1, mu = 0, s2 = 0.1, 
                      s2_MH = 1e-3, t_thres = 1000, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 4, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_IV_chain_", t, "_ZZ.rds"))
  
  ### DM-ZIDM
  start_time <- Sys.time()
  mod <- DMZIDM(iter = 5000, Kmax = 10, z = as.matrix(dat[[t]]$dat), 
                beta_mat = matrix(0, nrow = 10, ncol = 150), 
                ci_init = rep(0, 150), theta = 1, launch_iter = 30, 
                r0c = 1, r1c = 4, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_IV_chain_", t, "_DZ.rds"))
  
  ### DM-sDM
  start_time <- Sys.time()
  mod <- DMDM(iter = 5000, Kmax = 10, z = as.matrix(dat[[t]]$dat), 
              beta_mat = matrix(0, nrow = 10, ncol = 150), 
              ci_init = rep(0, 150), theta = 1e-10, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_IV_chain_", t, "_DsD.rds"))
  
}
stopImplicitCluster()

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
                         tree = read.tree(text = paste("(", paste(1:150, collapse = ", "), ");")), 
                         init_c = rep(1, 150), mcmc_iter = 5000)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = t(clus_result$post_c))
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Simulation Study/simu_data_case_IV_chain_", t, "_DMDP.rds"))
  
}
stopImplicitCluster()

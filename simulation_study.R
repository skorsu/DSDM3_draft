library(salso)
library(mcclust)
library(gridExtra)
library(foreach)
library(doParallel)
library(sparseMbClust)

### Import the function
im_path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/"
mb_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/"
path <- NULL
if(dir.exists(mb_path)){
  path <- mb_path
} else {
  path <- im_path
}

source(paste0(path, "data_sim_DM.R"))

registerDoParallel(detectCores() - 3)
data_sim <- foreach(t = 1:30) %dopar% {
  set.seed(2*t)
  signal <- sample(1:20, 5)
  pat_mat <- diag(20)[signal, ]
  sim_dat <- simDM(n = 200, pattern = pat_mat, xi_conc = 100, pi_gm = c(0.85), 
                   pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                   theta = 0.01)
  sim_dat_sum <- simDM_sum(sim_dat)$summary_zero
  return(list(data = sim_dat, zero = sim_dat_sum))
}
stopImplicitCluster()

zero_sum_mat <- matrix(NA, ncol = 3, nrow = 30)
for(i in 1:30){
  zero_sum_mat[i, ] <- data_sim[[i]]$zero
}

quantile(zero_sum_mat[, 1], probs = c(0.025, 0.975)) ## Proportion of zero
quantile(zero_sum_mat[, 2], probs = c(0.025, 0.975)) ## P(at risk|zero)
quantile(zero_sum_mat[, 3], probs = c(0.025, 0.975)) ## P(structure|zero)

### ZIDM-ZIDM, DM-ZIDM, and Shi
ss <- Sys.time()
registerDoParallel(detectCores() - 3)
shi_zidm <- foreach(t = 1:2) %dopar% {
  set.seed(3*t)
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  zz_mod <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$data$z, 
                      theta_vec = rep(1, 10), launch_iter = 5, MH_var = 10, 
                      mu = 0, s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                      print_iter = 2500)
  zz_time <- difftime(Sys.time(), start_time, units = "secs")
  zz_assign <- as.numeric(salso(zz_mod$assign[-c(1:5000), ], maxNClusters = 10))
  
  ### DM-ZIDM
  start_time <- Sys.time()
  dz_mod <- DM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$data$z, 
                    theta_vec = rep(1, 10), launch_iter = 5, MH_var = 10, mu = 0, 
                    s2 = 1, r0c = 1, r1c = 1, print_iter = 2500)
  dz_time <- difftime(Sys.time(), start_time, units = "secs")
  dz_assign <- as.numeric(salso(dz_mod$assign[-c(1:5000), ], maxNClusters = 10))
  
  ### Shi
  start_time <- Sys.time()
  shi_mod <- PerformClustering(t(data_sim[[t]]$data$z), "DP", w = 1, 
                               beta1 = 1000, beta2 = 1, totaliter = 10000, 
                               burnin = 5000, thin = 1)
  shi_time <- difftime(Sys.time(), start_time, units = "secs")
  shi_assign <- as.numeric(salso(shi_mod$crec, maxNClusters = 10))

  return(list(zz_assign = zz_assign, dz_assign = dz_assign, shi_assign = shi_assign,
              zz_time = zz_time, dz_time = dz_time, shi_time = shi_mod))
  
}
stopImplicitCluster()
paste0("Next: Perform DM-DM ----- Spend ", 
       difftime(Sys.time(), ss), " mins already.")

### DM-DM
registerDoParallel(detectCores() - 3)
dm_dm_result <- foreach(t = 1:2) %:%
  foreach(k = 2:5) %do% {
    set.seed((k * t) + 12)
    start_time <- Sys.time()
    dd_mod <- DM_DM(iter = 10000, K_max = k, z = data_sim[[t]]$data$z,
                    theta_vec = rep(1, k), MH_var = 10, mu = 0, s2 = 1,
                    print_iter = 2500)
    dd_time <- difftime(Sys.time(), start_time, units = "secs")
    dd_assign <- as.numeric(salso(dd_mod[-c(1:5000), ], maxNClusters = 10))
    return(list(t = t, k = k, dd_assign = dd_assign, dd_time = dd_time))
  }
stopImplicitCluster()






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

### Data Simulation
pat_mat1 <- diag(20)[1:5, ]
set.seed(72)
sim_dat <- simDM(n = 100, pattern = pat_mat1, xi_conc = 10, pi_gm = c(0.95), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 750, 
                 theta = 0.01)
simDM_sum(sim_dat)
log_at_risk(z = sim_dat$z, gamma_mat = sim_dat$gamma)


# ------------------------------------------------------------------------------

### Check: log_gamma_ijk function
### Check for only z_ijk = 0

index_mat <- cbind(which(sim_list$z[, 1:4] == 0) - floor(which(sim_list$z[, 1:4] == 0)/100) * 100,
                   ceiling(which(sim_list$z[, 1:4] == 0)/100))

result_calc_mat <- matrix(NA, nrow = 11, ncol = 9)

for(i in 1:11){
  
  ### Calculate by hand
  gijk <- sim_list$gamma[index_mat[i, 1], index_mat[i, 2]]
  gi <- sim_list$gamma[index_mat[i, 1], 1:4]
  zi <- sim_list$z[index_mat[i, 1], 1:4]
  ci <- sim_list$ci[index_mat[i, 1]] + 1
  xi <- exp(sim_list$beta[ci, 1:4])
  
  hand_calc_1 <- lbeta(1 + gijk, 1 + (1 - gijk)) + lgamma(sum(gi * xi)) - lgamma(sum(zi + (gi*xi)))
  
  ### By using the function
  jj <- index_mat[i, 2] - 1
  func_calc_1 <- log_g_ijk(j = jj, zi = sim_list$z[index_mat[i, 1], ], gi = sim_list$gamma[index_mat[i, 1], ], 
                         w = c(rep(1, 4), rep(0, 6)), beta_k = sim_list$beta[ci, ],
                         r0g = 1, r1g = 1)
  
  ### Adjusted the gijk = (1 - gijk)
  gijk <- 1 - gijk
  gi[index_mat[i, 2]] <- gijk
  hand_calc_0 <- lbeta(1 + gijk, 1 + (1 - gijk)) + lgamma(sum(gi * xi)) - lgamma(sum(zi + (gi*xi)))
  
  ### By using the function
  gi_adj <- sim_list$gamma[index_mat[i, 1], ]
  gi_adj[index_mat[i, 2]] <- gijk
  func_calc_0 <- log_g_ijk(j = jj, zi = sim_list$z[index_mat[i, 1], ], gi = gi_adj, 
                           w = c(rep(1, 4), rep(0, 6)), beta_k = sim_list$beta[ci, ],
                           r0g = 1, r1g = 1)
  
  result_calc_mat[i, ] <- c(index_mat[i, 1], index_mat[i, 2], ci - 1, 1 - gijk, 
                            xi[index_mat[i, 2]], hand_calc_0, func_calc_0, hand_calc_1, func_calc_1)
  
}

result_mat_1 <- data.frame(result_calc_mat)
xtable(result_mat_1, digits = c(rep(0, 5), rep(5, 5)))

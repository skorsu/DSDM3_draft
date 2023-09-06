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
pat_mat1 <- diag(10)[1:5, ]
set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 10, pi_gm = c(0.75), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 50, z_sum_U = 100, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)

xx <- sample(c(2, 3), 200, replace = TRUE)
bb <- matrix(0, nrow = 5, ncol = 10)
bb[2:3, ] <- rnorm(20, 0, sqrt(0.01))

tt <- sm(K_max = 5, z = sim_dat$z, xx - 1,
         gamma_mat = sim_dat$gamma, beta_mat = bb,
         tau_vec = c(0, 1, 1, 0, 0), theta_vec = rep(1, 5), launch_iter = 10, 
         mu = 0, s2 = 0.01, r0c = 1, r1c = 1)

tt$expand_ind

log_proposal(clus_after = tt$proposed_assign, clus_before = rep(1, 200), 
             z = sim_dat$z, gamma_mat = sim_dat$gamma, 
             beta_mat = tt$proposed_beta, S = tt$S, 
             clus_sm = tt$samp_clus)

log_proposal(clus_after = rep(1, 200), clus_before = tt$proposed_assign, 
             z = sim_dat$z, gamma_mat = sim_dat$gamma, 
             beta_mat = tt$proposed_beta, S = tt$S, 
             clus_sm = tt$samp_clus)


table(tt$proposed_assign)

tt$expand_ind

table(proposed = tt$proposed_assign, launch = tt$launch_assign)

cbind(tt$launch_tau, tt$proposed_tau)
tt$launch_beta
tt$proposed_beta


table(proposed = tt$proposed_assign)

table(`launch` = tt$launch_assign)
table(`launch` = tt$launch_assign, init = xx)
table(sim_dat$ci - 1, xx)
table(sim_dat$ci - 1, tt$launch_assign)
tt$samp_ind
tt$samp_clus



(sim_dat$ci - 1)[which((sim_dat$ci - 1) %in% c(4, 2))]

table(sim_dat$ci - 1)


table(x - 1)

colSums(sim_dat$z[which(sim_dat$ci == 2), ])/sum(sim_dat$z[which(sim_dat$ci == 2), ])
start_time <- Sys.time()
result <- beta_ar_update(K = 5, iter = 10000, z = sim_dat$z, clus_assign = sim_dat$ci - 1, 
                         r0g = 1, r1g = 1, mu = 0, s2 = 1000, s2_MH = 0.01)
Sys.time() - start_time
k <- 2
xi_mcmc <- exp(t(result$beta[k, , -c(1:5000)]))/rowSums(exp(t(result$beta[k, , -c(1:5000)])))
plot(xi_mcmc[, 2], type = "l")


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

library(foreach)
library(doParallel)
library(salso)
library(mclustcomp)
library(xtable)

### Import function
path_mb <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data"
path_pc <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data"

tryCatch({
  message("Try MB")
  source(paste0(path_mb, "/new_data_sim.R"))
  }, error=function(cond){
    message("Try PC")
    source(paste0(path_pc, "/new_data_sim.R"))
  })

### Function
report_table <- function(mat){
  mm <- round(colMeans(mat), 4)
  ss <- round(apply(mat, 2, sd), 4)
  paste0(mm, " (", ss, ")")
}

### Create a plot
set.seed(12)
sim_plot(n = 100, J = 10, K = 5, scenario = 1, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_plot(n = 100, J = 10, K = 5, scenario = 2, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_plot(n = 100, J = 10, K = 5, scenario = 3, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_plot(n = 100, J = 10, K = 5, scenario = 4, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)

### Run the models
registerDoParallel(detectCores() - 3)
result <- foreach(t = 1:10) %dopar% {
  
  ### Matrix for storing result
  final_re <- matrix(NA, ncol = 4, nrow = 3)
  colnames(final_re) <- c("Time", "AdjRand", "Jaccard", "VI")
  rownames(final_re) <- c("Our", "DMZ", "DMDM")
  
  ### Generate the data
  set.seed(2023 + t)
  sim_list <- dat_sim(n = 100, J = 10, K = 5, scenario = 4, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
  
  ### Our model
  start_time <- Sys.time()
  KK <- 10
  mod <- full_func_at_risk(iter = 20000, K_max = KK, z = sim_list$z,
                           w = rep(1, KK), theta = rep(1, KK), launch_iter = 10,
                           MH_var = 0.01, s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1)
  tot_time <- as.numeric(difftime(Sys.time(), start_time, "secs"))
  ci_result <- as.numeric(salso(t(mod$ci)[-(1:5000), ], maxNClusters = KK)) 
  
  final_re[1, ] <- c(tot_time, mclustcomp(sim_list$ci, ci_result)[c(1, 5, 22), 2])
  
  ### DM-ZIDM
  start_time <- Sys.time()
  mod <- full_func(iter = 20000, K_max = KK, z = sim_list$z,
                   gamma_mat = matrix(1, ncol = 10, nrow = 100), 
                   w = rep(1, KK), theta = rep(1, KK), launch_iter = 10,
                   MH_var = 0.01, s2 = 1, r0c = 1, r1c = 1)
  tot_time <- as.numeric(difftime(Sys.time(), start_time, "secs"))
  ci_result <- as.numeric(salso(t(mod$ci)[-(1:5000), ], maxNClusters = KK)) 
  
  final_re[2, ] <- c(tot_time, mclustcomp(sim_list$ci, ci_result)[c(1, 5, 22), 2])
  
  ### DM-DM
  start_time <- Sys.time()
  mod <- clus_no_SM(iter = 20000, K_max = KK, z = sim_list$z,
                    gamma_mat = matrix(1, ncol = 10, nrow = 100), 
                    w = rep(1, KK), MH_var = 0.01, s2 = 1, theta_vec = rep(1, KK))
  tot_time <- as.numeric(difftime(Sys.time(), start_time, "secs"))
  ci_result <- as.numeric(salso(mod$ci[-(1:5000), ], maxNClusters = KK)) 
  
  final_re[3, ] <- c(tot_time, mclustcomp(sim_list$ci, ci_result)[c(1, 5, 22), 2])
  
  return(final_re)
  
  }
stopImplicitCluster()

### Summarize
our_mod <- matrix(NA, ncol = 4, nrow = 10)
dmz_mod <- matrix(NA, ncol = 4, nrow = 10)
dmdm_mod <- matrix(NA, ncol = 4, nrow = 10)

for(i in 1:10){
  our_mod[i, ] <- result[[i]]["Our", ]
  dmz_mod[i, ] <- result[[i]]["DMZ", ]
  dmdm_mod[i, ] <- result[[i]]["DMDM", ]
}

xtable(t(sapply(list(our_mod, dmz_mod, dmdm_mod), report_table)))

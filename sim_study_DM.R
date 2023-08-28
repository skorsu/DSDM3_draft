### Required Libraries
library(foreach)
library(doParallel)

### Import the function
im_path <- NULL
mb_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/"
path <- NULL
if(dir.exists(mb_path)){
  path <- mb_path
} else {
  path <- im_path
}

source(paste0(path, "data_sim_DM.R"))

registerDoParallel(detectCores() - 1)
pat_mat <- diag(5)
result <- foreach(t = 1:50, .combine = "rbind") %dopar% {
  set.seed(t + 21)
  sim_dat <- simDM(n = 100, pattern = pat_mat, xi_conc = 10, pi_gm = c(0.95), 
                   pi_c = c(1, 1, 1, 1, 1), z_sum = 2500, theta = 0.01)
  summary_data <- simDM_sum(sim_dat)
  return(summary_data$summary_zero)
}
stopImplicitCluster()

result

set.seed(2)
pat_mat <- diag(3)
pat_mat <- diag(5)
pat_mat[3, 4] <- 1
pat_mat[4, 5] <- 1
pat_mat[5, 3] <- 1

simDM(n = 100, pattern = pat_mat, xi_conc = 10, pi_gm = c(0.95), 
      pi_c = c(1, 1, 1, 1, 1), z_sum = 2500, theta = 0.01) %>% simDM_sum() %>% .$summary_zero


sim_dat <- simDM(n = 100, pattern = pat_mat, xi_conc = 10, pi_gm = c(1), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum = 2500, theta = 0.01)
summary_data <- simDM_sum(sim_dat)
grid.arrange(grobs = summary_data$plot)
sim_dat$z

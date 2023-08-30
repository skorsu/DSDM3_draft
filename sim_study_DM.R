### Required Libraries
library(foreach)
library(doParallel)

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

### Generate 50 datasets to see the ratio of zero and proportion of at-risk zero
pat_mat <- diag(10)[1:5, ] ### Scenario 1
pat_mat <- rbind(c(1, rep(0, 9)), c(0, 1, rep(0, 8)), c(rep(0, 7), 1, 1, 0),
                 c(rep(0, 7), 0, 1, 1), c(rep(0, 7), 1, 0, 1)) ### Scenario 2
pat_mat <- rbind(c(rep(1, 10), rep(0, 20)), c(rep(0, 10), rep(1, 10), rep(0, 10)), 
                 c(rep(0, 20), rep(1, 10))) ### Scenario 3
pat_mat <- rbind(c(rep(1, 12), rep(0, 18)), c(rep(0, 9), rep(1, 12), rep(0, 9)), 
                 c(rep(0, 18), rep(1, 12))) ### Scenario 4

registerDoParallel(detectCores() - 1)
result <- foreach(t = 1:20, .combine = "rbind") %dopar% {
  set.seed(t + 83)
  sim_dat <- simDM(n = 100, pattern = pat_mat, xi_conc = 10, pi_gm = c(0.75), 
                   pi_c = c(1, 1, 1, 1, 1), z_sum_L = 250, z_sum_U = 500, 
                   theta = 0.01)
  summary_data <- simDM_sum(sim_dat)
  return(summary_data$summary_zero)
}
stopImplicitCluster()

mean(result[, 1] * (100 * ncol(pat_mat)))
sd(result[, 1] * (100 * ncol(pat_mat)))

apply(result, 2, mean)
apply(result, 2, sd)

### Create a plot
set.seed(72)
sim_dat <- simDM(n = 100, pattern = pat_mat, xi_conc = 10, pi_gm = c(0.75), 
                 pi_c = c(1, 1, 1), z_sum_L = 250, z_sum_U = 500, 
                 theta = 0.01)
summary_data <- simDM_sum(sim_dat)
grid.arrange(grobs = summary_data$plot)
sim_dat$z


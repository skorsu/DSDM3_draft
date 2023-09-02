### Required Libraries
library(foreach)
library(doParallel)
library(salso)

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

### Scenario 1: Strong signal w/ no overlap
pat_mat1 <- diag(20)[1:5, ]
set.seed(72)
sim_dat <- simDM(n = 100, pattern = pat_mat1, xi_conc = 10, pi_gm = c(0.95), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 750, 
                 theta = 0.01)
test <- simDM_sum(sim_dat)
sum(apply(sim_dat$z, 1, function(x){sum(x == 0)}) == 9)
grid.arrange(grobs = test$plot)
start_time <- Sys.time()
result <- full_func_at_risk(iter = 25000, K_max = 10, z = sim_dat$z,
                            w = c(rep(1, 20), rep(0, 0)), theta = rep(1, 10), launch_iter = 10,
                            MH_var = 0.1, s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1)
Sys.time() - start_time
table(sim_dat$ci, as.numeric(salso(t(result$ci)[-(1:5000), ], maxNClusters = 10)))

plot(apply(result$ci, 2, function(x){length(unique(x))}))
summary(result$logA)

mean(apply(sim_dat$z, 1, function(x){sum(x == 0)}))


plot(apply(result$ci, 2, function(x){length(unique(x))}))
summary(result$logA)
sort(result$logA, decreasing = TRUE)[1:10]

library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_VI.rds"))

### Function
bal_tab <- function(x, r = 4){
  paste0(round(mean(x, na.rm = TRUE), r), " (", round(sd(x, na.rm = TRUE), r), ")")
}

unique_clus <- function(x){
  length(unique(x))
}

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Import the data
dat_case_V <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_V.rds"))
dat <- dat_case_V[[1]]$dat ### Data Used

### Play with K_max
K_max_List <- c(5, 5, 5, 50, 50, 50)
chainIndex <- rep(1:3, 2)
  
### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
foreach(t = 1:6) %dopar% {
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = K_max_List[t], nbeta_split = 50, 
                      z = as.matrix(dat), atrisk_init = matrix(1, nrow = 50, ncol = 250), 
                      beta_init = matrix(0, nrow = K_max_List[t], ncol = 250), 
                      ci_init = rep(0, 50), theta = 1, mu = 0, s2 = 0.1, 
                      s2_MH = 1e-3, t_thres = 1000, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Sensitivity/sensitivity_K_max_", K_max_List[t],
          "_chain_", chainIndex[t], "_.rds"))
  
}
stopImplicitCluster()

# 
# ### Try with Kmax
# mod <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 50, 
#                     z = dat, atrisk_init = matrix(1, nrow = 50, ncol = 250), 
#                     beta_init = matrix(0, nrow = 10, ncol = 250), 
#                     ci_init = rep(0, 50), theta = 1, mu = 0, s2 = 0.1, 
#                     s2_MH = 1e-3, t_thres = 1000, launch_iter = 30, 
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# 
# apply(mod$ci_result, 1, function(x){length(unique(x))}) %>% plot(type = "l")
# as.numeric(salso(mod$ci_result)) %>% table(dat_case_V[[1]]$clus)


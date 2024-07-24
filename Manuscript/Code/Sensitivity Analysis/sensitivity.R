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

### Hyperparameters
hyperSet <- matrix(NA, ncol = 9, nrow = 19)
colnames(hyperSet) <- c("Kmax", "n_eta", "theta", "s2", "s2MH", "ag", "bg", "ac", "bc")
hyperSet[1, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[2, ] <- c(5, 50, 1, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[3, ] <- c(20, 50, 1, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[4, ] <- c(10, 10, 1, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[5, ] <- c(10, 100, 1, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[6, ] <- c(10, 50, 0.1, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[7, ] <- c(10, 50, 10, 0.1, 1e-3, 1, 1, 1, 1)
hyperSet[8, ] <- c(10, 50, 1, 1, 1e-3, 1, 1, 1, 1)
hyperSet[9, ] <- c(10, 50, 1, 0.01, 1e-3, 1, 1, 1, 1)
hyperSet[10, ] <- c(10, 50, 1, 0.1, 0.1, 1, 1, 1, 1)
hyperSet[11, ] <- c(10, 50, 1, 0.1, 1e-5, 1, 1, 1, 1)
hyperSet[12, ] <- c(10, 50, 1, 0.1, 1e-3, 4, 1, 1, 1)
hyperSet[13, ] <- c(10, 50, 1, 0.1, 1e-3, 9, 1, 1, 1)
hyperSet[14, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 4, 1, 1)
hyperSet[15, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 9, 1, 1)
hyperSet[16, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 1, 4, 1)
hyperSet[17, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 1, 9, 1)
hyperSet[18, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 1, 1, 4)
hyperSet[19, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 1, 1, 9)

### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
foreach(t = 1:19) %dopar% {
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = hyperSet[t, "Kmax"], nbeta_split = hyperSet[t, "n_eta"], 
                      z = as.matrix(dat), atrisk_init = matrix(1, nrow = 50, ncol = 250), 
                      beta_init = matrix(0, nrow = hyperSet[t, "Kmax"], ncol = 250), 
                      ci_init = rep(0, 50), theta = hyperSet[t, "theta"], mu = 0, 
                      s2 = hyperSet[t, "s2"], s2_MH = hyperSet[t, "s2MH"], t_thres = 1000, launch_iter = 30, 
                      r0g = hyperSet[t, "ag"], r1g = hyperSet[t, "bg"], 
                      r0c = hyperSet[t, "ac"], r1c = hyperSet[t, "bc"], thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Sensitivity/sensitivity_case_", (t-1), "_.rds"))
  
}
stopImplicitCluster()

### Play with K_max
# K_max_List <- c(5, 5, 5, 50, 50, 50)
# chainIndex <- rep(1:3, 2)

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


library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

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
dat_case <- readRDS(paste0(path, "Data/Simulation Study/simu_data_case_III.rds"))
dat <- dat_case[[1]]$dat ### Data Used

### Hyperparameters
hyperSet <- matrix(NA, ncol = 9, nrow = 15)
colnames(hyperSet) <- c("Kmax", "n_eta", "theta", "s2", "s2MH", "bg", "bc")
hyperSet[1, ] <- c(10, 20, 1, 0.1, 1e-3, 1, 4) ### Default

### Kmax
hyperSet[2, ] <- c(5, 20, 1, 0.1, 1e-3, 1, 4)
hyperSet[3, ] <- c(25, 20, 1, 0.1, 1e-3, 1, 4)

### n_eta
hyperSet[4, ] <- c(10, 10, 1, 0.1, 1e-3, 1, 4)
hyperSet[5, ] <- c(10, 50, 1, 0.1, 1e-3, 1, 4)

### theta
hyperSet[6, ] <- c(10, 20, 0.1, 0.1, 1e-3, 1, 4)
hyperSet[7, ] <- c(10, 20, 10, 0.1, 1e-3, 1, 4)

### s2
hyperSet[8, ] <- c(10, 20, 1, 0.01, 1e-3, 1, 4)
hyperSet[9, ] <- c(10, 20, 1, 1, 1e-3, 1, 4)

### s2MH
hyperSet[10, ] <- c(10, 20, 1, 0.1, 1e-1, 1, 4)
hyperSet[11, ] <- c(10, 20, 1, 0.1, 1e-5, 1, 4)

### bg
hyperSet[12, ] <- c(10, 20, 1, 0.1, 1e-3, 4, 4)
hyperSet[13, ] <- c(10, 20, 1, 0.1, 1e-3, 9, 4)

### bc
hyperSet[14, ] <- c(10, 20, 1, 0.1, 1e-3, 1, 1)
hyperSet[15, ] <- c(10, 20, 1, 0.1, 1e-3, 1, 9)


### Run the model
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
foreach(t = 1:15) %dopar% {
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = hyperSet[t, "Kmax"], nbeta_split = hyperSet[t, "n_eta"], 
                      z = as.matrix(dat), atrisk_init = matrix(1, nrow = 50, ncol = 250), 
                      beta_init = matrix(0, nrow = hyperSet[t, "Kmax"], ncol = 250), 
                      ci_init = rep(0, 50), theta = hyperSet[t, "theta"], mu = 0, 
                      s2 = hyperSet[t, "s2"], s2_MH = hyperSet[t, "s2MH"], t_thres = 1000, launch_iter = 30, 
                      r0g = 1, r1g = hyperSet[t, "bg"], 
                      r0c = 1, r1c = hyperSet[t, "bc"], thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(mod = mod, time = tot_time), 
          file = paste0(path, "Result/Sensitivity/sensitivity_scenario_III_case_", (t-1), "_.rds"))
  
}
stopImplicitCluster()


rm(list = ls())

### Required Libraries
library(salso)
library(mclustcomp)
library(foreach)
library(doParallel)

### Import the data
path <- "/Users/kevinkvp/Desktop/result_local/"
case_name <- "strong_sig_no_overlap_zi"
dat <- readRDS(paste0(path, "data_", case_name, ".RData"))
result <- readRDS(paste0(path, "result_", case_name, ".RData"))

### Run the DM-sDM
registerDoParallel(5)
dmsdm <- foreach(t = 1:30) %dopar% {
  stime <- Sys.time()
  mod <- DM_DM(iter = 10000, K_max = 10, z = dat$dat[[t]]$data$z,
               theta_vec = rep(0.001, 10), MH_var = 1, mu = 0, s2 = 1,
               print_iter = 2500)
  dd_time <- difftime(Sys.time(), stime, "secs")
  return(list(assign = mod, time = dd_time))
}
stopImplicitCluster()

dmsdm[[1]]$time

### DM-DM: Calculate the measure quantity (Jaccard, AdjRand, and VI)
dd_quan <- matrix(NA, ncol = 4, nrow = 30)
dsd_quan <- matrix(NA, ncol = 4, nrow = 30)
for(i in 1:30){
  
  dd_quan[i, ] <- c(result$dd_result[[i]][[4]]$dd_time, 
                    mclustcomp(dat$dat[[i]]$data$ci, 
                               as.numeric(salso(result$dd_result[[i]][[4]]$dd_mod[-(1:5000), ], 
                                                maxNClusters = 10)))[c(5, 1, 22), 2])
  dsd_quan[i, ] <- c(dmsdm[[i]]$time, 
                     mclustcomp(dat$dat[[i]]$data$ci, 
                                as.numeric(salso(dmsdm[[1]]$assign)))[c(5, 1, 22), 22])
  
  
}

colMeans(dd_quan)
apply(dd_quan, 2, sd)



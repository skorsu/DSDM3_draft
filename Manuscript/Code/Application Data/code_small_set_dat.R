### Load libraries: ------------------------------------------------------------
library(tidyverse)
library(foreach)
library(doParallel)

### Change the settings: -------------------------------------------------------
path <- "~/Desktop/Github - Repo/ClusterZI/Manuscript/"
datpath <- paste0(path, "Data/Application Data/Cleaned Data/Species/")
resultpath <- paste0(path, "Result/Application Data/")

### Hyperparameters: -----------------------------------------------------------
datname <- rep(c("crc_xiang", "crc_zhao", "mhe_zhang"), 2)
xiSplit <- sort(rep(c(10, 25), 3))

### Start the code -------------------------------------------------------------
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
foreach(t = 1:length(datname)) %dopar% {
  start_time <- Sys.time()
  dat <- readRDS(paste0(datpath, datname[t], "_cleaned_species.rds"))
  mod <- mod_adaptive(iter = 10000, Kmax = 10, nbeta_split = xiSplit[t],
                      z = as.matrix(dat),
                      atrisk_init = matrix(1, nrow = dim(dat)[1], ncol = dim(dat)[2]),
                      beta_init = matrix(0, nrow = 10, ncol = dim(dat)[2]),
                      ci_init = rep(0, dim(dat)[1]),
                      theta = 1e-3, mu = 0, s2 = 1e-3,
                      s2_MH = 1e-3, t_thres = 1000, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 10)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time), file = paste0(resultpath, "result_", datname[t], "_xiSplit_", xiSplit[t], "_theta_1e3_ADP_1000_iter_10000.rds"))
}
stopImplicitCluster()

################################################################################









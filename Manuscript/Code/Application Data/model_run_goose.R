### Load libraries
library(tidyverse)
library(foreach)
library(doParallel)
library(ClusterZI)

### Change the settings --------------------------------------------------------
path <- "~/Desktop/Github - Repo/ClusterZI/Manuscript/Data/Application Data/Cleaned Data/Species/"
nCores <- 5

### Start the code -------------------------------------------------------------
datname <- rep(c("cdi_vincent_v3v5", "ob_ross", "edd_singh", "cdi_schubert", 
                 "ibd_alm", "ibd_huttenhower", "ob_gordon_2008_v2", "crc_zeller"), 2)
xiChange <- sort(rep(c(10, 50), 8))

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(nCores)
foreach(t = 1:length(datname)) %dopar% {
  start_time <- Sys.time()
  dat <- readRDS(paste0(data_path, datname[t], "_cleaned_species.rds"))
  mod <- mod_adaptive(iter = 2, Kmax = 10, nbeta_split = xiChange[t],
                      z = as.matrix(dat),
                      atrisk_init = matrix(1, nrow = dim(dat)[1], ncol = dim(dat)[2]),
                      beta_init = matrix(0, nrow = 10, ncol = dim(dat)[2]),
                      ci_init = rep(0, dim(dat)[1]),
                      theta = 1e-3, mu = 0, s2 = 1e-3,
                      s2_MH = 1e-3, t_thres = 1000, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 10)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time), file = paste0(result_path, "result_", datname[t], "_betaSplit_", xiChange[t], "_theta_1e3_ADP_1000_iter_10000.rds"))
}
stopImplicitCluster()

################################################################################









### Load libraries
library(foreach)
library(doParallel)
library(ClusterZI)

### Change the settings --------------------------------------------------------
path <- "~/Desktop/Github Repo/ClusterZI/Manuscript/Data/Application Data/Cleaned Data/Species/"
datname <- c("cdi_vincent_v3v5", "ob_ross", "edd_singh", "cdi_schubert", 
             "ibd_alm", "ibd_huttenhower", "ob_gordon_2008_v2", "crc_zeller")

paste0(path, datname, "_cleaned_species.rds") %>%
  file.exists()

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
foreach(t = 1:length(datname)) %dopar% {
  dumDat <- matrix(t, ncol = t, nrow = t)
  saveRDS(dumDat, file = paste0(path, "TEST_", datname[t], ".rds"))
}
stopImplicitCluster()



################################################################################









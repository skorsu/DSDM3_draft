### Load libraries
library(foreach)
library(doParallel)
library(salso)
library(ggplot2)

### Change the settings --------------------------------------------------------
case_name <- "diffindex_3_K_5"
nData <- 20

### Import the data and result -------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

dat <- readRDS(paste0(path, "ClusterZI/Manuscript/Data/", case_name, "_simDat.RData"))
result_path <- paste0(path, "ClusterZI/Manuscript/Result/Simulation Study/",
                      case_name, paste0(c("_ZZ", "_DZ", "_DD", "_DsD"), ".RData"))
registerDoParallel(5)
result_list <- foreach(t = 1:length(result_path)) %dopar% {
 readRDS(result_path[[t]]) 
}
stopImplicitCluster()
names(result_list) <- c("ZZ", "DZ", "DD", "DsD")

### Analyze the result ---------------------------------------------------------
#### Plot MCMC

table(salso(result_list$ZZ[[1]]$result$ci_result), dat[[1]]$clus)






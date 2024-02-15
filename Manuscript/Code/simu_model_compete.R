### Load libraries
library(tidyverse)
library(foreach)
library(doParallel)
library(sparseMbClust)
library(ecodist)
library(ape)
library(DTMM)
library(salso)

### Change the settings --------------------------------------------------------
case_name <- "diffindex_3_K_5"
nData <- 20

### Import the data ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

dat <- readRDS(paste0(path, "ClusterZI/Manuscript/Data/", case_name, "_simDat.RData"))
dat_J <- ncol(dat[[1]]$dat) 
dat_n <- nrow(dat[[1]]$dat) 

save_path <- paste0(path, "ClusterZI/Manuscript/Result/Simulation Study/")

### Adjusted the DTMM function -------------------------------------------------
adjDTMM <- function(Y, tree, tau_vec = 10 ^ seq(-1, 4, 0.5), nu_vec = 1, 
                    theta_vec = seq(0.01, 0.99, 0.08), init_c, 
                    beta = "default", mcmc_iter = 2000){
  
  if(beta == "default"){
    temp = runif(1)
    alpha = temp/(1 - temp)
  }
  
  tree <- reorder(tree, order = "postorder")
  edge <- apply(tree$edge, 2, rev)
  p = dim(Y)[2]
  gamma_sample = rbinom(p - 1, 1, 0.5)
  
  HDTMcpp(Y, edge, tau_vec, nu_vec, theta_vec,
          init_c, gamma_sample, alpha, mcmc_iter, select = FALSE)
}


### Run the models -------------------------------------------------------------
#### DTMM: with no structures
set.seed(3124, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
resultDTMM <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- adjDTMM(dat[[t]]$dat, 
                         tree = read.tree(text = paste("(", paste(1:50, collapse = ", "), ");")), 
                         init_c = rep(1, 100), mcmc_iter = 100000)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = t(clus_result$post_c[, (1:1000) * 100]))
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(resultDTMM, paste0(save_path, case_name, "_DTMM_no_structure.RData"))

################################################################################

### Reconstruct the phylogenetic trees based on the dataset
### https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html#objectives
tre <- nj(dist(t(dat[[1]]$dat)))
tre <- ladderize(tre)
plot(tre, cex = 0.6)
plot(as.vector(dist(t(dat[[1]]$dat))), as.vector(as.dist(cophenetic(tre))))
cor(as.vector(dist(t(dat[[1]]$dat))), as.vector(as.dist(cophenetic(tre))))^2









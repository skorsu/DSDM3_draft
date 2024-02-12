library(Rcpp)
library(salso)
library(foreach)
library(doParallel)
library(mclustcomp)
library(cluster)
library(ecodist)
library(ggplot2)
library(plotrix)
library(latex2exp)
library(sparseMbClust)
library(tidyverse)
library(pheatmap)
library(mixtools)
library(coda.base)

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
# sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

### Set of the hyperparameter --------------------------------------------------
hyperParam <- list(c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 5, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 20, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 1, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 0.1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 10, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 0.1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 10, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 0.1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 5, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 20, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 0.1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 10, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 0.1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 10, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 4, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 9, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 4),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 9))

### Import the data and result -------------------------------------------------
save_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/sensitivity/"
# save_path <- "/Users/kevin-imac/Desktop/sensitivity_0119/"
case_name <- "diffindex_3_K_5"
dat <- readRDS(paste0(save_path, case_name, "_simDat.RData")) ## Data
nHyperSet <- length(hyperParam)

registerDoParallel(5)
result <- foreach(t = 1:nHyperSet) %dopar% {
  file_case <- paste(paste0(names(hyperParam[[t]]), "_", hyperParam[[t]]), 
                     collapse = "_")
  readRDS(paste0(save_path, case_name, "_", file_case, "_MB.RData"))
}
stopImplicitCluster()

### Example of the data --------------------------------------------------------
pheatmap(dat[[1]]$dat[sort(dat[[1]]$clus, index.return = TRUE)$ix, ], 
         display_numbers = F, color = colorRampPalette(c('white','lightblue3'))(100), 
         cluster_rows = F, cluster_cols = F)

### Analyze the result ---------------------------------------------------------
meanSD <- function(x, dplace = 3){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

uniqueClus <- function(x){
  length(unique(x))
}

nData <- length(dat)
mclustcomp(as.numeric(salso(result[[1]]$result[[1]]$result)), dat[[1]]$clus)
ssVI <- sapply(1:nData, function(x){salso(result[[1]]$result[[x]]$result[-c(1:500), ])})


actual_clus <- sapply(1:nData, function(x){dat[[x]]$clus})

for(i in 1:length(result)){
  
  print(noquote(" ================================================================= "))
  
  print(result[[i]]$hyper)
  print(noquote(" "))
  
  ### Computational Time
  comp_time <- sapply(1:nData, function(x){result[[i]]$result[[x]]$time}) %>% meanSD(2)
  print(noquote(paste0("Time: ", comp_time)))
  
  ### Cluster: MCMC
  clusMCMC <- sapply(1:nData,
                     function(x){mean(apply(result[[i]]$result[[x]]$result, 1, uniqueClus))}) %>%
    meanSD()
  print(noquote(paste0("Avg Cluster in MCMC: ", clusMCMC)))
  
  ### Performance
  ### Loss: VI
  salso_VI <- sapply(1:nData, function(x){salso(result[[i]]$result[[x]]$result[-c(1:500), ])})
  VI_result <- suppressWarnings(sapply(1:nData, 
                                       function(x){mclustcomp(salso_VI[, x], actual_clus[, x], type = c("adjrand", "jaccard", "vi"))[, 2]}) %>%
                                  apply(1, meanSD))
  VI_clus <- apply(salso_VI, 2, uniqueClus) %>% meanSD()
  
  ### Loss: binder
  salso_binder <- sapply(1:nData, function(x){salso(result[[i]]$result[[x]]$result[-c(1:500), ], loss = "binder")})
  salso_result <- suppressWarnings(sapply(1:nData, 
                                          function(x){mclustcomp(salso_binder[, x], actual_clus[, x], type = c("adjrand", "jaccard", "vi"))[, 2]}) %>%
                                     apply(1, meanSD))
  binder_clus <- apply(salso_binder, 2, uniqueClus) %>% meanSD()
  
  rbind(VI = c(VI_result, VI_clus), salso = c(salso_result, binder_clus)) %>%
    `colnames<-`(c("Adj. Rand. Index", "Jaccard", "VI", "Active Clusters")) %>%
    print()
  
}

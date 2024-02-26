### Load libraries
library(foreach)
library(doParallel)
library(salso)
library(tidyverse)
library(xtable)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(mclustcomp)
library(Rcpp)
library(cluster)
library(ecodist)
library(coda.base)

### User-defined Functions -----------------------------------------------------
meanSD <- function(x, dplace = 5){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

uniqueClus <- function(clus_assign){
  length(unique(clus_assign))
}

### Change the settings --------------------------------------------------------
case_name <- "diffindex_3_K_2"
nData <- 20

### Import the data and result -------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

sourceCpp(paste0(path, "ClusterZI/src/clusterZI.cpp"))

dat <- readRDS(paste0(path, "ClusterZI/Manuscript/Data/", case_name, "_simDat.RData"))
result_path <- paste0(path, "ClusterZI/Manuscript/Result/Simulation Study/",
                      case_name, paste0(c("_ZZ", "_DZ", "_DD", "_DsD", "_DTMM_no_structure"), ".RData"))
registerDoParallel(5)
result_list <- foreach(t = 1:length(result_path)) %dopar% {
 readRDS(result_path[[t]]) 
}
stopImplicitCluster()
modName <- c("ZZ", "DZ", "DD", "DsD", "DTMM_no_structure")
names(result_list) <- modName

### Analyze the result ---------------------------------------------------------
#### Plot MCMC
##### ZIDM-ZIDM
mcmcZZ <- sapply(1:nData, function(x){apply(result_list$ZZ[[x]]$result$ci_result, 1, uniqueClus)}) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = value)) +
  geom_line(aes(color = factor(Var2))) + 
  labs(x = "Iteration (Thinning every 1000 iterations from 100,000 iterations)",
       y = "Active Clusters",
       title = "Active clusters via MCMC for the ZIDM-ZIDM model") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(2, 20, 2), limits = c(0, 20))

##### DM-ZIDM
mcmcDZ <- sapply(1:nData, function(x){apply(result_list$DZ[[x]]$result$ci_result, 1, uniqueClus)}) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = value)) +
  geom_line(aes(color = factor(Var2))) + 
  labs(x = "Iteration (Thinning every 1000 iterations from 100,000 iterations)",
       y = "Active Clusters",
       title = "Active clusters via MCMC for the DM-ZIDM model") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(2, 20, 2), limits = c(0, 20))

##### DTMM
mcmcDTMM <- sapply(1:nData, function(x){apply(result_list$DTMM_no_structure[[x]]$result, 1, uniqueClus)}) %>%
  melt() %>%
  ggplot(aes(x = Var1, y = value)) +
  geom_line(aes(color = factor(Var2))) + 
  labs(x = "Iteration (Thinning every 1000 iterations from 100,000 iterations)",
       y = "Active Clusters",
       title = "Active clusters via MCMC for the DTMM model") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(2, 20, 2), limits = c(0, 20))

grid.arrange(mcmcZZ, mcmcDZ, mcmcDTMM)

#### Computation Time
timeMat <- matrix(NA, ncol = length(modName), nrow = nData)
colnames(timeMat) <- modName
timeMat[, -3] <- sapply(c("ZZ", "DZ", "DsD", "DTMM_no_structure"), function(x){sapply(1:nData, function(y){result_list[[x]][[y]]$time})})
timeMat[, 3] <- colSums(sapply(1:nData, function(x){sapply(1:9, function(y){result_list$DD[[x]][[y]]$time})}))

#### Active Clusters via MCMC
actveMCMC <- rep(NA, 5)
actveMCMC[1] <- meanSD(colMeans(sapply(1:nData, function(x){apply(result_list$ZZ[[x]]$result$ci_result, 1, uniqueClus)})), 2)
actveMCMC[2] <- meanSD(colMeans(sapply(1:nData, function(x){apply(result_list$DZ[[x]]$result$ci_result, 1, uniqueClus)})), 2)
actveMCMC[5] <- meanSD(colMeans(sapply(1:nData, function(x){apply(result_list$DTMM_no_structure[[x]]$result, 1, uniqueClus)})), 2)

t(data.frame(apply(timeMat/60, 2, meanSD, dplace = 2), actveMCMC)) %>%
  `rownames<-`(c("Computational Time (mins)", "Active Cluster via MCMC")) %>%
  `colnames<-`(c("ZIDM-ZIDM", "DM-ZIDM", "DM-DM", "DM-sDM", "DTMM")) %>%
  xtable()

#### Post-processing
aClus <- sapply(1:nData, function(x){dat[[x]]$clus})
loss_type <- c("VI", "binder")

##### ZIDM-ZIDM
clusZZ <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$ZZ[[x]]$result$ci_result[-(1:500), ], loss = loss_type[a]))})})})
activeZZ <- sapply(1:2, function(x){apply(clusZZ[[x]], 2, uniqueClus)})
ariZZ <- sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusZZ[[a]][, x], aClus[, x], type = "adjrand")[2])})})
apply(activeZZ, 2, meanSD, dplace = 2) ### Active Cluster
apply(ariZZ, 2, meanSD, dplace = 4)

##### DM-ZIDM
clusDZ <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$DZ[[x]]$result$ci_result[-(1:500), ], loss = loss_type[a]))})})})
activeDZ <- sapply(1:2, function(x){apply(clusDZ[[x]], 2, uniqueClus)})
ariDZ <- sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusDZ[[a]][, x], aClus[, x], type = "adjrand")[2])})})
apply(activeDZ, 2, meanSD, dplace = 2) ### Active Cluster
apply(ariDZ, 2, meanSD, dplace = 4)

##### DM-sDM
clusDsD <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$DsD[[x]]$result[-(1:500), ], loss = loss_type[a]))})})})
activeDSD <- sapply(1:2, function(x){apply(clusDsD[[x]], 2, uniqueClus)})
ariDsD <- sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusDsD[[a]][, x], aClus[, x], type = "adjrand")[2])})})
apply(activeDSD, 2, meanSD, dplace = 2) ### Active Cluster
apply(ariDsD, 2, meanSD, dplace = 4)

#### DTMM
clusDTMM <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$DTMM_no_structure[[x]]$result[-(1:500), ], loss = loss_type[a]))})})})
activeDTMM <- sapply(1:2, function(x){apply(clusDTMM[[x]], 2, uniqueClus)})
ariDTMM <- sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusDTMM[[a]][, x], aClus[, x], type = "adjrand")[2])})})
apply(activeDTMM, 2, meanSD, dplace = 2) ### Active Cluster
apply(ariDTMM, 2, meanSD, dplace = 4)

#### Create the boxplot for the adjusted Rand index
fullmod_name <- c("ZIDM-ZIDM", "DM-ZIDM", "DM-sDM", "DTMM")
ari_list <- list(ariZZ, ariDZ, ariDsD, ariDTMM)

registerDoParallel(5)
boxplotDat <- foreach(t = 1:4, .combine = "rbind") %dopar% {
  rbind(data.frame(adjRI = ari_list[[t]][, 1], mod = fullmod_name[t], lossF = loss_type[1]),
        data.frame(adjRI = ari_list[[t]][, 2], mod = fullmod_name[t], lossF = loss_type[2]))
}
stopImplicitCluster()

boxplotDat$mod <- factor(boxplotDat$mod, levels = fullmod_name)
boxplotDat$lossF <- factor(boxplotDat$lossF, levels = loss_type)

ggplot(boxplotDat, aes(x = mod, y = adjRI)) +
  geom_boxplot() +
  facet_grid(. ~ lossF) +
  theme_bw() +
  labs(x = "Model", y = "Adjusted Rand Index", 
       title = "The boxplot of the adjusted Rand index for each models and loss functions.")
### DM-DM
#### Calculate the log-likelihood for each k
registerDoParallel(5)
bicCalc <- foreach(t = 1:nData) %:%
  foreach(k = 1:9) %dopar% {
    clusK <- sapply(1:9, function(x){salso(result_list[["DD"]][[t]][[x]]$result[-(1:500), ])})
    ll <- logmar(z = dat[[t]]$dat, atrisk = matrix(1, ncol = 50, nrow = 50), 
                 beta_mat = matrix(0, ncol = 50, nrow = 10))
    (-2 * sum(ll[cbind(1:50, clusK[, k])])) + (((k + 1) * 50) * log(50))
  }
stopImplicitCluster()

optClusIndex <- apply(sapply(1:20, function(x){bicCalc[[x]]}), 2, which.min)

registerDoParallel(5)
DDresult <- foreach(t = 1:nData, .combine = "rbind") %dopar% {
  clusVI <- as.numeric(salso(result_list[["DD"]][[t]][[optClusIndex[t]]]$result[-(1:500), ]))
  clusBD <- as.numeric(salso(result_list[["DD"]][[t]][[optClusIndex[t]]]$result[-(1:500), ],
                             loss = "binder"))
  c(as.numeric(mclustcomp(clusVI, aClus[, t], type = "adjrand")[2]),
    as.numeric(mclustcomp(clusBD, aClus[, t], type = "adjrand")[2]))
}
stopImplicitCluster()

### Distance-based Methods -----------------------------------------------------
#### PAM: Bray-Curtis
registerDoParallel(5)
pamBC <- foreach(t = 1:nData) %dopar% {
  silBC <- sapply(2:10, function(x){mean(silhouette(pam(bcdist(dat[[t]]$dat), x)$clustering, bcdist(dat[[t]]$dat))[, 3])})
  K <- which.max(silBC) + 1
  list(K = K, clus = pam(bcdist(dat[[t]]$dat), K)$clustering) 
}
stopImplicitCluster()
as.numeric(mclustcomp(pamBC[[1]]$clus, aClus[, 1], type = "adjrand")[2])



pam(dist(dat[[1]]$dat, "aitchison"), 2)


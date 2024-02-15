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
case_name <- "diffindex_3_K_5"
nData <- 20

### Import the data and result -------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

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
  scale_y_continuous(breaks = seq(5, 30, 5), limits = c(0, 30))

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
  scale_y_continuous(breaks = seq(5, 30, 5), limits = c(0, 30))

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
  scale_y_continuous(breaks = seq(5, 30, 5), limits = c(0, 30))

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
clusZZ <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$ZZ[[x]]$result$ci_result, loss = loss_type[a]))})})})
apply(sapply(1:2, function(x){apply(clusZZ[[x]], 2, uniqueClus)}), 2, meanSD, dplace = 2) ### Active Cluster
apply(sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusZZ[[a]][, x], aClus[, x], type = "adjrand")[2])})}), 2, meanSD, dplace = 4)

##### DM-ZIDM
clusDZ <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$DZ[[x]]$result$ci_result, loss = loss_type[a]))})})})
apply(sapply(1:2, function(x){apply(clusDZ[[x]], 2, uniqueClus)}), 2, meanSD, dplace = 2) ### Active Cluster
apply(sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusDZ[[a]][, x], aClus[, x], type = "adjrand")[2])})}), 2, meanSD, dplace = 4)

##### DM-sDM
clusDsD <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$DsD[[x]]$result, loss = loss_type[a]))})})})
apply(sapply(1:2, function(x){apply(clusDsD[[x]], 2, uniqueClus)}), 2, meanSD, dplace = 2) ### Active Cluster
apply(sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusDsD[[a]][, x], aClus[, x], type = "adjrand")[2])})}), 2, meanSD, dplace = 4)

#### DTMM
clusDsD <- suppressWarnings({lapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(salso(result_list$DsD[[x]]$result, loss = loss_type[a]))})})})
apply(sapply(1:2, function(x){apply(clusDsD[[x]], 2, uniqueClus)}), 2, meanSD, dplace = 2) ### Active Cluster
apply(sapply(1:2, function(a){sapply(1:nData, function(x){as.numeric(mclustcomp(clusDsD[[a]][, x], aClus[, x], type = "adjrand")[2])})}), 2, meanSD, dplace = 4)

table(salso(result_list$DTMM_no_structure[[1]]$result), aClus[, 1])
table(salso(t(result_list$DZ[[1]]$result$ci_result)), aClus[, 1])

dim(result_list$DTMM_no_structure[[1]]$result)







# devtools::install_github("YushuShi/MicrobiomeCluster")

library(ape)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(salso)
library(ecodist)
library(foreach)
library(doParallel)
library(cluster)
library(ggplot2)
library(ecodist)
library(coda.base)

### User-defined functions: ----------------------------------------------------
uniqueClus <- function(x){
  length(unique(x))
}

### Global Objects: ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

### Dinh 2015, HIV: ------------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/hiv_dinh_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/hiv_dinh_results.tar"))

#### Metadata
metData <- read.table(paste0(path, "hiv_dinh_results/hiv_dinh.metadata.txt"), sep = "\t", header = TRUE)
metData$Sample_Name_s
dim(metData)

#### OTU Table
otuTab <- read.table(paste0(path, "hiv_dinh_results/RDP/hiv_dinh.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
dim(otuTab)

##### Check the samples: All samples in OTU has the metadata.
str_extract(rownames(otuTab), "[^X]+") %in% metData$Sample_Name_s
metData$Sample_Name_s %in% str_extract(rownames(otuTab), "[^X]+") 

##### Filter OTU table
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)

##### Demographic
table(metData$DiseaseState)
table(metData$Alcohol_s)
table(metData$Education_s)
table(metData$sex_s)

##### PAM
data.frame(k = 2:10,
           Euclidean = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(otuTab), y))[, 3])}),
           BrayCurtis = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(bcdist(otuTab), y))[, 3])}),
           Aitchison = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(otuTab + 1e-20, method = "aitchison"), y))[, 3])})) %>%
  pivot_longer(cols = !k) %>%
  ggplot(aes(x = k, y = value)) +
  geom_line() +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_bw() +
  labs(y = "Average Silhouette", x = "Number of clusters")

##### Cluster: PAM with Euclidean
table(pam(dist(otuTab), 7)$clustering, metData$DiseaseState)
table(pam(dist(otuTab), 7)$clustering, metData$Alcohol_s)
table(pam(dist(otuTab), 7)$clustering, metData$Education_s)
table(pam(dist(otuTab), 7)$clustering, metData$sex_s)

##### Cluster: PAM with Euclidean
table(pam(bcdist(otuTab), 7)$clustering, metData$DiseaseState)
table(pam(bcdist(otuTab), 7)$clustering, metData$Alcohol_s)
table(pam(bcdist(otuTab), 7)$clustering, metData$Education_s)
table(pam(bcdist(otuTab), 7)$clustering, metData$sex_s)

##### Import the result form our model
##### Note: Need to check convergence and try random starting points
result <- readRDS(paste0(path, "Manuscript/Result/Application Data/hiv_dinh_result_split_10_theta_1e2.rds"))
plot(apply(result[[1]]$mod$ci_result, 1, uniqueClus), type = "l")

set.seed(1)
salso(result[[1]]$mod$ci_result[-(1:200), ], loss = "binder")

clusIndex <- as.numeric(salso(result[[1]]$mod$ci_result[-(1:200), ], loss = "binder"))
table(clusIndex, metData$DiseaseState)
table(clusIndex, metData$Alcohol_s)
table(clusIndex, metData$sex_s)

table(clusIndex, paste0(metData$sex_s, ": ", metData$is_gay_s))

data.frame(clus = paste0("Cluster ", clusIndex), val = metData$BMI_s) %>%
  ggplot(aes(x = clus, y = val, group = clus)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "BMI")

data.frame(clus = paste0("Cluster ", clusIndex), val = metData$age_s) %>%
  ggplot(aes(x = clus, y = val, group = clus)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Age")

data.frame(clus = paste0("Cluster ", clusIndex), val = metData$Fat_s) %>%
  ggplot(aes(x = clus, y = val, group = clus)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Fat")


data.frame(clusIndex, val = metData$Fat_s) %>%
  group_by(clusIndex) %>%
  summarise(mean(val), sd(val))

table(pam(dist(otuTab), 7)$clustering, metData$Alcohol_s)
table(clusIndex, metData$Education_s)

### Try: Random Starting Point (Theta = 1)
### Note: Dihn 6/27
set.seed(1)
clusInit <- matrix(0, ncol = 6, nrow = 36)
clusInit[, 3] <- sample(0:2, size = 36, replace = TRUE)
clusInit[, 4] <- sample(0:2, size = 36, replace = TRUE)
clusInit[, 5] <- 0:35
clusInit[, 6] <- 0:35

KmaxVec <- c(10, 10, 20, 20, 36, 36)

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
start_time_GLOBAL <- Sys.time()

result_split_10_theta_1_rn <- foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = KmaxVec[t], nbeta_split = 10,
                      z = as.matrix(otuTab),
                      atrisk_init = matrix(1, nrow = 36, ncol = 1024),
                      beta_init = matrix(0, nrow = KmaxVec[t], ncol = 1024),
                      ci_init = clusInit[, t],
                      theta = 1, mu = 0, s2 = 1e-3,
                      s2_MH = 1e-3, t_thres = 2500, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 5)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  return(list(time = comp_time, mod = mod))
}
stopImplicitCluster()
difftime(Sys.time(), start_time_GLOBAL)


# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(5)
# start_time_GLOBAL <- Sys.time()
# s2Mat <- data.frame(expand.grid(c(1e-3, 1e-4, 1e-5), c(1e-3, 1e-4, 1e-5)))
# 
# result_split_10_r0g_4 <- foreach(t = 1:9) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 3000, Kmax = 10, nbeta_split = 10, 
#                       z = as.matrix(otuTab), 
#                       atrisk_init = matrix(1, nrow = 36, ncol = 1024), 
#                       beta_init = matrix(0, nrow = 10, ncol = 1024), 
#                       ci_init = rep(0, 36), 
#                       theta = 1, mu = 0, s2 = s2Mat[t, 1], 
#                       s2_MH = s2Mat[t, 2], t_thres = 1500, launch_iter = 30, 
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 3)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   return(list(time = comp_time, mod = mod))
# }
# stopImplicitCluster()
# difftime(Sys.time(), start_time_GLOBAL)
# 
# saveRDS(result_split_10_r0g_4, file = paste0(path, "Manuscript/Data/Application Data/hiv_dinh_result_split_10_r0g_4.rds"))
# 
# # set.seed(1)
# # startTime <- Sys.time()
# # modResult <- mod_adaptive(iter = 1000, Kmax = 10, nbeta_split = 10, 
# #                           z = as.matrix(otuTab), 
# #                           atrisk_init = matrix(1, nrow = 36, ncol = 1024), 
# #                           beta_init = matrix(0, nrow = 10, ncol = 1024), 
# #                           ci_init = rep(0, 36), theta = 1, mu = 0, s2 = 1e-3, 
# #                           s2_MH = 1e-5, t_thres = 500, launch_iter = 30, 
# #                           r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 5)
# # totTime <- difftime(Sys.time(), startTime)
# 
# 
# 
# 
# plot(apply(modResult$ci_result, 1, activeClus), type = "l")
# 
# 
# #clus_assign <- data.frame(ID = str_extract(rownames(otuTab), "[^X]+"),
# #                          cluster = as.numeric(salso(modResult$ci_result, loss = "binder")))
# 
# ## This is the result when using 1e-3 for both s2.
# table(metData$sex_s, clus_assign$cluster)




### Papa 2012, IBD: ------------------------------------------------------------
# # file.exists(paste0(path, "Manuscript/Data/Application Data/ibd_alm_results.tar"))
# untar(paste0(path, "Manuscript/Data/Application Data/ibd_alm_results.tar"))
# info <- read.table(paste0(path, "ibd_alm_results/ibd_alm.metadata.txt"), sep = "\t", header = TRUE) %>%
#   dplyr::select(-X)
# info <- data.frame(ID = str_replace(str_extract(info$sample, pattern = "[:digit:]+\\-[:alpha:]"), "\\-", ""),
#                    info %>% dplyr::select(-sample))
# 
# table(info$DiseaseState)
# table(info$ibd)
# table(info$ibd, info$DiseaseState)
# 
# taxa <- read.table(paste0(path, "ibd_alm_results/RDP/ibd_alm.otu_table.100.denovo.rdp_assigned"))
# taxa <- t(taxa)
# taxa <- taxa[, -which(colMeans(taxa > 0) < 0.1)]
# 
# dim(taxa)
# 
# ### K-means
# plot(sapply(2:10, function(x){sqrt(sum(kmeans(dist(taxa), x)$withinss))}), type = "b")
# kmean_clus <- kmeans(dist(taxa), 2)$cluster
# kmean_info <- data.frame(ID = str_extract(names(kmean_clus), "[:digit:]+[:alpha:]"),
#                          clus = kmean_clus) %>%
#   inner_join(info)
# 
# table(kmean_info$clus, kmean_info$ibd)
# table(kmean_info$clus, kmean_info$gender)
# table(kmean_info$clus, kmean_info$DiseaseState)
# 
# ### K-mean with BC dist
# plot(sapply(2:10, function(x){sqrt(sum(kmeans(bcdist(taxa), x)$withinss))}), type = "b")
# kmean_BC_clus <- kmeans(bcdist(taxa), 2)$cluster
# kmean_BC_info <- data.frame(ID = str_extract(names(kmean_clus), "[:digit:]+[:alpha:]"),
#                             clus = kmean_BC_clus) %>%
#   inner_join(info)
# 
# table(kmean_BC_info$clus, kmean_BC_info$ibd)
# table(kmean_BC_info$clus, kmean_BC_info$gender)
# table(kmean_BC_info$clus, kmean_BC_info$DiseaseState)
# table(kmean_info$clus, kmean_info$abx)
# 
# ### Our model
# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(5)
# start_time_GLOBAL <- Sys.time()
# s2Mat <- data.frame(expand.grid(c(1e-3, 1e-4, 1e-5), c(1e-3, 1e-4, 1e-5)))
# 
# result_split_10 <- foreach(t = 1:9) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 10, 
#                       z = as.matrix(taxa), 
#                       atrisk_init = matrix(1, nrow = 91, ncol = 616), 
#                       beta_init = matrix(0, nrow = 10, ncol = 616), 
#                       ci_init = rep(0, 91), 
#                       theta = 1, mu = 0, s2 = s2Mat[t, 1], 
#                       s2_MH = s2Mat[t, 2], t_thres = 2500, launch_iter = 30, 
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 25)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   return(list(time = comp_time, mod = mod))
# }
# stopImplicitCluster()
# difftime(Sys.time(), start_time_GLOBAL)
# 
# result_split_10[[6]]$time/3600
# 
# salso(result_split_10[[9]]$mod$ci_result[-(1:100), ])
# 
# ourMod_clus <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 100, 
#                             z = as.matrix(taxa), 
#                             atrisk_init = matrix(1, nrow = 91, ncol = 616), 
#                             beta_init = matrix(0, nrow = 10, ncol = 616), 
#                             ci_init = rep(0, 91), 
#                             theta = 1, mu = 0, s2 = 1e-3, 
#                             s2_MH = 1e-5, t_thres = 2500, launch_iter = 30, 
#                             r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 25)
# 
# apply(ourMod_clus$ci_result, 1, function(x){length(unique(x))}) %>%
#   plot(type = "l")
# 
# 
# 
# 
# table(ourMod_clus$ci_result[100, ], 
#       ourMod_clus$ci_result[102, ])
# 
# salso(ourMod_clus$ci_result)
# 
# 
# 
# # salso(resultSmits$ci_result, loss = "binder") %>% table()
# 
# colSums(resultSmits$MH_accept == 1)/colSums(resultSmits$MH_accept != -1)
# 
# 
# 
# 
# plot(apply(resultSmits$ci_result, 1, function(x){length(unique(x))}), type = "l")
# 
# table(summary(salso(resultSmits$ci_result, loss = "binder"))$estimate)
# 
# ourMod <- inner_join(data.frame(ID = str_extract(rownames(taxa), pattern = "[:digit:]+[:alpha:]"),
#                                 cluster = as.numeric(salso(resultSmits$ci_result, loss = "binder"))),
#                      data.frame(ID = str_replace(str_extract(info$sample, pattern = "[:digit:]+\\-[:alpha:]"), "\\-", ""),
#                                 info %>% dplyr::select(-sample)))
# 
# table(ourMod$cluster, ourMod$DiseaseState)
# 
# str_replace(str_extract(info$sample, pattern = "[:digit:]+\\-[:alpha:]"), "\\-", "")
# 
# data.frame(ID = rownames(taxa), 
#            cluster = as.numeric(salso(resultSmits$ci_result, loss = "binder")))
# 
# rownames(taxa)
# 
# rownames(info)
# 
# info$sample
# 
# table(info$gender)
# table(info$DiseaseState)
# 
# table(resultSmits$ci_result[99, ], resultSmits$ci_result[100, ])
# 
# ### Smits Dataset: -------------------------------------------------------------
# load(paste0(path, "Manuscript/Data/Application Data/Smits/Smits.RData")) 
# datSmits <- t(Smits$otutable)
# dim(datSmits)
# 
# mean(colMeans(datSmits > 0) < 0.1) ### Non-zero proportion
# datSmitsADJ <- datSmits[, -which(colMeans(datSmits > 0) < 0.9)] ### Filter: not less than 10% non-zero columns
# 
# dim(datSmitsADJ)
# 
# set.seed(1)
# resultSmits <- mod_adaptive(iter = 1000, Kmax = 5, nbeta_split = 10, 
#                             z = as.matrix(datSmitsADJ), 
#                             atrisk_init = matrix(1, nrow = 259, ncol = 47), 
#                             beta_init = matrix(0, nrow = 5, ncol = 47), 
#                             ci_init = rep(0, 259), 
#                             theta = 1, mu = 0, s2 = 1e-5, 
#                             s2_MH = 1e-3, t_thres = 500, launch_iter = 30, 
#                             r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# apply(resultSmits$ci_result, 1, function(x){length(unique(x))})
# 
# 
# 
# table(Smits$sampleinfo, resultSmits$ci_result[12, ])
# 
# table(Smits$sampleinfo, salso(resultSmits$ci_result, loss = "binder"))
# 
# data.frame(LD = colMeans(datSmitsADJ[which(Smits$sampleinfo == "Late Dry"), ]),
#            ED = colMeans(datSmitsADJ[-which(Smits$sampleinfo == "Late Dry"), ]))
# 
# names(Smits$sampleinfo) == rownames(datSmitsADJ)
# 
# 
# table(Smits$sampleinfo)
# datSmits <- datSmits[, -which(colMeans(datSmits == 0) == 1)]
# 
# 
# mean(datSmits == 0) ### Proportion of zero
# colMeans(datSmits == 0) %>%
#   hist()
# 
# prop_zero <- seq(0, 1, 0.01)
# zeroOVA <- rep(NA, length(prop_zero))
# propCol <- rep(NA, length(prop_zero))
# 
# for(i in 1:length(prop_zero)){
#   colInt <- colMeans(datSmits == 0) < prop_zero[i]
#   propCol[i] <- mean(colInt)
#   zeroOVA[i] <- mean(datSmits[, colInt] == 0)
# }
# 
# data.frame(prop_zero, zeroOVA, propCol) %>%
#   pivot_longer(!prop_zero) %>%
#   dplyr::mutate(name = factor(name, levels = c("zeroOVA", "propCol"),
#                               labels = c("Overall zero", "Columns included"))) %>%
#   ggplot(aes(x = prop_zero, y = value, color = name)) +
#   geom_line() +
#   theme_bw() +
#   labs(x = "The maximum proportion of zero for each columns allowed into the analysis",
#        y = "Proportion",
#        title = "Smits: Taxa selection") +
#   theme(legend.position = "bottom", legend.title = element_blank())
# 
# propCol[which(prop_zero == 0.9)] * dim(datSmits)
# 
# relaTab <- datSmits[, colMeans(datSmits == 0) < 0.75]/rowSums(datSmits[, colMeans(datSmits == 0) < 0.75])
# data.frame(x = colnames(relaTab), y = relaTab[1, ], index = 1) %>%
#   `rownames<-`(NULL) %>%
#   ggplot(aes(x = index, y = y, fill = x)) +
#   geom_bar(position = "fill", stat = "identity")
# 
# relaTab[, which(colMeans(relaTab) > 0.01)] %>%
#   data.frame() %>%
#   dplyr::mutate(Other = 1 - rowSums(relaTab[, which(colMeans(relaTab) > 0.01)]),
#                 Sample = rownames(relaTab[, which(colMeans(relaTab) > 0.01)])) %>%
#   pivot_longer(!Sample) %>%
#   ggplot(aes(x = factor(Sample), y = value, fill = factor(name))) +
#   geom_bar(position = "fill", stat = "identity")
# 
# 1 - rowSums(relaTab[, which(colMeans(relaTab) > 0.001)])
# 
# hist(colMeans(relaTab))
# 
# which()
# 
# propCol[prop_zero == 0.75] * 12000
# table(kmeans(dist(datSmits[, colMeans(datSmits == 0) < 0.75]), 2)$cluster, Smits$sampleinfo)
# 
# data.frame(x = names(datSmits[1, colMeans(datSmits == 0) < 0.75]),
#            y = datSmits[1, colMeans(datSmits == 0) < 0.75]) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_point()
# 
# colMeans(datSmits == 0) < 0.5
# sum(colMeans(datSmits == 0) < 0.75)
# datSmits[, as.numeric(which(colMeans(datSmits == 0) < 0.75))] %>% View()
# 
# # t1 <- read.csv(paste0(path, "Manuscript/Data/Application Data/Smits/aan4834_table_s4.csv"))
# # t2 <- read.csv(paste0(path, "Manuscript/Data/Application Data/Smits/aan4834_table_s3.csv"))
# 
# tab <- read.csv(paste0(path, "Manuscript/Data/Application Data/Smits/aan4834_table_s1.csv"))
# tab$X.SampleID
# 
# str_extract(names(Smits$sampleinfo), "^[:alpha:]+") %>%
#   as.factor() %>%
#   table()
# 
# names(Smits$sampleinfo) %in% t1$Sample
# 
# colnames(t1)
# colnames(t2)
# 
# MicrobiomeCluster::convCounts2Abun()
# 
# # aan4834_table_s4.csv
# 
# Smits$sampleinfo
# 
# Smits$otutable
# 

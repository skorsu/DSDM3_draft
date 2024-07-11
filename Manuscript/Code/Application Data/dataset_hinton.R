# Required Libraries
library(selbal) # Dataset
library(tidyverse)
library(foreach)
library(doParallel)
library(salso)
library(ggplot2)
library(mclustcomp)
library(latex2exp)
library(RColorBrewer)
library(viridis)
library(ggpubr)

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

# Path
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}
resultpath <- paste0(path, "Result/selbal_HIV/")
# file.exists(resultpath)

# HIV dataset
datHIV <- HIV
sexHIV <- factor(datHIV[, 61], labels = c("non-MSM", "MSM"))
statusHIV <- factor(datHIV[, 62], labels = c("Healthy", "HIV-patient"))
otuHIV <- datHIV[, 1:60]

# Sequencing Dept
median(rowSums(otuHIV))
min(rowSums(otuHIV))
max(rowSums(otuHIV))

which(rowSums(otuHIV) == min(rowSums(otuHIV)))
otuHIV[30, ]

ggplot(data.frame(x = rowSums(otuHIV)), aes(x = x)) +
  geom_histogram() +
  labs(x = "Total Read Count", title = "HIV dataset: Distribution of the Total Read Count") +
  theme_bw() +
  theme(axis.title.y = element_blank())

# Shannon Diversity and Simpson Index
shannon_d <- sapply(1:155, function(x){
  pi <- otuHIV[x, otuHIV[x, ] != 0]/sum(otuHIV[x, otuHIV[x, ] != 0])
  -sum(pi * log(pi))})

simpson_d <- sapply(1:155, function(x){
  pi <- otuHIV[x, otuHIV[x, ] != 0]/sum(otuHIV[x, otuHIV[x, ] != 0])
  1 - sum(pi^2)})

statusSex <- paste0(statusHIV, " \n ", sexHIV)

numericLonger <- data.frame(trc = rowSums(otuHIV), shannon_d, simpson_d) %>%
  mutate(ID = 1:155, statusSex) %>%
  pivot_longer(!c(ID, statusSex)) 

numericLonger$name <- factor(numericLonger$name, labels = c("Shannon", "Simpson", "Total Read Count"))

ggplot(numericLonger, aes(x = statusSex, y = value)) +
  geom_boxplot() +
  facet_wrap(. ~ name, scales = "free_y")

table(statusSex)

data.frame(trc = rowSums(otuHIV), shannon_d, simpson_d) %>%
  mutate(ID = 1:155, statusSex) %>%
  group_by(statusSex) %>%
  summarise(n(), median(trc), min(trc), max(trc), median(shannon_d), min(shannon_d), 
            max(shannon_d), median(simpson_d), min(simpson_d), 
            max(simpson_d)) %>% View()

data.frame(trc = rowSums(otuHIV), shannon_d, simpson_d) %>%
  mutate(ID = 1:155, statusSex) %>%
  summarise(n(), median(trc), min(trc), max(trc), median(shannon_d), min(shannon_d), 
            max(shannon_d), median(simpson_d), min(simpson_d), 
            max(simpson_d)) %>% View()

which.min(shannon_d)
which.min(simpson_d)
which.max(rowSums(otuHIV))
which.max(shannon_d)
which.max(simpson_d)

# Demographic
data.frame(table(statusHIV, sexHIV)) %>%
  ggplot(aes(x = sexHIV, y = Freq)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(. ~ statusHIV) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.9), vjust = -0.5) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(title = "Demographic: HIV dataset from selbal packages",
       y = "Frequency")

## Adjust model based on meeting (7/10): ---------------------------------------
set.seed(1)
ciInit <- matrix(0, nrow = 155, ncol = 12)
ciInit[, 4] <- sample(0:2, 155, replace = TRUE)
ciInit[, 5] <- sample(0:2, 155, replace = TRUE)
ciInit[, 6] <- sample(0:2, 155, replace = TRUE)
ciInit[, 7] <- sample(0:4, 155, replace = TRUE)
ciInit[, 8] <- sample(0:4, 155, replace = TRUE)
ciInit[, 9] <- sample(0:4, 155, replace = TRUE)
ciInit[, 10] <- sample(0:19, 155, replace = TRUE)
ciInit[, 11] <- sample(0:19, 155, replace = TRUE)
ciInit[, 12] <- sample(0:19, 155, replace = TRUE)

xiInitDum <- lapply(4:12, function(y){sapply(0:max(ciInit[, y]), function(x){
  p <- colSums(otuHIV[which(ciInit[, y] == x), ])/sum(otuHIV[which(ciInit[, y] == x), ])
  ifelse(is.infinite(log(p/(1-p))), -20, log(p/(1-p)))
}) %>% t()
})

xiInit <- vector("list", 12)
xiInit[[1]] <- matrix(0, nrow = 20, ncol = 60)
xiInit[[2]] <- matrix(0, nrow = 20, ncol = 60)
xiInit[[3]] <- matrix(0, nrow = 20, ncol = 60)

xiInit[[4]] <- rbind(xiInitDum[[1]], matrix(0, nrow = 17, ncol = 60))
xiInit[[5]] <- rbind(xiInitDum[[2]], matrix(0, nrow = 17, ncol = 60))
xiInit[[6]] <- rbind(xiInitDum[[3]], matrix(0, nrow = 17, ncol = 60))

xiInit[[7]] <- rbind(xiInitDum[[4]], matrix(0, nrow = 15, ncol = 60))
xiInit[[8]] <- rbind(xiInitDum[[5]], matrix(0, nrow = 15, ncol = 60))
xiInit[[9]] <- rbind(xiInitDum[[6]], matrix(0, nrow = 15, ncol = 60))

xiInit[[10]] <- xiInitDum[[7]]
xiInit[[11]] <- xiInitDum[[8]]
xiInit[[12]] <- xiInitDum[[9]]

resultName <- c(paste0("result_selbal_HIV_chain_", 1:3, "_init_oneClus_JUL10_fixed.rds"),
                paste0("result_selbal_HIV_chain_", 1:3, "_init_3clus_JUL10_fixed.rds"),
                paste0("result_selbal_HIV_chain_", 1:3, "_init_5clus_JUL10_fixed.rds"),
                paste0("result_selbal_HIV_chain_", 1:3, "_init_20clus_JUL10_fixed.rds"))

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:12) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 20, nbeta_split = 5,
                      z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60),
                      beta_init = as.matrix(xiInit[[t]]),
                      ci_init = ciInit[, t],
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
                      t_thres = 2500, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

# ### Default set of Hyperparameters: 1e-3
# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# globalTime <- Sys.time()
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5,
#                       z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60),
#                       beta_init = matrix(0, nrow = 10, ncol = 60),
#                       ci_init = rep(0, 155),
#                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                       t_thres = 2500, launch_iter = 30,
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_oneClus_defaultHyper.rds"))
# }
# stopImplicitCluster()
# difftime(Sys.time(), globalTime)
# 
# ### Different starting point: 1e-3 - PART 1
# set.seed(1)
# ciInit <- matrix(NA, nrow = 155, ncol = 6)
# ciInit[, 1] <- sample(0:4, 155, replace = TRUE)
# ciInit[, 2] <- sample(0:4, 155, replace = TRUE)
# ciInit[, 3] <- sample(0:4, 155, replace = TRUE)
# ciInit[, 4] <- sample(0:19, 155, replace = TRUE)
# ciInit[, 5] <- sample(0:19, 155, replace = TRUE)
# ciInit[, 6] <- sample(0:19, 155, replace = TRUE)
# 
# KmaxVec <- c(20, 20, 20, 50, 50, 50)
# 
# xiInit <- lapply(1:6, function(y){sapply(0:max(ciInit[, y]), function(x){
#   colSums(otuHIV[which(ciInit[, y] == x), ])/sum(otuHIV[which(ciInit[, y] == x), ])
# }) %>% t()
# })
# 
# xiInit[[1]] <- rbind(xiInit[[1]], matrix(0, nrow = 15, ncol = 60))
# xiInit[[2]] <- rbind(xiInit[[2]], matrix(0, nrow = 15, ncol = 60))
# xiInit[[3]] <- rbind(xiInit[[3]], matrix(0, nrow = 15, ncol = 60))
# xiInit[[4]] <- rbind(xiInit[[4]], matrix(0, nrow = 30, ncol = 60))
# xiInit[[5]] <- rbind(xiInit[[5]], matrix(0, nrow = 30, ncol = 60))
# xiInit[[6]] <- rbind(xiInit[[6]], matrix(0, nrow = 30, ncol = 60))
# 
# resultName <- c(paste0("result_selbal_HIV_chain_", 1:3, "_init_5clus_Kmax_20_defaultHyper.rds"),
#                 paste0("result_selbal_HIV_chain_", 1:3, "_init_20clus_Kmax_50_defaultHyper.rds"))
# 
# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# globalTime <- Sys.time()
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 25000, Kmax = KmaxVec[t], nbeta_split = 5,
#                       z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60),
#                       beta_init = as.matrix(xiInit[[t]]),
#                       ci_init = ciInit[, t],
#                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                       t_thres = 2500, launch_iter = 30,
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
# }
# stopImplicitCluster()
# difftime(Sys.time(), globalTime)
# 
# ### Different starting point: 1e-3 - PART 2
# set.seed(1)
# ciInit <- matrix(NA, nrow = 155, ncol = 3)
# ciInit[, 1] <- sample(0:2, 155, replace = TRUE)
# ciInit[, 2] <- sample(0:2, 155, replace = TRUE)
# ciInit[, 3] <- sample(0:2, 155, replace = TRUE)
# 
# xiInit <- lapply(1:3, function(y){sapply(0:max(ciInit[, y]), function(x){
#   colSums(otuHIV[which(ciInit[, y] == x), ])/sum(otuHIV[which(ciInit[, y] == x), ])
# }) %>% t()
# })
# 
# xiInit[[1]] <- rbind(xiInit[[1]], matrix(0, nrow = 12, ncol = 60))
# xiInit[[2]] <- rbind(xiInit[[2]], matrix(0, nrow = 12, ncol = 60))
# xiInit[[3]] <- rbind(xiInit[[3]], matrix(0, nrow = 12, ncol = 60))
# 
# resultName <- paste0("result_selbal_HIV_chain_", 1:3, "_init_3clus_Kmax_15_defaultHyper.rds")
# 
# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(3)
# globalTime <- Sys.time()
# foreach(t = 1:3) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 25000, Kmax = 15, nbeta_split = 5,
#                       z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60),
#                       beta_init = as.matrix(xiInit[[t]]),
#                       ci_init = ciInit[, t],
#                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                       t_thres = 2500, launch_iter = 30,
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
# }
# stopImplicitCluster()
# difftime(Sys.time(), globalTime)

### Post Analysis: -------------------------------------------------------------
#### Read the result
# resultFilename <- c(paste0(resultpath, "result_selbal_HIV_chain_", 1:6, "_init_oneClus_defaultHyper.rds"),
#                     paste0(resultpath, "result_selbal_HIV_chain_", 1:6, "_init_oneClus_s2MH_1en5.rds"))

resultFilename <- c(paste0(resultpath, "result_selbal_HIV_chain_", 1:3, "_init_oneClus_JUL10_fixed.rds"),
                    paste0(resultpath, "result_selbal_HIV_chain_", 1:3, "_init_3clus_JUL10_fixed.rds"),
                    paste0(resultpath, "result_selbal_HIV_chain_", 1:3, "_init_5clus_JUL10_fixed.rds"),
                    paste0(resultpath, "result_selbal_HIV_chain_", 1:3, "_init_20clus_JUL10_fixed.rds"))

### Computational time
registerDoParallel(6)
compTime <- foreach(t = 1:12, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  as.numeric(result$time)
}
stopImplicitCluster()

mean(compTime/3600)
sd(compTime/3600)

### Active Cluster
registerDoParallel(6)
activeClusMat <- foreach(t = 1:12, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  apply(result$mod$ci_result, 1, uniqueClus)
}
stopImplicitCluster()

activeClusMatPlot <- activeClusMat %>%
  as.data.frame() %>%
  mutate(iter = 1:25000) %>%
  pivot_longer(!iter)

activeClusMatPlot$name <- factor(activeClusMatPlot$name, 
                                 levels = paste0("result.", 1:12), labels = paste0("Chain ", 1:12))

ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusteres",
       color = "MCMC Chain") +
  guides(color = guide_legend(ncol = 12))

### Check the convergence of xi
data.frame(colMeans(otuHIV), apply(otuHIV, 2, var))

registerDoParallel(6)
xiFirst <- foreach(t = 1:12) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$beta_result[, 1, ] %>% t() %>%
    as.data.frame() %>%
    mutate(Iteration = 1:25000) %>%
    pivot_longer(!Iteration) %>%
    transmute(Iteration, xi = value,
              Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
}
stopImplicitCluster()

xiFirstLong <- lapply(1:12, function(x){data.frame(xiFirst[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows()

xiFirstLong$chain <- factor(xiFirstLong$chain, levels = paste0("Chain ", 1:12))

ggplot(xiFirstLong, aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k1}$")), x = "Iteration", y = TeX("\\xi"))

registerDoParallel(6)
xiThird <- foreach(t = 1:12) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$beta_result[, 3, ] %>% t() %>%
    as.data.frame() %>%
    mutate(Iteration = 1:25000) %>%
    pivot_longer(!Iteration) %>%
    transmute(Iteration, xi = value,
              Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
}
stopImplicitCluster()

xiThirdLong <- lapply(1:12, function(x){data.frame(xiThird[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows()

xiThirdLong$chain <- factor(xiThirdLong$chain, levels = paste0("Chain ", 1:12))

ggplot(xiThirdLong, aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k3}$")), x = "Iteration", y = TeX("\\xi"))

### Check the acceptance rate
# KmaxVec <- c(10, 10, 10, 15, 15, 15, 20, 20, 20, 50, 50, 50)
# registerDoParallel(6)
# xiAcceptRXN <- foreach(t = 1:12) %dopar% {
#   result <- readRDS(resultFilename[t])
#   sapply(1:KmaxVec[t], function(x){sum(result$mod$MH_accept[, x] == 1)/sum(result$mod$MH_accept[, x] != -1)})
# }
# stopImplicitCluster()
# 
# registerDoParallel(6)
# xiAccept <- foreach(t = 1:12) %dopar% {
#   result <- readRDS(resultFilename[t])
#   sapply(1:KmaxVec[t], function(x){sum(result$mod$MH_accept[, x] != -1)})
# }
# stopImplicitCluster()
# 
# colorHM <- matrix(NA, nrow = 50, ncol = 12)
# labelHM <- matrix(NA, nrow = 50, ncol = 12)
# 
# for(i in 1:12){
#   
#   for(j in 1:KmaxVec[i]){
#     
#     colorHM[j, i] <- xiAcceptRXN[[i]][j]
#     labelHM[j, i] <- paste0(round(xiAcceptRXN[[i]][j], 2), " (", xiAccept[[i]][j], ")")
#     
#   }
#   
# }
# 
# labelHM[labelHM == "NaN(0)"] <- NA
# 
# xiAcceptLab <- as.data.frame(colorHM) %>%
#   mutate(Component = paste0("Component ", 1:50)) %>%
#   pivot_longer(!Component, values_to = "AcceptRate") %>%
#   inner_join(as.data.frame(labelHM) %>%
#                mutate(Component = paste0("Component ", 1:50)) %>%
#                pivot_longer(!Component, values_to = "RateLabel"))
#   
# xiAcceptLab$name <- factor(xiAcceptLab$name, levels = paste0("V", 1:12), labels = paste0("Chain ", 1:12))
# xiAcceptLab$RateLabel[xiAcceptLab$RateLabel == "NaN (0)"] <- NA
# 
# xiAcceptLab %>%
#   ggplot(aes(x = name, y = Component, fill = AcceptRate)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red") +
#   geom_text(aes(label = RateLabel), color = "black", size = 2) +
#   theme_bw() +
#   theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
#         legend.position = "bottom") +
#   labs(x = "MCMC Chain", y = "Component", fill = "Acceptance Rate",
#        title = TeX(paste0("Acceptance Rate of the Adpative MH for ", "$\\textbf{\\xi}_{k}$")))

### Cluster Assignment - Individual
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
clusSALSO <- foreach(t = 1:12, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  as.numeric(salso(result$mod$ci_result[-(1:5000), ]))
}
stopImplicitCluster()

#### Individual - Demographic
salsoDemo <- lapply(1:12, function(x){data.frame(table(clusSALSO[, x], statusHIV, sexHIV)) %>%
    mutate(status_sex = paste0(statusHIV, ": ", sexHIV),
           Chain = paste0("Chain ", x), 
           Cluster = paste0("Cluster ", Var1))}) %>%
  bind_rows() 

salsoDemo$Cluster <- factor(salsoDemo$Cluster, levels = paste0("Cluster ", 1:max(clusSALSO)))
salsoDemo$Chain <- factor(salsoDemo$Chain, levels = paste0("Chain ", 1:12))

salsoDemo %>%
  ggplot(aes(x = Cluster, y = Freq, fill = status_sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.5), vjust = -0.15, size = 2) + 
  scale_fill_manual(values = c("springgreen3", "springgreen4", "coral1", "coral3")) +
  facet_wrap(. ~ Chain, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4)) +
  labs(y = "Frequency", fill = "Demographic", title = "Demographic by Cluster for every MCMC chain.")

#### Individual - Relative Abundance
taxaName <- data.frame(c = str_remove(str_extract(colnames(otuHIV), "c_[:alpha:]+"), "c_"),
                       o = str_remove(str_extract(colnames(otuHIV), "o_[:alpha:]+"), "o_"),
                       f = str_remove(str_extract(colnames(otuHIV), "f_[:alpha:]+"), "f_"),
                       g = str_remove(str_extract(colnames(otuHIV), "g_[:alpha:]+"), "g_"))

highTaxa <- sapply(1:3, function(x){
(otuHIV[which(clusSALSO[, 1] == x), ]/rowSums(otuHIV[which(clusSALSO[, 1] == x), ])) %>%
  colMeans() %>%
  sort(decreasing = TRUE, index.return = TRUE) %>%
  .$ix %>%
  .[1:10]})

otuHIVvisual <- otuHIV/rowSums(otuHIV)
highTaxaIndex <- union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])
highTaxaIndexLabel <- ifelse(taxaName$g %in% c("unclassified", "Incertae"), 
       ifelse(is.na(taxaName$f), 
              ifelse(is.na(taxaName$o), taxaName$c, taxaName$o), taxaName$f), taxaName$g) 
highTaxaIndexLabel[12] <- "Lachnospiraceae - Incertae"
highTaxaIndexLabel <- highTaxaIndexLabel %>% .[union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])]
highTaxaIndexLabel[is.na(highTaxaIndexLabel)] <- "Unclassified"

otuRela <- matrix(NA, ncol = length(highTaxaIndex) + 1, nrow = 155)
for(i in 1:length(highTaxaIndex)){
  otuRela[, i] <- otuHIVvisual[, highTaxaIndex[i]]
}
otuRela[, length(highTaxaIndex) + 1] <- 1 - rowSums(otuRela[, 1:length(highTaxaIndex)])
colnames(otuRela) <- c(highTaxaIndexLabel, "Others")


otuRelaPlot <- otuRela %>%
  as.data.frame() %>% 
  mutate(ID = str_extract(rownames(otuHIV), "[:digit:]+"), cluster = paste0("Cluster ", clusSALSO[, 1])) %>%
  pivot_longer(!c(ID, cluster))

otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(highTaxaIndexLabel, "Others"))
# otuRelaPlot$ID <- factor(otuRelaPlot$ID, levels = c(rownames(otuHIVvisual)[which(clusSALSO[, 1] == 1)][sort(otuHIVvisual[which(clusSALSO[, 1] == 1), 3], decreasing = TRUE, index.return = TRUE)$ix],
#                                                     rownames(otuHIVvisual)[which(clusSALSO[, 1] == 2)][sort(otuHIVvisual[which(clusSALSO[, 1] == 2), 6], decreasing = TRUE, index.return = TRUE)$ix],
#                                                     rownames(otuHIVvisual)[which(clusSALSO[, 1] == 3)][sort(otuHIVvisual[which(clusSALSO[, 1] == 3), 1], decreasing = TRUE, index.return = TRUE)$ix]))
# 
# rownames(otuHIVvisual)[which(clusSALSO[, 1] == 1)][sort(otuHIVvisual[which(clusSALSO[, 1] == 1), 3], decreasing = TRUE, index.return = TRUE)$ix]
# rownames(otuHIVvisual)[which(clusSALSO[, 1] == 2)][sort(otuHIVvisual[which(clusSALSO[, 1] == 2), 6], decreasing = TRUE, index.return = TRUE)$ix]
# rownames(otuHIVvisual)[which(clusSALSO[, 1] == 3)][sort(otuHIVvisual[which(clusSALSO[, 1] == 3), 1], decreasing = TRUE, index.return = TRUE)$ix]

ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(. ~ cluster, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(values = c(turbo(n = 16), "gray90")) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill = "Taxa", y = "Relative Abundance",
       title = paste0("Relative Abundance for each cluster"))

# relaPlot <- lapply(1:12, function(x){
#   otuRelaPlot <- otuRela %>%
#     as.data.frame() %>% 
#     mutate(ID = str_extract(rownames(otuHIV), "[:digit:]+"), cluster = paste0("Cluster ", clusSALSO[, x])) %>%
#     pivot_longer(!c(ID, cluster))
#   otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(highTaxaIndexLabel, "Others"))
#   ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
#     geom_bar(position = "stack", stat = "identity") +
#     facet_grid(. ~ cluster, scales = "free_x") +
#     theme_bw() +
#     scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
#     scale_y_continuous(labels = scales::percent) +
#     theme(axis.text.x = element_text(angle = 90)) +
#     theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#     guides(fill = guide_legend(nrow = 2)) +
#     labs(fill = "Taxa", y = "Relative Abundance",
#          title = paste0("Chain ", x, " - Relative Abundance for each cluster"))
# })
# 
# ggarrange(plotlist = relaPlot, common.legend = TRUE, legend = "bottom")

### Combine all chains
registerDoParallel(6)
MCMCCombine <- foreach(t = 1:12, .combine = rbind) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$ci_result[15001:25000, ]
}
stopImplicitCluster()

set.seed(1)
clusComb <- as.numeric(salso(MCMCCombine))

#### Demographic - Combined Chain
salsoCombDemo <- data.frame(table(clusComb, statusHIV, sexHIV)) %>%
  mutate(status_sex = paste0(statusHIV, ": ", sexHIV), 
         Cluster = paste0("Cluster ", clusComb))
salsoCombDemo$Cluster <- factor(salsoCombDemo$Cluster, levels = paste0("Cluster ", 1:3))

salsoCombDemo %>%
  ggplot(aes(x = Cluster, y = Freq, fill = status_sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.5), vjust = -0.5) + 
  scale_fill_manual(values = c("springgreen3", "springgreen4", "coral1", "coral3")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 4)) +
  labs(y = "Frequency", fill = "Demographic", title = "Demographic by Cluster for the combined MCMC chain.")

#### Relative Abundance - Combined Chain
otuRelaCombPlot <- otuRela %>%
  as.data.frame() %>% 
  mutate(ID = str_extract(rownames(otuHIV), "[:digit:]+"), cluster = paste0("Cluster ", clusComb)) %>%
  pivot_longer(!c(ID, cluster))

otuRelaCombPlot$name <- factor(otuRelaCombPlot$name, levels = c(highTaxaIndexLabel, "Others"))

ggplot(otuRelaCombPlot, aes(x = ID, y = value, fill = name)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(. ~ cluster, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(values = c(turbo(n = 16), "gray90")) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill = "Taxa", y = "Relative Abundance",
       title = paste0("Relative Abundance for each cluster of the combined MCMC chain"))

data.frame(trc = rowSums(otuHIV), shannon_d, simpson_d) %>%
  mutate(ID = 1:155, clusComb) %>%
  group_by(clusComb) %>%
  summarise(n(), median(trc), min(trc), max(trc), median(shannon_d), min(shannon_d), 
            max(shannon_d), median(simpson_d), min(simpson_d), 
            max(simpson_d)) %>% View()

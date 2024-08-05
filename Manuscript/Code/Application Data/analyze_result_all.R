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

### Dataset: -------------------------------------------------------------------
#### HIV dataset
datHIV <- HIV
sexHIV <- factor(datHIV[, 61], labels = c("non-MSM", "MSM"))
statusHIV <- factor(datHIV[, 62], labels = c("Healthy", "HIV-patient"))
otuHIV <- datHIV[, 1:60]

mean(otuHIV == 0)

#### EDD
##### Metadata
metadata <- read.delim(paste0(path, "Data/Application Data/singh/edd_singh.metadata.txt"))

##### OTU Tables
otuEDD_genus <- readRDS(paste0(path, "Data/Application Data/singh/clean_singh.rds"))
otuEDD_species <- readRDS(paste0(path, "Data/Application Data/singh/clean_singh_species.rds"))
sum(rownames(otuEDD_genus) %in% rownames(otuEDD_species))
sum(rownames(otuEDD_species) %in% rownames(otuEDD_genus))
metadata <- metadata[metadata$SampleID %in% rownames(otuEDD_genus), ]

### Basic Summary for each dataset (Based on Shi 2022 Paper): ------------------
#### HIV Dataset
dim(otuHIV)
sexstatusHIV <- paste0(statusHIV, ": ", sexHIV)
table(sexstatusHIV)

seqDepthHIV <- rowSums(otuHIV)
shannonHIV <- sapply(1:155, function(x){
  p <- otuHIV[x, which(otuHIV[x, ] != 0)]/sum(otuHIV[x, which(otuHIV[x, ] != 0)])
  -sum(p * log(p))})

data.frame(sexstatusHIV, seqDepthHIV, shannonHIV) %>%
  summarise(mean(seqDepthHIV), sd(seqDepthHIV), mean(shannonHIV), sd(shannonHIV)) %>%
  round(4)

data.frame(sexstatusHIV, seqDepthHIV, shannonHIV) %>%
  group_by(sexstatusHIV) %>%
  summarise(round(mean(seqDepthHIV), 2), 
            round(sd(seqDepthHIV), 2), 
            round(mean(shannonHIV), 4), 
            round(sd(shannonHIV), 4)) %>%
  as.data.frame()

### Analysis: HIV Dataset ------------------------------------------------------
resultFilename <- c(paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1:2, "_init_oneClus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_3clus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_5clus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_20clus_JUL10_fixed.rds"))

### Computational time
# registerDoParallel(6)
# compTime <- foreach(t = 1:12, .combine = cbind) %dopar% {
#   result <- readRDS(resultFilename[t])
#   as.numeric(result$time)
# }
# stopImplicitCluster()
# 
# mean(compTime/3600)
# sd(compTime/3600)

### Active Cluster
registerDoParallel(5)
activeClusMat <- foreach(t = 1:5, .combine = cbind) %dopar% {
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
  labs(x = "Iteration", y = "Number of the active clusters", color = "MCMC Chain") +
  guides(color = guide_legend(ncol = 5)) +
  scale_y_continuous(breaks = seq(0, 12, 1))

### Check the convergence of xi
data.frame(colMeans(otuHIV), apply(otuHIV, 2, var))

registerDoParallel(5)
xiFirst <- foreach(t = 1:5) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$beta_result[, 1, ] %>% t() %>%
    as.data.frame() %>%
    mutate(Iteration = 1:25000) %>%
    pivot_longer(!Iteration) %>%
    transmute(Iteration, xi = value,
              Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
}
stopImplicitCluster()

xiFirstLong <- lapply(1:5, function(x){data.frame(xiFirst[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows()

xiFirstLong$chain <- factor(xiFirstLong$chain, levels = paste0("Chain ", 1:5))

ggplot(xiFirstLong, aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Iteration", y = TeX("$\\xi_{k1}$"))

registerDoParallel(5)
xiThird <- foreach(t = 1:5) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$beta_result[, 3, ] %>% t() %>%
    as.data.frame() %>%
    mutate(Iteration = 1:25000) %>%
    pivot_longer(!Iteration) %>%
    transmute(Iteration, xi = value,
              Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
}
stopImplicitCluster()

xiThirdLong <- lapply(1:5, function(x){data.frame(xiThird[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows()

xiThirdLong$chain <- factor(xiThirdLong$chain, levels = paste0("Chain ", 1:12))

ggplot(xiThirdLong, aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Iteration", y = TeX("$\\xi_{k3}$"))

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
registerDoParallel(5)
clusSALSO <- foreach(t = 1:5, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  as.numeric(salso(result$mod$ci_result[-(1:5000), ]))
}
stopImplicitCluster()

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
registerDoParallel(5)
MCMCCombine <- foreach(t = 1:5, .combine = rbind) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$ci_result[15001:25000, ]
}
stopImplicitCluster()

set.seed(1)
clusComb <- as.numeric(salso(MCMCCombine))

table(clusComb)

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
  labs(y = "Frequency", fill = " ")

#### Relative Abundance - Combined Chain
otuRelaCombPlot <- otuRela %>%
  as.data.frame() %>% 
  mutate(ID = str_extract(rownames(otuHIV), "[:digit:]+"), cluster = paste0("Cluster ", clusComb)) %>%
  pivot_longer(!c(ID, cluster))

otuRelaCombPlot$name <- factor(otuRelaCombPlot$name, levels = c(highTaxaIndexLabel, "Others"))

ggplot(otuRelaCombPlot, aes(x = ID, y = value, fill = name)) +
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.25) +
  facet_wrap(. ~ cluster, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(values = c(turbo(n = 16), "gray90")) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill = "Taxa", y = "Relative Abundance")

### Mao
lapply(1:3, function(y){
  
  meanOTU_OVA <- as.numeric(colMeans(otuHIV))
  meanOTU_j <- as.numeric(colMeans(otuHIV[which(clusComb == y), ]))
  meanOTU_not_j <- as.numeric(colMeans(otuHIV[which(clusComb != y), ]))
  first_bottom <- (otuHIV[which(clusComb == y), ] - meanOTU_j)^2  
  second_bottom <- (otuHIV[which(clusComb != y), ] - meanOTU_not_j)^2
  top <- (clusComb[y] * ((meanOTU_j - meanOTU_OVA) ^ 2)) + ((155 - clusComb[y]) * ((meanOTU_not_j - meanOTU_OVA) ^ 2))
  
  data.frame(sapply(1:60, function(x){top[x]/(sum(first_bottom[, x]) + sum(second_bottom[, x]))}))
  
}) %>%
  bind_cols() %>%
  `colnames<-`(paste0("Cluster ", 1:3)) %>%
  mutate(Taxa = colnames(otuHIV)) %>%
  pivot_longer(!Taxa) %>%
  ggplot(aes(x = Taxa, y = name, fill = value)) +
  geom_tile()




meanOTU_OVA <- as.numeric(colMeans(otuHIV))
meanOTU_j <- as.numeric(colMeans(otuHIV[which(clusComb == 1), ]))
meanOTU_not_j <- as.numeric(colMeans(otuHIV[which(clusComb != 1), ]))
first_bottom <- (otuHIV[which(clusComb == 1), ] - meanOTU_j)^2  
second_bottom <- (otuHIV[which(clusComb != 1), ] - meanOTU_not_j)^2
top <- (clusComb[1] * ((meanOTU_j - meanOTU_OVA) ^ 2)) + ((155 - clusComb[1]) * ((meanOTU_not_j - meanOTU_OVA) ^ 2))

sapply(1:60, function(x){top[x]/(sum(first_bottom[, x]) + sum(second_bottom[, x]))})

dim(otuHIV)

as.numeric(colMeans(otuHIV[which(clusComb == 1), ]) - meanOTU)^2
as.numeric(colMeans(otuHIV[which(clusComb != 1), ]) - meanOTU)^2



### Richness and Shannon
data.frame(richness = rowSums(otuHIV > 0), clusComb) %>%
  group_by(clusComb) %>%
  summarise(mean(richness), sd(richness))

relaHIV <- otuHIV/rowSums(otuHIV)
data.frame(shannon = apply(relaHIV, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))}), clusComb) %>%
  group_by(clusComb) %>%
  summarise(mean(shannon), sd(shannon))

data.frame(Cluster = paste0("Cluster ", clusComb), 
           Richness = rowSums(otuHIV > 0), 
           Shannon = apply(relaHIV, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})) %>%
  pivot_longer(!Cluster) %>%
  ggplot(aes(y = value, x = Cluster)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(. ~ name, scales = "free_y") +
  labs(y = "") +
  theme_bw()

kruskal.test(as.numeric(rowSums(otuHIV > 0)) ~ clusComb)
pairwise.wilcox.test(as.numeric(rowSums(otuHIV > 0)), clusComb, p.adjust.method = "BH")

kruskal.test(as.numeric(apply(relaHIV, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})) ~ clusComb)
pairwise.wilcox.test(as.numeric(apply(relaHIV, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})), clusComb, p.adjust.method = "BH")





# ### Singh - Genus --------------------------------------------------------------
# resultFilename <- c(paste0(path, "Result/singh/result_cleaned_singh_chain_", 1:3, "_init_oneClus_s2_1en1_s2MH_1en3.rds"),
#                     paste0(path, "Result/singh/result_cleaned_singh_chain_", c(2, 3), "_init_3clus_s2_1en1_s2MH_1en3.rds"))
# 
# ### Active Cluster
# registerDoParallel(2)
# activeClusMat <- foreach(t = 1:5, .combine = cbind) %dopar% {
#   result <- readRDS(resultFilename[t])
#   apply(result$mod$ci_result, 1, uniqueClus)
# }
# stopImplicitCluster()
# 
# activeClusMatPlot <- activeClusMat %>%
#   as.data.frame() %>%
#   mutate(iter = 1:25000) %>%
#   pivot_longer(!iter)
# 
# activeClusMatPlot$name <- factor(activeClusMatPlot$name, 
#                                  levels = paste0("result.", 1:12), labels = paste0("Chain ", 1:12))
# 
# ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
#   geom_line() +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   labs(title = "Singh (Genus Level) - Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusters",
#        color = "MCMC Chain") +
#   guides(color = guide_legend(ncol = 5)) +
#   scale_y_continuous(breaks = seq(0, 12, 1))
# 
# ### Check the convergence of xi
# data.frame(colMeans(otuHIV), apply(otuHIV, 2, var))
# 
# registerDoParallel(2)
# xiFirst <- foreach(t = 1:5) %dopar% {
#   result <- readRDS(resultFilename[t])
#   result$mod$beta_result[, 9, ] %>% t() %>%
#     as.data.frame() %>%
#     mutate(Iteration = 1:25000) %>%
#     pivot_longer(!Iteration) %>%
#     transmute(Iteration, xi = value,
#               Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
# }
# stopImplicitCluster()
# 
# xiFirstLong <- lapply(1:5, function(x){data.frame(xiFirst[[x]], chain = paste0("Chain ", x))}) %>%
#   bind_rows()
# 
# xiFirstLong$chain <- factor(xiFirstLong$chain, levels = paste0("Chain ", 1:5))
# 
# ggplot(xiFirstLong, aes(x = Iteration, y = xi, color = Cluster)) +
#   geom_line() +
#   facet_wrap(. ~ chain) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(title = TeX(paste0("Singh (Genus Level) - Trace plot for ", "$\\xi_{k9}$ ", "in every MCMC chains")), x = "Iteration", y = TeX("\\xi"))
# 
# registerDoParallel(5)
# xiThird <- foreach(t = 1:5) %dopar% {
#   result <- readRDS(resultFilename[t])
#   result$mod$beta_result[, 3, ] %>% t() %>%
#     as.data.frame() %>%
#     mutate(Iteration = 1:25000) %>%
#     pivot_longer(!Iteration) %>%
#     transmute(Iteration, xi = value,
#               Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
# }
# stopImplicitCluster()
# 
# xiThirdLong <- lapply(1:5, function(x){data.frame(xiThird[[x]], chain = paste0("Chain ", x))}) %>%
#   bind_rows()
# 
# xiThirdLong$chain <- factor(xiThirdLong$chain, levels = paste0("Chain ", 1:12))
# 
# ggplot(xiThirdLong, aes(x = Iteration, y = xi, color = Cluster)) +
#   geom_line() +
#   facet_wrap(. ~ chain) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(title = TeX(paste0("Singh (Genus Level) - Trace plot for ", "$\\xi_{k3}$ ", "in every MCMC chains")), x = "Iteration", y = TeX("\\xi"))
# 
# 
# #### Combined Chain
# set.seed(2, kind = "L'Ecuyer-CMRG")
# registerDoParallel(1)
# globalTime <- Sys.time()
# combineClus <- foreach(t = 1:5, .combine = "rbind") %dopar% {
#   result <- readRDS(resultFilename[t])
#   result$mod$ci_result[15001:25000, ]
# }
# stopImplicitCluster()
# 
# set.seed(1)
# combSalso <- as.numeric(salso(combineClus))
# table(combSalso)
# sumTab <- data.frame(SampleID = rownames(otuEDD_genus), clus = combSalso) %>%
#   inner_join(metadata)
# 
# data.frame(table(Cluster = paste0("Cluster ", sumTab$clus), status = factor(sumTab$Status))) %>%
#   group_by(Cluster) %>%
#   mutate(Percent = Freq/sum(Freq)) %>%
#   ggplot(aes(x = Cluster, y = Percent, fill = status)) +
#   geom_bar(stat = "identity", position = "stack", width = 0.5) +
#   geom_text(aes(label = Freq), position = position_stack(vjust = 0.5)) + 
#   scale_fill_manual(values = c("springgreen3", "salmon")) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   guides(color = guide_legend(ncol = 4)) +
#   labs(y = "Frequency", fill = " ", 
#        title = "Singh et al. (Genus) - EDD status by Cluster for the combined MCMC chain.")
# 
# ### Relative Taxa Plot
# taxaName <- data.frame(k = str_remove(str_extract(colnames(otuEDD_genus), "k__[:alpha:]+"), "k__"),
#                        p = str_remove(str_extract(colnames(otuEDD_genus), "p__[:alpha:]+"), "p__"),
#                        c = str_remove(str_extract(colnames(otuEDD_genus), "c__[:alpha:]+"), "c__"),
#                        o = str_remove(str_extract(colnames(otuEDD_genus), "o__[:alpha:]+"), "o__"),
#                        f = str_remove(str_extract(colnames(otuEDD_genus), "f__[:alpha:]+"), "f__"),
#                        g = str_remove(str_extract(colnames(otuEDD_genus), "g__[:alpha:]+"), "g__"))
# 
# highTaxa <- sapply(1:2, function(x){
#   (otuEDD_genus[which(combSalso == x), ]/rowSums(otuEDD_genus[which(combSalso == x), ])) %>%
#     colMeans() %>%
#     sort(decreasing = TRUE, index.return = TRUE) %>%
#     .$ix %>%
#     .[1:10]})
# 
# otu_visual <- otuEDD_genus/rowSums(otuEDD_genus)
# # highTaxaIndex <- union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])
# highTaxaIndex <- union(highTaxa[, 1], highTaxa[, 2])
# highTaxaIndexLabel <- ifelse(is.na(taxaName$g[highTaxaIndex]), taxaName$f[highTaxaIndex], taxaName$g[highTaxaIndex])
# 
# otuRela <- matrix(NA, ncol = length(highTaxaIndex) + 1, nrow = 303)
# for(i in 1:length(highTaxaIndex)){
#   otuRela[, i] <- otu_visual[, highTaxaIndex[i]]
# }
# otuRela[, length(highTaxaIndex) + 1] <- 1 - rowSums(otuRela[, 1:length(highTaxaIndex)])
# colnames(otuRela) <- c(highTaxaIndexLabel, "Others")
# 
# otuRelaPlot <- otuRela %>%
#   as.data.frame() %>% 
#   mutate(ID = str_extract(rownames(otuEDD_genus), "[:digit:]+"), cluster = paste0("Cluster ", combSalso)) %>%
#   pivot_longer(!c(ID, cluster))
# 
# sortIndex <- data.frame(taxaName[highTaxaIndex, ], index = highTaxaIndex) %>% arrange(p) %>% .$index
# sortIndexLabel <- ifelse(is.na(taxaName$g[sortIndex]), taxaName$f[sortIndex], taxaName$g[sortIndex])
# otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(sortIndexLabel, "Others"))
# bacteria <- taxaName[sortIndex, ]$p
# table(bacteria)
# 
# # Define color palettes
# green_shades <- c("#98FB98", "#90EE90", "#00FA9A", "#00FF7F", "#2E8B57")
# blue_shades <- c("#E6E6FA", "#D8BFD8", "#DDA0DD", "#DA70D6", "#BA55D3")
# red_shades <- c("#FFFFE0", "#FFFACD", "#FAFAD2", "#FFEFD5", "#FFE4B5", "#FFD700")
# 
# # Function to assign colors
# assign_color <- function(x, index) {
#   if (x == "Bacteroidetes") {
#     return(green_shades[(index %% length(green_shades)) + 1])
#   } else if (x == "Proteobacteria") {
#     return(blue_shades[(index %% length(blue_shades)) + 1])
#   } else if (x == "Firmicutes") {
#     return(red_shades[(index %% length(red_shades)) + 1])
#   } else {
#     return("black") 
#   }
# }
# 
# colors <- sapply(seq_along(bacteria), function(i) assign_color(bacteria[i], i))
# 
# ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
#   geom_bar(position = "stack", stat = "identity") +
#   facet_wrap(. ~ cluster, scales = "free_x") +
#   theme_bw() +
#   scale_fill_manual(values = c(colors, "gray90")) + 
#   scale_y_continuous(labels = scales::percent) +
#   theme(axis.text.x = element_text(angle = 90, size = 3)) +
#   theme(legend.position = "bottom") +
#   guides(fill = guide_legend(nrow = 2)) +
#   labs(fill = "Taxa", y = "Relative Abundance",
#        title = paste0("Singh et al. (Genus) - Relative Abundance for each cluster"))













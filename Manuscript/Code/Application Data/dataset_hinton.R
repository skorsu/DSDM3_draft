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

# ### Default set of Hyperparameters
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
# ### Lower the s2_MH to 1e-5
# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# globalTime <- Sys.time()
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5, 
#                       z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60), 
#                       beta_init = matrix(0, nrow = 10, ncol = 60), 
#                       ci_init = rep(0, 155), 
#                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-5, 
#                       t_thres = 2500, launch_iter = 30, 
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_oneClus_s2MH_1en5.rds"))
# }
# stopImplicitCluster()
# difftime(Sys.time(), globalTime)

### Post Analysis: -------------------------------------------------------------
#### Read the result
resultFilename <- c(paste0(resultpath, "result_selbal_HIV_chain_", 1:6, "_init_oneClus_defaultHyper.rds"),
                    paste0(resultpath, "result_selbal_HIV_chain_", 1:6, "_init_oneClus_s2MH_1en5.rds"))

### Computational time
registerDoParallel(6)
compTime <- foreach(t = 1:12, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  as.numeric(result$time)
}
stopImplicitCluster()

mean(compTime/3600)
sd(compTime/3600)

### Active Cluster - Combine both two hyperparameters
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

activeClusMatPlot$name <- factor(activeClusMatPlot$name, levels = paste0("result.", 1:12))

ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  scale_color_discrete(labels = unname(TeX(c(paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-3}$"), ": Chain ", 1:6),
                                             paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-5}$"), ": Chain ", 1:6))))) +
  theme(legend.position = "bottom") +
  labs(title = "Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusteres",
       color = "MCMC Chain") +
  guides(color = guide_legend(ncol = 6))

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

# labelPlotchain <- c(paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-3}$ \\n "), "Chain ", 1:6),
#                     paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-5}$"), ": Chain ", 1:6))

lapply(1:12, function(x){data.frame(xiFirst[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows() %>%
  ggplot(aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k1}$")), x = "Iteration", y = TeX("\\xi"),
       color = "MCMC Chain") +
  guides(color = guide_legend(ncol = 6))

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

lapply(1:12, function(x){data.frame(xiThird[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows() %>%
  ggplot(aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k3}$")), x = "Iteration", y = TeX("\\xi"),
       color = "MCMC Chain") +
  guides(color = guide_legend(ncol = 6))

### Check the acceptance rate
registerDoParallel(6)
xiAccept <- foreach(t = 1:12) %dopar% {
  result <- readRDS(resultFilename[t])
  sapply(1:10, function(x){sum(result$mod$MH_accept[, x] == 1)/sum(result$mod$MH_accept[, x] != -1)})
}
stopImplicitCluster()

xiAcceptLong <- do.call(rbind, xiAccept) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Component ", 1:10)) %>%
  mutate(Chain = paste0("Chain ", 1:12)) %>%
  pivot_longer(!Chain)

xiAcceptLong$name <- factor(xiAcceptLong$name, levels = paste0("Component ", 10:1))
xiAcceptLong$Chain <- factor(xiAcceptLong$Chain, levels = paste0("Chain ", 1:12))

xiAcceptLong %>%
  ggplot(aes(x = Chain, y = name, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  labs(x = "MCMC Chain", y = "Component", fill = "Acceptance Rate",
       title = TeX(paste0("Acceptance Rate of the Adpative MH for ", "$\\textbf{\\xi}_{k}$"))) +
  theme(legend.position = "bottom") +
  theme_bw()

### Cluster Assignment - Individual
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

salsoDemo$Cluster <- factor(salsoDemo$Cluster, levels = paste0("Cluster ", 1:5))
salsoDemo$Chain <- factor(salsoDemo$Chain, levels = paste0("Chain ", 1:12))

salsoDemo %>%
  ggplot(aes(x = Cluster, y = Freq, fill = status_sex)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.5), vjust = -0.5) + 
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
  scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill = "Taxa", y = "Relative Abundance",
       title = paste0("Relative Abundance for each cluster"))

relaPlot <- lapply(1:12, function(x){
  otuRelaPlot <- otuRela %>%
    as.data.frame() %>% 
    mutate(ID = str_extract(rownames(otuHIV), "[:digit:]+"), cluster = paste0("Cluster ", clusSALSO[, x])) %>%
    pivot_longer(!c(ID, cluster))
  otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(highTaxaIndexLabel, "Others"))
  ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
    geom_bar(position = "stack", stat = "identity") +
    facet_grid(. ~ cluster, scales = "free_x") +
    theme_bw() +
    scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2)) +
    labs(fill = "Taxa", y = "Relative Abundance",
         title = paste0("Chain ", x, " - Relative Abundance for each cluster"))
})

ggarrange(plotlist = relaPlot, common.legend = TRUE, legend = "bottom")

### Combine all chains
registerDoParallel(6)
MCMCCombine <- foreach(t = 1:12, .combine = rbind) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$ci_result[15001:25000, ]
}
stopImplicitCluster()

set.seed(1)
clusComb <- salso(MCMCCombine)

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
  scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(fill = "Taxa", y = "Relative Abundance",
       title = paste0("Relative Abundance for each cluster"))



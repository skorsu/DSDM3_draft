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

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

# Path
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}
resultpath <- paste0(path, "Result/selbal_crohn/")
# file.exists(resultpath)

# Crohn dataset
CrohnData <- selbal::Crohn
statusCrohn <- CrohnData[, 49]
otuCrohn <- CrohnData[, -49]

# Demographic
dim(otuCrohn)
table(statusCrohn)

### Default set of Hyperparameters
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5,
                      z = as.matrix(otuCrohn), atrisk_init = matrix(1, nrow = 975, ncol = 48),
                      beta_init = matrix(0, nrow = 10, ncol = 48),
                      ci_init = rep(0, 975),
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
                      t_thres = 2500, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_crohn_chain_", t, "_init_oneClus_defaultHyper.rds"))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

### Post Analysis: -------------------------------------------------------------

resultFilename <- c(paste0(resultpath, "result_selbal_crohn_chain_", 1:6, "_init_oneClus_defaultHyper.rds"))

### Computational time
registerDoParallel(3)
compTime <- foreach(t = 1:6, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  as.numeric(result$time)
}
stopImplicitCluster()

mean(compTime/3600)
sd(compTime/3600)

### Active Cluster
registerDoParallel(3)
activeClusMat <- foreach(t = 1:6, .combine = cbind) %dopar% {
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
  guides(color = guide_legend(ncol = 6))

registerDoParallel(3)
xiFirst <- foreach(t = 1:6) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$beta_result[, 1, ] %>% t() %>%
    as.data.frame() %>%
    mutate(Iteration = 1:25000) %>%
    pivot_longer(!Iteration) %>%
    transmute(Iteration, xi = value,
              Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
}
stopImplicitCluster()

xiFirstLong <- lapply(1:6, function(x){data.frame(xiFirst[[x]], chain = paste0("Chain ", x))}) %>%
  bind_rows()

xiFirstLong$chain <- factor(xiFirstLong$chain, levels = paste0("Chain ", 1:6))

ggplot(xiFirstLong, aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ chain) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k1}$")), x = "Iteration", y = TeX("\\xi"))

result <- readRDS(resultFilename[2])
clusSALSO <- as.numeric(salso(result$mod$ci_result[-(1:5000), ]))
table(clusSALSO, statusCrohn)

colnames(otuCrohn)
taxaName <- data.frame(c = str_remove(str_extract(colnames(otuCrohn), "c__[:alpha:]+"), "c__"),
                       o = str_remove(str_extract(colnames(otuCrohn), "o__[:alpha:]+"), "o__"),
                       f = str_remove(str_extract(colnames(otuCrohn), "f__[:alpha:]+"), "f__"),
                       g = str_remove(str_extract(colnames(otuCrohn), "g__[:alpha:]+"), "g__"))

View(taxaName)
highTaxa <- sapply(1:9, function(x){
  (otuCrohn[which(clusSALSO == x), ]/rowSums(otuCrohn[which(clusSALSO == x), ])) %>%
    colMeans() %>%
    sort(decreasing = TRUE, index.return = TRUE) %>%
    .$ix %>%
    .[1:10]})

otuHIVvisual <- otuCrohn/otuCrohn(otuHIV)
highTaxaIndex <- union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])
highTaxaIndexLabel <- ifelse(taxaName$g %in% c("unclassified", "Incertae"), 
                             ifelse(is.na(taxaName$f), 
                                    ifelse(is.na(taxaName$o), taxaName$c, taxaName$o), taxaName$f), taxaName$g) 
highTaxaIndexLabel[highTaxaIndex]

# #### Read the result
# resultFilename <- c(paste0(resultpath, "result_selbal_HIV_chain_", 1:6, "_init_oneClus_defaultHyper.rds"),
#                     paste0(resultpath, "result_selbal_HIV_chain_", 1:6, "_init_oneClus_s2MH_1en5.rds"))
# 
# ### Active Cluster - Combine both two hyperparameters
# registerDoParallel(6)
# activeClusMat <- foreach(t = 1:12, .combine = cbind) %dopar% {
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
# activeClusMatPlot$name <- factor(activeClusMatPlot$name, levels = paste0("result.", 1:12))
# 
# ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
#   geom_line() +
#   theme_bw() +
#   scale_color_discrete(labels = unname(TeX(c(paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-3}$"), ": Chain ", 1:6),
#                                              paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-5}$"), ": Chain ", 1:6))))) +
#   theme(legend.position = "bottom") +
#   labs(title = "Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusteres",
#        color = "MCMC Chain")
# 
# ### Check the convergence of xi
# result <- readRDS(paste0(resultpath, "result_selbal_HIV_chain_", 5, "_init_oneClus_s2MH_1en5.rds"))
# result$mod$beta_result[, 1, ] %>% t() %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:25000) %>%
#   pivot_longer(!Iteration) %>%
#   transmute(Iteration, xi = value,
#             Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+"))) %>%
#   ggplot(aes(x = Iteration, y = xi, color = Cluster)) +
#   geom_line() +
#   theme_bw()
# 
# ### Check the acceptance rate (xi overall)
# sapply(1:10, function(x){sum(result$mod$MH_accept[, x] == 1)/sum(result$mod$MH_accept[, x] != -1)})
# 
# sapply(1:10, function(x){
#   
#   c(sum(result$mod$MH_accept[, x] != -1), ### Total Active
#     sum(result$mod$MH_accept[, x] == 1), ### Total Accept
#     sum(result$mod$MH_accept[1:2500, x] != -1), ### Active Before Adaptive
#     sum(result$mod$MH_accept[1:2500, x] == 1)) ### Accept Before Adaptive
#   
# }) %>% t()
# 
# ### Cluster Assignment - Individual
# registerDoParallel(6)
# clusSALSO <- foreach(t = 1:12, .combine = cbind) %dopar% {
#   result <- readRDS(resultFilename[t])
#   as.numeric(salso(result$mod$ci_result[-(1:5000), ]))
# }
# stopImplicitCluster()
# 
# expand.grid(1:12, 1:12) %>%
#   apply(1, function(x){mclustcomp(clusSALSO[, x[1]], clusSALSO[, x[2]])[1, 2]}) %>%
#   matrix(ncol = 12)
# 
# table(clusSALSO[, 1], statusHIV)
# table(clusSALSO[, 1], sexHIV)
# 
# taxaName <- data.frame(c = str_remove(str_extract(colnames(otuHIV), "c_[:alpha:]+"), "c_"),
#                        o = str_remove(str_extract(colnames(otuHIV), "o_[:alpha:]+"), "o_"),
#                        f = str_remove(str_extract(colnames(otuHIV), "f_[:alpha:]+"), "f_"),
#                        g = str_remove(str_extract(colnames(otuHIV), "g_[:alpha:]+"), "g_"))
# 
# highTaxa <- sapply(1:3, function(x){
# (otuHIV[which(clusSALSO[, 1] == x), ]/rowSums(otuHIV[which(clusSALSO[, 1] == x), ])) %>%
#   colMeans() %>%
#   sort(decreasing = TRUE, index.return = TRUE) %>%
#   .$ix %>%
#   .[1:10]})
# 
# otuHIVvisual <- otuHIV/rowSums(otuHIV)
# highTaxaIndex <- union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])
# highTaxaIndexLabel <- ifelse(taxaName$g %in% c("unclassified", "Incertae"), 
#        ifelse(is.na(taxaName$f), 
#               ifelse(is.na(taxaName$o), taxaName$c, taxaName$o), taxaName$f), taxaName$g) 
# highTaxaIndexLabel[12] <- "Lachnospiraceae - Incertae"
# highTaxaIndexLabel <- highTaxaIndexLabel %>% .[union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])]
# highTaxaIndexLabel[is.na(highTaxaIndexLabel)] <- "Unclassified"
# 
# otuRela <- matrix(NA, ncol = length(highTaxaIndex) + 1, nrow = 155)
# for(i in 1:length(highTaxaIndex)){
#   otuRela[, i] <- otuHIVvisual[, highTaxaIndex[i]]
# }
# otuRela[, length(highTaxaIndex) + 1] <- 1 - rowSums(otuRela[, 1:length(highTaxaIndex)])
# colnames(otuRela) <- c(highTaxaIndexLabel, "Others")
# 
# otuRelaPlot <- otuRela %>%
#   as.data.frame() %>% 
#   mutate(ID = rownames(otuHIV), cluster = paste0("Cluster ", clusSALSO[, 1])) %>%
#   pivot_longer(!c(ID, cluster))
# 
# otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(highTaxaIndexLabel, "Others"))
# otuRelaPlot$ID <- factor(otuRelaPlot$ID, levels = c(rownames(otuHIVvisual)[which(clusSALSO[, 1] == 1)][sort(otuHIVvisual[which(clusSALSO[, 1] == 1), 3], decreasing = TRUE, index.return = TRUE)$ix],
#                                                     rownames(otuHIVvisual)[which(clusSALSO[, 1] == 2)][sort(otuHIVvisual[which(clusSALSO[, 1] == 2), 6], decreasing = TRUE, index.return = TRUE)$ix],
#                                                     rownames(otuHIVvisual)[which(clusSALSO[, 1] == 3)][sort(otuHIVvisual[which(clusSALSO[, 1] == 3), 1], decreasing = TRUE, index.return = TRUE)$ix]))
# 
# rownames(otuHIVvisual)[which(clusSALSO[, 1] == 1)][sort(otuHIVvisual[which(clusSALSO[, 1] == 1), 3], decreasing = TRUE, index.return = TRUE)$ix]
# rownames(otuHIVvisual)[which(clusSALSO[, 1] == 2)][sort(otuHIVvisual[which(clusSALSO[, 1] == 2), 6], decreasing = TRUE, index.return = TRUE)$ix]
# rownames(otuHIVvisual)[which(clusSALSO[, 1] == 3)][sort(otuHIVvisual[which(clusSALSO[, 1] == 3), 1], decreasing = TRUE, index.return = TRUE)$ix]
# 
# ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
#   geom_bar(position = "stack", stat = "identity") +
#   facet_wrap(. ~ cluster, scales = "free_x") +
#   theme_bw() +
#   scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
#   scale_y_continuous(labels = scales::percent) +
#   theme(axis.text.x = element_text(angle = 90))
# 
# table(clusSALSO[, 1], statusHIV)
# 
# ### Combine all chains
# registerDoParallel(6)
# MCMCCombine <- foreach(t = 1:12, .combine = rbind) %dopar% {
#   result <- readRDS(resultFilename[t])
#   result$mod$ci_result[seq(15000, 25000, 500), ]
# }
# stopImplicitCluster()
# 
# clusComb <- salso(MCMCCombine)
# table(clusComb)
# table(clusSALSO[, 7], statusHIV)
# table(clusComb, statusHIV)

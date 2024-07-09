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
resultpath <- paste0(path, "Result/selbal_HIV/")
# file.exists(resultpath)

# Crohn dataset
sCD14Data <- selbal::sCD14
status_sCD14 <- sCD14Data[, 61]
otu_sCD14 <- sCD14Data[, -61]

# otu_sCD14[, -which(colMeans(otu_sCD14 > 0) < 0.1)]

# Demographic
# boxplot(status_sCD14)
# dim(otu_sCD14)

# Run the model
# start_time <- Sys.time()
# mod <- mod_adaptive(iter = 500, Kmax = 10, nbeta_split = 5, 
#                     z = as.matrix(otu_sCD14), atrisk_init = matrix(1, nrow = 151, ncol = 60), 
#                     beta_init = matrix(0, nrow = 10, ncol = 60), 
#                     ci_init = rep(0, 151), 
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3, 
#                     t_thres = 50, launch_iter = 30, 
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# comp_time <- difftime(Sys.time(), start_time, units = "secs")
# apply(mod$ci_result, 1, uniqueClus) %>% plot(type = "l")
# 
# clusSALSO <- salso(mod$ci_result) %>% as.numeric()
# data.frame(cluster = clusSALSO, sCD14 = status_sCD14) %>%
#   group_by(cluster) %>%
#   summarize(mean(sCD14), sd(sCD14), median(sCD14))
# 
# data.frame(cluster = clusSALSO, sCD14 = status_sCD14) %>%
#   ggplot(aes(x = cluster, y = sCD14, group = cluster)) +
#   geom_boxplot()

### Default set of Hyperparameters
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5, 
                      z = as.matrix(otu_sCD14), atrisk_init = matrix(1, nrow = 151, ncol = 60), 
                      beta_init = matrix(0, nrow = 10, ncol = 60), 
                      ci_init = rep(0, 151), 
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3, 
                      t_thres = 2500, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_sCD14_chain_", t, "_init_oneClus_defaultHyper.rds"))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

### Lower the s2_MH to 1e-5
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5, 
                      z = as.matrix(otu_sCD14), atrisk_init = matrix(1, nrow = 151, ncol = 60), 
                      beta_init = matrix(0, nrow = 10, ncol = 60), 
                      ci_init = rep(0, 151), 
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-5, 
                      t_thres = 2500, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_sCD14_chain_", t, "_init_oneClus_s2MH_1en5.rds"))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

### Post Analysis: -------------------------------------------------------------
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

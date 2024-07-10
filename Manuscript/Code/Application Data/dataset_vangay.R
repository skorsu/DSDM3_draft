# Required Libraries
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
datapath <- paste0(path, "Data/Application Data/vangay_data/")
resultpath <- paste0(path, "Result/vangay_data/")
# file.exists(resultpath)

dat <- read.delim(paste0(datapath, "ravel/refseq/taxatable.txt"))
dat <- dat[-which(rowMeans(dat[, -1] > 0) < 0.1), ]
taxaName <- dat[, 1]
dat <- t(dat)
dat <- dat[-1, ]
mode(dat) <- "numeric"
colnames(dat) <- taxaName

# Demographic
demoBH <- read.delim(paste0(datapath, "ravel/task-black-hispanic.txt"))
demoWB <- read.delim(paste0(datapath, "ravel/task-white-black.txt"))
colnames(demoBH)[1] <- "ID"
colnames(demoWB)[1] <- "ID"
demoWBH <- full_join(demoBH, demoWB)

demoN <- read.delim(paste0(datapath, "ravel/task-nugent-category.txt"))
colnames(demoN)[1] <- "ID"
colnames(demoN)[2] <- "Nugent"
demoWBHN <- full_join(demoWBH, demoN)
colnames(demoWBHN) <- c("ID", "Race", "Nugent")

dat <- dat[which(rownames(dat) %in% intersect(rownames(dat), demoWBHN$ID)), ]
dim(dat)

# mod <- mod_adaptive(iter = 2500, Kmax = 10, nbeta_split = 5,
#                     z = as.matrix(dat), atrisk_init = matrix(1, nrow = 394, ncol = 56),
#                     beta_init = matrix(0, nrow = 10, ncol = 56),
#                     ci_init = rep(0, 394),
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                     t_thres = 1000, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# 
# sapply(1:10, function(x){sum(mod$MH_accept[, x] == 1)/sum(mod$MH_accept[, x] != -1)})
# 
# plot(apply(mod$ci_result, 1, uniqueClus), type = "l")
# clusZZ <- as.numeric(salso(mod$ci_result[-(1:1000), ]))
# 
# demoClus <- data.frame(ID = rownames(dat), clus = clusZZ) %>%
#   inner_join(demoWBHN)
# 
# table(demoClus$Nugent, demoClus$clus)
# table(demoClus$Var, demoClus$clus)
# 
# rownames(dat)

### Default set of Hyperparameters
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5,
                      z = as.matrix(dat), atrisk_init = matrix(1, nrow = 375, ncol = 56),
                      beta_init = matrix(0, nrow = 10, ncol = 56),
                      ci_init = rep(0, 375),
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
                      t_thres = 2500, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_ravel_chain_", t, "_init_oneClus_defaultHyper.rds"))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

### Different starting point: 1e-3 - PART 1
set.seed(1)
ciInit <- matrix(NA, nrow = 375, ncol = 6)
ciInit[, 1] <- sample(0:4, 375, replace = TRUE)
ciInit[, 2] <- sample(0:4, 375, replace = TRUE)
ciInit[, 3] <- sample(0:4, 375, replace = TRUE)
ciInit[, 4] <- sample(0:19, 375, replace = TRUE)
ciInit[, 5] <- sample(0:19, 375, replace = TRUE)
ciInit[, 6] <- sample(0:19, 375, replace = TRUE)

KmaxVec <- c(20, 20, 20, 50, 50, 50)

xiInit <- lapply(1:6, function(y){sapply(0:max(ciInit[, y]), function(x){
  colSums(dat[which(ciInit[, y] == x), ])/sum(dat[which(ciInit[, y] == x), ])
}) %>% t()
})

xiInit[[1]] <- rbind(xiInit[[1]], matrix(0, nrow = 15, ncol = 56))
xiInit[[2]] <- rbind(xiInit[[2]], matrix(0, nrow = 15, ncol = 56))
xiInit[[3]] <- rbind(xiInit[[3]], matrix(0, nrow = 15, ncol = 56))
xiInit[[4]] <- rbind(xiInit[[4]], matrix(0, nrow = 30, ncol = 56))
xiInit[[5]] <- rbind(xiInit[[5]], matrix(0, nrow = 30, ncol = 56))
xiInit[[6]] <- rbind(xiInit[[6]], matrix(0, nrow = 30, ncol = 56))

resultName <- c(paste0("result_ravel_chain_", 1:3, "_init_5clus_Kmax_20_defaultHyper.rds"),
                paste0("result_ravel_chain_", 1:3, "_init_20clus_Kmax_50_defaultHyper.rds"))

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = KmaxVec[t], nbeta_split = 5,
                      z = as.matrix(dat), atrisk_init = matrix(1, nrow = 375, ncol = 56),
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

### Different starting point: 1e-3 - PART 2
set.seed(1)
ciInit <- matrix(NA, nrow = 375, ncol = 3)
ciInit[, 1] <- sample(0:2, 375, replace = TRUE)
ciInit[, 2] <- sample(0:2, 375, replace = TRUE)
ciInit[, 3] <- sample(0:2, 375, replace = TRUE)

xiInit <- lapply(1:3, function(y){sapply(0:max(ciInit[, y]), function(x){
  colSums(otuHIV[which(ciInit[, y] == x), ])/sum(otuHIV[which(ciInit[, y] == x), ])
}) %>% t()
})

xiInit[[1]] <- rbind(xiInit[[1]], matrix(0, nrow = 12, ncol = 56))
xiInit[[2]] <- rbind(xiInit[[2]], matrix(0, nrow = 12, ncol = 56))
xiInit[[3]] <- rbind(xiInit[[3]], matrix(0, nrow = 12, ncol = 56))

resultName <- paste0("result_ravel_chain_", 1:3, "_init_3clus_Kmax_15_defaultHyper.rds")

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(3)
globalTime <- Sys.time()
foreach(t = 1:3) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 15, nbeta_split = 5,
                      z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 375, ncol = 56),
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

# ### Post Analysis: -------------------------------------------------------------
# resultFilename <- c(paste0(resultpath, "result_selbal_crohn_chain_", 1:6, "_init_oneClus_defaultHyper.rds"),
#                     paste0(resultpath, "result_selbal_crohn_chain_", 1:6, "_init_oneClus_s2MH_1en5.rds"))
# 
# ### Active Cluster - Combine both two hyperparameters
# registerDoParallel(2)
# activeClusMat <- foreach(t = 1:6, .combine = cbind) %dopar% {
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
# activeClusMatPlot$name <- factor(activeClusMatPlot$name, levels = paste0("result.", 1:6))
# 
# ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
#   geom_line() +
#   theme_bw() +
#   scale_color_discrete(labels = unname(TeX(c(paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-3}$"), ": Chain ", 1:6),
#                                              paste0(rep("$\\sigma^{2}_{MH} = 1 \\times 10^{-5}$"), ": Chain ", 1:6))))) +
#   theme(legend.position = "bottom") +
#   labs(title = "Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusteres",
#        color = "MCMC Chain") +
#   guides(color = guide_legend(ncol = 6))
# 
# result <- readRDS(resultFilename[1])
# set.seed(1)
# clusSALSO <- salso(result$mod$ci_result[seq(5000, 25000, 100), ])
# clusSALSO2 <- salso(result$mod$ci_result[seq(5000, 25000, 100), ])
# 
# table(clusSALSO, statusCrohn)
# 
# View(otuCrohn)
# sapply(1:10, function(x){sum(result$mod$MH_accept[, x] == 1)/sum(result$mod$MH_accept[, x] != -1)})
# 
# t(result$mod$beta_result[, 3, ]) %>%
#   as.data.frame() %>%
#   mutate(Iter = 1:25000) %>%
#   pivot_longer(!Iter) %>%
#   ggplot(aes(x = Iter, y = value, color = name)) +
#   geom_line()
# 
# ### Dummy: ----
# mod <- mod_adaptive(iter = 2500, Kmax = 10, nbeta_split = 5,
#                     z = as.matrix(otuCrohn), atrisk_init = matrix(1, nrow = 975, ncol = 48),
#                     beta_init = matrix(0, nrow = 10, ncol = 48),
#                     ci_init = rep(0, 975),
#                     theta = 1, mu = 0, s2 = 10, s2_MH = 1e-3,
#                     t_thres = 1000, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# 
# 
# sapply(1:10, function(x){sum(mod$MH_accept[, x] == 1)/sum(mod$MH_accept[, x] != -1)})
# t(mod$beta_result[, 3, ]) %>%
#   as.data.frame() %>%
#   mutate(Iter = 1:2500) %>%
#   pivot_longer(!Iter) %>%
#   ggplot(aes(x = Iter, y = value, color = name)) +
#   geom_line()






### Load libraries: ------------------------------------------------------------
library(foreach)
library(doParallel)
library(stringr)
library(tidyverse)
library(salso)
library(ggplot2)
library(coda)
library(latex2exp)

### User-defined Functions: ----------------------------------------------------
meanSD <- function(x, dplace = 3){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

uniqueClus <- function(x){
  length(unique(x))
}

### Import the data: -----------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

month_analysis <- "6m"
caseName <- c("One Cluster", "Singletons", "30 Clusters", "60 Clusters")

r1 <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_oneClus.RData")) 
r2 <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_singleton.RData"))
r3 <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_init30.RData"))
r4 <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_init60.RData"))

listData <- list(r1 = r1, r2 = r2, r3 = r3, r4 = r4)

### Active Clusters: -----------------------------------------------------------
lapply(1:4, function(y){sapply(1:5, function(x){apply(listData[[y]][[x]]$result$ci_result, 1, uniqueClus)}) %>%
    as.data.frame() %>%
    `colnames<-`(paste0("Chain ", 1:5)) %>%
    mutate(Iteration = 1:15000, case = caseName[y])}) %>%
  bind_rows(.id = NULL) %>%
  pivot_longer(!(c(Iteration, case)), names_to = "Chain", values_to = "Cluster") %>%
  mutate(caseName = paste0(case, ": ", Chain)) %>%
  ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
  geom_line() +
  theme_bw() +
  facet_wrap(. ~ factor(case, levels = c("One Cluster", "Singletons", "30 Clusters", "60 Clusters")),
             scales = "free_y") +
  theme(legend.position = "bottom") +
  labs(y = "Active Cluster", 
       title = paste0(str_extract(month_analysis, "^[:digit:]+"), "-Month: Active Clusters via MCMC chains")) +
  scale_y_continuous(limits = c(0, 45), breaks = c(seq(2, 10, 2), seq(15, 45, 5)))

###: ---------------------------------------------------------------------------

### Active Clusters
sapply(1:5, function(x){apply(testResult1[[x]]$result$ci_result, 1, uniqueClus)}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Iteration = 1:15000) %>%
  pivot_longer(!Iteration) %>%
  ggplot(aes(x = Iteration, y = value)) +
  geom_line() +
  theme_bw() +
  facet_wrap(. ~ name) +
  ylim(1, 10) + 
  labs(title = TeX("Active Clusters: Adaptive Metropolis-Hasting for $\\beta$ (6-Month: One Cluster)"))

sapply(1:5, function(x){apply(testResult2[[x]]$result$ci_result, 1, uniqueClus)}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Iteration = 1:15000) %>%
  pivot_longer(!Iteration) %>%
  ggplot(aes(x = Iteration, y = value)) +
  geom_line() +
  theme_bw() +
  facet_wrap(. ~ name) +
  ylim(1, 10) + 
  labs(title = TeX("Active Clusters: Adaptive Metropolis-Hasting for $\\beta$ (8-Month: One Cluster)"))

sapply(1:5, function(x){apply(testResult3[[x]]$result$ci_result, 1, uniqueClus)}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Iteration = 1:15000) %>%
  pivot_longer(!Iteration) %>%
  ggplot(aes(x = Iteration, y = value)) +
  geom_line() +
  theme_bw() +
  facet_wrap(. ~ name) +
  ylim(1, 10) + 
  labs(title = TeX("Active Clusters: Adaptive Metropolis-Hasting for $\\beta$ (12-Month: One Cluster)"))

# sapply(1:5, function(x){apply(testResult2[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   `colnames<-`(paste0("Chain ", 1:5)) %>%
#   mutate(Iteration = 1:15000) %>%
#   pivot_longer(!Iteration) %>%
#   ggplot(aes(x = Iteration, y = value)) +
#   geom_line() +
#   theme_bw() +
#   facet_wrap(. ~ name) +
#   ylim(1, 90) + 
#   geom_hline(yintercept = 2, linetype = "dotted", color = "blue", size = 0.25) +
#   labs(title = TeX("Active Clusters: Adaptive Metropolis-Hasting for $\\beta$ (Singleton)"))
# 
# sapply(1:5, function(x){apply(testResult3[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   `colnames<-`(paste0("Chain ", 1:5)) %>%
#   mutate(Iteration = 1:15000) %>%
#   pivot_longer(!Iteration) %>%
#   ggplot(aes(x = Iteration, y = value)) +
#   geom_line() +
#   theme_bw() +
#   facet_wrap(. ~ name) +
#   ylim(1, 30) + 
#   geom_hline(yintercept = 2, linetype = "dotted", color = "blue", size = 0.25) +
#   labs(title = TeX("Active Clusters: Adaptive Metropolis-Hasting for $\\beta$ (30 Clusters)"))
# 
# sapply(1:5, function(x){apply(testResult4[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   `colnames<-`(paste0("Chain ", 1:5)) %>%
#   mutate(Iteration = 1:15000) %>%
#   pivot_longer(!Iteration) %>%
#   ggplot(aes(x = Iteration, y = value)) +
#   geom_line() +
#   theme_bw() +
#   facet_wrap(. ~ name) +
#   ylim(1, 60) + 
#   geom_hline(yintercept = 2, linetype = "dotted", color = "blue", size = 0.25) +
#   labs(title = TeX("Active Clusters: Adaptive Metropolis-Hasting for $\\beta$ (60 Clusters)"))

### Acceptance Rate for beta
#### One Cluster
AcceptRate <- lapply(1:5, function(y){sapply(1:10, function(x){mean(testResult1[[y]]$result$MH_accept[which(testResult1[[y]]$result$MH_accept[, x] != -1), x])})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:10)) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  round(3)

nUpdate <- lapply(1:5, function(y){sapply(1:10, function(x){sum(testResult1[[y]]$result$MH_accept[, x] != -1)})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:10)) %>%
  `colnames<-`(paste0("Chain ", 1:5))

AcceptLonger <- data.frame(AcceptRate) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:10)) %>%
  pivot_longer(!(Cluster))

nLonger <- sapply(1:5, function(k){paste0(AcceptRate[, k], " (", nUpdate[, k], ")")}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:10)) %>%
  pivot_longer(!(Cluster))

ggplot() +
  geom_tile(data = AcceptLonger, aes(x = name, y = Cluster, fill = value)) +
  geom_text(data = nLonger, aes(x = name, y = Cluster, label = value), color = "white") +
  labs(x = "Chain", y = "Cluster", title = TeX("Acceptance Rate when updating $\\beta$ via Adaptive MH (One Cluster)")) +
  theme_minimal()

#### Singleton
AcceptRate <- lapply(1:5, function(y){sapply(1:90, function(x){mean(testResult2[[y]]$result$MH_accept[which(testResult2[[y]]$result$MH_accept[, x] != -1), x])})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:90)) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  round(3)

nUpdate <- lapply(1:5, function(y){sapply(1:90, function(x){sum(testResult2[[y]]$result$MH_accept[, x] != -1)})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:90)) %>%
  `colnames<-`(paste0("Chain ", 1:5))

AcceptLonger <- data.frame(AcceptRate) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:90)) %>%
  pivot_longer(!(Cluster))

nLonger <- sapply(1:5, function(k){paste0(AcceptRate[, k], " (", nUpdate[, k], ")")}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:90)) %>%
  pivot_longer(!(Cluster))

ggplot() +
  geom_tile(data = AcceptLonger, aes(x = name, y = Cluster, fill = value)) +
  geom_text(data = nLonger, aes(x = name, y = Cluster, label = value), color = "white") +
  labs(x = "Chain", y = "Cluster", title = TeX("Acceptance Rate when updating $\\beta$ via Adaptive MH (Singleton)")) +
  theme_minimal()

#### 30 Clusters
AcceptRate <- lapply(1:5, function(y){sapply(1:30, function(x){mean(testResult3[[y]]$result$MH_accept[which(testResult3[[y]]$result$MH_accept[, x] != -1), x])})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:30)) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  round(3)

nUpdate <- lapply(1:5, function(y){sapply(1:30, function(x){sum(testResult3[[y]]$result$MH_accept[, x] != -1)})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:30)) %>%
  `colnames<-`(paste0("Chain ", 1:5))

AcceptLonger <- data.frame(AcceptRate) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:30)) %>%
  pivot_longer(!(Cluster))

nLonger <- sapply(1:5, function(k){paste0(AcceptRate[, k], " (", nUpdate[, k], ")")}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:30)) %>%
  pivot_longer(!(Cluster))

ggplot() +
  geom_tile(data = AcceptLonger, aes(x = name, y = Cluster, fill = value)) +
  geom_text(data = nLonger, aes(x = name, y = Cluster, label = value), color = "white") +
  labs(x = "Chain", y = "Cluster", title = TeX("Acceptance Rate when updating $\\beta$ via Adaptive MH (30 Clusters)")) +
  theme_minimal()

#### 60 Clusters
AcceptRate <- lapply(1:5, function(y){sapply(1:60, function(x){mean(testResult4[[y]]$result$MH_accept[which(testResult4[[y]]$result$MH_accept[, x] != -1), x])})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:60)) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  round(3)

nUpdate <- lapply(1:5, function(y){sapply(1:60, function(x){sum(testResult4[[y]]$result$MH_accept[, x] != -1)})}) %>%
  as.data.frame() %>%
  bind_rows(.id = NULL) %>%
  `rownames<-`(paste0("Cluster ", 1:60)) %>%
  `colnames<-`(paste0("Chain ", 1:5))

AcceptLonger <- data.frame(AcceptRate) %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:60)) %>%
  pivot_longer(!(Cluster))

nLonger <- sapply(1:5, function(k){paste0(AcceptRate[, k], " (", nUpdate[, k], ")")}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:5)) %>%
  mutate(Cluster = paste0("Cluster ", 1:60)) %>%
  pivot_longer(!(Cluster))

ggplot() +
  geom_tile(data = AcceptLonger, aes(x = name, y = Cluster, fill = value)) +
  geom_text(data = nLonger, aes(x = name, y = Cluster, label = value), color = "white") +
  labs(x = "Chain", y = "Cluster", title = TeX("Acceptance Rate when updating $\\beta$ via Adaptive MH (60 Clusters)")) +
  theme_minimal()

### Trace plot for beta
lapply(1:5, function(y){sapply(1:10, function(x){sum(testResult1[[y]]$result$MH_accept[, x] != -1)})})

aC <- which(sapply(1:10, function(x){sum(testResult1[[4]]$result$MH_accept[, x] != -1)}) >= 10000)
lapply(aC, function(x){t(testResult1[[4]]$result$beta_result[x, , ]) %>%
    as.data.frame() %>%
    `colnames<-`(paste0("b", 1:38)) %>%
    mutate(Iteration = 1:15000, Cluster = paste0("Cluster ", x))}) %>%
  bind_rows(.id = NULL) %>%
  pivot_longer(!(c(Iteration, Cluster))) %>%
  ggplot(aes(x = Iteration, y = value, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ name) +
  labs(title = TeX("Trace plot for $\\beta$ via Adpative MH (One Cluster - Chain 4)"))


aC <- which(sapply(1:90, function(x){sum(testResult2[[1]]$result$MH_accept[, x] != -1)}) >= 10000)
lapply(aC, function(x){t(testResult2[[1]]$result$beta_result[x, , ]) %>%
    as.data.frame() %>%
    `colnames<-`(paste0("b", 1:38)) %>%
    mutate(Iteration = 1:15000, Cluster = paste0("Cluster ", x))}) %>%
  bind_rows(.id = NULL) %>%
  pivot_longer(!(c(Iteration, Cluster))) %>%
  ggplot(aes(x = Iteration, y = value, color = Cluster)) +
  geom_line() +
  facet_wrap(. ~ name) +
  labs(title = TeX("Trace plot for $\\beta$ via Adpative MH (Actual Singleton - Chain 1)"))

### ACF plot
# sapply(1:38, function(x){acf(as.numeric(testResult1[[1]]$result$beta_result[8, x, ]), lag.max = 100, plot = FALSE)[[1]]}) %>%
#   `colnames<-`(paste0("b", 1:38)) %>%
#   as.data.frame() %>%
#   mutate(Lags = 0:100) %>%
#   pivot_longer(!Lags) %>%
#   ggplot(aes(x = Lags, y = value, color = name)) +
#   geom_line() +
#   ylim(-1, 1)

### Geweke
sapply(1:90, function(y){sapply(1:38, function(x){geweke.diag(mcmc(testResult1[[5]]$result$beta_result[y, x, ], start = 5000))$z}) %>%
    abs() %>%
    pnorm(lower.tail = FALSE) * 2}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Cluster ", 1:90)) %>%
  mutate(beta = paste0("b", 1:38)) %>%
  pivot_longer(!beta) %>%
  mutate(significant = (value < 0.05)) %>%
  ggplot(aes(x = name, y = beta, fill = significant)) +
  geom_tile() +
  labs(title = "Geweke Statistics (Chain 5) -- MH update (Significant = not converged)")

### salso
lapply(1:5, function(x){table(salso(testResult1[[x]]$result$ci_result[-(1:5000), ]))})
lapply(1:5, function(x){data.frame(testResult1[[x]]$result$ci_result[seq(5000, 15000, 100), ])}) %>%
  bind_rows(.id = NULL) %>%
  as.matrix() %>%
  salso() %>%
  table()

lapply(1:5, function(x){table(salso(testResult2[[x]]$result$ci_result[-(1:5000), ]))})
lapply(1:5, function(x){data.frame(testResult2[[x]]$result$ci_result[seq(5000, 15000, 100), ])}) %>%
  bind_rows(.id = NULL) %>%
  as.matrix() %>%
  salso() %>%
  table()

lapply(1:5, function(x){table(salso(testResult3[[x]]$result$ci_result[-(1:5000), ]))})
lapply(1:5, function(x){data.frame(testResult3[[x]]$result$ci_result[seq(5000, 15000, 100), ])}) %>%
  bind_rows(.id = NULL) %>%
  as.matrix() %>%
  salso() %>%
  table()

lapply(1:5, function(x){table(salso(testResult4[[x]]$result$ci_result[-(1:5000), ]))})
lapply(1:5, function(x){data.frame(testResult4[[x]]$result$ci_result[seq(5000, 15000, 100), ])}) %>%
  bind_rows(.id = NULL) %>%
  as.matrix() %>%
  salso() %>%
  table()

###
lapply(1:5, function(x){data.frame(testResult4[[x]]$result$ci_result[seq(5000, 15000, 100), ])}) %>%
  bind_rows(.id = NULL) %>%
  as.matrix() %>%
  salso() %>%
  as.numeric() %>%
  table(dat06[, 1])

### Run the models -------------------------------------------------------------
#### Set of the hyperparameters with a longer MCMC chains
# setHyper <- list(c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 10, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 1, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 10, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 0.1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 10, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 0.1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 0.1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 4, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 9, r1g = 1, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 4, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 9, r0c = 1, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 4, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 9, r1c = 1),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 4),
#                  c(nbeta_split = 5, theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 9))

### Run with different starting point with less thinning
#### One Clusters
set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
testResult1 <- foreach(t = 1:5) %dopar% {
  start_time <- Sys.time()
  clus_result <- mod_adaptive(iter = 15000, Kmax = 10, nbeta_split = 1,
                              z = as.matrix(dat06[, -c(1:5)]), atrisk_init = matrix(1, ncol = 38, nrow = 90),
                              beta_init = matrix(0, ncol = 38, nrow = 10),
                              ci_init = rep(0, 90), theta = 1, mu = 0,
                              s2 = 1, s2_MH = 1e-3, t_thres = 5000, launch_iter = 30, 
                              r0g = 4, r1g = 1, r0c = 1, r1c = 1,
                              thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(testResult1, paste0(path, "Manuscript/Result/microbiome_result_6m_oneClus.RData"))
rm(testResult1)

#### Singleton
set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
testResult2 <- foreach(t = 1:5) %dopar% {
  start_time <- Sys.time()
  clus_result <- mod_adaptive(iter = 15000, Kmax = 90, nbeta_split = 1,
                              z = as.matrix(dat06[, -c(1:5)]), atrisk_init = matrix(1, ncol = 38, nrow = 90),
                              beta_init = as.matrix(dat06[, -c(1:5)])/rowSums(as.matrix(dat06[, -c(1:5)])),
                              ci_init = 0:89, theta = 1, mu = 0,
                              s2 = 1, s2_MH = 1e-3, t_thres = 5000, launch_iter = 30, 
                              r0g = 4, r1g = 1, r0c = 1, r1c = 1,
                              thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(testResult2, paste0(path, "Manuscript/Result/microbiome_result_6m_singleton.RData"))
rm(testResult2)

### 1/3 Clusters = Start with 30 clusters
set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
testResult3 <- foreach(t = 1:5) %dopar% {
  
  beta_mat_init <- matrix(NA, nrow = 30, ncol = 38)
  ci_init_30 <- sample(rep(1:30, 3), size = 90, replace = FALSE)
  for(i in 1:30){
    beta_mat_init[i, ] <- as.numeric(colSums(dat06[which(ci_init_30 == i), -(1:5)])/sum(dat06[which(ci_init_30 == i), -(1:5)]))
  }
  
  start_time <- Sys.time()
  clus_result <- mod_adaptive(iter = 15000, Kmax = 30, nbeta_split = 1,
                              z = as.matrix(dat06[, -c(1:5)]), atrisk_init = matrix(1, ncol = 38, nrow = 90),
                              beta_init = beta_mat_init,
                              ci_init = ci_init_30 - 1, theta = 1, mu = 0,
                              s2 = 1, s2_MH = 1e-3, t_thres = 5000, launch_iter = 30, 
                              r0g = 4, r1g = 1, r0c = 1, r1c = 1,
                              thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(testResult3, paste0(path, "Manuscript/Result/microbiome_result_6m_init30.RData"))
rm(testResult3)

### 1/1.5 Clusters = Start with 60 clusters
set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
testResult4 <- foreach(t = 1:5) %dopar% {
  
  beta_mat_init <- matrix(NA, nrow = 60, ncol = 38)
  ci_init_60 <- sample(c(1:60, 1:30), size = 90, replace = FALSE)
  for(i in 1:60){
    beta_mat_init[i, ] <- as.numeric(colSums(dat06[which(ci_init_60 == i), -(1:5)])/sum(dat06[which(ci_init_60 == i), -(1:5)]))
  }
  
  start_time <- Sys.time()
  clus_result <- mod_adaptive(iter = 15000, Kmax = 60, nbeta_split = 1,
                              z = as.matrix(dat06[, -c(1:5)]), atrisk_init = matrix(1, ncol = 38, nrow = 90),
                              beta_init = beta_mat_init,
                              ci_init = ci_init_60 - 1, theta = 1, mu = 0,
                              s2 = 1, s2_MH = 1e-3, t_thres = 5000, launch_iter = 30, 
                              r0g = 4, r1g = 1, r0c = 1, r1c = 1,
                              thin = 1)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(testResult4, paste0(path, "Manuscript/Result/microbiome_result_6m_init60.RData"))
rm(testResult4)

### Try running case 3+10 with shorter iterations ------------------------------
##### matrix(0, ncol = 38, nrow = 10)
##### as.matrix(dat06[, -c(1:5)])/rowSums(as.matrix(dat06[, -c(1:5)]))
# set.seed(1415, kind = "L'Ecuyer-CMRG")
# start_ova <- Sys.time()
# registerDoParallel(7)
# testResult <- foreach(t = 1:7) %dopar% {
#     start_time <- Sys.time()
#     clus_result <- mod(iter = 100000, Kmax = 10, nbeta_split = 1,
#                        z = as.matrix(dat08[, -c(1:5)]), atrisk_init = matrix(1, ncol = 38, nrow = 90),
#                        beta_init = matrix(0, ncol = 38, nrow = 10),
#                        ci_init = rep(0, 90), theta = 1, mu = 0,
#                        s2 = 10, s2_MH = 1, launch_iter = 1, r0g = 4, r1g = 1, r0c = 1, r1c = 1,
#                        thin = 100)
#     tot_time <- difftime(Sys.time(), start_time, units = "secs")
#     list(time = tot_time, result = clus_result)
# 
#   }
# stopImplicitCluster()
# difftime(Sys.time(), start_ova)
# 
# 
# saveRDS(testResult, paste0(path, "Manuscript/Result/microbiome_result_at_risk_8m_nb_1_s2_10_s2H_1_r0g_4_oneclus.RData"))

### Import the result
# m6s1 <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_at_risk_6m_nb_1_s2_1_s2H_1_r0g_4_singleton.RData"))
# sapply(1:7, function(x){mean(m6s1[[x]]$result$sm_accept)}) ### This is over 100,000 iterations
# sapply(1:7, function(x){apply(m6s1[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:1000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "Singleton: nB = 1, r0g = 4, s2 = 1")
# 
# as.numeric(colMeans(as.matrix(dat06[, -(1:5)])))
# meanSD(colMeans(as.matrix(dat06[, -(1:5)])))
# meanSD(log(colMeans(as.matrix(dat06[, -(1:5)]))))
# 
# m6s10one <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_at_risk_6m_nb_1_s2_10_s2H_1_r0g_4_oneclus.RData")) 
# sapply(1:7, function(x){mean(m6s10one[[x]]$result$sm_accept)}) ### This is over 100,000 iterations
# sapply(1:7, function(x){apply(m6s10one[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:1000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "One Cluster: nB = 1, r0g = 4, s2 = 1")
# 
# m6s10sin <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_at_risk_6m_nb_1_s2_10_s2H_1_r0g_4_singleton.RData")) 
# sapply(1:7, function(x){mean(m6s10sin[[x]]$result$sm_accept)}) ### This is over 100,000 iterations
# sapply(1:7, function(x){apply(m6s10sin[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:1000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "Singleton: nB = 1, r0g = 4, s2 = 1")
# 
# m8s10one <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_at_risk_8m_nb_1_s2_10_s2H_1_r0g_4_oneclus.RData")) 
# sapply(1:7, function(x){mean(m8s10one[[x]]$result$sm_accept)}) ### This is over 100,000 iterations
# sapply(1:7, function(x){apply(m8s10one[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:1000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "(8-Month) One cluster: nB = 1, r0g = 4, s2 = 1")

## microbiome_result_at_risk_8m_nb_1_s2_10_s2H_1_r0g_4_oneclus

### Import the result ----------------------------------------------------------
# result6m <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_at_risk_6m_hyper.RData")) 
# lapply(1:17, function(y){sapply(1:3, function(x){apply(result6m[[y]][[x]]$result$ci_result, 1, uniqueClus)}) %>%
#     as.data.frame() %>%
#     `colnames<-`(paste0("Chain ", 1:3)) %>%
#     mutate(Iteration = 1:1000, Case = paste0("Case ", y))}) %>%
#   bind_rows(.id = NULL) %>%
#   pivot_longer(!c(Iteration, Case), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   facet_wrap(factor(Case, levels = paste0("Case ", 1:17)) ~ .) +
#   theme_bw()

### Try the best set of hyperparameters with suitable hyperparameters: ---------
### ZIDM-ZIDM 
# set.seed(1415, kind = "L'Ecuyer-CMRG") 
# start_ova <- Sys.time()
# registerDoParallel(7) 
# resultZZ <- foreach(t = 1:20) %dopar% {
#     start_time <- Sys.time()
#     clus_result <- mod(iter = 500000, Kmax = 10, nbeta_split = 5,
#                        z = as.matrix(dat06[, -c(1:5)]), atrisk_init = matrix(1, ncol = 38, nrow = 90),
#                        beta_init = matrix(0, ncol = 38, nrow = 10),
#                        ci_init = rep(0, 90), theta = 1, mu = 0,
#                        s2 = 10, s2_MH = 1, launch_iter = 1, r0g = 4, r1g = 1, r0c = 1, r1c = 1,
#                        thin = 500)
#     tot_time <- difftime(Sys.time(), start_time, units = "secs")
#     list(time = tot_time, result = clus_result)
# 
#   }
# stopImplicitCluster() 
# difftime(Sys.time(), start_ova)
# saveRDS(resultZZ, paste0(path, "Manuscript/Result/microbiome_result_at_risk_6m_s2_1_s2H_1_r0g_4.RData"))

# resultZZ <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_at_risk_6m_s2_1_s2H_1_r0g_4.RData"))
# sapply(1:20, function(x){apply(resultZZ[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:1000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw()
# 
# sapply(1:20, function(x){salso(resultZZ[[x]]$result$ci_result[-(1:500), ])}) %>%
#   apply(2, uniqueClus) %>%
#   table()


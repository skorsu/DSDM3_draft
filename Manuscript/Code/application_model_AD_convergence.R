### Load libraries: ------------------------------------------------------------
library(foreach)
library(doParallel)
library(stringr)
library(tidyverse)
library(salso)
library(ggplot2)
library(coda)
library(latex2exp)
library(xtable)
library(mclustcomp)

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

month_analysis <- "8m"
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

### Final Cluster Assignment: --------------------------------------------------
salIndUngroup <- lapply(1:4, function(y){
  sapply(1:5, function(x){as.numeric(salso(listData[[y]][[x]]$result$ci_result[-(1:5000), ]))})})

salIndUngroup <- cbind(salIndUngroup[[1]], salIndUngroup[[2]], 
                       salIndUngroup[[3]], salIndUngroup[[4]])

salInd <- salIndUngroup %>%
  apply(2, table) %>% lapply(`length<-`, max(lengths(salIndUngroup %>%
                                                       apply(2, table)))) %>%
  lapply(as.numeric) %>%
  do.call(what = rbind) %>%
  as.data.frame()
  
colnames(salInd) <- paste0("Cluster ", 1:ncol(salInd))

#### salso Cluster for each chains
salInd %>%
  mutate(case = c(rep(caseName[1], 5), rep(caseName[2], 5), 
                  rep(caseName[3], 5), rep(caseName[4], 5)),
         chain = paste0("Chain ", rep(1:5, 4))) %>%
  pivot_longer(!(c(case, chain)), values_to = "size", names_to = "cluster") %>%
  ggplot(aes(x = chain, y = size, fill = cluster, label = size)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(. ~ factor(case, levels = c("One Cluster", "Singletons", "30 Clusters", "60 Clusters"))) +
  geom_text(size = 4, color = "white", position = position_stack(vjust = 0.5)) +
  labs(x = "Chain", y = "Cluster Size", fill = "Cluster",
       title = paste0(str_extract(month_analysis, "^[:digit:]+"), "-Month: Post-MCMC Cluster size with salso")) +
  theme_bw() +
  theme(legend.position = "bottom")

#### Adj. Rand Index
allComb <- expand.grid(1:20, 1:20)
caseChain <- paste0(c(rep(caseName[1], 5), rep(caseName[2], 5), 
                      rep(caseName[3], 5), rep(caseName[4], 5)), " \n ",
                    paste0("Chain ", rep(1:5, 4)))

data.frame(c1 = caseChain[allComb[, 1]], c2 = caseChain[allComb[, 2]],
           val = sapply(1:400, function(x){mclustcomp(salIndUngroup[, allComb[x, 1]], salIndUngroup[, allComb[x, 2]])[1, 2]})) %>%
  ggplot(aes(x = c1, y = c2, fill = val)) +
  geom_tile() +
  scale_fill_gradient(name = "ARI", low = "pink", high = "lightgreen", limit = c(0, 1)) +
  geom_text(aes(label = round(val, 2))) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = paste0(str_extract(month_analysis, "^[:digit:]+"), "-Month: Post-MCMC Adjusted Rand Index with salso"))




### Acceptance Rate: -----------------------------------------------------------
lapply(1:5, function(y){sapply(1:90, function(x){if(x > dim(r1[[y]]$result$MH_accept)[2]){return(c(NA, NA))} else{
  indexUpdate <- r1[[y]]$result$MH_accept[, x] != -1
  acceptVec <- r1[[y]]$result$MH_accept[indexUpdate, x]
  c(length(acceptVec), mean(acceptVec))}
  }) %>% t() %>% as.data.frame() %>%
  `colnames<-`(c("nPropose", "acceptRate")) %>%
  mutate(Cluster = paste0("Cluster ", 1:90), case = caseName[1], chain = paste0("Chain ", y))}) 






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



###: ---------------------------------------------------------------------------



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





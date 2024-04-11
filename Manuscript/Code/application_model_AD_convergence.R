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

month_analysis <- "12m"
caseName <- c("One Cluster", "Singletons", "30 Clusters", "60 Clusters")

listData <- list(r1 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_oneClus.RData")), 
                 r2 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_singleton.RData")), 
                 r3 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_init30.RData")), 
                 r4 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_init60.RData")))

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
lapply(1:4, function(z){lapply(1:5, function(y){sapply(1:90, function(x){if(x > dim(listData[[z]][[y]]$result$MH_accept)[2]){return(c(NA, NA))} else{
  indexUpdate <- listData[[z]][[y]]$result$MH_accept[, x] != -1
  acceptVec <- listData[[z]][[y]]$result$MH_accept[indexUpdate, x]
  c(length(acceptVec), mean(acceptVec))}
}) %>% t() %>% as.data.frame() %>%
    `colnames<-`(c("nPropose", "acceptRate")) %>%
    mutate(Cluster = paste0("Cluster ", 1:90), case = caseName[z], chain = paste0("Chain ", y))}) %>%
    bind_rows(.id = NULL)}) %>%
  bind_rows(.id = NULL) %>%
  mutate(caseChain = paste0(case, " \n ", chain)) %>%
  ggplot(aes(x = caseChain, y = factor(Cluster, levels = paste0("Cluster ", 90:1)), fill = nPropose)) +
  geom_tile() +
  geom_text(aes(label = round(acceptRate, 2)), color = "white", size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(fill = "Number of Propose",
       title = paste0(str_extract(month_analysis, "^[:digit:]+"), "-Month: Acceptance Rate for the cluster concentation."))

### Beta Convergence: ----------------------------------------------------------
intVar <- c(4, 8, 10, 16, 27, 36, 38)
lapply(1:4, function(z){lapply(1:5, function(y){lapply(which(colSums(listData[[z]][[y]]$result$MH_accept != -1) >= 10000),
                                                       function(x){data.frame(t(listData[[z]][[y]]$result$beta_result[x, intVar, ])) %>%
                                                           `colnames<-`(paste0("b", intVar)) %>%
                                                           mutate(iter = 1:15000, case = caseName[z], chain = paste0("Chain ", y), cluster = paste0("Cluster ", x)) %>%
                                                           pivot_longer(!(c(iter, case, chain, cluster)))}) %>%
    bind_rows(.id = NULL)}) %>%
    bind_rows(.id = NULL)}) %>%
  bind_rows(.id = NULL) %>%
  mutate(caseChainClus = paste0(case, " \n ", chain, " - ", cluster)) %>%
  ggplot(aes(x = iter, y = value, color = factor(name, levels = paste0("b", intVar)))) +
  geom_line() +
  facet_wrap(. ~ caseChainClus) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1)) +
  labs(color = "beta",
       title = paste0(str_extract(month_analysis, "^[:digit:]+"), "-Month: Traceplot for the beta (Active Cluster only)"))

### Save the result for further analysis: --------------------------------------
burn_in <- 5000
thin <- 100
mcmcIter <- seq(burn_in + thin, 15000, thin)

indexRun <- expand.grid(1:4, 1:5, mcmcIter)
indexRun <- indexRun %>% arrange(Var1, Var2, Var3)

resultList <- list(ci_result = lapply(1:4, function(y){lapply(1:5, function(x){listData[[y]][[x]]$result$ci_result[mcmcIter, ] %>%
    as.data.frame()}) %>%
    bind_rows(.id = NULL)}) %>%
      bind_rows(.id = NULL) %>%
      as.matrix(),
    atrisk_result = lapply(1:2000, function(x){listData[[indexRun[x, 1]]][[indexRun[x, 2]]]$result$atrisk_result[, , indexRun[x, 3]]}) %>%
      unlist() %>%
      array(c(90, 38, 2000)),
    beta_result = lapply(1:2000, function(x){listData[[indexRun[x, 1]]][[indexRun[x, 2]]]$result$beta_result[, , indexRun[x, 3]]}))

saveRDS(resultList,
        paste0(path, "Manuscript/Result/microbiome_result_", month_analysis, "_combined.RData"))
###: ---------------------------------------------------------------------------

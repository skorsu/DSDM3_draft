### Load libraries -------------------------------------------------------------
library(Rcpp)
library(foreach)
library(doParallel)
library(stringr)
library(tidyverse)
library(salso)
library(ggplot2)

### Import the data ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

datpath <- "/Users/kevin-imac/Desktop/Annika/"
if(! file.exists(datpath)){
  datpath <- "/Users/kevinkvp/Desktop/Annika/Application/"
}

### Data: 6 and 8 Months
ni68 <- read.csv(paste0(datpath, "Data/Nicaragua_6mo_8mo_genus.csv"))
ml68 <- read.csv(paste0(datpath, "Data/Mali_6mo_8mo_genus.csv"))

### Data: 12 Months
ni12 <- read.csv(paste0(datpath, "Data/Nicaragua_12mo_Metadata_csv.csv"))
ml12 <- read.csv(paste0(datpath, "Data/Mali_12mo_Metadata_csv.csv"))

### Result
result6mo <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_6m.RData"))
result8mo <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_8m.RData"))
result12mo <- readRDS(paste0(path, "Manuscript/Result/microbiome_result_12m.RData"))

### User-defined Functions -----------------------------------------------------
uniqueClus <- function(x){
  length(unique(x))
}

### MCMC: Convergence ----------------------------------------------------------
#### Trace plot
trPlot <- function(datList, monthLab){
  plotTitle <- paste0(monthLab, ": Active Clusters via MCMC with different initial cluster assignment")
  
  lapply(1:8, function(x){data.frame(Cluster = apply(datList[[x]]$MCMC$result$ci_result, 1, uniqueClus),
                                     Chain = paste0(datList[[x]]$init, ": Chain ", (x%%2) + 1),
                                     Iteration = 1:10000)}) %>%
    bind_rows(.id = NULL) %>%
    ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
    geom_line() +
    theme_bw() +
    scale_y_continuous(name = "Active Clusters", limits = c(1, 10), breaks = seq(1, 10)) +
    labs(title = plotTitle) +
    theme(legend.position = "bottom", legend.box = "horizontal")
}

trPlot(result6mo, "6 Month")
trPlot(result8mo, "8 Month")
trPlot(result12mo, "12 Month")

# ### salso: combine all chains --------------------------------------------------
# rbind(result8mo[[2]]$MCMC$result$ci_result[seq(8000, 10000, 10), ],
#       result8mo[[3]]$MCMC$result$ci_result[seq(8000, 10000, 10), ],
#       result8mo[[5]]$MCMC$result$ci_result[seq(8000, 10000, 10), ],
#       result8mo[[7]]$MCMC$result$ci_result[seq(8000, 10000, 10), ]) %>%
#   salso()
# 
# salso(result8mo[[8]]$MCMC$result$ci_result[seq(6000, 10000, 100), ]) %>% table()
# 
# result12mo <- vector("list", length = 8)
# result12mo[[1]] <- list(MCMC = result01[[4]], init = "One Cluster")
# result12mo[[2]] <- list(MCMC = result01[[5]], init = "One Cluster")
# result12mo[[3]] <- list(MCMC = result90[[2]], init = "Singleton")
# result12mo[[4]] <- list(MCMC = result90[[6]], init = "Singleton")
# result12mo[[5]] <- list(MCMC = result30[[1]], init = "30 Clusters")
# result12mo[[6]] <- list(MCMC = result30[[7]], init = "30 Clusters")
# result12mo[[7]] <- list(MCMC = result60[[2]], init = "60 Clusters")
# result12mo[[8]] <- list(MCMC = result60[[5]], init = "60 Clusters")
# saveRDS(result12mo, paste0(path, "Manuscript/Result/microbiome_result_12m.RData"))
# # 
# # 
# # 
# # 
# sapply(1:7, function(x){apply(result01[[x]]$result$ci_result, 1, uniqueClus)})[10000, ]
# sapply(c(4, 5), function(x){apply(result01[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:10000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "One Cluster: nB = 1, r0g = 4, s2 = 1")
# 
# sapply(1:7, function(x){apply(result90[[x]]$result$ci_result, 1, uniqueClus)})[10000, ]
# sapply(c(2, 6), function(x){apply(result90[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:10000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "One Cluster: nB = 1, r0g = 4, s2 = 1")
# 
# sapply(1:7, function(x){apply(result30[[x]]$result$ci_result, 1, uniqueClus)})[10000, ]
# sapply(c(1, 7), function(x){apply(result30[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:10000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "One Cluster: nB = 1, r0g = 4, s2 = 1")
# 
# sapply(1:7, function(x){apply(result60[[x]]$result$ci_result, 1, uniqueClus)})[10000, ]
# sapply(c(2, 5), function(x){apply(result60[[x]]$result$ci_result, 1, uniqueClus)}) %>%
#   as.data.frame() %>%
#   mutate(Iteration = 1:10000) %>%
#   pivot_longer(!c(Iteration), names_to = "Chain", values_to = "Cluster") %>%
#   ggplot(aes(x = Iteration, y = Cluster, color = Chain)) +
#   geom_line() +
#   theme_bw() +
#   labs(title = "One Cluster: nB = 1, r0g = 4, s2 = 1")





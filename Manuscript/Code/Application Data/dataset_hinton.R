# Required Libraries
library(selbal) # Dataset
library(tidyverse)
library(foreach)
library(doParallel)
library(salso)
library(ggplot2)
library(mclustcomp)

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

# Path
path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
resultpath <- paste0(path, "Result/selbal_HIV/")
# file.exists(resultpath)

# HIV dataset
datHIV <- HIV
statusHIV <- datHIV[, 62]
otuHIV <- datHIV[, 1:60]

# s2 Intuition: ----------------------------------------------------------------
nTaxa <- 10
n_iter <- 1000
s2List <- c(1e-3, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)

s2Intui <- lapply(1:length(s2List), function(y){
  sapply(1:n_iter, function(x){
    set.seed(x + 1)
    xiDum <- rnorm(nTaxa, 0, sqrt(s2List[y]))
    relativeDum <- exp(xiDum)/sum(exp(xiDum))
    var(exp(xiDum)/sum(exp(xiDum)))
  }) %>% as.data.frame()
}) %>%
  bind_cols() %>%
  mutate(Iter = 1:n_iter) %>%
  pivot_longer(!Iter)

s2Intui$name <- factor(s2Intui$name, levels = paste0("....", 1:length(s2List)),
                       labels = paste0("s2 = ", s2List))

s2Intui %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot()

nTaxaList <- c(10, 25, 50, 100, 250, 500)
s2List <- c(1e-3, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)
nS2 <- expand.grid(n = nTaxaList, s2 = s2List)
n_iter <- 1000

s2Intui <- lapply(1:nrow(nS2), function(y){
  
  lapply(1:n_iter, function(x){
    set.seed(x)
    xiDum <- rnorm(nS2[y, 1], 0, sqrt(nS2[y, 2]))
    relativeDum <- exp(xiDum)/sum(exp(xiDum))
    c(Iter = x, n = nS2[y, 1], s2 = nS2[y, 2], var = var(relativeDum))
  }) %>%
    bind_rows()
  
}) %>%
  bind_rows() %>%
  transmute(n = paste0("n = ", n), s2 = paste0("s2 = ", s2), var)

s2Intui$n <- factor(s2Intui$n, levels = paste0("n = ", nTaxaList))
s2Intui$s2 <- factor(s2Intui$s2, levels = paste0("s2 = ", s2List))

ggplot(s2Intui, aes(x = s2, y = var)) +
  geom_boxplot() +
  facet_wrap(. ~ n, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 270))

### Dummy: ---------------------------------------------------------------------
nTaxa <- 50
s2List <- c(1e-3, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)
set.seed(1)
x <- rnorm(nTaxa, 0, sqrt(1))

data.frame(j = 1:50, p = exp(x)/sum(exp(x))) %>%
  ggplot(aes(x = factor(j), y = p)) +
  geom_bar(stat = "identity")

otuHIV

### n_xi intuition: ------------------------------------------------------------
data.frame(x = paste0("OTU ", 1:60),
           actual_y = as.numeric(colSums(otuHIV)/sum(otuHIV))) %>%
  ggplot(aes(x = x, y = actual_y)) +
  geom_bar(stat = "identity")

x <- rnorm(60, 0, sqrt(1))
estimated_Y <- exp(x)/sum(exp(x))

var(as.numeric(colSums(otuHIV)/sum(otuHIV)))

data.frame(x = paste0("OTU ", 1:60), estimated_Y) %>%
  ggplot(aes(x = x, y = estimated_Y)) +
  geom_bar(stat = "identity")


n <- 60
s2 <- 1


# (otuHIV/rowSums(otuHIV)) %>%
#   colMeans() %>%
#   as.numeric() %>%
#   var()
# 
# (otuHIV/rowSums(otuHIV)) %>%
#   colMeans() %>%
#   as.numeric() %>%
#   sort(decreasing = TRUE)
# 
# otuHIV %>%
#   colMeans() %>%
#   as.numeric() %>%
#   var()

### Default set of Hyperparameters
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5, 
                      z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60), 
                      beta_init = matrix(0, nrow = 10, ncol = 60), 
                      ci_init = rep(0, 155), 
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3, 
                      t_thres = 2500, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_oneClus_defaultHyper.rds"))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

### Default set of Hyperparameters - Singleton
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 155, nbeta_split = 5, 
                      z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60), 
                      beta_init = as.matrix(otuHIV/rowSums(otuHIV)), 
                      ci_init = 0:154, 
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3, 
                      t_thres = 2500, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_singleton_defaultHyper.rds"))
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
                      z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60), 
                      beta_init = matrix(0, nrow = 10, ncol = 60), 
                      ci_init = rep(0, 155), 
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-5, 
                      t_thres = 2500, launch_iter = 30, 
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_oneClus_s2MH_1en5.rds"))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

### Post Analysis: -------------------------------------------------------------
# result <- readRDS(paste0(resultpath, "result_selbal_HIV_chain_", 1, "_init_oneClus_defaultHyper.rds"))

### Active Cluster - Default Case
registerDoParallel(6)
activeClusMat <- foreach(t = 1:6, .combine = cbind) %dopar% {
  result <- readRDS(paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_oneClus_defaultHyper.rds"))
  apply(result$mod$ci_result, 1, uniqueClus)
}
stopImplicitCluster()

activeClusMat %>%
  as.data.frame() %>%
  mutate(iter = 1:25000) %>%
  pivot_longer(!iter) %>%
  ggplot(aes(x = iter, y = value, color = name)) +
  geom_line() +
  theme_bw()

### Active Cluster - s2MH to 1e-5
registerDoParallel(6)
activeClusMat <- foreach(t = 1:6, .combine = cbind) %dopar% {
  result <- readRDS(paste0(resultpath, "result_selbal_HIV_chain_", t, "_init_oneClus_s2MH_1en5.rds"))
  apply(result$mod$ci_result, 1, uniqueClus)
}
stopImplicitCluster()

activeClusMat %>%
  as.data.frame() %>%
  mutate(iter = 1:25000) %>%
  pivot_longer(!iter) %>%
  ggplot(aes(x = iter, y = value, color = name)) +
  geom_line() +
  theme_bw()

### Check the convergence of xi
result <- readRDS(paste0(resultpath, "result_selbal_HIV_chain_", 5, "_init_oneClus_s2MH_1en5.rds"))
result$mod$beta_result[, 1, ] %>% t() %>%
  as.data.frame() %>%
  mutate(Iteration = 1:25000) %>%
  pivot_longer(!Iteration) %>%
  transmute(Iteration, xi = value,
            Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+"))) %>%
  ggplot(aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  theme_bw()

### Check the acceptance rate (xi overall)
sapply(1:10, function(x){sum(result$mod$MH_accept[, x] == 1)/sum(result$mod$MH_accept[, x] != -1)})

sapply(1:10, function(x){
  
  c(sum(result$mod$MH_accept[, x] != -1), ### Total Active
    sum(result$mod$MH_accept[, x] == 1), ### Total Accept
    sum(result$mod$MH_accept[1:2500, x] != -1), ### Active Before Adaptive
    sum(result$mod$MH_accept[1:2500, x] == 1)) ### Accept Before Adaptive
  
}) %>% t()

### Dummy: ---------------------------------------------------------------------
result <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 60, 
                       z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60), 
                       beta_init = rbind(log(as.numeric(colSums(otuHIV)/sum(otuHIV))),
                                         matrix(0, nrow = 9, ncol = 60)), 
                       ci_init = rep(0, 155), 
                       theta = 1, mu = 0, s2 = 2.5, s2_MH = 1, 
                       t_thres = 2000, launch_iter = 30, 
                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)

xiAccept <- sapply(1:10, function(x){
  
  c(sum(result$MH_accept[, x] != -1), ### Total Active
    sum(result$MH_accept[, x] == 1), ### Total Accept
    sum(result$MH_accept[1:2500, x] != -1), ### Active Before Adaptive
    sum(result$MH_accept[1:2500, x] == 1)) ### Accept Before Adaptive
  
}) %>% t()

(xiAccept[, 2]/xiAccept[, 1])[which(xiAccept[, 1] >= 2500)]
(xiAccept[, 4]/xiAccept[, 3])[which(xiAccept[, 1] >= 2500)]

result$beta_result[, 1, ] %>% t() %>%
  as.data.frame() %>%
  mutate(Iteration = 1:5000) %>%
  pivot_longer(!Iteration) %>%
  transmute(Iteration, xi = value,
            Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+"))) %>%
  ggplot(aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  theme_bw()

apply(result$ci_result, 1, uniqueClus) %>% plot(type = "l")

as.numeric(salso(result$ci_result[-(1:2000), ])) %>%
  mclustcomp(kmeans(otuHIV, 2)$cluster) %>%
  .[1, ]
  
  
  
  table(statusHIV)

table(kmeans(otuHIV, 2)$cluster, statusHIV)

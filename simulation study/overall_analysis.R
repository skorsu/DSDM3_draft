### Required Library
library(tidyverse)
library(salso)
library(mcclust)
library(mclustcomp)

### User defined function: -----------------------------------------------------
ms <- function(ms_mat, r = 4){
  paste0(round(colMeans(ms_mat), r), " (", 
         round(apply(ms_mat, 2, sd), 4), ")")
}

ctofac <- function(x){
  as.numeric(factor(x))
}

n_unique <- function(x){
  length(unique(x))
}

### Import the data: -----------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/result_local/"
case <- "strong_sig_with_overlap"
dat <- readRDS(paste0(path, "data_", case, ".RData"))
result_z <- readRDS(paste0(path, "result_", case, ".RData"))
result_sparse_DM <- readRDS(paste0(path, "result_sparseDM_", case, ".RData"))
result_shi_DM <- readRDS(paste0(path, "fix_shi_result_", case, ".RData"))
result_shi_MFM <- readRDS(paste0(path, "MFM_", case, ".RData"))

### Time: ----------------------------------------------------------------------
time_mat <- matrix(NA, nrow = 30, ncol = 6) ### ZZ, DZ, DP, MFM, DD, DsD
for(i in 1:30){
  time_mat[i, 1] <- result_z$zz_result[[i]]$zz_time
  time_mat[i, 2] <- result_z$dz_result[[i]]$dz_time
  ## time_mat[i, 3] <- result_shi_DM$fix_shi_result[[i]]$shi_time
  time_mat[i, 3] <- result_z$shi_result[[i]]$shi_time
  time_mat[i, 4] <- result_shi_MFM$list_MFM[[i]]$MFM_time
  time_mat[i, 5] <- result_z$dd_result[[i]][[4]]$dd_time
  time_mat[i, 6] <- result_sparse_DM$result_sDM[[i]]$time
}

ms(time_mat)

### VI: ------------------------------------------------------------------------

#### ZIDM-ZIDM
zz_VI <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  zz_VI[, i] <- result_z$zz_result[[i]]$zz_mod$assign[-(1:5000), ] %>%
    salso(maxNClusters = 10, loss = "VI") %>%
    as.numeric()
}

#### DM-ZIDM
dz_VI <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dz_VI[, i] <- result_z$dz_result[[i]]$dz_mod$assign[-(1:5000), ] %>%
    salso(maxNClusters = 10, loss = "VI") %>%
    as.numeric()
}

#### Shi (DP)
shi_dp_VI <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  # shi_dp_VI[, i] <- result_shi_DM$fix_shi_result[[i]]$shi_mod$crec %>%
  #   salso(maxNClusters = 10, loss = "VI") %>%
  #   as.numeric()
  shi_dp_VI[, i] <- result_z$shi_result[[i]]$shi_mod$crec %>%
    salso(maxNClusters = 10, loss = "VI") %>%
    as.numeric()
}

#### Shi (MFM)
shi_mfm_VI <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  shi_mfm_VI[, i] <- result_shi_MFM$list_MFM[[i]]$MFM_mod$crec %>%
    salso(maxNClusters = 10, loss = "VI") %>%
    as.numeric()
}

#### DM-DM
dd_VI <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dd_VI[, i] <- result_z$dd_result[[i]][[4]]$dd_mod[-(1:5000), ] %>%
    salso(maxNClusters = 10, loss = "VI") %>%
    as.numeric()
}

#### DM-sDM
dsd_VI <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dsd_VI[, i] <- result_sparse_DM$result_sDM[[i]]$assign[-(1:5000), ] %>%
    salso(maxNClusters = 10, loss = "VI") %>%
    as.numeric()
}


### Binder: --------------------------------------------------------------------

#### ZIDM-ZIDM
zz_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  zz_binder[, i] <- t(apply(result_z$zz_result[[i]]$zz_mod$assign[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder(max.k = 10) %>%
    .$cl
}

#### DM-ZIDM
dz_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dz_binder[, i] <- t(apply(result_z$dz_result[[i]]$dz_mod$assign[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder(max.k = 10) %>%
    .$cl
}

#### Shi (DP)
shi_dp_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  # shi_dp_binder[, i] <- t(apply(result_shi_DM$fix_shi_result[[i]]$shi_mod$crec, 1, ctofac)) %>%
  #   comp.psm() %>% 
  #   minbinder(max.k = 10) %>%
  #   .$cl
  shi_dp_binder[, i] <- t(apply(result_z$shi_result[[i]]$shi_mod$crec, 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder(max.k = 10) %>%
    .$cl
}

#### Shi (MFM)
shi_mfm_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  shi_mfm_binder[, i] <- t(apply(result_shi_MFM$list_MFM[[i]]$MFM_mod$crec, 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder(max.k = 10) %>%
    .$cl
}

#### DM-DM
dd_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dd_binder[, i] <- t(apply(result_z$dd_result[[i]][[4]]$dd_mod[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder(max.k = 10) %>%
    .$cl
}

#### DM-sDM
dsd_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dsd_binder[, i] <- t(apply(result_sparse_DM$result_sDM[[i]]$assign[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder(max.k = 10) %>%
    .$cl
}

### Analyze: -------------------------------------------------------------------
#### VI
zz_vi_sum <- matrix(NA, nrow = 30, ncol = 4)
dz_vi_sum <- matrix(NA, nrow = 30, ncol = 4)
dp_vi_sum <- matrix(NA, nrow = 30, ncol = 4)
mfm_vi_sum <- matrix(NA, nrow = 30, ncol = 4)
dd_vi_sum <- matrix(NA, nrow = 30, ncol = 3)
dsd_vi_sum <- matrix(NA, nrow = 30, ncol = 3)

for(i in 1:30){
  actual_ci <- dat$dat[[i]]$data$ci
  
  ### ZIDM-ZIDM
  zz_vi_sum[i, 1:3] <- mclustcomp(zz_VI[, i], actual_ci)[c(5, 1, 22), 2]
  zz_vi_sum[i, 4] <- n_unique(zz_VI[, i])
  
  ### DM-ZIDM
  dz_vi_sum[i, 1:3] <- mclustcomp(dz_VI[, i], actual_ci)[c(5, 1, 22), 2]
  dz_vi_sum[i, 4] <- n_unique(dz_VI[, i])
  
  ### Shi (DP)
  dp_vi_sum[i, 1:3] <- mclustcomp(shi_dp_VI[, i], actual_ci)[c(5, 1, 22), 2]
  dp_vi_sum[i, 4] <- n_unique(shi_dp_VI[, i])
  
  ### Shi (MFM)
  mfm_vi_sum[i, 1:3] <- mclustcomp(shi_mfm_VI[, i], actual_ci)[c(5, 1, 22), 2]
  mfm_vi_sum[i, 4] <- n_unique(shi_mfm_VI[, i])
  
  ### DM-DM
  dd_vi_sum[i, 1:3] <- mclustcomp(dd_VI[, i], actual_ci)[c(5, 1, 22), 2]
  
  ### DM-sDM
  dsd_vi_sum[i, 1:3] <- mclustcomp(dsd_VI[, i], actual_ci)[c(5, 1, 22), 2]
  
}

ms(zz_vi_sum)
ms(dz_vi_sum)
ms(dp_vi_sum)
ms(mfm_vi_sum)
ms(dd_vi_sum)
ms(dsd_vi_sum)

#### Binder
zz_bd_sum <- matrix(NA, nrow = 30, ncol = 4)
dz_bd_sum <- matrix(NA, nrow = 30, ncol = 4)
dp_bd_sum <- matrix(NA, nrow = 30, ncol = 4)
mfm_bd_sum <- matrix(NA, nrow = 30, ncol = 4)
dd_bd_sum <- matrix(NA, nrow = 30, ncol = 3)
dsd_bd_sum <- matrix(NA, nrow = 30, ncol = 3)

for(i in 1:30){
  actual_ci <- dat$dat[[i]]$data$ci
  
  ### ZIDM-ZIDM
  zz_bd_sum[i, 1:3] <- mclustcomp(zz_binder[, i], actual_ci)[c(5, 1, 22), 2]
  zz_bd_sum[i, 4] <- n_unique(zz_binder[, i])
  
  ### DM-ZIDM
  dz_bd_sum[i, 1:3] <- mclustcomp(dz_binder[, i], actual_ci)[c(5, 1, 22), 2]
  dz_bd_sum[i, 4] <- n_unique(dz_binder[, i])
  
  ### Shi (DP)
  dp_bd_sum[i, 1:3] <- mclustcomp(shi_dp_binder[, i], actual_ci)[c(5, 1, 22), 2]
  dp_bd_sum[i, 4] <- n_unique(shi_dp_binder[, i])
  
  ### Shi (MFM)
  mfm_bd_sum[i, 1:3] <- mclustcomp(shi_mfm_binder[, i], actual_ci)[c(5, 1, 22), 2]
  mfm_bd_sum[i, 4] <- n_unique(shi_mfm_binder[, i])
  
  ### DM-DM
  dd_bd_sum[i, 1:3] <- mclustcomp(dd_binder[, i], actual_ci)[c(5, 1, 22), 2]
  
  ### DM-sDM
  dsd_bd_sum[i, 1:3] <- mclustcomp(dsd_binder[, i], actual_ci)[c(5, 1, 22), 2]
  
}

ms(zz_bd_sum)
ms(dz_bd_sum)
ms(dp_bd_sum)
ms(mfm_bd_sum)
ms(dd_bd_sum)
ms(dsd_bd_sum)



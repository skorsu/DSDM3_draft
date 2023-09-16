### Required Library
library(salso)
library(mcclust)
library(mclustcomp)

### Import the data
path <- "/Users/kevinkvp/Desktop/result_local/"
case <- "strong_sig_no_overlap_zi"
dat <- readRDS(paste0(path, "data_", case, ".RData"))
result_z <- readRDS(paste0(path, "result_", case, ".RData"))
result_sparse_DM <- readRDS(paste0(path, "result_sparseDM_", case, ".RData"))
result_shi_DM <- readRDS(paste0(path, "fix_shi_result_", case, ".RData"))
result_shi_MFM <- readRDS(paste0(path, "MFM_", case, ".RData"))

### VI

### Binder: --------------------------------------------------------------------
ctofac <- function(x){
  as.numeric(factor(x))
}

#### ZIDM-ZIDM
zz_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  zz_binder[, i] <- t(apply(result_z$zz_result[[i]]$zz_mod$assign[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder() %>%
    .$cl
}

#### DM-ZIDM
dz_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dz_binder[, i] <- t(apply(result_z$dz_result[[i]]$dz_mod$assign[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder() %>%
    .$cl
}

#### Shi (DP)
shi_dp_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  shi_dp_binder[, i] <- t(apply(result_shi_DM$fix_shi_result[[i]]$shi_mod$crec, 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder() %>%
    .$cl
}

#### Shi (MFM)
shi_mfm_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  shi_mfm_binder[, i] <- t(apply(result_shi_MFM$list_MFM[[i]]$MFM_mod$crec, 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder() %>%
    .$cl
}

#### DM-DM
dd_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dd_binder[, i] <- t(apply(result_z$dd_result[[i]][[4]]$dd_mod[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder() %>%
    .$cl
}

#### DM-sDM
dsd_binder <- matrix(NA, nrow = 200, ncol = 30)
for(i in 1:30){
  dsd_binder[, i] <- t(apply(result_sparse_DM$result_sDM[[i]]$assign[-(1:5000), ], 1, ctofac)) %>%
    comp.psm() %>% 
    minbinder() %>%
    .$cl
}

### Time: ----------------------------------------------------------------------
ms <- function(ms_mat, r = 4){
  paste0(round(colMeans(ms_mat), r), " (", 
         round(apply(ms_mat, 2, sd), 4), ")")
}

time_mat <- matrix(NA, nrow = 30, ncol = 6) ### ZZ, DZ, DP, MFM, DD, DsD
for(i in 1:30){
  time_mat[i, 1] <- result_z$zz_result[[i]]$zz_time
  time_mat[i, 2] <- result_z$dz_result[[i]]$dz_time
  time_mat[i, 3] <- result_shi_DM$fix_shi_result[[i]]$shi_time
  time_mat[i, 4] <- 60 * result_shi_MFM$list_MFM[[i]]$MFM_time
  time_mat[i, 5] <- result_z$dd_result[[i]][[4]]$dd_time
  time_mat[i, 6] <- result_sparse_DM$result_sDM[[i]]$time
}

ms(time_mat)



table(dz_binder[, 1], zz_binder[, 1])
shi_dp_binder[, 1]


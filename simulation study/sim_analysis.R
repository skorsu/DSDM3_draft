### Required Libraries
library(mclustcomp)
library(xtable)
library(gridExtra)

### Function: Get the mean and SD
msd <- function(mc_mat, time_vec, r = 4){
  mm <- round(colMeans(cbind(mc_mat, time_vec)), r)
  ss <- round(apply(cbind(mc_mat, time_vec), 2, sd), r)
  paste0(mm, " (", ss, ")")
}

### Import the result from goose
im_path <- "/Users/kevin-imac/Desktop/"
mb_path <- "/Users/kevinkvp/Desktop/"
scenario_name <- "strong_sig_with_overlap"

path <- NULL
if(dir.exists(mb_path)){
  path <- mb_path
} else {
  path <- im_path
}

### Import the function for simulated data
source(paste0(path, "Github - Repo/ClusterZI/data/data_sim_DM.R"))

#### Data
pat_mat <- diag(20)[1:5, ]
pat_mat[3, 3] <- 0
pat_mat[4, 4] <- 0
pat_mat[5, 5] <- 0
pat_mat[3, c(10, 12)] <- 1
pat_mat[4, c(11, 12)] <- 1
pat_mat[5, c(10, 11)] <- 1
lplot <- simDM(n = 200, pattern = pat_mat, xi_conc = 10, pi_gm = 1, 
      pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, theta = 0.01) %>%
  simDM_sum(col_plot = "lightpink") %>%
  .$plot

grid.arrange(grobs = lplot)

### Normal Case 
result_list <- readRDS(paste0(path, scenario_name, ".RData"))
repli <- 30

### Number of zero
zero_mat <- matrix(NA, nrow = repli, ncol = 3)
for(i in 1:repli){
  zero_mat[i, ] <- result_list$actual_ci[[i]]$zero
}

colMeans(zero_mat)
apply(zero_mat, 2, sd) ## Proportion of zero, at-risk, structure zero

quantile(zero_mat[, 1])

### ZIDM-ZIDM and DM-ZIDM
zz_mat <- matrix(NA, ncol = 3, nrow = repli)
dz_mat <- matrix(NA, ncol = 3, nrow = repli)

for(i in 1:repli){
  zz_mat_dum <- mclustcomp(result_list$zz_salso[, i], 
                           result_list$actual_ci[[i]]$actual_ci)
  zz_mat[i, ] <- zz_mat_dum[c(1, 5, 22), 2]
  
  dz_mat_dum <- mclustcomp(result_list$dz_salso[, i], 
                           result_list$actual_ci[[i]]$actual_ci)
  dz_mat[i, ] <- dz_mat_dum[c(1, 5, 22), 2]
}

zz_normal <- msd(zz_mat, result_list$time_mat[, 1])
dz_normal <- msd(dz_mat, result_list$time_mat[, 2])

#### Zero-Inflated Case
scenario_name <- "strong_sig_with_overlap_zi"
result_list <- readRDS(paste0(path, scenario_name, ".RData"))

### Number of zero
zero_mat <- matrix(NA, nrow = repli, ncol = 3)
for(i in 1:repli){
  zero_mat[i, ] <- result_list$actual_ci[[i]]$zero
}

colMeans(zero_mat)
apply(zero_mat, 2, sd) ## Proportion of zero, at-risk, structure zero

quantile(zero_mat[, 1])

### ZIDM-ZIDM and DM-ZIDM
zz_mat <- matrix(NA, ncol = 3, nrow = repli)
dz_mat <- matrix(NA, ncol = 3, nrow = repli)

for(i in 1:repli){
  zz_mat_dum <- mclustcomp(result_list$zz_salso[, i], 
                           result_list$actual_ci[[i]]$actual_ci)
  zz_mat[i, ] <- zz_mat_dum[c(1, 5, 22), 2]
  
  dz_mat_dum <- mclustcomp(result_list$dz_salso[, i], 
                           result_list$actual_ci[[i]]$actual_ci)
  dz_mat[i, ] <- dz_mat_dum[c(1, 5, 22), 2]
}

zz_ZI <- msd(zz_mat, result_list$time_mat[, 1])
dz_ZI <- msd(dz_mat, result_list$time_mat[, 2])

zz_normal
dz_normal
zz_ZI
dz_ZI































### Required Libraries: --------------------------------------------------------
library(foreach)
library(doParallel)
library(salso)
library(mclustcomp)
library(tidyverse)

### Import the result: ---------------------------------------------------------
path <- "/Users/kevinkvp/Desktop/"
case_name <- "weak_sig_with_overlap_zi"

dat <- readRDS(paste0(path, "result_local/data_", case_name, ".Rdata"))
result <- readRDS(paste0(path, "result_local/result_", case_name, ".Rdata"))

### Functions: -----------------------------------------------------------------
### Function: Summarize the zero proportion
zero_summary <- function(data_list){
  zero_mat <- matrix(NA, ncol = 3, nrow = 30)
  for(i in 1:30){
    zero_mat[i, ] <- dat$dat[[i]]$zero
  }
  zero_mat
}

### Function: Retrive the cluster assignment for each methods (applying salso)
result_assign <- function(result_list){
  
  ### ZIDM-ZIDM, DM-ZIDM, and Shi
  zz_assign <- matrix(NA, ncol = 30, nrow = 200)
  dz_assign <- matrix(NA, ncol = 30, nrow = 200)
  shi_assign <- matrix(NA, ncol = 30, nrow = 200)
  
  for(i in 1:30){
    zz_assign[, i] <- as.numeric(salso(result_list$zz_result[[i]]$zz_mod$assign[-c(1:5000), ],
                                       maxNClusters = 10))
    dz_assign[, i] <- as.numeric(salso(result_list$dz_result[[i]]$dz_mod$assign[-c(1:5000), ],
                                       maxNClusters = 10))
    shi_assign[, i] <- as.numeric(salso(result_list$shi_result[[i]]$shi_mod$crec,
                                        maxNClusters = 10))
  }
  
  ### DM-DM
  registerDoParallel(5)
  dd_assign <- foreach(t = 1:30) %:%
    foreach(k = 1:9, .combine = cbind) %dopar% {
      
      salso(result$dd_result[[t]][[k]]$dd_mod[-c(1:5000), ])
      
    }
  stopImplicitCluster()
  
  list(zz_assign = zz_assign, dz_assign = dz_assign, 
       shi_assign = shi_assign, dd_assign = dd_assign)
  
}

### Function: Retrieve the actual cluster assignment
actual_ci <- function(list_data){
  actual_clus <- matrix(NA, nrow = 200, ncol = 30)
  for(i in 1:30){
    actual_clus[, i] <- dat$dat[[i]]$data$ci
  }
  actual_clus
}

### Function: Preprocessing for the DM-DM method
dd_procedure <- function(actual_clus, dmdm_clus){
  
  ## DM-DM: First, find the highest jaccard index for each replicated data
  dd_jac <- foreach(t = 1:30, .combine = rbind) %dopar% {
    apply(dmdm_clus[[t]], 2, 
          function(x, y){ as.numeric(mclustcomp(x, y, types = "jaccard")[, 2])}, 
          y = actual_clus[, t])
    }
  stopImplicitCluster()
  preferred_dd <- as.numeric(apply(dd_jac, 1, which.max))
  
  ### DM-DM: Second, create a matrix of cluster assignment based on the highest
  ### Jaccard Index
  dd_assign_prefer <- matrix(NA, ncol = 30, nrow = 200)
  for(i in 1:30){
    dd_assign_prefer[, i] <- dmdm_clus[[i]][, preferred_dd[i]]
  }
  
  list(dd_assign = dd_assign_prefer, preferred_dd = preferred_dd)

}

### Function: Calculate the Jaccard, Adjusted Rand, and VI
clus_summarize <- function(actual_clus, zz_clus, dz_clus, dd_clus, shi_clus){
  
  clus_array <- array(NA, dim = c(3, 4, 30))
  
  for(i in 1:30){
    clus_assign <- cbind(zz_clus[, i], dz_clus[, i], dd_clus[, i], shi_clus[, i])
    clus_array[, , i] <- apply(clus_assign, 2, 
                               function(x, y){mclustcomp(x, y)[c(5, 1, 22), 2]}, 
                               y = actual_clus[, i])
  }
  
  clus_array
  
}

### Function: Create a time matrix
time_mat <- function(result_list, dd_prefer){
  time_m <- matrix(NA, ncol = 4, nrow = 30)
  for(i in 1:30){
    time_m[i, ] <- c(result_list$zz_result[[i]]$zz_time,
                     result_list$dz_result[[i]]$dz_time,
                     result_list$dd_result[[i]][[dd_prefer[i]]]$dd_time,
                     result_list$shi_result[[i]]$shi_time)
  }
  time_m 
}

### Function: Balanced Quantities
ms <- function(mean_vec, sd_vec, r = 4){
  mm <- round(mean_vec, r)
  ss <- round(sd_vec, r)
  paste0(mm, " (", ss, ")")
}

### Function: Summarize the array matrics
bal_tab <- function(sum_array, time_m){
  time_mean <- colMeans(time_m)
  time_sd <- apply(time_m, 2, sd)
  
  final_result <- matrix(NA, ncol = 4, nrow = 4) ## Row = Method, Col = Time & Metric
  final_result[, 1] <- ms(time_mean, time_sd)
  
  for(i in 1:4){
    final_result[i, 2:4] <- ms(rowMeans(sum_array[, i, ]), 
                               apply(sum_array[, i, ], 1, sd))
  }
  
  colnames(final_result) <- c("Time", "Jaccard", "AdjRand", "VI")
  rownames(final_result) <- c("ZIDM-ZIDM", "DM-ZIDM", "DM-DM", "Shi")
  final_result
  
}

### Analyze: -------------------------------------------------------------------
assign <- result_assign(result)
actual <- actual_ci(dat)
dd_clus <- dd_procedure(actual, assign$dd_assign)
clus_sum <- clus_summarize(actual, assign$zz_assign, assign$dz_assign, 
                           dd_clus$dd_assign, assign$shi_assign)
time_mm <- time_mat(result, dd_clus$preferred_dd)
bal_tab(clus_sum, time_mm)

zero_s <- zero_summary(dat)
ms(colMeans(zero_s), apply(zero_s, 2, sd))

### Analyze (Only Shi): --------------------------------------------------------
path <- "/Users/kevinkvp/Desktop/"
case_name <- "weak_sig_with_overlap_zi"

dat <- readRDS(paste0(path, "result_local/data_", case_name, ".Rdata"))
result_shi <- readRDS(paste0(path, "result_local/fix_shi_result_", case_name, ".Rdata"))

registerDoParallel(5)
shi_assign <- foreach(t = 1:30, .combine = cbind) %dopar% {
  salso(result_shi$fix_shi_result[[t]]$shi_mod$crec, maxNClusters = 10)
}
stopImplicitCluster()

registerDoParallel(5)
shi_time <- foreach(t = 1:30, .combine = cbind) %dopar% {
  result_shi$fix_shi_result[[t]]$shi_time
}
stopImplicitCluster()

quan <- matrix(NA, ncol = 3, nrow = 30)
for(i in 1:30){
  quan[i, ] <- mclustcomp(dat$dat[[i]]$data$ci, shi_assign[, i])[c(5, 1, 22), 2] 
}

ms(c(mean(shi_time), colMeans(quan)), c(sd(shi_time), apply(quan, 2, sd)), r = 4)


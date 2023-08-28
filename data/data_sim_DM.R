### Required Libraries
library(tidyverse)
library(ggplot2)
library(dirmult)

### Function: simulate the zero-inflated DM.
simDM <- function(n, pattern, xi_conc, pi_gm, pi_c, z_sum, theta){
  
  ### The number of cluster
  K <- nrow(pattern)
  J <- ncol(pattern)
  
  ### Cluster assignment
  ci <- sample(1:K, n, replace = TRUE, prob = pi_c/sum(pi_c))
  
  ### xi matrix
  xi <- pattern * xi_conc
  xi[xi == 0] <- 1
  
  ### pi_gamma matrix
  pi_gm_mat <- matrix(pi_gm, ncol = J, nrow = K)
  
  ### gm_matrix and data
  gm_mat <- matrix(NA, ncol = J, nrow = n)
  z_mat <- matrix(NA, ncol = J, nrow = n)
  
  for(i in 1:n){
    gm_mat[i, ] <- rbinom(J, 1, pi_gm_mat[ci[i], ])
    z_mat[i, ] <- simPop(J = 1, n = z_sum, 
                         pi = (gm_mat[i, ] * xi[ci[i], ])/sum(gm_mat[i, ] * xi[ci[i], ]),
                         theta = theta)$data
  }
  
  list(K = K, J = J, xi = xi, ci = ci, gamma = gm_mat, z = z_mat, z_sum = z_sum)
  
}

### Function: Summarize the simulated dataset
simDM_sum <- function(simDM_list){
  
  ### Proportion of the zero
  ind_zero <- which(simDM_list$z == 0, arr.ind = TRUE)
  print(paste0("Number of zero: ", nrow(ind_zero)))
  print(paste0("Number of at-risk zero: ", sum(simDM_list$gamma[ind_zero] == 1)))
  print(paste0("Number of structure zero: ", sum(simDM_list$gamma[ind_zero] == 0)))
  
  sum_list <- c(nrow(ind_zero)/prod(dim(simDM_list$z)), 
                sum(simDM_list$gamma[ind_zero] == 1)/nrow(ind_zero),
                sum(simDM_list$gamma[ind_zero] == 0)/nrow(ind_zero))
  names(sum_list) <- c("P(zero)", "P(at risk|zero)", "P(structure|zero)")
  print(round(sum_list, 5))
  
  ### Plot
  dat_k <- simDM_list$z[which(simDM_list$ci == 1), ]
  title_text <- paste0("Cluster ", 1)
  data.frame(x = paste0("V", str_pad(1:simDM_list$J, ceiling(log10(simDM_list$J)), pad = "0")),
             y = colMeans(dat_k)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_bar(stat = "identity", fill = "mediumaquamarine") +
    theme_bw() +
    # ylim(0, 0.75 * simDM_list$z_sum) +
    labs(x = "Variables", y = "Count", title = title_text)

}

set.seed(1)
pat_matrix <- as.matrix(rbind(c(rep(1, 10), rep(0, 20)),
                              c(rep(0, 10), rep(1, 10), rep(0, 10)),
                              c(rep(0, 20), rep(1, 10))))
test <- simDM(n = 100, pattern = pat_matrix, xi_conc = 10, pi_gm = c(0.95, 0.9, 0.75), 
              pi_c = c(0.2, 0.5, 0.8), z_sum = 3000, theta = 0.01)
simDM_sum(test)

test$z[test$ci == 1, ]
test$z

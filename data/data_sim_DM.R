### Required Libraries
library(tidyverse)
library(ggplot2)
library(dirmult)
library(gridExtra)

### Function: simulate the zero-inflated DM.
simDM <- function(n, pattern, xi_conc, pi_gm, pi_c, z_sum_L, z_sum_U, theta){
  
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
    gm_vec <- rbinom(J, 1, pi_gm_mat[ci[i], ])
    gm_mat[i, ] <- as.numeric(gm_vec | pattern[ci[i], ])
    z_sum <- round(runif(1, z_sum_L, z_sum_U))
    z_mat[i, ] <- simPop(J = 1, n = z_sum, 
                         pi = (gm_mat[i, ] * xi[ci[i], ])/sum(gm_mat[i, ] * xi[ci[i], ]),
                         theta = theta)$data
  }
  
  list(K = K, J = J, xi = xi, ci = ci, gamma = gm_mat, z = z_mat)
  
}

### Function: Summarize the simulated dataset
simDM_sum <- function(simDM_list, col_plot = "mediumaquamarine"){
  
  ### Proportion of the zero
  ind_zero <- which(simDM_list$z == 0, arr.ind = TRUE)

  sum_list <- c(nrow(ind_zero)/prod(dim(simDM_list$z)), 
                sum(simDM_list$gamma[ind_zero] == 1)/nrow(ind_zero),
                sum(simDM_list$gamma[ind_zero] == 0)/nrow(ind_zero))
  names(sum_list) <- c("P(zero)", "P(at risk|zero)", "P(structure|zero)")
  
  ### Plot
  ylim_u <- data.frame(z = simDM_list$z, ci = simDM_list$ci) %>%
    group_by(ci) %>%
    summarise_all("mean") %>% 
    max() %>%
    plyr::round_any(10, f = ceiling)
  
  list_plot <- vector(mode = "list", length =  simDM_list$K)
  for(k in 1:simDM_list$K){
    dat_k <- simDM_list$z[which(simDM_list$ci == k), ]
    title_text <- paste0("Cluster ", k)
    list_plot[[k]] <- data.frame(x = paste0("V", str_pad(1:simDM_list$J, ceiling(log10(simDM_list$J)) + 1, pad = "0")),
                                 y = colMeans(dat_k)) %>%
      ggplot(aes(x = x, y = y)) +
      geom_bar(stat = "identity", fill = col_plot) +
      scale_x_discrete(guide = guide_axis(angle = 90))  + 
      theme_bw() +
      ylim(0, ylim_u) +
      labs(x = "Variables", y = "Count", title = title_text)
    
  }
  
  list(summary_zero = sum_list, plot = list_plot)

}

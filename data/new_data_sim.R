library(ggplot2)
library(gridExtra)
library(tidyverse)

### Function: Simulating the data
dat_sim <- function(n, J, K, scenario, xi_conc, pi_gm, sum_z){
  
  if(! scenario %in% 1:4){
    warning("the scenario does not existed. use the first scenario instead.")
    scenario <- 1
  }
  
  ### Generate the cluster assignment
  ci <- sample(1:K, n, replace = TRUE)
  
  ### Create an important matrix
  ind_peak <- sample(1:J, 3, replace = FALSE)
  if(scenario == 1){
    ind_peak <- sample(1:J, 5, replace = FALSE)
  } else if(scenario == 2){
    ind_peak <- NULL
  }
  
  mix_peak <- NULL
  mix_index <- NULL
  if(scenario == 2){
    mix_peak <- sample(1:J, replace = FALSE)
  } else if(scenario == 3){
    mix_index <- sample((1:J)[! 1:J %in% ind_peak], 3, replace = FALSE)
    t1 <- sample(mix_index, 2, replace = FALSE)
    t2 <- sample(mix_index, 2, replace = FALSE)
    while(sum(t2 %in% t1) == 2){
      t2 <- sample(mix_index, 2, replace = FALSE)
    }
    mix_peak <- c(t1, t2)
  } else if(scenario == 4){
    mix_index <- sample(ind_peak, 3, replace = FALSE)
    t1 <- sample(mix_index, 2, replace = FALSE)
    t2 <- sample(mix_index, 2, replace = FALSE)
    while(sum(t2 %in% t1) == 2){
      t2 <- sample(mix_index, 2, replace = FALSE)
    }
    mix_peak <- c(t1, t2)
  }
  
  ind_mat <- diag(J)[ind_peak, ]
  mix_mat <- diag(J)[mix_peak, ]
  imp_mat <- rbind(ind_mat, rowsum(mix_mat, as.integer(gl(nrow(mix_mat), 2, nrow(mix_mat)))))
  imp_mat <- matrix(imp_mat, ncol = J)
  
  ### Create a xi matrix
  xi <- imp_mat
  xi[xi == 0] <- xi_conc
  
  ### Create an at-risk indicator matrix and simulate the data
  gm <- matrix(NA, nrow = n, ncol = J)
  z <- matrix(NA, nrow = n, ncol = J)
  pi_gamma <- matrix(pi_gm, ncol = J, nrow = K)
  
  for(i in 1:n){
    gm[i, ] <- as.numeric(rbinom(J, 1, pi_gamma[ci[i], ]) | imp_mat[ci[i], ])
    z[i, ] <- rmultinom(1, sum_z, gm[i, ] * xi[ci[i], ])
  }
  
  list(ci = ci, clus_var = imp_mat, xi = xi, gamma = gm, z = z)

}

### Function: Plot
sim_plot <- function(n, J, K, scenario, xi_conc, pi_gm, sum_z){
  
  if(! scenario %in% 1:4){
    warning("the scenario does not existed. use the first scenario instead.")
    scenario <- 1
  }
  
  ### Simulate the data
  sim_list <- dat_sim(n, J, K, scenario, xi_conc, pi_gm, sum_z)
  
  ### Create a plot
  clus_plot <- vector(mode = "list", length = K)
  for(k in 1:K){
    title_text <- paste0("Cluster ", k)
    clus_plot[[k]] <- data.frame(x = paste0("X", str_pad(1:10, 2, pad = "0")), 
                                 y = colMeans(sim_list$z[which(sim_list$ci == k), ])) %>%
      ggplot(aes(x = x, y = y)) +
      geom_bar(stat = "identity", fill = "mediumaquamarine") +
      theme_bw() +
      ylim(0, 0.75 * sum_z) +
      labs(x = "Variables", y = "Count", title = title_text)
  }
  
  grid.arrange(grobs = clus_plot, layout_matrix = matrix(1:6, ncol = 2))
  
}

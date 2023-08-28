### Required Libraries
library(tidyverse)
library(dirmult)

### Pattern Matrix
pat_matrix <- as.matrix(rbind(c(rep(1, 10), rep(0, 20)),
                              c(rep(0, 10), rep(1, 10), rep(0, 10)),
                              c(rep(0, 20), rep(1, 10))))

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
  
  list(xi = xi, ci = ci, gamma = gm_mat, z = z_mat)
  
}

test <- simDM(n = 100, pattern = pat_matrix, xi_conc = 10, pi_gm = c(0.95, 0.9, 0.75), 
              pi_c = c(0.2, 0.5, 0.8), z_sum = 200, theta = 0.01)

test$z[test$ci == 1, ]
test$z

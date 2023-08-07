data_sim <- function(n, K, J_imp, pi_gm_mat, xi_scale, sum_zi){
  
  J <- ncol(pi_gm_mat)
  
  ### Generate the cluster assignment
  ci <- sample(1:K, n, replace = TRUE)
  
  ### Generate the at-risk indicator variables
  gm <- matrix(sapply(pi_gm_mat[ci, ], function(x){rbinom(1, 1, x)}), nrow = n)
  for(i in 1:n){
    gm[i, ci[i]] <- 1
  }
  
  ### Generate xi matrix
  xi_imp <- diag(J_imp)[1:K, ] * xi_scale
  xi_imp[xi_imp == 0] <- xi_scale/10
  xi_mat <- cbind(xi_imp, matrix(xi_scale/100, ncol = J - J_imp, nrow = K))
  
  gm_xi <- gm * xi_mat[ci, ]
  
  z <- matrix(NA, ncol = J, nrow = n)
  
  for(i in 1:n){
    
    ### Normalize the alpha vector
    norm_alp <- gm_xi[i, ]/sum(gm_xi[i, ])
    z[i, ] <- rmultinom(1, sum_zi, norm_alp)

  }

  list(ci = ci - 1, gamma = gm, beta = log(xi_mat), z = z)
  
}
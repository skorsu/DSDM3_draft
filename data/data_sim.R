data_sim <- function(n, K, J_imp, pi_gm_mat, xi_scale, sum_zi_imp, sum_zi_unimp){
  
  J <- ncol(pi_gm_mat)
  
  ### Generate the cluster assignment
  ci <- sample(1:K, n, replace = TRUE)
  
  ### Generate the at-risk indicator variables
  gm <- matrix(sapply(pi_gm_mat[ci, ], function(x){rbinom(1, 1, x)}), 
                  nrow = n)
  for(i in 1:n){
    gm[i, ci[i]] <- 1
  }
  
  ### Generate xi matrix
  xi_imp <- xi_scale * diag(J_imp)[1:K, ]
  xi_imp[xi_imp == 0] <- xi_scale/2
  cbind(xi_imp, matrix(runif(K * (J - J_imp)), nrow = K))
  xi_mat <- matrix(runif(K * J), nrow = K, ncol = J)
  diag(xi_mat) <- xi_scale
  
  gm_xi <- gm * xi_mat[ci, ]
  
  z <- matrix(NA, ncol = J, nrow = n)
  
  for(i in 1:100){
    
    gm_xi_imp_i <- gm_xi[i, 1:J_imp]/sum(gm_xi[i, 1:J_imp])
    gm_xi_unimp_i <- gm_xi[i, -(1:J_imp)]
    
    if(sum(gm_xi_unimp_i) != 0){
      gm_xi_unimp_i <- gm_xi_unimp_i/sum(gm_xi_unimp_i)
    }
    
    ### zi for the important variables
    zi_active <- rmultinom(1, sum_zi_imp, gm_xi_imp_i)
    zi_inactive <- gm_xi_unimp_i
    
    if(sum(zi_inactive) != 0){
      zi_inactive <- rmultinom(1, sum_zi_unimp, gm_xi_unimp_i)
    }
    
    z[i, ] <- c(zi_active, zi_inactive)
    
  }
  
  list(ci = ci - 1, gamma = gm, beta = log(xi_mat), z = z)
  
}


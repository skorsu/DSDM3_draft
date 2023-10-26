library(Rcpp)
library(dirmult)
library(salso)
library(foreach)
library(doParallel)
library(mclustcomp)

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")

### Simulate the data (Easiest Pattern)
data_sim <- function(n, pat_mat, pi_gamma, xi_conc, xi_non_conc, sum_z){
  K <- nrow(pat_mat)
  ci <- sample(1:K, n, replace = TRUE)
  
  pmat <- (pat_mat[ci, ] * xi_conc) + xi_non_conc
  at_risk_mat <- matrix(rbinom(n * ncol(pat_mat), 1, pi_gamma),
                        nrow = n, ncol = ncol(pat_mat))
  prob_mat <- matrix(NA, nrow = n, ncol = ncol(pat_mat))
  z <- matrix(NA, nrow = n, ncol = ncol(pat_mat))
  for(i in 1:n){
    prob_mat[i, ] <- as.numeric(at_risk_mat[i, ]) * pmat[i, ]
    z[i, ] <- simPop(J = 1, K = ncol(pat_mat), n = sum_z, 
                     pi = prob_mat[i, ]/sum(prob_mat[i, ]),
                     theta = 0.001)$data
  }

  return(list(ci = ci, at_risk_mat = at_risk_mat, prob_mat = prob_mat,
              z = z))
  
}

### Simulate the data
registerDoParallel(5)
data_sim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

### First
registerDoParallel(5)
result_clus <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  test_result <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$z,
                           theta_vec = rep(1, 10), launch_iter = 5,
                           MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                           r0c = 1, r1c = 1, print_iter = 2000)
  time_diff <- difftime(Sys.time(), start_time)
  return(list(result = test_result$assign, time_diff = time_diff))
}
stopImplicitCluster()

adj_rand <- rep(NA, 15)
for(i in 1:15){
  adj_rand[i] <- mclustcomp(as.numeric(salso(result_clus[[i]]$result)), 
                            data_sim[[i]]$ci)[1, 2]
}



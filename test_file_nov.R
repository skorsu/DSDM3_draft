library(Rcpp)
library(dirmult)
library(salso)
library(foreach)
library(doParallel)
library(mclustcomp)
library(cluster)
library(ecodist)
library(ggplot2)
library(plotrix)
library(latex2exp)
library(sparseMbClust)

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
# sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

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

### Shi's model
registerDoParallel(5)
result_shi <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  PerformClustering(otutable = t(data_sim[[t]]$z), ClusteringMethod = "DP", 
                    c = rep(0, 50), totaliter = 10000, burnin = 0, thin = 1)
}
stopImplicitCluster()

sapply(1:15, function(x){apply(result_shi[[x]]$crec, 1, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, lwd = 0.5, 
          main = "Active Cluster (Shi's model)", ylab = "# Active Cluster",
          xlab = "Iteration")

salso_result <- sapply(1:15, 
       function(x){as.numeric(salso(result_shi[[x]]$crec[-(1:5000), ]))})

adj_rand <- sapply(1:15,
                   function(x){mclustcomp(salso_result[, x], data_sim[[x]]$ci)[1, 2]})
mean(adj_rand)
sd(adj_rand)

### Our model ==================================================================
#### Assume that we initialize the cluster assignment to be a singleton
#### beta = 0 all i,j,k

margi_i <- log_marginal(zi = data_sim[[1]]$z[1, ], gmi = rep(1, 20), beta_k = rep(1, 20))
assign_prob <- log_sum_exp(log(c(1, rep(1 + 1, 49))) + margi_i)
sample(1:50, size = 1, prob = assign_prob)

rep()




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

### OK: When sampling the beta, the beta should be cluster specific.
### Need to: Sampling the beta
### Number along the way
### Try taking out the marginal out
### Collapse the cluster

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
# sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

### Function: Simulating the data
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

### Log marginal for each cluster
### Easy Case
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

beta_mat <- matrix(1, ncol = 20, nrow = 5)
test <- logmar_k(z = datsim[[1]]$z, atrisk = datsim[[1]]$at_risk_mat, beta_k = diag(20)[1, ])
test_all <- logmar(z = datsim[[1]]$z, atrisk = datsim[[1]]$at_risk_mat, beta_mat = diag(20))
test_all[, 1] == test
lgamma(2500 + 1)

ub <- update_beta(z = datsim[[1]]$z, atrisk = datsim[[1]]$at_risk_mat, 
                  beta_old = diag(20)[1:5, ], ci = datsim[[1]]$ci - 1, 
                  mu = 0, s2 = 1, s2_MH = 1)

ci_test <- sample(c(8, 5, 1), size = 50, replace = TRUE)

test <- update_ci(Kmax = 10, z = datsim[[1]]$z, atrisk = datsim[[1]]$at_risk_mat, 
                  beta_mat = diag(20)[1:10, ], ci_old = ci_test, 
                  theta = 1, mu = 0, s2 = 1)
test$split_index
test$ci_realloc[test$samp_ind + 1, ]
datsim[[1]]$z[test$ci_realloc == 2, ]


msamp <- matrix(NA, ncol = 10000, nrow = 4)
for(i in 1:10000){
  msamp[, i] <- rmultinom_1(probs_arma = c(0.1, 0.2, 0.3, 0.4))
}
rowMeans(msamp)


### Update 11/14 ===============================================================
#### Try Reallocate without SM and Marginal ====================================

### Easy Case
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

set.seed(1)
betmat <- matrix(runif(200), ncol = 20)
ci_rand <- sample(c(1, 3), size = 50, replace = TRUE)
test_beta <- array(NA, dim = c(3, 20, 10000))
for(i in 1:10000){
  test_beta[, , i] <- update_beta_lmar(z = datsim[[1]]$z, atrisk = datsim[[1]]$at_risk_mat, 
                         beta_mat = matrix(1, ncol = 20, nrow = 3), ci = datsim[[1]]$ci - 1, mu = 0, s2 = 1,
                         s2_MH = 1)
}

test <- realloc_beta(iter = 10000, Kmax = 10, z = datsim[[1]]$z, 
                     atrisk = datsim[[1]]$at_risk_mat, 
                     init_beta = betmat, init_ci = rep(0, 50),
                     theta = 1, mu = 0, s2 = 1, s2_MH = 1)

salso(test$ci_rec)
test$beta_rec[, , 5000]


plot(apply(test$ci_rec, 1, function(x){length(unique(x))}), type = "l")

matplot(exp(t(test$beta_rec[3, , ]))/rowSums(exp(t(test$beta_rec[3, , ]))), 
        type = "l")



sum(test_beta[1, 1, ] != 1)/10000






















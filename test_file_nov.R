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
### Take out the SM, and look at the log(observation) only.
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

registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  test <- realloc_n(Kmax = 50, iter = 10000, z = datsim[[t]]$z, 
                    init_assign = 0:49, theta_vec = rep(1, 50))
}
stopImplicitCluster() 

sapply(1:15, function(x){apply(result[[x]]$assign_result, 2, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, main = "Active Cluster: Easy Case", 
          xlab = "Iteration", ylab = "# Active Cluster", lwd = 3)

### More Difficult Case
patmat <- matrix(0, nrow = 3, ncol = 20)
patmat[1, 1] <- 1
patmat[2, c(1, 2)] <- 1
patmat[3, 1:3] <- 1

registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t + 10)
  test <- realloc_n(Kmax = 50, iter = 10000, z = datsim[[t]]$z, 
                    init_assign = 0:49, theta_vec = rep(1, 50))
}
stopImplicitCluster() 

sapply(1:15, function(x){apply(result[[x]]$assign_result, 2, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, main = "Active Cluster: Difficult Case", 
          xlab = "Iteration", ylab = "# Active Cluster", lwd = 3)

result[[1]]$process[, ,50]

#### Log Marginal ==============================================================

### Data Simulation: Easy Case
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

test <- logmar(z = datsim[[1]]$z, atrisk = matrix(1, ncol = 20, nrow = 50), 
               beta_mat = matrix(1, ncol = 20, nrow = 3)) ### Kmax = 3 with J = 20

rowSums(lgamma(datsim[[1]]$z + 1))

set.seed(1)
dattest <- matrix(round(runif(100)), ncol = 20, nrow = 5)

test <- logmar(z = dattest, atrisk = matrix(1, ncol = 20, nrow = 5), 
               beta_mat = matrix(1, ncol = 20, nrow = 7))
test
dim(test)

### Shi's Result ===============================================================

### Easiest Dataset
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

registerDoParallel(5)
result_shi <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  PerformClustering(otutable = t(datsim[[t]]$z), ClusteringMethod = "DP", 
                    c = rep(0, 50), totaliter = 10000, burnin = 0, thin = 1)
}
stopImplicitCluster()

sapply(1:15, function(x){apply(result_shi[[x]]$crec, 1, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, lwd = 0.5, 
          main = "Active Cluster (Shi's model) -- Easy", ylab = "# Active Cluster",
          xlab = "Iteration")

salso_result <- sapply(1:15, 
       function(x){as.numeric(salso(result_shi[[x]]$crec[-(1:5000), ], loss = binder()))})

adj_rand <- sapply(1:15,
                   function(x){mclustcomp(salso_result[, x], datsim[[x]]$ci)[1, 2]})
mean(adj_rand)
sd(adj_rand)

matplot(t(diag(20)[1:3, ]), type = "p", lty = 1, lwd = 0.75,
        main = "Signal: Easy", xlab = "Variable Index", ylab = "Signal")


### Compare to our model with random beta ~ Unif(0,1)
set.seed(1)
betmat <- matrix(runif(1000, 0, 1), ncol = 20, nrow = 50)
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_sm_nobeta(Kmax = 50, iter = 10000, z = datsim[[t]]$z, 
                              init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                              beta_mat = betmat, tau_vec = rep(1, 50), 
                              theta_vec = rep(1, 50), r0c = 1, r1c = 1, launch_iter = 5)
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ], loss = binder()))
  adj_rand <- mclustcomp(assign, datsim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()
sapply(1:15, function(x){result_salso[[x]]$adj_rand})

### More difficult Dataset (The one that we discussed during the last meeting)
patmat <- matrix(0, nrow = 3, ncol = 20)
patmat[1, 1] <- 1
patmat[2, c(1, 2)] <- 1
patmat[3, 1:3] <- 1

registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

matplot(t(patmat), type = "p", lty = 1, lwd = 0.75,
        main = "Signal: Difficult", xlab = "Variable Index", ylab = "Signal")

registerDoParallel(5)
result_shi <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  PerformClustering(otutable = t(datsim[[t]]$z), ClusteringMethod = "DP", 
                    c = rep(0, 50), totaliter = 10000, burnin = 0, thin = 1)
}
stopImplicitCluster()

sapply(1:15, function(x){apply(result_shi[[x]]$crec, 1, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, lwd = 0.5, 
          main = "Active Cluster (Shi's Model) -- Difficult", ylab = "# Active Cluster",
          xlab = "Iteration")

salso_result <- sapply(1:15, 
                       function(x){as.numeric(salso(result_shi[[x]]$crec[-(1:5000), ], loss = binder()))})

adj_rand <- sapply(1:15,
                   function(x){mclustcomp(salso_result[, x], datsim[[x]]$ci)[1, 2]})
mean(adj_rand)
sd(adj_rand)

### Compare to our model with random beta ~ Unif(0,1)
set.seed(1)
betmat <- matrix(runif(1000, 0, 1), ncol = 20, nrow = 50)
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_sm_nobeta(Kmax = 50, iter = 10000, z = datsim[[t]]$z, 
                              init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                              beta_mat = betmat, tau_vec = rep(1, 50), 
                              theta_vec = rep(1, 50), r0c = 1, r1c = 1, launch_iter = 5)
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ], loss = binder()))
  adj_rand <- mclustcomp(assign, datsim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()
mean(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))
sd(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))

### Label Switching ============================================================
table(labswitch(z = data_sim[[1]]$z, ci = data_sim[[1]]$ci), data_sim[[1]]$ci)

### Exp
test1 <- rep(NA, 1000)
for(i in 1:1000){
  test1[i] <- sample(0:3, size = 1, prob = log_sum_exp(c(0, log(4), 0, 0)))
}
table(test1)

### Our model ==================================================================
#### Since we know that random noise work, but is it because we are lucky to have that?
#### I will try many set of random beta and see.

### Easiest Dataset
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

registerDoParallel(5)
result <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  start <- Sys.time()
  betmat <- matrix(runif(1000, 0, 1), ncol = 20, nrow = 50)
  assign <- realloc_sm_nobeta(Kmax = 50, iter = 10000, z = datsim[[1]]$z, 
                              init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                              beta_mat = betmat, tau_vec = rep(1, 50), 
                              theta_vec = rep(1, 50), r0c = 1, r1c = 1, launch_iter = 5)
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, datsim[[1]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()
sapply(1:20, function(x){result_salso[[x]]$adj_rand})

table(result_salso[[1]]$assign, datsim[[1]]$ci)
apply(result[[1]]$assign$assign_result, 2, function(x){length(unique(x))}) |>
  plot(type = "l")

table(result_salso[[2]]$assign, datsim[[1]]$ci)
apply(result[[2]]$assign$assign_result, 2, function(x){length(unique(x))}) |>
  plot(type = "l")

table(result_salso[[3]]$assign, datsim[[1]]$ci)
apply(result[[3]]$assign$assign_result, 2, function(x){length(unique(x))}) |>
  plot(type = "l")

table(sm = result[[1]]$assign$sm_status, accept = result[[1]]$assign$sm_ac)
table(sm = result[[3]]$assign$sm_status, accept = result[[3]]$assign$sm_ac)


mean(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))
sd(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))



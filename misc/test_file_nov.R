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

# sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

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

### Max Cluster = 50
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(3 * (t + 10) - 1)
  test <- realloc_n(Kmax = 50, iter = 10000, z = datsim[[t]]$z, 
                    init_assign = rep(0, 50), theta_vec = rep(1, 50))
}
stopImplicitCluster() 

sapply(1:15, function(x){apply(result[[x]]$assign_result, 2, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, main = "Active Cluster: Easy Case (Kmax = 50)", 
          xlab = "Iteration", ylab = "# Active Cluster", lwd = 3,
          ylim = c(1, 50))

salso_clus <- sapply(1:15,
                     function(x){salso(t(result[[x]]$assign_result[, -(1:5000)]))})
sd(apply(salso_clus, 2, function(x){length(unique(x))}))


### Max Cluster = 10
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(3 * (t + 10) - 1)
  test <- realloc_n(Kmax = 10, iter = 10000, z = datsim[[t]]$z, 
                    init_assign = rep(0, 50), theta_vec = rep(1, 10))
}
stopImplicitCluster() 

sapply(1:15, function(x){apply(result[[x]]$assign_result, 2, 
                               function(y){length(unique(y))})}) |>
  matplot(type = "l", lty = 1, main = "Active Cluster: Easy Case (Kmax = 10)", 
          xlab = "Iteration", ylab = "# Active Cluster", lwd = 3,
          ylim = c(1, 50))

sapply(1:15,
       function(x){salso(t(result[[x]]$assign_result[, -(1:5000)]))})

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

#### Try Reallocate without SM (include the marginal) ==========================

### Data Simulation: Easy Case
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

### Test: all beta are same -- expected to see the same result as the previous case
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(3 * (t + 10) - 1)
  realloc_lmar(iter = 10000, Kmax = 50, z = datsim[[t]]$z, 
               atrisk = datsim[[t]]$at_risk_mat, 
               beta_mat = matrix(1, ncol = 20, nrow = 50), 
               init_ci = 0:49, theta = 1)
}
stopImplicitCluster() 

nactive <- sapply(1:15,
                  function(x){apply(result[[x]]$assign_result, 1, 
                                    function(y){length(unique(y))})})
matplot(nactive, type = "l", lty = 1, main = "Active Cluster: Easy Case", 
        xlab = "Iteration", ylab = "# Active Cluster", lwd = 3,
        ylim = c(1, 50))
abline(h = 35, col = "red")
abline(h = 34, col = "blue")

salso_assign <- sapply(1:15,
                       function(x){as.numeric(salso(result[[x]]$assign_result[-(1:5000), ]))})
mean(apply(salso_assign, 2, function(x){length(unique(x))}))
sd(apply(salso_assign, 2, function(x){length(unique(x))}))

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
  set.seed(3 * (t + 10) - 1)
  realloc_lmar(iter = 10000, Kmax = 50, z = datsim[[t]]$z, 
               atrisk = datsim[[t]]$at_risk_mat, 
               beta_mat = matrix(1, ncol = 20, nrow = 50), 
               init_ci = 0:49, theta = 1)
}
stopImplicitCluster() 

nactive <- sapply(1:15,
                  function(x){apply(result[[x]]$assign_result, 1, 
                                    function(y){length(unique(y))})})
matplot(nactive, type = "l", lty = 1, main = "Active Cluster: Difficult Case", 
        xlab = "Iteration", ylab = "# Active Cluster", lwd = 3,
        ylim = c(1, 50))

salso_assign <- sapply(1:15,
                       function(x){as.numeric(salso(result[[x]]$assign_result[-(1:5000), ]))})
mean(apply(salso_assign, 2, function(x){length(unique(x))}))
sd(apply(salso_assign, 2, function(x){length(unique(x))}))


### LogSumExp explanation
ni <- 1:50 ## Number of count in each clusters (Let's say we have Kmax = 50)
theta <- 1 ## theta
LSE_mat <- matrix(NA, ncol = 1000, nrow = 50)
for(i in 1:1000){
  lmar <- runif(1, -100, 100) ## Random number for the log marginal
  ## Calculate the reallocation probability by applying LogSumExp trick
  LSE_mat[, i] <- log_sum_exp(log(ni + theta) + lmar) 
}

### Check that all rows (same number of the observation in the cluster)
### has the same reallocation probability
length(unique(as.list(round(LSE_mat, 13)))) == 50

### Test: different beta
### Use the same dataset, but change the beta matrix

### Data Simulation: Easy Case
registerDoParallel(5)
datsim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

registerDoParallel(5)
result <- foreach(t = 1:30) %dopar% {
  set.seed(t)
  bmat <- matrix(runif(1000), ncol = 20, nrow = 50)
  mod_result <- realloc_lmar(iter = 10000, Kmax = 50, z = datsim[[1]]$z, 
                             atrisk = datsim[[1]]$at_risk_mat, 
                             beta_mat = bmat, init_ci = 0:49, theta = 1)
  return(list(bmat = bmat, mod_result = mod_result))
}
stopImplicitCluster()

adj_rand <- sapply(1:30,
                   function(x){mclustcomp(datsim[[1]]$ci,
                                          as.numeric(salso(result[[x]]$mod_result$assign_result[-(1:5000), ])))[1, 2]})
adj_rand

### Why Iteration 1 works perfect?
#### Look at the log marginal
table(actual = datsim[[1]]$ci, 
      highest_marginal = apply(result[[1]]$mod_result$log_marginal, 1, which.max))
#### By using the marginal, we can also detect the cluster correctly.

log_sum_exp(sort(result[[1]]$mod_result$log_marginal[1, ], decreasing = TRUE))

-218.9284 + log(50)

### Maximuize log marginal
lmar_max <- function(zi, gmi, bet_k){
  as.numeric(logmar(z = matrix(zi, ncol = 20), atrisk = matrix(gmi, ncol = 20), 
                    beta_mat = matrix(bet_k, ncol = 20)))
}
lmar_max(datsim[[1]]$z[1, ], rep(1, 20), c(50, rep(0, 19)))

### par maximize
optim_result <- optim(rep(0, 20), lmar_max, z = datsim[[1]]$z[1, ], 
                      gmi = rep(1, 20),
                      control = list(fnscale = -1), hessian = FALSE)
exp(optim_result$par[1:3])/sum(exp(optim_result$par[1:3]))
exp(result[[1]]$bmat[7, 1:3])/sum(exp(result[[1]]$bmat[7, 1:3]))

bb <- result[[1]]$bmat
result[[1]]$mod_result$log_marginal[datsim[[1]]$ci == 1, ]
bb[7, ]

result[[1]]$mod_result$log_marginal[datsim[[1]]$ci == 2, ]
bb[47, ]

bb[c(12, 41), ]


#### Consider the beta for the first three columns

imp_index <- apply(exp(bb[, 1:3])/rowSums(exp(bb[, 1:3])), 1, which.max)




lty_pattern <- rep(3, 17)
lty_pattern[3] <- 1
lwd_pattern <- rep(0.25, 17)
lwd_pattern[3] <- 1

matplot(t(exp(bb[imp_index == 1, ])/rowSums(exp(bb[imp_index == 1, ]))),
        type = "l", lty = lty_pattern, lwd = lwd_pattern)

which(which(imp_index == 2) == 46)

lty_pattern <- rep(3, 14)
lty_pattern[c(10, 12)] <- 1
lwd_pattern <- rep(0.25, 14)
lwd_pattern[10] <- 1

matplot(t(exp(bb[imp_index == 2, ])/rowSums(exp(bb[imp_index == 2, ]))),
        type = "l", lty = lty_pattern, lwd = lwd_pattern)


table(apply(result[[1]]$mod_result$log_marginal, 1, which.max), datsim[[1]]$ci)
table(apply(result[[2]]$mod_result$log_marginal, 1, which.max), datsim[[1]]$ci)

lty_pattern <- rep(3, 50)
lty_pattern[44] <- 1
lwd_pattern <- rep(0.25, 50)
lwd_pattern[44] <- 1

matplot(t(exp(result[[2]]$bmat)), type = "l", lty = lty_pattern, lwd = lwd_pattern)

matplot(t(exp(result[[2]]$bmat)/rowSums(exp(result[[2]]$bmat))), type = "l",
        lty = lty_pattern, lwd = lwd_pattern)

bb <- result[[2]]$bmat
bb_ratio <- exp(bb)/rowSums(bb)
matplot(t(bb_ratio[, 1:3]/rowSums(bb_ratio[, 1:3])), type = "l")
apply(result[[5]]$mod_result$log_marginal, 1, which.max)
table(apply(result[[5]]$mod_result$log_marginal, 1, which.max), datsim[[1]]$ci)
table(salso(result[[5]]$mod_result$assign_result[-(1:5000), ]), datsim[[1]]$ci)

bb <- result[[1]]$bmat
bb_ratio <- exp(bb)/rowSums(bb)
bb_ratio[, 1:3]/rowSums(bb_ratio[, 1:3])
lty_control <- rep(3, 50)
lty_control[c(7, 12, 41, 44)] <- 1
lwd_control <- rep(0.25, 50)
lwd_control[c(7, 12, 41, 44)] <- 1
matplot(t(bb_ratio[, 1:3]/rowSums(bb_ratio[, 1:3])), type = "l",
        lty = lty_control, lwd = lwd_control)
apply(result[[1]]$mod_result$log_marginal, 1, which.max)
## c(7, 12, 41, 44)
sort(exp(bb)[, 3]/rowSums(bb), index.return = TRUE, decreasing = TRUE)

apply((exp(bb)/rowSums(bb))[c(12, 21), ],1, which.max)

apply(bb, 1, which.max)
which(apply(exp(bb)/rowSums(bb), 1, which.max) == 1)
bb_ratio <- apply((exp(bb)/rowSums(bb))[c(7, 18), ], 1, 
                  function(x){sort(x, decreasing = TRUE)})
apply(bb_ratio, 2, function(x){(x[1] - x[2])/x[2]})

which(apply(exp(bb)/rowSums(bb), 1, which.max) == 12)
bb_ratio <- apply((exp(bb)/rowSums(bb))[c(27, 46), ], 1, 
                  function(x){sort(x, decreasing = TRUE)})
apply(bb_ratio, 2, function(x){(x[1] - x[2])/x[2]})


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



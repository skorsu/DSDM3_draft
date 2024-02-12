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

patmat <- matrix(0, ncol = 4, nrow = 3)
patmat[1, 1] <- 1
patmat[2, c(1, 2)] <- 1
patmat[3, c(1, 2, 3)] <- 1

### Simulate the data
registerDoParallel(5)
data_sim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

### The way that we simulate the data, pretty same for all clusters. Try make it more different.

colSums(data_sim[[1]]$z[data_sim[[1]]$ci == 1, ])/sum(data_sim[[1]]$z[data_sim[[1]]$ci == 1, ])
colSums(data_sim[[1]]$z[data_sim[[1]]$ci == 2, ])/sum(data_sim[[1]]$z[data_sim[[1]]$ci == 2, ])
colSums(data_sim[[1]]$z[data_sim[[1]]$ci == 3, ])/sum(data_sim[[1]]$z[data_sim[[1]]$ci == 3, ])

### Label switching should fix it. If it does not, find the explanation why adding random noise work.
### Choose the column with the higest vvariance with the raw data
### Check why adding noise is work?
### Why using same prior does not work? -- Likelihood?
### Check the likelihood
### Try with Shi model (Priority)

### 11/02: ---------------------------------------------------------------------
#### Try beta: -----------------------------------------------------------------
set.seed(1)
betmat <- matrix(0, ncol = 4, nrow = 50)
assign <- realloc_sm_nobeta(Kmax = 50, iter = 10000, z = data_sim[[1]]$z, 
                            init_assign = 0:49, gamma_mat = matrix(1, ncol = 4, nrow = 50), 
                            beta_mat = betmat, tau_vec = rep(1, 50), 
                            theta_vec = rep(15, 50), r0c = 1, r1c = 1, launch_iter = 10)

table(data_sim[[1]]$ci, as.numeric(salso(t(assign$assign_result[, -(1:5000)]))))

### Cluster 37
cc <- 6
n0 <- rep(NA, 10000)
pp <- matrix(NA, ncol = 20, nrow = 10000)
for(i in 1:10000){
  l0 <- which(assign$assign_result[, i] == cc)
  if(length(l0) == 1){
    pp[i, ] <- data_sim[[1]]$z[l0, ]/2500
  } else if(length(l0) > 1) {
    z0 <- data_sim[[1]]$z[(which(assign$assign_result[, i] == cc)), ]
    n0[i] <- dim(z0)[1]
    pp[i, ] <- colSums(z0)/sum(z0)
  } else {
    pp[i, ] <- rep(0, 20)
  }
}

matplot(pp, lty = 1, type = "l")

n44 <- rep(NA, 10000)
pp <- matrix(NA, ncol = 20, nrow = 10000)
for(i in 1:10000){
  l44 <- which(assign$assign_result[, i] == 44)
  if(length(l44) == 1){
    pp[i, ] <- data_sim[[1]]$z[l44, ]/2500
  } else if(length(l44) > 1) {
    z44 <- data_sim[[1]]$z[(which(assign$assign_result[, i] == 44)), ]
    n44[i] <- dim(z44)[1]
    pp[i, ] <- colSums(z44)/sum(z44)
  } else {
    pp[i, ] <- rep(0, 20)
  }
}

matplot(pp, lty = 1, type = "l")


plot(apply(assign$assign_result, 2, function(x){length(unique(x))}), type = "l")
table(data_sim[[1]]$ci, as.numeric(salso(t(assign$assign_result[, -(1:5000)]))))

#### During the meeting, I set the betamat to be the same for all dataset. "set.seed(1)" outside foreach
### Maybe we are lucky about the betmat, as when I try different randomization. It is not good.
set.seed(1)
betmat <- matrix(runif(1000), ncol = 20, nrow = 50)

registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_sm_nobeta(Kmax = 50, iter = 10000, z = data_sim[[t]]$z, 
                                 init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                                 beta_mat = matrix(runif(1000), ncol = 20, nrow = 50), tau_vec = rep(1, 50), 
                                 theta_vec = rep(1, 50), r0c = 1, r1c = 1, launch_iter = 5)
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Measure Time
comp_time <- sapply(1:15, function(x){as.numeric(result[[x]]$runtime)})
mean(comp_time)
sd(comp_time)

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = "SM",
        ylim = c(1, 50))
abline(h = 3, lty = "dotted")

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()
sapply(1:15, function(x){result_salso[[x]]$adj_rand})

set.seed(3)
bb <- matrix(runif(1000), ncol = 20, nrow = 50)
matplot(t(bb), lty = 1, lwd = 0.5, type = "l")

### Before 10/31: --------------------------------------------------------------
#### Try: Fix beta to some constant, same cluster ------------------------------

registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_test(Kmax = 10, iter = 10000, z = data_sim[[t]]$z, 
                         init_assign = rep(0, 50), gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                         beta_mat = matrix(0, ncol = 20, nrow = 10), 
                         tau_vec = rep(1, 10), theta_vec = rep(1, 10))
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = TeX("$\\beta_{jk} = 0 \\forall jk, c_{i} = 0 \\forall i$"))

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()

sapply(1:15, function(x){result_salso[[x]]$assign})
sapply(1:15, function(x){result_salso[[x]]$adj_rand})

table(result[[1]]$assign$assign_result[, 5001],
      result[[1]]$assign$assign_result[, 9001])

#### Try: Fix beta to some constant, singleton ---------------------------------
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_test(Kmax = 50, iter = 20000, z = data_sim[[t]]$z, 
                         init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                         beta_mat = matrix(0, ncol = 20, nrow = 50), 
                         tau_vec = rep(1, 10), theta_vec = rep(1, 10))
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = TeX("$\\beta_{jk} = 0 \\forall jk, c_{i} = i$"))

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()

sapply(1:15, function(x){result_salso[[x]]$assign})
sapply(1:15, function(x){result_salso[[x]]$adj_rand})
mean(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))
sd(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))

apply(sapply(1:15, function(x){result_salso[[x]]$assign}), 2, 
      function(y){length(unique(y))})

table(sapply(1:15, function(x){result_salso[[x]]$assign})[, 1],
      data_sim[[1]]$ci)

#### Try: beta = relative count, same cluster ----------------------------------
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_test(Kmax = 50, iter = 10000, z = data_sim[[t]]$z, 
                         init_assign = rep(0, 50), gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                         beta_mat = data_sim[[t]]$z/2500, 
                         tau_vec = rep(1, 50), theta_vec = rep(1, 50))
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = TeX("$\\beta_{jk} = relative, c_{i} = 0 \\forall i$"),
        ylim = c(1, 50))

ci1 <- as.numeric(salso(t(result[[1]]$assign$assign_result[, -(1:5000)])))
data_sim[[1]]$z[data_sim[[1]]$ci == 1 & ci1 == 1, ]/2500
data_sim[[1]]$z[data_sim[[1]]$ci == 1 & ci1 == 3, ]/2500
table(ci1, data_sim[[1]]$ci)

### Thought: Group by the relative count (strong signal) and the location of zero.
### Need to be exactly the same in terms of the location of zero.

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()

sapply(1:15, function(x){result_salso[[x]]$assign})
mean(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))
sd(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))

#### Try: beta = relative count, singleton -------------------------------------
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_test(Kmax = 50, iter = 10000, z = data_sim[[t]]$z, 
                         init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                         beta_mat = data_sim[[t]]$z/2500, 
                         tau_vec = rep(1, 50), theta_vec = rep(1, 50))
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = TeX("$\\beta_{jk} = relative, c_{i} = singleton$"),
        ylim = c(1, 50))

ci1 <- as.numeric(salso(t(result[[1]]$assign$assign_result[, -(1:5000)])))
data_sim[[1]]$z[data_sim[[1]]$ci == 1 & ci1 == 1, ]/2500
data_sim[[1]]$z[data_sim[[1]]$ci == 1 & ci1 == 3, ]/2500
table(ci1, data_sim[[1]]$ci)

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()

sapply(1:15, function(x){result_salso[[x]]$assign})
mean(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))
sd(sapply(1:15, function(x){result_salso[[x]]$adj_rand}))

### Investigate
set.seed(12)
c1 <- as.numeric(salso(t(result[[1]]$assign$assign_result[, -(1:5000)])))
table(actual = data_sim[[1]]$ci, model = c1)
z1 <- data_sim[[1]]$z[data_sim[[1]]$ci == 1, ]
c1 <- c1[data_sim[[1]]$ci == 1]

z1[order(c1), ]

rownames(z1) <- paste0("index: ", which(data_sim[[1]]$ci == 1), " | clus: ", c1)
loc_zero <- ifelse(z1 == 0, 1, 0)
heatmap(loc_zero[order(c1), ], Rowv = NA, Colv = NA) 


#### Try: beta is a pattern matrix, same cluster -------------------------------
par(mfrow = c(2, 2))
for(x in c(0.1, 0.25, 0.5, 1)){
  sigmat <- x * cbind(expand.grid(c(0, 1), c(0, 1), c(0, 1)), 
                          matrix(0, nrow = 8, ncol = 17))
  rel <- exp(sigmat)/rowSums(exp(sigmat))
  matplot(t(rel[c(2, 4, 8), ]), type = "b", ylim = c(0, 0.15), 
          xlab = "Variable", ylab = "Relative exp(beta)",
          main = paste0("s = ", x))
}

conct <- 0.1
sigmat <- conct * cbind(expand.grid(c(0, 1), c(0, 1), c(0, 1)), 
                        matrix(0, nrow = 8, ncol = 17))
rel <- exp(sigmat[2, ])/sum(exp(sigmat[2, ]))
plot(as.numeric(rel), ylim = c(0, 1), xlab = "Variable", 
     ylab = "Relative exp(beta)")
abline(h = rel[2], col = "red")

registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  assign <- realloc_test(Kmax = 8, iter = 10000, z = data_sim[[t]]$z, 
                         init_assign = rep(0, 50), gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                         beta_mat = as.matrix(sigmat), 
                         tau_vec = rep(1, 8), theta_vec = rep(1, 8))
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = TeX("$\\beta_{jk} = 0 \\forall jk, c_{i} = 0 \\forall i$"),
        ylim = c(1, 10))

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()

sapply(1:15, function(x){result_salso[[x]]$assign})
sapply(1:15, function(x){result_salso[[x]]$adj_rand})

#### Include Split-Merge: ------------------------------------------------------
#### Extension to (Try: beta = relative count)
registerDoParallel(5)
result <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start <- Sys.time()
  betmat <- matrix(runif(1000, 1, 3), ncol = 20, nrow = 50)
  assign <- realloc_sm_nobeta(Kmax = 50, iter = 10000, z = data_sim[[t]]$z, 
                              init_assign = 0:49, gamma_mat = matrix(1, ncol = 20, nrow = 50), 
                              beta_mat = betmat, tau_vec = rep(1, 50), 
                              theta_vec = rep(1, 50), r0c = 1, r1c = 1, launch_iter = 5)
  runtime <- difftime(Sys.time(), start)
  return(list(assign = assign, runtime = runtime))
}
stopImplicitCluster()

### Measure Time
comp_time <- sapply(1:15, function(x){as.numeric(result[[x]]$runtime)})
mean(comp_time)
sd(comp_time)

### Active Clusters
aclus_mat <- sapply(1:15, function(x){apply(result[[x]]$assign$assign_result, 2, 
                                            function(y){length(unique(y))})})
matplot(aclus_mat, type = "l", lty = 1, lwd = 0.5, xlab = "Iteration",
        ylab = "Active Clusters", col = 1:10,
        main = "SM",
        ylim = c(1, 50))
abline(h = 3, lty = "dotted")

### Cluster Assignment Performance
registerDoParallel(5)
result_salso <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  assign <- as.numeric(salso(t(result[[t]]$assign$assign_result)[-(1:5000), ]))
  adj_rand <- mclustcomp(assign, data_sim[[t]]$ci)[1, 2]
  return(list(assign = assign, adj_rand = adj_rand))
}
stopImplicitCluster()
sapply(1:15, function(x){result_salso[[x]]$adj_rand})

### Try: relabel the stuff (sum of j = 1 for all observations in each clusters)
### ()
### If it not, find the reason why noise help.

### Before 10/26: --------------------------------------------------------------
### Try pam with BrayCurtis
pamBC <- matrix(NA, nrow = 50, ncol = 15)
pamBC_adjR <- rep(NA, 15)
for(i in 1:15){
  pam_bc <- lapply(2:10, 
                   function(x, dat_mat){pam(bcdist(dat_mat), x)$silinfo$avg.width}, 
                   dat_mat = data_sim[[i]]$z) |>
    unlist() |>
    which.max() + 1
  pamBC[, i] <- pam(bcdist(data_sim[[i]]$z), pam_bc)$clustering
  pamBC_adjR[i] <- mclustcomp(pamBC[, i], data_sim[[i]]$ci)[1, 2]
}

### Our Model
registerDoParallel(5)
result_clus <- foreach(t = 1:5) %dopar% {
  set.seed(3)
  start_time <- Sys.time()
  test_result <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$z,
                           theta_vec = rep(1, 10), launch_iter = 10,
                           MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                           r0c = 1, r1c = 1, print_iter = 2000)
  time_diff <- difftime(Sys.time(), start_time)
  return(list(result = test_result$assign, time_diff = time_diff))
}
stopImplicitCluster()

### Comment

mod_result <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[1]]$z,
                        theta_vec = rep(1, 10), launch_iter = 10,
                        MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                        r0c = 1, r1c = 1, print_iter = 2000)
table(salso(mod_result$assign[-(1:5000), ]), data_sim[[1]]$ci)

adj_rand <- rep(NA, 5)
for(i in 1:5){
  adj_rand[i] <- mclustcomp(as.numeric(salso(result_clus[[i]]$result[-(1:5000), ])), 
                            data_sim[[i]]$ci)[1, 2]
}
mean(adj_rand)

nactive <- lapply(1:15,
       function(x){apply(result_clus[[x]]$result, 1, function(y){length(unique(y))})}) |>
  unlist() |>
  matrix(ncol = 15, byrow = FALSE)

par(mfrow = c(2, 2))
matplot(nactive[, which(adj_rand == 1)], type = "l", lty = 1, 
        ylim = c(1, 10), main = "Perfect Case", ylab = "# Active Cluster")
matplot(nactive[, which(adj_rand == 0)], type = "l", lty = 1,
        ylim = c(1, 10), main = "Worst Case", ylab = "# Active Cluster")
matplot(nactive[, which(0 < adj_rand & adj_rand < 1)], type = "l", lty = 1,
        ylim = c(1, 10), main = "Not-Good Case", ylab = "# Active Cluster")

### Update only beta
for(t in 1:15){
  set.seed(t)
  tt <- beta_mat_update(K = 3, iter = 10000, z = data_sim[[t]]$z, 
                        clus_assign = data_sim[[t]]$ci - 1, 
                        gm = data_sim[[t]]$at_risk_mat, 
                        mu = 0, s2 = 1, s2_MH = 1)
  
  path_pic <- "/Users/kevin-imac/Desktop/pic/"
  jpeg(paste0(path_pic, 'beta_rep', t, '.jpg'), width = 1532, height = 931,
       quality = 90)
  par(mfrow = c(1, 3))
  for(i in 1:3){
    plot_main <- paste0("Rep: ", t, " - Cluster: ", i)
    est_relative <- exp(t(tt[i, , ]))/rowSums(exp(t(tt[i, , ])))
    matplot(est_relative, type = "l", lty = 1, lwd = 1.25, ylab = "Relative Count", 
            main = plot_main)
    abline(h = colMeans(data_sim[[t]]$z[which(data_sim[[t]]$ci == i), ]/2500),
           lty = "dashed")
  }
  dev.off()
  
}

### Update only at-risk indicator
atrisk_q <- matrix(NA, nrow = 15, ncol = 3)
for(i in 1:15){
  app_beta <- rbind(colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == 1), ]/2500),
                    colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == 2), ]/2500),
                    colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == 3), ]/2500))
  tt <- atrisk_update(K = 3, iter = 10000, z = data_sim[[i]]$z, 
                      clus_assign = data_sim[[i]]$ci - 1, beta_mat = app_beta, 
                      r0g = 1, r1g = 1)
  zero_loc <- which(data_sim[[i]]$z == 0)
  atrisk_zero <- matrix(NA, nrow = 5000, ncol = length(zero_loc))
  for(t in 1:5000){
    atrisk_zero[t, ] <- as.numeric(tt[, , t + 5000])[zero_loc]
  }
  
  atrisk_q[i, ] <- quantile(apply(atrisk_zero, 2, mean), c(0.5, 0.025, 0.975))
}

colnames(atrisk_q) <- c("y", "low", "up")
atrisk_q <- data.frame(r = 1:15, atrisk_q)
plotCI(x = atrisk_q$r, y = atrisk_q$y, li = atrisk_q$low, ui = atrisk_q$up,
       ylab = "Rate", xlab = "Replicated Data Set")

### Update the cluster assignment

### Update both at-risk and beta
ab_list <- vector(mode = "list", length = 15)
for(i in 1:15){
  ab_list[[i]] <- beta_ar_update(K = 3, iter = 10000, z = data_sim[[i]]$z, 
                                 clus_assign = data_sim[[i]]$ci - 1, r0g = 1, r1g = 1, 
                                 mu = 0, s2 = 1, s2_MH = 1)
}


ttt <- beta_ar_update(K = 3, iter = 10000, z = data_sim[[1]]$z, 
               clus_assign = data_sim[[1]]$ci - 1, r0g = 1, r1g = 1, 
               mu = 0, s2 = 1, s2_MH = 1)
exp_r <- exp(t(ttt$beta[1, ,]))/rowSums(exp(t(ttt$beta[1, ,])))
matplot(exp_r, type = "l", lty = 1, lwd = 1.25)
abline(h = colMeans(data_sim[[1]]$z[which(data_sim[[1]]$ci == 1), ]/2500),
       lty = "dashed")


i <- 11
par(mfrow = c(1, 3))
for(k in 1:3){
  plot_main <- paste0("Rep: ", i, " | Cluster: ", k)
  rel_exp <- exp(t(ab_list[[i]]$beta[k, , ]))/rowSums(exp(t(ab_list[[i]]$beta[k, , ])))
  matplot(rel_exp, type = "l", lty = 1, lwd = 1.25, main = plot_main)
  abline(h = colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == k), ]/2500),
         lty = "dashed")
}





### Our Model (with 20,000 iterations)
registerDoParallel(5)
result_clus1 <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  test_result <- ZIDM_ZIDM(iter = 20000, K_max = 10, z = data_sim[[t]]$z,
                           theta_vec = rep(1, 10), launch_iter = 5,
                           MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                           r0c = 1, r1c = 1, print_iter = 2000)
  time_diff <- difftime(Sys.time(), start_time)
  return(list(result = test_result$assign, time_diff = time_diff))
}
stopImplicitCluster()

adj_rand1 <- rep(NA, 15)
for(i in 1:15){
  adj_rand1[i] <- mclustcomp(as.numeric(salso(result_clus1[[i]]$result[-(1:5000), ])), 
                             data_sim[[i]]$ci)[1, 2]
}

      
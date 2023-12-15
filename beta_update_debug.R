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
library(tidyverse)

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
              z = z, xi_mat = pmat))
  
}

## Try to update the beta ======================================================
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  patmat <- matrix(0, ncol = 50, nrow = 3)
  patmat[1, 1] <- 1
  patmat[2, 1:2] <- 1
  patmat[3, 1:3] <- 1
  dat <- data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
                  xi_conc = 10, xi_non_conc = 0.1, sum_z = 2500)
  list(dat = dat, patmat = patmat)
}
stopImplicitCluster()

registerDoParallel(5)
result <- foreach(t = 1:20) %dopar% {
  
  ### S1: beta_init = datsim[[t]]$patmat, ci_init = datsim[[t]]$dat$ci - 1
  ### S2: beta_init = datsim[[t]]$patmat, ci_init = sample(0:Kmax, 50, replace = TRUE)
  ### S3: Kmax = 50, beta_init = datsim[[t]]$dat$z/rowSums(datsim[[t]]$dat$z), ci_init = 0:49
  ### S4: Kmax = 10, beta_init = matrix(0, ncol = 50, nrow = 10), ci_init = rep(0, 50)
  
  set.seed(t)
  debug_brs(iter = 10000, Kmax = 10, nbeta_split = 5,
            z = datsim[[t]]$dat$z, 
            atrisk = datsim[[t]]$dat$at_risk_mat, 
            beta_init = matrix(0, ncol = 50, nrow = 10), 
            ci_init = rep(0, 50),
            theta = 10, mu = 0, s2 = 1, s2_MH = 10, launch_iter = 10, 
            r0c = 1, r1c = 1)
  
}
stopImplicitCluster()

### Plot
sapply(1:20, 
       function(x){apply(result[[x]]$ci_result, 1, function(y){length(unique(y))})}) |>
  matplot(type = "l", ylim = c(1, 10), ylab = "Active Clusters",
          xlab = "Iteration", main = "One cluster - theta = 10 and s2MH = 10 with n_beta = 5",
          lwd = 1.5)

### Individual Replicated Data

par(mfrow = c(2, 2))

for(i in 1:20){
  title_plot <- paste0("Case 6d: Replicate # ", i)
  adjRandindex <- mclustcomp(as.numeric(salso(result[[i]]$ci_result[-(1:5000), ])),
                             datsim[[i]]$dat$ci)[1, 2]
  nsplit <- sum(result[[i]]$sm_status)
  naccept_split <- sum(result[[i]]$sm_accept[result[[i]]$sm_status == 1])
    
  apply(result[[i]]$ci_result, 1, function(x){length(unique(x))}) %>%
  plot(type = "l", ylim = c(1, 10), main = title_plot, 
       xlab = "Iteration", ylab = "Number of the active cluster",
       sub = paste0("Adj Rand index = ", round(adjRandindex, 3), 
                    " -- P(Accept|Split) = ", naccept_split, "/", nsplit))
}

### Acceptance Probability
accept_prob <- function(result_list){
  overall_accept <- mean(result_list$sm_accept)
  merge_accept <- mean(result_list$sm_accept[result_list$sm_status == 0])
  split_accept <- mean(result_list$sm_accept[result_list$sm_status == 1])
  c(overall_accept, merge_accept, split_accept)
}

accept_prob(result[[1]])

rowMeans(sapply(1:20, function(x){accept_prob(result[[x]])}))
apply(sapply(1:20, function(x){accept_prob(result[[x]])}), 1, sd)

### Adjusted Rand Index
set.seed(3)
adj_rand <- sapply(1:20,
                   function(x){mclustcomp(as.numeric(salso(result[[x]]$ci_result[-(1:5000), ])),
                                          datsim[[x]]$dat$ci)[1, 2]})
mean(adj_rand)
sd(adj_rand)

table(as.numeric(salso(result[[5]]$ci_result[-(1:5000), ])), datsim[[5]]$dat$ci)

## Now: Include the at-risk zero part ==========================================
#### Note: Need to find the new way to simulated the data
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  patmat <- matrix(0, ncol = 50, nrow = 3)
  patmat[1, 1] <- 1
  patmat[2, 1:2] <- 1
  patmat[3, 1:3] <- 1
  dat <- data_sim(n = 50, pat_mat = patmat, pi_gamma = 0.95,
                  xi_conc = 10, xi_non_conc = 0.1, sum_z = 250)
  list(dat = dat, patmat = patmat)
}
stopImplicitCluster()

log_atrisk(zi = datsim[[1]]$dat$z[1, ], atrisk_i = datsim[[1]]$dat$at_risk_mat[1, ], 
           beta_k = datsim[[1]]$patmat[datsim[[1]]$dat$ci[1], ],
           r0g = 3, r1g = 2)
which(datsim[[1]]$dat$z[1, ] == 0)
datsim[[1]]$dat$z
datsim[[1]]$dat$ci



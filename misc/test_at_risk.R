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
library(pheatmap)

# sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

### Function: Simulating the data ----------------------------------------------
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

### Test DM-DM: Can it go down if we have a obvious pattern ====================
set.seed(34)
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

registerDoParallel(5)
result <- foreach(t = 1:20) %dopar% {
  DMDM(iter = 10000, Kmax = 10, z = datsim[[t]]$z, 
       beta_mat = matrix(0, ncol = 20, nrow = 10), ci_init = rep(0, 50), 
       theta = 1, thin = 100)
}
stopImplicitCluster()

resultCluster <- sapply(1:20, function(x){salso(result[[x]][-(1:50), ])})
table(apply(resultCluster, 2, function(x){length(unique(x))}))

### Try sparse
registerDoParallel(5)
resultSparse <- foreach(t = 1:20) %dopar% {
  DMDM(iter = 10000, Kmax = 10, z = datsim[[t]]$z, 
       beta_mat = matrix(0, ncol = 20, nrow = 10), ci_init = rep(0, 50), 
       theta = 0.01, thin = 10)
}
stopImplicitCluster()

resultSparseCluster <- sapply(1:20, function(x){salso(resultSparse[[x]][-(1:50), ])})
table(apply(resultSparseCluster, 2, function(x){length(unique(x))}))

### Test DM-x ==================================================================
set.seed(34)
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

### Test: DM-ZIDM
test <- DMZIDM(iter = 10000, Kmax = 50, z = datsim[[1]]$z, 
               beta_mat = matrix(0, ncol = 20, nrow = 50), 
               ci_init = rep(0, 50), 
               theta = 1, launch_iter = 10, r0c = 1, r1c = 9, 
               thin = 10)
plot(apply(test$ci_result, 1, function(x){length(unique(x))}))
table(salso(test$ci_result), datsim[[1]]$ci)

samp_obs <- sample(0:49, 2) ## all of them are form cluster 0
samp_clus <- sample(1:4, 1)
ci_launch <- rep(0, 50)
ci_launch[-(samp_obs)] <- sample(c(samp_clus, 0), 48, replace = TRUE)
S_set <- 0:49
S_set <- S_set[-(samp_obs + 1)]
table(ci_launch)
test_launch <- DMZIDM_launch_mcmc(z = datsim[[1]]$z, beta_mat = bm, ci_old = ci_launch,
                                  launch_iter = 10, S = S_set, 
                                  samp_clus = c(samp_clus, 0), nk = c(19, 31, 0, 0, 0))
realloc_new <- DMZIDM_realloc(z = datsim[[1]]$z, beta_mat = bm, 
                              ci_old = test_launch$ci_result, theta = 1)
table(realloc_new, datsim[[1]]$ci)
table(realloc_new)



which(datsim[[1]]$ci - 1 == 1)

i <- 17
datpart <- lfactorial(sum(datsim[[1]]$z[i, ])) - sum(lfactorial(datsim[[1]]$z[i, ]))
k <- 2
ck <- which(datsim[[1]]$ci[-i] == k)
zk_sum <- colSums(datsim[[1]]$z[-i, ][ck, ])
sum(lgamma(zk_sum + datsim[[1]]$z[i, ] + bm[k, ])) - 
  lgamma(sum(zk_sum + datsim[[1]]$z[i, ] + bm[k, ])) - 
  (sum(lgamma(zk_sum + bm[k, ])) - lgamma(sum(zk_sum + bm[k, ]))) + datpart
dm_data(Kmax = 5, zi = datsim[[1]]$z[i, ], z_not_i = datsim[[1]]$z[-i, ],
        ci_not_i = datsim[[1]]$ci[-i] - 1, beta_mat = bm)


test <- DMZIDM_realloc(z = datsim[[1]]$z, beta_mat = matrix(1, nrow = 5, ncol = 20), 
                       ci_old = sample(0:4, size = 50, replace = TRUE), theta = 1)

table(test, datsim[[1]]$ci)

test <- DMDM(iter = 5000, Kmax = 3, z = datsim[[1]]$z, 
             beta_mat = diag(20)[1:3, ] + 0.01, ci_init = datsim[[1]]$ci - 1,
             theta = 1, thin = 1000)

test <- mod(iter = 100000, Kmax = 3, nbeta_split = 5, 
            z = datsim[[1]]$z, 
            atrisk_init = matrix(1, ncol = 20, nrow = 50), 
            beta_init = diag(20)[1:3, ], 
            ci_init = datsim[[1]]$ci - 1, 
            theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, 
            r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1000)
table(salso(test$ci_result[-(1:20), ]), datsim[[1]]$ci)


dim(test)

salso(test[-(1:25), ])

table(as.numeric(salso(test[-(1:2000), ])), datsim[[1]]$ci)


### Data Simulation followed Shi's paper ---------------------------------------
data_sim_shi <- function(N, J, pi_gamma, z_case, aPhi = 1, bPhi = 9,
                         aLambda = 2, bLambda = 8, zsum){
  
  ### (a) Get the "marginal" probability
  pPhi <- rbeta(J/2, aPhi, bPhi)
  pLambda <- rbeta(J/2, aLambda, bLambda)
  
  ### (b) Shi's adjustment
  pPhi_A <- (1 - (z_case/5)) * pPhi
  pLambda_A <- ((sum(pLambda) + ((z_case/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  pPhi_B <- (1 + (z_case/5)) * pPhi
  pLambda_B <- ((sum(pLambda) - ((z_case/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  
  ### (c) At-risk indicator
  gamma_mat <- matrix(rbinom(N * J, 1, pi_gamma), nrow = N, ncol = J)
  
  ### (d) Calculate the probability matrix for each observations and normalize
  prob_mat <- rbind(t(matrix(rep(c(pPhi_A, pLambda_A), N/2), ncol = N/2)),
                    t(matrix(rep(c(pPhi_B, pLambda_B), N/2), ncol = N/2)))
  pg <- prob_mat * gamma_mat
  
  ### (f) Generate the Dirichlet random variables
  dat <- matrix(NA, ncol = J, nrow = N)
  for(i in 1:N){
    prob <- rdirichlet(1, 200 * pg[i, ]/sum(pg[i, ]))
    dat[i, ] <- rmultinom(1, zsum, prob)
  }
  
  dat
  
}

### Moderate Case of Shi's ------------------------------------------------------

#### Simulate the data
registerDoParallel(5)
datsim <- foreach(t = 1:5) %dopar% {
  set.seed(t)
  cbind(data_sim_shi(N = 50, J = 40, pi_gamma = 1, z_case = 1, 
               aPhi = 4, bPhi = 1, aLambda = 4, bLambda = 1, 
               zsum = 20000),
        data_sim_shi(N = 50, J = 10, pi_gamma = 1, z_case = 4, 
               aPhi = 4, bPhi = 1, aLambda = 4, bLambda = 1,
               zsum = 10000))
  
}
stopImplicitCluster()

#### Example of the data
pheatmap(datsim[[1]], display_numbers = F, 
         color = colorRampPalette(c('white','red'))(100), 
         cluster_rows = F, cluster_cols = F)

rbind(c(rep(0, 25), rep(1, 25)), 
      c(rep(1, 25), rep(0, 25)))

rbind(matrix(c(rep(0, 25), rep(1, 25)), ncol = 50, nrow = 25, byrow = T), 
      matrix(c(rep(1, 25), rep(0, 25)), ncol = 50, nrow = 25, byrow = T))

#### Run the model
registerDoParallel(5)
result <- foreach(t = 1:5) %dopar% {
  
  ## Truth for beta: beta_init = rbind(c(rep(0, 25), rep(1, 25)), c(rep(1, 25), rep(0, 25)))
  ## S1: Holding the truth (/) ci_init = c(rep(0, 25), rep(1, 25))
  ## S2: Random cluster assignment (/) ci_init = sample(c(0, 1), size = 50, replace = TRUE)
  ## S3: Singleton (X): beta_init = datsim[[t]]/15000, ci_init = 0:49, Kmax = 50
  ## For S3, need to adjust the hyperparameters
  ## rbind(matrix(c(rep(0, 25), rep(1, 25)), ncol = 50, nrow = 25, byrow = T), 
  ####### matrix(c(rep(1, 25), rep(0, 25)), ncol = 50, nrow = 25, byrow = T))
  set.seed(t)
  mod(iter = 30000, Kmax = 50, nbeta_split = 10, z = datsim[[t]], 
      atrisk_init = matrix(1, ncol = 50, nrow = 50), 
      beta_init = datsim[[t]]/30000, 
      ci_init = 0:49,
      theta = 1, mu = 0, s2 = 1, s2_MH = 1, launch_iter = 10, 
      r0g = 1, r1g = 1, r0c = 1, r1c = 1)
  
}
stopImplicitCluster()

### Plot
sapply(1:20, 
       function(x){apply(result[[x]]$ci_result, 1, function(y){length(unique(y))})}) |>
  matplot(type = "l", ylim = c(1, 5), ylab = "Active Clusters",
          xlab = "Iteration", main = " ",
          lwd = 1.5)

table(factor(result[[1]]$sm_status, levels = c(0, 1), labels = c("Merge", "Split")),
      factor(result[[1]]$sm_accept, levels = c(0, 1), labels = c("Reject", "Accept")))


### Adjusted Rand Index
set.seed(3)
adj_rand <- sapply(1:5,
                   function(x){mclustcomp(as.numeric(salso(result[[x]]$ci_result[-(1:15000), ])),
                                          c(rep(0, 25), rep(1, 25)))[1, 2]})
mean(adj_rand)
sd(adj_rand)





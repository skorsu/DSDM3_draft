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
              z = z, xi_mat = pmat))
  
}

### Data Simulation ------------------------------------------------------------
#### Case 1
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  patmat <- diag(10)[1:5, ]
  dat <- data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
                  xi_conc = 10, xi_non_conc = 0.1, sum_z = 1000)
  list(dat = dat, patmat = patmat)
}
stopImplicitCluster()

#### Case 1a
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  patmat <- diag(10)[sample(1:10, size = 5), ]
  dat <- data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
                  xi_conc = 10, xi_non_conc = 0.1, sum_z = 1000)
  list(dat = dat, patmat = patmat)
}
stopImplicitCluster()

#### Case 2 (reduced case: # variable = 5; n = 200 with J = 10)
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  patmat <- matrix(0, ncol = 20, nrow = 3)
  patmat[1, 1] <- 1
  patmat[2, 1:2] <- 1
  patmat[3, 1:3] <- 1
  dat <- data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
                  xi_conc = 10, xi_non_conc = 0.1, sum_z = 2500)
  list(dat = dat, patmat = patmat)
}
stopImplicitCluster()

#### Case 2a (reduced case: # variable = 5)
registerDoParallel(5)
datsim <- foreach(t = 1:20) %dopar% {
  set.seed(t)
  patmat <- matrix(0, ncol = 10, nrow = 3)
  patmat[1, sample(1:10, 1, replace = FALSE)] <- 1
  patmat[2, sample(1:10, 2, replace = FALSE)] <- 1
  patmat[3, sample(1:10, 3, replace = FALSE)] <- 1
  dat <- data_sim(n = 50, pat_mat = patmat, pi_gamma = 1,
                  xi_conc = 10, xi_non_conc = 0.1, sum_z = 2500)
  list(dat = dat, patmat = patmat)
}
stopImplicitCluster()

### Data Plot ------------------------------------------------------------------
datplot <- data.frame(datsim[[1]]$dat$z, datsim[[1]]$dat$ci) %>%
  group_by(datsim[[1]]$dat$ci) %>%
  summarise(across(everything(), mean)) %>%
  select(- c("datsim..1...dat.ci")) %>%
  rename(Cluster = "datsim[[1]]$dat$ci") %>%
  pivot_longer(-"Cluster", names_to = "variable")

datplot <- mutate(datplot, 
                  Variable = paste0("X", str_pad(str_match(datplot$variable, "[:digit:]+"), 3, pad = "0")))

datplot %>%
  ggplot(aes(x = Variable, y = Cluster, fill = value)) + 
  geom_tile() +
  theme_bw()

### Reallocation ---------------------------------------------------------------
### Run the model

### patmat for case 2 when all incorrect.
### patmat <- matrix(0, ncol = 20, nrow = 20)
### patmat[lower.tri(matrix(1, ncol = 20, nrow = 20), diag = TRUE)] <- 1

### patmat for case 2a when all incorrect.
### patmat <- datsim[[1]]$patmat
### while(nrow(patmat) < 20){
  
###   pm <- rbind(patmat, matrix(0, ncol = 20, nrow = 20 - nrow(patmat)))
###   for(k in (nrow(patmat) + 1):20){
###     nimp <- sample(1:3, 1)
###     pm[k, sample(1:20, size = nimp, replace = FALSE)] <- 1
###   }
###   patmat <- as.matrix(distinct(as.data.frame(pm)))
  
### }

### patmat[sample(1:20, size = 20, replace = FALSE), ]

registerDoParallel(5)
result <- foreach(t = 1:20) %dopar% {
  
  ### ci_init = datsim[[t]]$dat$ci - 1
  ### sample(0:4, size = 50, replace = TRUE)
  ### datsim[[t]]$patmat
  
  patmat <- datsim[[t]]$patmat
  while(nrow(patmat) < 20){
    
    pm <- rbind(patmat, matrix(0, ncol = 20, nrow = 20 - nrow(patmat)))
    for(k in (nrow(patmat) + 1):20){
      nimp <- sample(1:3, 1)
      pm[k, sample(1:20, size = nimp, replace = FALSE)] <- 1
    }
    patmat <- as.matrix(distinct(as.data.frame(pm)))
    
  }
  
  patmat[sample(1:20, size = 20, replace = FALSE), ]
  
  debug_r(iter = 1000, Kmax = 20, z = datsim[[t]]$dat$z, 
          atrisk = datsim[[t]]$dat$at_risk_mat, 
          beta_fixed = patmat, 
          ci_init = sample(0:19, size = 50, replace = TRUE), theta = 1)
}
stopImplicitCluster()

### Analyze the result
sapply(1:20, 
       function(x){apply(result[[x]], 1, function(y){length(unique(y))})}) |>
  matplot(type = "l", ylim = c(1, 10), ylab = "Active Clusters",
          xlab = "Iteration", main = "Case 2a: Incorrect Kmax, incorrect ci")

### Adjusted Rand Index
adj_rand <- sapply(1:20,
                   function(x){mclustcomp(as.numeric(salso(result[[x]][-(1:500), ])),
                                          datsim[[x]]$dat$ci)[1, 2]})
mean(adj_rand)
sd(adj_rand)

### Reallocation and Beta update -----------------------------------------------

### Run the model
registerDoParallel(5)
result <- foreach(t = 1:20) %dopar% {
  
  ### S1: beta_init = datsim[[t]]$patmat, ci_init = datsim[[t]]$dat$ci - 1
  ### S2: beta_init = datsim[[t]]$patmat, ci_init = sample(0:Kmax, 50, replace = TRUE)
  ### S3: Kmax = 50, beta_init = datsim[[t]]$dat$z/rowSums(datsim[[t]]$dat$z), ci_init = 0:49
  ### S4: Kmax = 50, beta_init = matrix(0, ncol = 10, nrow = 50), ci_init = rep(0, 50)
  ### S5: Kmax = 10, beta_init = matrix(0, ncol = 10, nrow = 10), ci_init = rep(0, 50)
  
  set.seed(t)
  debug_rb(iter = 10000, Kmax = 3, z = datsim[[t]]$dat$z, 
           atrisk = datsim[[t]]$dat$at_risk_mat, 
           beta_init = datsim[[t]]$patmat, 
           ci_init = datsim[[t]]$dat$ci - 1, 
           theta = 10, mu = 0, s2 = 1, s2_MH = 1)
  
}
stopImplicitCluster()

### Analyze the result
sapply(1:20, 
       function(x){apply(result[[x]]$ci_result, 1, function(y){length(unique(y))})}) |>
  matplot(type = "l", ylim = c(1, 10), ylab = "Active Clusters",
          xlab = "Iteration", main = "Case 2 (20 Variables) - S1 with theta = 10")

### Adjusted Rand Index
set.seed(3)
adj_rand <- sapply(1:20,
                   function(x){mclustcomp(as.numeric(salso(result[[x]]$ci_result[-(1:5000), ])),
                                          datsim[[x]]$dat$ci)[1, 2]})
mean(adj_rand)
sd(adj_rand)

for(i in 1:20){
  print(table(as.numeric(salso(result[[i]]$ci_result[-(1:5000), ])), datsim[[i]]$dat$ci))
}

###
which(adj_rand != 1)
colMeans(datsim[[2]]$dat$z[datsim[[2]]$dat$ci == 1, ]/2500)
colMeans(datsim[[2]]$dat$z[datsim[[2]]$dat$ci == 2, ]/2500)
colMeans(datsim[[2]]$dat$z[datsim[[2]]$dat$ci == 3, ]/2500)

### Beta: Accept-Reject
active_clus <- apply(result[[2]]$ci_result, 1, function(x){length(unique(x))})
which(active_clus == 3)
result[[2]]$ci_result[84, ]
exp(result[[2]]$beta_result[, , 86])/rowSums(exp(result[[2]]$beta_result[, , 86]))

apply(result[[2]]$ac_beta, 2, function(x){length(unique(x))})
table(factor(result[[2]]$ac_beta[, 3], levels = c(-1, 0, 1), labels = c("NA", "Reject", "Accept")))

matplot(exp(t(result[[2]]$beta_result[3, , ]))/rowSums(exp(t(result[[2]]$beta_result[3, , ]))),
        type = "l")

factor(result[[2]]$ac_beta[, 1], labels = c("NA", "Reject", "Accept"))
t(result[[1]]$beta_result[1, , ]) %>%
  head(10)
result[[1]]$ci_result[10, ]

datsim[[1]]$dat$z




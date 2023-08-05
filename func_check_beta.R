library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)
library(latex2exp)

### Import the external function
# source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(59)
sim_list <- data_sim(n = 100, K = 2, J_imp = 3, 
                     pi_gm_mat = matrix(c(1), ncol = 5, nrow = 2),
                     xi_scale = 5, sum_zi = 100)

table(sim_list$ci)
exp(sim_list$beta)

### Check: log_w function
### go through every active cluster
w <- c(1, 1, 1, 0, 0)


### Hand Calculation
beta_hand <- function(z, gm, bet, betk, s2){
  sum(dnorm(bet, mean = 0, sd = sqrt(s2), log = TRUE))
  
  xik <- exp(t(matrix(rep(betk, dim(z)[1]), nrow = length(betk))))
  gxk <- gm * xik
  z_gxk <- z + gxk
  
  sum(dnorm(bet, mean = 0, sd = sqrt(s2), log = TRUE)) + 
    sum(lgamma(apply(gxk, 1, sum))) - sum(lgamma(apply(z_gxk, 1, sum))) +
    sum(lgamma(z_gxk)) - sum(lgamma(gxk))
  
  #print(gxk)
  #print(z_gxk)
  
}

k <- 1
zk <- sim_list$z[which(sim_list$ci == k), which(w == 1)]
gmk <- sim_list$gamma[which(sim_list$ci == k), which(w == 1)]
betk <- sim_list$beta[(k+1), which(w == 1)]

beta_hand(z = zk, gm = gmk, bet = sim_list$beta[(k+1), ],
          betk = betk, s2 = 1)

log_beta_k(k = k, z = sim_list$z, clus_assign = sim_list$ci, 
           gamma_mat = sim_list$gamma, w = w, beta_k = sim_list$beta[(k+1), ], 
           s2 = 1)


### Sampler
set.seed(124)
beta_result <- debug_beta(iter = 10000, K = 2, z = sim_list$z, 
                          clus_assign = sim_list$ci, gamma_mat = sim_list$gamma, 
                          w = c(1, 1, 1, 0, 0), MH_var = 0.01, s2 = 1)

### Cluster 0
sim_list$z[which(sim_list$ci == 0), ]/100
pop_x1 <- apply(sim_list$z[which(sim_list$ci == 0), ], 2, sum)/sum(sim_list$z[which(sim_list$ci == 0), ])

#### MCMC
mcmc_b1 <- exp(t(beta_result[1, , 5001:10000]))
mcmc_x1 <- apply(mcmc_b1, 1, function(x){x/sum(x)}) %>% t()
apply(mcmc_x1, 2, sum)/sum(mcmc_x1)

par(mfrow = c(2, 3))
plot(mcmc_x1[, 1], type = "l", ylim = c(0, 1), 
     main = TeX("$\\xi_{10}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x1[1], b = 0, col = "red")
plot(mcmc_x1[, 2], type = "l", ylim = c(0, 1),
     main = TeX("$\\xi_{20}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x1[2], b = 0, col = "red")
plot(mcmc_x1[, 3], type = "l", ylim = c(0, 1), 
     main = TeX("$\\xi_{30}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x1[3], b = 0, col = "red")
plot(mcmc_x1[, 4], type = "l", ylim = c(0, 1),
     main = TeX("$\\xi_{40}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x1[4], b = 0, col = "red")
plot(mcmc_x1[, 5], type = "l", ylim = c(0, 1), 
     main = TeX("$\\xi_{50}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x1[5], b = 0, col = "red")

### Cluster 1
sim_list$z[which(sim_list$ci == 1), ]/100
pop_x2 <- apply(sim_list$z[which(sim_list$ci == 1), ], 2, sum)/sum(sim_list$z[which(sim_list$ci == 1), ])

#### MCMC
mcmc_b2 <- exp(t(beta_result[2, , 5001:10000]))
mcmc_x2 <- apply(mcmc_b2, 1, function(x){x/sum(x)}) %>% t()
apply(mcmc_x2, 2, sum)/sum(mcmc_x2)

par(mfrow = c(2, 3))
plot(mcmc_x2[, 1], type = "l", ylim = c(0, 1), 
     main = TeX("$\\xi_{11}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x2[1], b = 0, col = "red")
plot(mcmc_x2[, 2], type = "l", ylim = c(0, 1),
     main = TeX("$\\xi_{21}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x2[2], b = 0, col = "red")
plot(mcmc_x2[, 3], type = "l", ylim = c(0, 1), 
     main = TeX("$\\xi_{31}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x2[3], b = 0, col = "red")
plot(mcmc_x2[, 4], type = "l", ylim = c(0, 1),
     main = TeX("$\\xi_{41}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x2[4], b = 0, col = "red")
plot(mcmc_x2[, 5], type = "l", ylim = c(0, 1), 
     main = TeX("$\\xi_{51}$"), xlab = "Iteration", ylab = TeX("\\xi"))
abline(a = pop_x2[5], b = 0, col = "red")





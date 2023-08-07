library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)
library(latex2exp)

### Import the external function
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(184)
sim_list <- data_sim(n = 100, K = 3, J_imp = 5, 
                     pi_gm_mat = matrix(c(1), ncol = 7, nrow = 3),
                     xi_scale = 10, sum_zi = 100)

table(sim_list$ci)
exp(sim_list$beta)

10/14.2
1/14.2
0.1/14.2

### Check beta
set.seed(12345)
beta_result <- debug_beta(iter = 10000, K = 3, z = sim_list$z, 
                          clus_assign = sim_list$ci, gamma_mat = sim_list$gamma, 
                          w = c(1, 1, 1, 1, 1, 0, 0), MH_var = 0.01, s2 = 1)

### Cluster 0
b0_mcmc <- t(beta_result[1, , 5001:10000])
sum_x0 <- apply(exp(b0_mcmc), 1, sum)
norm_x0 <- data.frame(iter = 1:5000, exp(b0_mcmc)/sum_x0)
melt_norm_x0 <- melt(norm_x0 ,  id.vars = 'iter', variable.name = 'Variable')
ggplot(melt_norm_x0, aes(x = iter, y = value)) + 
  geom_line(aes(colour = Variable)) +
  theme_bw()

### Cluster Reallocate (without SM)
### Update both cluster assignment and beta
### keep gamma and xi to be fixed
table(sim_list$ci)
realloc_no_sm(K_max = 3, z = sim_list$z, clus_assign = sim_list$ci, 
              gamma_mat = sim_list$gamma, w = c(rep(1, 5), 0, 0), 
              beta_mat = sim_list$beta, theta = c(1, 1, 1))
t <- realloc_no_sm(K_max = 3, z = sim_list$z, clus_assign = sim_list$ci, 
                   gamma_mat = sim_list$gamma, w = c(rep(1, 5), 0, 0), 
                   beta_mat = beta_result[, , 10000], theta = c(1, 1, 1))
sim_list$ci

t$log_postpred[1, ]
sim_list$z[1, ]

lgamma(sum(exp(beta_result[3, 1:5, 10000]))) - 
  lgamma(sum(sim_list$z[1, 1:5] + exp(beta_result[3, 1:5, 10000]))) - 
  sum(lgamma(exp(beta_result[3, 1:5, 10000]))) + 
  sum(lgamma(sim_list$z[1, 1:5] + exp(beta_result[3, 1:5, 10000])))

library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)

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

ww <- matrix(NA, ncol = 5, nrow = 1000)
for(i in 1:1000){
  tt <- update_w(z = sim_list$z, clus_assign = sim_list$ci, gamma_mat = sim_list$gamma,
                 w = c(rep(1, 3), rep(0, 2)), beta_mat = sim_list$beta, launch_iter_w = 10,
                 r0w = 1, r1w = 1)
  ww[i, ] <- tt
}

mean(ac)
colMeans(ww)

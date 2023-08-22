library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)
library(latex2exp)
library(salso)
library(sparseMbClust)

### Import the external function
# source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Meeting: 08-17-2023
### n = 100; J = 10 (all variables are important); z = 2500; try 2 pi_g (1 and 0.75)

### Data Simulation
set.seed(2189)
sim_list <- data_sim(n = 100, K = 3, J_imp = 5, 
                     pi_gm_mat = matrix(c(0.5), ncol = 10, nrow = 5),
                     xi_scale = 10, sum_zi = 7500)

sum(sim_list$z == 0)

sim_list$z[which(sim_list$z == 0, arr.ind = TRUE)]
sim_list$gamma[which(sim_list$z == 0, arr.ind = TRUE)]

### Plot
dat <- cbind(index = 1:100, sim_list$z)
colnames(dat)[2:11] <- str_pad(1:10, 3, pad = "0")

t2 <- data.frame(dat)
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)

inner_join(data.frame(ci = factor(sim_list$ci), index = 1:100), 
           tt2, by = "index") %>%
  ggplot(aes(x = variable, y = value)) +
  geom_line(aes(color = ci, group = interaction(index, ci))) +
  theme_bw()

### Model: without at-risk part
KK <- 10
start_time <- Sys.time()
result_1 <- full_func(iter = 20000, K_max = KK, z = sim_list$z,
                      gamma_mat = matrix(1, nrow = dim(sim_list$z)[1], ncol = dim(sim_list$z)[2]), 
                      w = c(rep(1, 5), rep(0, 5)), theta = rep(1, KK), launch_iter = 10,
                      MH_var = 0.01, s2 = 1, r0c = 1, r1c = 1)
Sys.time() - start_time 

ci1 <- as.numeric(salso(t(result_1$ci)[-c(1:5000), ], maxNClusters = KK))
table(actual = sim_list$ci + 1, model = ci1)

### Model: with at-risk part
KK <- 10
start_time <- Sys.time()
result_2 <- full_func_at_risk(iter = 20000, K_max = KK, z = sim_list$z,
                              w = c(rep(1, 5), rep(0, 5)), theta = rep(1, KK), 
                              launch_iter = 10, MH_var = 0.01, s2 = 1, r0g = 1, 
                              r1g = 1, r0c = 1, r1c = 1)

Sys.time() - start_time
ci2 <- as.numeric(salso(t(result_2$ci)[-c(1:5000), ], maxNClusters = KK))
table(actual = sim_list$ci + 1, model = ci2)

plot(apply(result_2$ci, 2, function(x){length(unique(x))}), type = "l")

### Model: update all
KK <- 10
start_time <- Sys.time()
result_3 <- full_func_all(iter = 20000, K_max = KK, z = sim_list$z,
                          theta = rep(1, KK), launch_iter_w = 10, 
                          launch_iter = 10, MH_var = 0.01, s2 = 1, 
                          r0g = 1, r1g = 1, r0w = 1, r1w = 1, r0c = 1, r1c = 1)
Sys.time() - start_time 
ci3 <- as.numeric(salso(t(result_3$ci)[-c(1:5000), ], maxNClusters = KK))
table(actual = sim_list$ci + 1, model = ci3)

par(mfrow = c(2, 2))
plot(apply(result_1$ci, 2, function(x){length(unique(x))}), type = "l",
     main = "Number of the active clusters (beta, ci)", xlab = "Iteration", 
     ylab = "Active Clusters", ylim = c(1, 10))

plot(apply(result_2$ci, 2, function(x){length(unique(x))}), type = "l",
     main = "Number of the active clusters (gamma, beta, ci)", xlab = "Iteration", 
     ylab = "Active Clusters", ylim = c(1, 10))

plot(apply(result_3$ci, 2, function(x){length(unique(x))}), type = "l",
     main = "Number of the active clusters (all)", xlab = "Iteration", 
     ylab = "Active Clusters", ylim = c(1, 10))

scenario3 <- list(data = sim_list, 
                  result_no_at_risk = result_1, 
                  result_at_risk = result_2, 
                  result_all = result_3)

save(scenario3, file = "/Users/kevinkvp/Desktop/s3.RData")
rm(list = ls())

### Try Shi
load("/Users/kevinkvp/Desktop/s1.RData")

set.seed(1)
start_time <- Sys.time()
result_shi <- PerformClustering(t(scenario1$data$z), "DP", 
                                w = 1, totaliter = 20000, burnin = 5000, thin = 1)
Sys.time() - start_time 

salso(result_shi$crec)
table(salso(result_shi$crec, maxNClusters = 10), scenario3$data$ci + 1)

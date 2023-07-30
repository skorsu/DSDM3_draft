library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)

### Import the external function
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(12)
sim_list <- data_sim(n = 100, K = 3, J_imp = 4, 
                     pi_gm_mat = matrix(c(1), ncol = 10, nrow = 3),
                     xi_scale = 5, sum_zi_imp = 80, sum_zi_unimp = 20)

table(sim_list$ci)

### Plot
dat <- cbind(index = 1:100, sim_list$z)
colnames(dat)[2:11] <- str_pad(1:10, 3, pad = "0")
t2 <- data.frame(dat[which(sim_list$ci == 0), ])
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)
p0 <- ggplot(tt2, aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "", title = "Simulated Data: First Cluster")
t2 <- data.frame(dat[which(sim_list$ci == 1), ])
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)
p1 <- ggplot(tt2, aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "", title = "Simulated Data: Second Cluster")
t2 <- data.frame(dat[which(sim_list$ci == 2), ])
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)
p2 <- ggplot(tt2, aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "", title = "Simulated Data: Third Cluster")
grid.arrange(p0, p1, p2)

### Experiment 1: r0g = r1g = 1
set.seed(12)
start_time <- Sys.time()
test_gamma <- debug_gamma(iter = 10000, K = 3, z = sim_list$z, clus_assign = sim_list$ci - 1, 
                          w = c(rep(1, 4), rep(0, 6)), beta_mat = sim_list$beta, r0g = 1, r1g = 1)
Sys.time() - start_time

### The function will update only the first four columns as they are the only 
### important variables.

sum(sim_list$z[, 1:4] == 0)

index_mat <- cbind(which(sim_list$z[, 1:4] == 0) - floor(which(sim_list$z[, 1:4] == 0)/100) * 100,
                   ceiling(which(sim_list$z[, 1:4] == 0)/100))

for(i in 1:sum(sim_list$z[, 1:4] == 0)){
  print(sim_list$z[index_mat[i, 1], index_mat[i, 2]])
}

gamma_result_1 <- matrix(NA, ncol = 5, nrow = sum(sim_list$z[, 1:4] == 0))
for(i in 1:sum(sim_list$z[, 1:4] == 0)){
  gamma_result_1[i, 1] <- index_mat[i, 1]
  gamma_result_1[i, 2] <- index_mat[i, 2]
  gamma_result_1[i, 3] <- sim_list$z[index_mat[i, 1], index_mat[i, 2]]
  gamma_result_1[i, 4] <- sim_list$gamma[index_mat[i, 1], index_mat[i, 2]]
  gamma_result_1[i, 5] <- mean(test_gamma[index_mat[i, 1], index_mat[i, 2], 5001:10000])
}

colnames(gamma_result_1) <- c("row", "column", "z", "actual", "mcmc")
gamma_result_1

### Experiment 2: r0g = 2, r1g = 3
set.seed(12)
start_time <- Sys.time()
test_gamma <- debug_gamma(iter = 10000, K = 3, z = sim_list$z, clus_assign = sim_list$ci - 1, 
                          w = c(rep(1, 4), rep(0, 6)), beta_mat = sim_list$beta, r0g = 2, r1g = 3)
Sys.time() - start_time

### The function will update only the first four columns as they are the only 
### important variables.

sum(sim_list$z[, 1:4] == 0)

index_mat <- cbind(which(sim_list$z[, 1:4] == 0) - floor(which(sim_list$z[, 1:4] == 0)/100) * 100,
                   ceiling(which(sim_list$z[, 1:4] == 0)/100))

for(i in 1:sum(sim_list$z[, 1:4] == 0)){
  print(sim_list$z[index_mat[i, 1], index_mat[i, 2]])
}

gamma_result_2 <- matrix(NA, ncol = 5, nrow = sum(sim_list$z[, 1:4] == 0))
for(i in 1:sum(sim_list$z[, 1:4] == 0)){
  gamma_result_2[i, 1] <- index_mat[i, 1]
  gamma_result_2[i, 2] <- index_mat[i, 2]
  gamma_result_2[i, 3] <- sim_list$z[index_mat[i, 1], index_mat[i, 2]]
  gamma_result_2[i, 4] <- sim_list$gamma[index_mat[i, 1], index_mat[i, 2]]
  gamma_result_2[i, 5] <- mean(test_gamma[index_mat[i, 1], index_mat[i, 2], 5001:10000])
}

colnames(gamma_result_2) <- c("row", "column", "z", "actual", "mcmc")
gamma_result_2

### Experiment 3: r0g = 4, r1g = 1
set.seed(12)
start_time <- Sys.time()
test_gamma <- debug_gamma(iter = 10000, K = 3, z = sim_list$z, clus_assign = sim_list$ci - 1, 
                          w = c(rep(1, 4), rep(0, 6)), beta_mat = sim_list$beta, r0g = 4, r1g = 1)
Sys.time() - start_time

### The function will update only the first four columns as they are the only 
### important variables.

sum(sim_list$z[, 1:4] == 0)

index_mat <- cbind(which(sim_list$z[, 1:4] == 0) - floor(which(sim_list$z[, 1:4] == 0)/100) * 100,
                   ceiling(which(sim_list$z[, 1:4] == 0)/100))

for(i in 1:sum(sim_list$z[, 1:4] == 0)){
  print(sim_list$z[index_mat[i, 1], index_mat[i, 2]])
}

gamma_result_3 <- matrix(NA, ncol = 5, nrow = sum(sim_list$z[, 1:4] == 0))
for(i in 1:sum(sim_list$z[, 1:4] == 0)){
  gamma_result_3[i, 1] <- index_mat[i, 1]
  gamma_result_3[i, 2] <- index_mat[i, 2]
  gamma_result_3[i, 3] <- sim_list$z[index_mat[i, 1], index_mat[i, 2]]
  gamma_result_3[i, 4] <- sim_list$gamma[index_mat[i, 1], index_mat[i, 2]]
  gamma_result_3[i, 5] <- mean(test_gamma[index_mat[i, 1], index_mat[i, 2], 5001:10000])
}

colnames(gamma_result_3) <- c("row", "column", "z", "actual", "mcmc")
gamma_result_3
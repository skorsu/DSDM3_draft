library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)

### Import the external function
# source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(12)
sim_list <- data_sim(n = 100, K = 3, J_imp = 5, 
                     pi_gm_mat = matrix(c(1), ncol = 10, nrow = 3),
                     xi_scale = 10, sum_zi = 100)

table(sim_list$ci)

sim_list$beta

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

### Check: log_w function
### Check for all variables (j = 1,2, ..., 10)

### Check: log_gamma_ijk function
### Check for only z_ijk = 0

index_mat <- cbind(which(sim_list$z[, 1:4] == 0) - floor(which(sim_list$z[, 1:4] == 0)/100) * 100,
                   ceiling(which(sim_list$z[, 1:4] == 0)/100))

result_calc_mat <- matrix(NA, nrow = 11, ncol = 9)

for(i in 1:11){
  
  ### Calculate by hand
  gijk <- sim_list$gamma[index_mat[i, 1], index_mat[i, 2]]
  gi <- sim_list$gamma[index_mat[i, 1], 1:4]
  zi <- sim_list$z[index_mat[i, 1], 1:4]
  ci <- sim_list$ci[index_mat[i, 1]] + 1
  xi <- exp(sim_list$beta[ci, 1:4])
  
  hand_calc_1 <- lbeta(1 + gijk, 1 + (1 - gijk)) + lgamma(sum(gi * xi)) - lgamma(sum(zi + (gi*xi)))
  
  ### By using the function
  jj <- index_mat[i, 2] - 1
  func_calc_1 <- log_g_ijk(j = jj, zi = sim_list$z[index_mat[i, 1], ], gi = sim_list$gamma[index_mat[i, 1], ], 
                         w = c(rep(1, 4), rep(0, 6)), beta_k = sim_list$beta[ci, ],
                         r0g = 1, r1g = 1)
  
  ### Adjusted the gijk = (1 - gijk)
  gijk <- 1 - gijk
  gi[index_mat[i, 2]] <- gijk
  hand_calc_0 <- lbeta(1 + gijk, 1 + (1 - gijk)) + lgamma(sum(gi * xi)) - lgamma(sum(zi + (gi*xi)))
  
  ### By using the function
  gi_adj <- sim_list$gamma[index_mat[i, 1], ]
  gi_adj[index_mat[i, 2]] <- gijk
  func_calc_0 <- log_g_ijk(j = jj, zi = sim_list$z[index_mat[i, 1], ], gi = gi_adj, 
                           w = c(rep(1, 4), rep(0, 6)), beta_k = sim_list$beta[ci, ],
                           r0g = 1, r1g = 1)
  
  result_calc_mat[i, ] <- c(index_mat[i, 1], index_mat[i, 2], ci - 1, 1 - gijk, 
                            xi[index_mat[i, 2]], hand_calc_0, func_calc_0, hand_calc_1, func_calc_1)
  
}

result_mat_1 <- data.frame(result_calc_mat)
xtable(result_mat_1, digits = c(rep(0, 5), rep(5, 5)))

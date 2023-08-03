library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)

### Import the external function
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(12)
sim_list <- data_sim(n = 100, K = 3, J_imp = 5, 
                     pi_gm_mat = matrix(c(0.95), ncol = 10, nrow = 3),
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

(sim_list$gamma * exp(sim_list$beta[sim_list$ci + 1, ]))[, 1]

log_w(z = sim_list$z, clus_assign = sim_list$ci, gamma_mat = sim_list$gamma, 
      w = c(rep(1, 9), rep(0, 1)), beta_mat = sim_list$beta, r0w = 1, r1w = 1) -
  log_w(z = sim_list$z, clus_assign = sim_list$ci, gamma_mat = sim_list$gamma, 
      w = c(rep(1, 10), rep(0, 0)), beta_mat = sim_list$beta, r0w = 1, r1w = 1)

exp(-2506.781)
lgamma(10)
lgamma(1 + c(rep(1, 5), rep(0, 5))) + lgamma(1 + (1 - c(rep(1, 5), rep(0, 5))))

-10 * lgamma(3)

### Check for only z_ijk = 0

sum(log(dnorm(sim_list$beta[1, ])))

update_w(z = sim_list$z, clus_assign = sim_list$ci, gamma_mat = sim_list$gamma,
         w = c(rep(1, 5), rep(0, 5)), beta_mat = sim_list$beta, r0w = 1, r1w = 1)

test_result <- debug_w(iter = 10000, K = 3, z = sim_list$z, clus_assign = sim_list$ci, 
                       gamma_mat = sim_list$gamma, beta_mat = sim_list$beta, r0w = 1, r1w = 1)

plot(apply(test_result, 1, sum), type = "l")

update_beta(z = sim_list$z, clus_assign = sim_list$ci, gamma_mat = sim_list$gamma,
            w = rep(1, 10), beta_mat = sim_list$beta, MH_var = 0.01, s2 = 1)


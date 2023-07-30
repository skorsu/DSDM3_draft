library(tidyverse)
library(reshape2)

### Import the external function
# source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(12)
sim_list <- data_sim(n = 100, K = 3, J_imp = 4, 
                     pi_gm_mat = matrix(c(0.5), ncol = 10, nrow = 3),
                     xi_scale = 5, sum_zi_imp = 80, sum_zi_unimp = 20)

### Check
table(sim_list$ci)

dat <- cbind(index = 1:100, sim_list$z)
colnames(dat)[2:11] <- str_pad(1:10, 3, pad = "0")
t2 <- data.frame(dat[which(sim_list$ci == 2), ])
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)
ggplot(tt2, aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "")

### Debugging: at-risk
set.seed(12)
start_time <- Sys.time()
test_gamma <- debug_gamma(iter = 10000, K = 3, z = sim_list$z, clus_assign = sim_list$ci - 1, 
                          w = c(rep(1, 4), rep(0, 6)), beta_mat = sim_list$beta, r0g = 1, r1g = 1)
Sys.time() - start_time

### The function will update only the first four columns as they are the only 
### important variables.

gamma_hat <- matrix(NA, ncol = 4, nrow = 100)
for(i in 1:100){
  gamma_hat[i, ] <- rowMeans(test_gamma[i, 1:4, 5001:1000])
}

gamma_hat[which(sim_list$ci == 1), ][1:4, ]
sim_list$gamma[which(sim_list$ci == 1), 1:4][1:4, ]
sim_list$z[which(sim_list$ci == 1), 1:4][1:4, ]

sim_list$z[, 1]
sim_list$z[, 3]

sim_list$z[1, 1:4]
sim_list$ci[1]
sim_list$gamma[1, 1:4] * exp(sim_list$beta[2, 1:4])
sim_list$z[1, 1:4] + (sim_list$gamma[1, 1:4] * exp(sim_list$beta[2, 1:4]))

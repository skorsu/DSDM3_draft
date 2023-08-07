library(tidyverse)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(xtable)
library(latex2exp)
library(salso)

### Import the external function
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Data Simulation
set.seed(1291)
sim_list <- data_sim(n = 250, K = 5, J_imp = 5, 
                     pi_gm_mat = matrix(c(1), ncol = 10, nrow = 5),
                     xi_scale = 10, sum_zi = 500)

### Observed data relative abundances
ORA <- matrix(NA, ncol = 10, nrow = 5)
for(k in 0:4){
  ORA[(k+1), ] <- colSums(sim_list$z[which(sim_list$ci == k), ])/sum(sim_list$z[which(sim_list$ci == k), ])
}

ORA

### Cluster Reallocate (without SM)
### Update both cluster assignment and beta while keep gamma and xi to be fixed

set.seed(1284)
start_time <- Sys.time()
test_result <- clus_no_SM(iter = 10000, K_max = 5, z = sim_list$z, 
                          gamma_mat = sim_list$gamma, w = c(rep(1, 5), rep(0, 5)),
                          MH_var = 0.01, s2 = 1, theta_vec = rep(1, 5))
Sys.time() - start_time ### 4.53407 seconds

### Cluster Assignment
mcmc_ci <- salso(test_result$ci[5001:10000, ], maxNClusters = 10)
table(mcmc_ci, sim_list$ci)

## Beta
plot_list <- list()

for(k in 1:5){
  ci_index <- paste0("cluster: ", k)
  b0_mcmc <- t(test_result$beta[k, , 5001:10000])
  sum_x0 <- apply(exp(b0_mcmc), 1, sum)
  norm_x0 <- data.frame(iter = 1:5000, exp(b0_mcmc)/sum_x0)
  melt_norm_x0 <- melt(norm_x0 ,  id.vars = 'iter', variable.name = 'Variable')
  plot_list[[k]] <- ggplot(melt_norm_x0, aes(x = iter, y = value)) + 
    geom_line(aes(colour = Variable)) +
    theme_bw() +
    labs(x = "Iteration", y = TeX("Normalized \\hat{\\xi}"),
         title = ci_index) +
    ylim(0, 1) +
    geom_hline(yintercept = ORA[k, ], linetype = "dashed", color = "grey") +
    theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          legend.position = "hide", plot.title = element_text(size = 20))
}

grid.arrange(grobs = plot_list)

### Save data and result
result_list <- list(data = sim_list, result = test_result)
save(result_list, file = "/Users/kevin-imac/Desktop/result_c2.Rdata")

rm(list = ls())

xtable(ORA, digits = c(1, rep(5, 10)))
xtable(table(mcmc_ci, sim_list$ci))

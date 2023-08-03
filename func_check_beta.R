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

int k, arma::mat z, arma::uvec clus_assign, 
arma::mat gamma_mat, arma::vec w, arma::vec beta_k, 
double s2

### Sampler
store_plot <- list()

set.seed(21408)
for(i in 1:10){
  w_samp <- debug_w(iter = 10000, K = 2, z = sim_list$z, clus_assign = sim_list$ci, 
                    gamma_mat = sim_list$gamma, beta_mat = sim_list$beta,
                    r0w = 1, r1w = 1)
  store_plot[[i]] <- ggplot(data.frame(x = 1:10000, y = apply(w_samp, 1, sum)), 
                            aes(x = x, y = y)) +
    geom_line() +
    theme_bw() +
    geom_vline(xintercept = 5000, color = "red", linetype = "dashed") +
    labs(title = "MCMC: Important variables", x = "Iteration", 
         y = "# Imp Variable")
  apply(w_samp[-(1:5000), ], 2, mean) %>% print()
}

grid.arrange(grobs = store_plot)
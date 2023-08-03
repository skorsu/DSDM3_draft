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

### Plot
dat <- cbind(index = 1:100, sim_list$z)
colnames(dat)[2:6] <- str_pad(1:10, 3, pad = "0")
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
# t2 <- data.frame(dat[which(sim_list$ci == 2), ])
# tt2 <- melt(t2, id.vars = "index")
# tt2$variable <- substr(tt2$variable, 2, 4)
# p2 <- ggplot(tt2, aes(x = variable, y = value)) + 
#   geom_line(aes(color = index, group = index)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
#         axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
#         legend.position = "hide", plot.title = element_text(size = 30)) +
#   labs(x = "Variable", y = "", title = "Simulated Data: Third Cluster")
grid.arrange(p0, p1)

### Check: log_w function
### go through every possible vector w => 2^5 = 32 - 1 vectors. (Eliminate all-0)
w_test <- expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1))
w_test <- w_test[-1, ]

w_test <- arrange(w_test, -Var1, -Var2, -Var3, -Var4, -Var5)

### Hand Calculation
w_hand <- function(z, ci, gm, w_vec, bet, r0, r1){
  J <- length(w_vec)
  t1 <- sum(lbeta(r0 + w_vec, r1 + (1 - w_vec)))
  xi_mat <- exp(bet[ci, which(w_vec == 1)])
  z_mat <- z[, which(w_vec == 1)]
  gm_mat <- gm[, which(w_vec == 1)]
  gx_mat <- gm_mat * xi_mat
  
  z_gx_sum <- z_mat + gx_mat
  gx_sum <- gx_mat
  if(!is.null(dim(z_mat))){
    z_gx_sum <- apply(z_mat + gx_mat, 1, sum)
    gx_sum <- apply(gx_mat, 1, sum)
  }
  
  t1 + sum(lgamma(gx_sum)) - sum(lgamma(z_gx_sum)) + 
    sum(lgamma(z_mat + gx_mat)) - sum(lgamma(gx_mat))

}

w_calc_result <- matrix(NA, ncol = 2, nrow = 31)
for(i in 1:31){
  w_int <- as.numeric(w_test[i, ])
  hand_calc <- w_hand(z = sim_list$z, ci = sim_list$ci + 1, gm = sim_list$gamma,
                      w_vec = w_int, bet = sim_list$beta, r0 = 1, r1 = 1)
  rcpp_calc <- log_w(z = sim_list$z, clus_assign = sim_list$ci, 
                     gamma_mat = sim_list$gamma, w = w_int, beta_mat = sim_list$beta, 
                     r0w = 1, r1w = 1)
  w_calc_result[i, ] <- c(hand_calc, rcpp_calc)
}

data.frame(cbind(w_test, w_calc_result)) %>%
  xtable(digits = c(rep(0, 6), rep(5, 2)))

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
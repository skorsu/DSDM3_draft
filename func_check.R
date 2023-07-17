rm(list = ls())
library(tidyverse)

### Simulate the data
set.seed(31)
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data_sim.R")
dat_test <- data_sim(N = 100, K = 4, J = 100, J_imp = 10, z_lb = 75, z_ub = 100, 
                     z_lb_ova = 50, z_ub_ova = 75, pi_gk = rep(0.75, 4),
                     pi_g_ova = 0.9)
dat_test$z

### Test: w_j
set.seed(3)
w <- rbinom(100, 1, 0.5)
beta_mat <- matrix(rnorm(5 * 100), nrow = 5)
w_mat <- t(matrix(rep(w, 100), nrow = 100))
gwx <- dat_test$gm * w_mat * exp(beta_mat[dat_test$ci + 1, ])

for(j in 1:100){
  w_rcpp <- log_w_j(j - 1, z = dat_test$z, clus_assign = dat_test$ci, 
                    gamma_mat = dat_test$gm, w = w, beta = beta_mat, 
                    b0w = 13, b1w = 20)
  w_base_r <- lbeta(13 + w[j], 20 + 1 - w[j]) + sum(lgamma(apply(gwx, 1, sum))) - 
    sum(lgamma(apply(dat_test$z + gwx, 1, sum))) +
    sum(apply(data.frame(dat_test$z[, j] + gwx[, j]), 1, function(x){ifelse(x == 0, 0, lgamma(x))})) -
    sum(apply(data.frame(gwx[, j]), 1, function(x){ifelse(x == 0, 0, lgamma(x))}))
  
  print(c(w_rcpp, w_base_r))
}

t <- update_w(z = dat_test$z, clus_assign = dat_test$ci, gamma_mat = dat_test$gm,
              w = c(rep(0, 99), 1), beta_mat = beta_mat, b0w = 10, b1w = 2)
table(t, w)

### Test: gamma_ijk
set.seed(3)
w <- rbinom(50, 1, 0.5)
g <- rbinom(50, 1, 0.5)
bet_k <- rnorm(50)

for(i in 1:50){
  ### Rcpp function
  rcpp <- log_g_ijk(j = (i-1), zi = dat_test$z[1, ], gi = g, w = w, beta_k = bet_k,
                    b0g = 6, b1g = 12)
  ### base R
  wgx_i <- w * g * exp(bet_k)
  z_gwx_i <- dat_test$z[1, ] + wgx_i
  base_R <- lbeta(6 + g[i], 12 + (1 - g[i])) - lbeta(6, 12) + lgamma(sum(wgx_i)) - 
    sum(lgamma(wgx_i[wgx_i != 0])) + sum(lgamma(z_gwx_i[z_gwx_i != 0])) - lgamma(sum(z_gwx_i))
  
  print(c(rcpp, base_R))
  
}


wgx_i <- w * g * exp(bet_k)
z_gwx_i <- dat_test$z[1, ] + wgx_i

lbeta(1 + g[1], 2 + (1 - g[1])) - lbeta(1 , 2) + lgamma(sum(wgx_i)) - 
  sum(lgamma(wgx_i[wgx_i != 0])) + sum(lgamma(z_gwx_i[z_gwx_i != 0])) - lgamma(sum(z_gwx_i))

set.seed(12)
w <- rbinom(50, 1, 0.5)
beta_mat <- matrix(rnorm(5 * 50), nrow = 5)
tt <- update_gamma(z = dat_test$z, clus_assign = dat_test$ci, gamma_mat = dat_test$gm, 
                   w = w, beta_mat = beta_mat, b0g = 1, b1g = 1)



which(dat_test$z[99, ] == 0) - 1

dat_test$z[100, ]
dat_test$gm[100, ]
w

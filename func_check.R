### Import the external function
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")

### Debugging: at-risk
#### Generate the data
#### 50 observations.
#### 3 clusters with 10 variables (all variables are active).

set.seed(12)
ci <- sort(rep(1:3, 20))
beta_mat <- matrix(c(1, 2, 3), nrow = 3, ncol = 10)
gamma_mat <- matrix(rbinom(10 * 60, 1, 0.5), nrow = 60, ncol = 10)
alpha_mat <- t(apply(gamma_mat * exp(beta_mat[ci, ]), 1, function(x){x/sum(x)}))
z <- matrix(NA, ncol = 10, nrow = 60)
for(i in 1:60){
  z[i, ] <- rmultinom(1, 10, alpha_mat[i, ])
}

set.seed(12)
start_time <- Sys.time()
test_gamma <- debug_gamma(iter = 10000, K = 3, z = z, clus_assign = ci - 1, 
                          w = rep(1, 10), beta_mat = beta_mat, r0g = 1, r1g = 1)
Sys.time() - start_time

gamma_hat <- matrix(NA, ncol = 10, nrow = 60)
for(i in 1:60){
  gamma_hat[i, ] <- rowMeans(test_gamma[i, , 5001:1000])
}

result_0726 <- list(z = z, gamma_mat = gamma_mat, gamma_hat = gamma_hat)
saveRDS(result_0726,
        file = "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/result_0726.rds")
result_0726$gamma_mat

mean(result_0726$gamma_hat[2, result_0726$gamma_hat[2, ] != 1])
sd(result_0726$gamma_hat[2, result_0726$gamma_hat[2, ] != 1])

mean(result_0726$gamma_hat[3, result_0726$gamma_hat[3, ] != 1])
sd(result_0726$gamma_hat[3, result_0726$gamma_hat[3, ] != 1])


gamma_hat[gamma_hat < 0.5]

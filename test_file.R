library(dirmult)
library(Rcpp)

rm(list = ls())
sourceCpp(file = "~/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")

set.seed(2)
zi_test <- simPop(J = 5, K = 10, n = 10, theta = 0.001)$data
gam_array <- array(matrix(1, nrow = 5, ncol = 10), dim = c(5, 10, 3))
gam_array[, , 1][zi_test == 0] <- rbinom(sum(zi_test == 0), 1, 0.25)
gam_array[, , 2][zi_test == 0] <- rbinom(sum(zi_test == 0), 1, 0.25)
gam_array[, , 3][zi_test == 0] <- rbinom(sum(zi_test == 0), 1, 0.25)
beta_mat <- matrix(rnorm(30), nrow = 3, ncol = 10)
ci_test <- c(2, 2, 3, 3, 2) - 1

s1 <- update_atrisk(ci = ci_test, z = zi_test, gamma_old = gam_array,
                    beta_mat = beta_mat, r0g = 1, r1g = 1)
gam_array[, , 2] == s1[, , 2]

tt <- update_beta(ci = ci_test, z = zi_test, gamma_cube = gam_array, 
                  beta_old = beta_mat, mu_beta = 0, s2_beta = 1,
                  s2_MH = 0.01)

tt
beta_mat

log_betak(k = 0, z = zi_test, gamma_mat = gam_array[, , 1], beta_k = beta_mat[2, ],
          ci = ci_test, mu_beta = 0, s2_beta = 10)

dnorm(beta_mat[2, ], 0, sqrt(10), log = TRUE)

exp(-32.29777)

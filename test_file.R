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

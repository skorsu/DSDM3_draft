library(tidyverse)

### Data Simulation
set.seed(1)
N <- 10
J <- 100
K <- 3
pi_g <- c(0.5, 0.3, 0.1)
pi_w <- 0.75
s2 <- 1

### cluster assignment
ci <- sample(1:K, N, replace = TRUE) 
### gamma for each observation based on its cluster (at-risk indicator)
gm_mat <- t(apply(matrix(pi_g[ci]), 1, rbinom, n = J, size = 1))
### w (important variable)
w <- rbinom(J, 1, pi_w)
### Generate xi = exp(beta) where beta ~ N(0, s2).
xi <- matrix(exp(rnorm(J * K, 0, sqrt(s2))), ncol = J)

### alpha matrix
gm_w <- sweep(gm_mat, MARGIN = 2, w, `*`)
w_xi <- sweep(xi[ci, ], MARGIN = 2, w, `*`)
alpha_mat <- matrix(rgamma(N * J, gm_w * w_xi, 1), ncol = J)

### normalize the alpha matrix
alpha_norm <- t(apply(alpha_mat, 1, function(x){x/sum(x)}))
rmultinom(1, round(runif(1, 500, 1000)), alpha_norm[1, ])



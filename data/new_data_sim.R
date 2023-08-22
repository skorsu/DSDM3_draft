library(tidyverse)

data_sim <- function(n, J, K, peak, prob_1, pi_gm, sum_z){
  
  ci <- sample(1:K, n, replace = TRUE)
  
  imp_base <- expand.grid(replicate(J, list(c(0, 1))))
  imp_clus <- NULL
  if(length(peak) == 1){
    imp_clus <- sample(which(rowSums(imp_base) %in% peak), K, replace = FALSE)
  } else {
    K_not_1 <- max(round((1 - prob_1) * K), 1)
    imp_clus <- c(sample(which(rowSums(imp_base) == peak[1]), K - K_not_1, replace = FALSE),
                  sample(which(rowSums(imp_base) == peak[2]), K_not_1, replace = FALSE))
  }
  
  imp_var_clus <- matrix(unlist(imp_base[imp_clus, ]), ncol = J)
  
  xi <- imp_base[imp_clus, ]
  xi[xi == 0] <- 0.1
  xi <- matrix(unlist(xi), ncol = J)
  
  pi_gm_mat <- matrix(pi_gm, nrow = K, ncol = J)
  
  zi <- matrix(NA, ncol = J, nrow = n)
  gm <- matrix(NA, ncol = J, nrow = n)
  
  for(i in 1:n){
    gm[i, ] <- as.numeric(rbinom(J, 1, pi_gm_mat[ci[i], ]) | imp_var_clus[ci[i], ])
    zi[i, ] <- rmultinom(1, sum_z, xi[ci[i], ] * gm[i, ])
  }
  
  list(ci = ci - 1, xi = xi, z = zi, gm = gm)
  
}

set.seed(1)
t <- data_sim(n = 100, J = 10, K = 5, peak = c(1, 2), prob_1 = 0.6, 
              sum_z = 2500, pi_gm = c(0.75, 0.85, 0.95, 1, 0.8))

data.frame(v = paste0("V", str_pad(1:10, 2, pad = "0")),
           m = colMeans(t$z[which(t$ci == 4), ])) %>%
  ggplot(aes(x = v, y = m)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylim(0, 2000)

matrix(unlist(t$xi), ncol = 10)
t$xi

as.numeric(t$xi)

t$xi[4, ]/1.09
t$xi[5, ]/2.08
rowSums(t$xi)

sapply(t$xi)

tt <- matrix(c(0.5, 0.5, 0.75, 0.75, 1), ncol = 10, nrow = 5)
xx <- rbinom(10, 1, tt[1, ])
xx
as.numeric(t$imp_clus[1, ] | xx)

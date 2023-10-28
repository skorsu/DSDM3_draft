library(Rcpp)
library(dirmult)
library(salso)
library(foreach)
library(doParallel)
library(mclustcomp)
library(cluster)
library(ecodist)
library(ggplot2)
library(plotrix)

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
# sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

### Simulate the data (Easiest Pattern)
data_sim <- function(n, pat_mat, pi_gamma, xi_conc, xi_non_conc, sum_z){
  K <- nrow(pat_mat)
  ci <- sample(1:K, n, replace = TRUE)
  
  pmat <- (pat_mat[ci, ] * xi_conc) + xi_non_conc
  at_risk_mat <- matrix(rbinom(n * ncol(pat_mat), 1, pi_gamma),
                        nrow = n, ncol = ncol(pat_mat))
  prob_mat <- matrix(NA, nrow = n, ncol = ncol(pat_mat))
  z <- matrix(NA, nrow = n, ncol = ncol(pat_mat))
  for(i in 1:n){
    prob_mat[i, ] <- as.numeric(at_risk_mat[i, ]) * pmat[i, ]
    z[i, ] <- simPop(J = 1, K = ncol(pat_mat), n = sum_z, 
                     pi = prob_mat[i, ]/sum(prob_mat[i, ]),
                     theta = 0.001)$data
  }

  return(list(ci = ci, at_risk_mat = at_risk_mat, prob_mat = prob_mat,
              z = z))
  
}

### Simulate the data
registerDoParallel(5)
data_sim <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  data_sim(n = 50, pat_mat = (diag(20)[1:3, ]), pi_gamma = 1,
           xi_conc = 10, xi_non_conc = 0.01, sum_z = 2500)
}
stopImplicitCluster()

### Try pam with BrayCurtis
pamBC <- matrix(NA, nrow = 50, ncol = 15)
pamBC_adjR <- rep(NA, 15)
for(i in 1:15){
  pam_bc <- lapply(2:10, 
                   function(x, dat_mat){pam(bcdist(dat_mat), x)$silinfo$avg.width}, 
                   dat_mat = data_sim[[i]]$z) |>
    unlist() |>
    which.max() + 1
  pamBC[, i] <- pam(bcdist(data_sim[[i]]$z), pam_bc)$clustering
  pamBC_adjR[i] <- mclustcomp(pamBC[, i], data_sim[[i]]$ci)[1, 2]
}

### Our Model
registerDoParallel(5)
result_clus <- foreach(t = 1:5) %dopar% {
  set.seed(3)
  start_time <- Sys.time()
  test_result <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[t]]$z,
                           theta_vec = rep(1, 10), launch_iter = 10,
                           MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                           r0c = 1, r1c = 1, print_iter = 2000)
  time_diff <- difftime(Sys.time(), start_time)
  return(list(result = test_result$assign, time_diff = time_diff))
}
stopImplicitCluster()

### Comment

mod_result <- ZIDM_ZIDM(iter = 10000, K_max = 10, z = data_sim[[1]]$z,
                        theta_vec = rep(1, 10), launch_iter = 10,
                        MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                        r0c = 1, r1c = 1, print_iter = 2000)
table(salso(mod_result$assign[-(1:5000), ]), data_sim[[1]]$ci)

adj_rand <- rep(NA, 5)
for(i in 1:5){
  adj_rand[i] <- mclustcomp(as.numeric(salso(result_clus[[i]]$result[-(1:5000), ])), 
                            data_sim[[i]]$ci)[1, 2]
}
mean(adj_rand)

nactive <- lapply(1:15,
       function(x){apply(result_clus[[x]]$result, 1, function(y){length(unique(y))})}) |>
  unlist() |>
  matrix(ncol = 15, byrow = FALSE)

par(mfrow = c(2, 2))
matplot(nactive[, which(adj_rand == 1)], type = "l", lty = 1, 
        ylim = c(1, 10), main = "Perfect Case", ylab = "# Active Cluster")
matplot(nactive[, which(adj_rand == 0)], type = "l", lty = 1,
        ylim = c(1, 10), main = "Worst Case", ylab = "# Active Cluster")
matplot(nactive[, which(0 < adj_rand & adj_rand < 1)], type = "l", lty = 1,
        ylim = c(1, 10), main = "Not-Good Case", ylab = "# Active Cluster")

### Update only beta
for(t in 1:15){
  set.seed(t)
  tt <- beta_mat_update(K = 3, iter = 10000, z = data_sim[[t]]$z, 
                        clus_assign = data_sim[[t]]$ci - 1, 
                        gm = data_sim[[t]]$at_risk_mat, 
                        mu = 0, s2 = 1, s2_MH = 1)
  
  path_pic <- "/Users/kevin-imac/Desktop/pic/"
  jpeg(paste0(path_pic, 'beta_rep', t, '.jpg'), width = 1532, height = 931,
       quality = 90)
  par(mfrow = c(1, 3))
  for(i in 1:3){
    plot_main <- paste0("Rep: ", t, " - Cluster: ", i)
    est_relative <- exp(t(tt[i, , ]))/rowSums(exp(t(tt[i, , ])))
    matplot(est_relative, type = "l", lty = 1, lwd = 1.25, ylab = "Relative Count", 
            main = plot_main)
    abline(h = colMeans(data_sim[[t]]$z[which(data_sim[[t]]$ci == i), ]/2500),
           lty = "dashed")
  }
  dev.off()
  
}

### Update only at-risk indicator
atrisk_q <- matrix(NA, nrow = 15, ncol = 3)
for(i in 1:15){
  app_beta <- rbind(colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == 1), ]/2500),
                    colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == 2), ]/2500),
                    colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == 3), ]/2500))
  tt <- atrisk_update(K = 3, iter = 10000, z = data_sim[[i]]$z, 
                      clus_assign = data_sim[[i]]$ci - 1, beta_mat = app_beta, 
                      r0g = 1, r1g = 1)
  zero_loc <- which(data_sim[[i]]$z == 0)
  atrisk_zero <- matrix(NA, nrow = 5000, ncol = length(zero_loc))
  for(t in 1:5000){
    atrisk_zero[t, ] <- as.numeric(tt[, , t + 5000])[zero_loc]
  }
  
  atrisk_q[i, ] <- quantile(apply(atrisk_zero, 2, mean), c(0.5, 0.025, 0.975))
}

colnames(atrisk_q) <- c("y", "low", "up")
atrisk_q <- data.frame(r = 1:15, atrisk_q)
plotCI(x = atrisk_q$r, y = atrisk_q$y, li = atrisk_q$low, ui = atrisk_q$up,
       ylab = "Rate", xlab = "Replicated Data Set")

### Update the cluster assignment

### Update both at-risk and beta
ab_list <- vector(mode = "list", length = 15)
for(i in 1:15){
  ab_list[[i]] <- beta_ar_update(K = 3, iter = 10000, z = data_sim[[i]]$z, 
                                 clus_assign = data_sim[[i]]$ci - 1, r0g = 1, r1g = 1, 
                                 mu = 0, s2 = 1, s2_MH = 1)
}


ttt <- beta_ar_update(K = 3, iter = 10000, z = data_sim[[1]]$z, 
               clus_assign = data_sim[[1]]$ci - 1, r0g = 1, r1g = 1, 
               mu = 0, s2 = 1, s2_MH = 1)
exp_r <- exp(t(ttt$beta[1, ,]))/rowSums(exp(t(ttt$beta[1, ,])))
matplot(exp_r, type = "l", lty = 1, lwd = 1.25)
abline(h = colMeans(data_sim[[1]]$z[which(data_sim[[1]]$ci == 1), ]/2500),
       lty = "dashed")


i <- 11
par(mfrow = c(1, 3))
for(k in 1:3){
  plot_main <- paste0("Rep: ", i, " | Cluster: ", k)
  rel_exp <- exp(t(ab_list[[i]]$beta[k, , ]))/rowSums(exp(t(ab_list[[i]]$beta[k, , ])))
  matplot(rel_exp, type = "l", lty = 1, lwd = 1.25, main = plot_main)
  abline(h = colMeans(data_sim[[i]]$z[which(data_sim[[i]]$ci == k), ]/2500),
         lty = "dashed")
}





### Our Model (with 20,000 iterations)
registerDoParallel(5)
result_clus1 <- foreach(t = 1:15) %dopar% {
  set.seed(t)
  start_time <- Sys.time()
  test_result <- ZIDM_ZIDM(iter = 20000, K_max = 10, z = data_sim[[t]]$z,
                           theta_vec = rep(1, 10), launch_iter = 5,
                           MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                           r0c = 1, r1c = 1, print_iter = 2000)
  time_diff <- difftime(Sys.time(), start_time)
  return(list(result = test_result$assign, time_diff = time_diff))
}
stopImplicitCluster()

adj_rand1 <- rep(NA, 15)
for(i in 1:15){
  adj_rand1[i] <- mclustcomp(as.numeric(salso(result_clus1[[i]]$result[-(1:5000), ])), 
                             data_sim[[i]]$ci)[1, 2]
}

      
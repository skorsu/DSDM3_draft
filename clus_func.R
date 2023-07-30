rm(list = ls())
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")

### Install packages from the Yushu shi
# install.packages("/Users/kevin-imac/Desktop/sparseMbClust_1.0.tar.gz", repos = NULL, type="source")
library(sparseMbClust)
library(salso)
library(reshape)
library(ggplot2)

data("children")
dim(children)

median(apply(children, 2, sum))

children <- as.matrix(children)
## PerformClustering(children, ClusteringMethod = "DP", totaliter = 2000, burnin = 1000, thin = 5)

### Try with the simulated dataset instead
set.seed(21)
dat_test <- data_sim(N = 100, K = 2, J = 10, J_imp = 2, pi_gk = rep(0.25, 5), 
                     pi_g_ova = 0.75, xi_conc = 5, U_imp = 100, U_unimp = 100,
                     shuffle = FALSE)
### Simulated Data
table(dat_test$ci)
dat_test$ci[99:100]
dat_test$z[99:100, ]
log(dat_test$xi)


### Our model
#### First Try
set.seed(135)
start_time <- Sys.time()
test_result <- clusterZI(K_max = 10, iter = 5000,
                         z = dat_test$z, theta = rep(1, 10), 
                         b0g = 1, b1g = 1, b0w = 1, b1w = 1, MH_var = 0.01, 
                         s2 = 1, launch_iter = 10, b0c = 1, b1c = 1, 
                         print_iter = 250)
Sys.time() - start_time
plot(apply(test_result$w, 1, sum), type = "l") ### Ends up at 6 important variables
first_imp_w <- which(apply(test_result$w[4000:5000, ], 2, mean) == 1)
first_imp_w
first_ci <- salso(test_result$assign[-c(1:3000), ])
table(model = first_ci, actual = dat_test$ci)

#### Second Try
set.seed(7941)
start_time <- Sys.time()
test_result <- clusterZI(K_max = 10, iter = 5000,
                         z = dat_test$z, theta = rep(1, 10), 
                         b0g = 1, b1g = 1, b0w = 1, b1w = 1, MH_var = 0.01, 
                         s2 = 1, launch_iter = 10, b0c = 1, b1c = 1, 
                         print_iter = 250)
Sys.time() - start_time
plot(apply(test_result$w, 1, sum), type = "l") ### Ends up at 6 important variables
second_imp_w <- which(apply(test_result$w[4000:5000, ], 2, mean) == 1)
cbind(first_imp_w, second_imp_w)
second_ci <- salso(test_result$assign[-c(1:3000), ])
table(model = second_ci, actual = dat_test$ci)

plot(apply(test_result$assign, 1, function(x){length(unique(x))}), type = "l")



### Yushu Shi
dat_z <- t(dat_test$z)
rownames(dat_z) <- paste0("var ", 1:100)
colnames(dat_z) <- paste0("obs", 1:100)
set.seed(1)
start_time <- Sys.time()
test_yss <- PerformClustering(dat_z, ClusteringMethod = "MFM", totaliter = 2000, burnin = 500, thin = 1)
Sys.time() - start_time
clus_yss <- as.numeric(salso(test_yss$crec, maxNClusters = 2))
table(clus_yss, dat_test$ci)

### Generate the data (Based on Yushu Shi)
set.seed(97)
data("children")
child_dat <- t(children)
otu_group <- sample(0:2, dim(child_dat)[2], replace = TRUE) ### 1 is Phi, 2 is lambda

table(otu_group)

### Marginal probability
mar_prob <- colSums(child_dat)/sum(colSums(child_dat))
mar_phi <- mar_prob[which(otu_group == 1)]
mar_lambda <- mar_prob[which(otu_group == 2)]

summary(mar_phi)
summary(mar_lambda)

z <- 1

### A
sum((1-(z/5)) * mar_phi)
sum(((sum(mar_lambda) + ((z/5) * sum(mar_phi)))/sum(mar_lambda)) * mar_lambda)





source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")
set.seed(21)
dat_test <- data_sim(N = 100, K = 4, J = 200, J_imp = 20, pi_gk = rep(0.25, 5), 
                     pi_g_ova = 0.75, xi_conc = 5, U_imp = 100, U_unimp = 100,
                     shuffle = FALSE)
start_time <- Sys.time()
test_result <- clusterZI(K_max = 10, iter = 5000,
                         z = dat_test$z, clus_assign = dat_test$ci, 
                         b0g = 1, b1g = 1, b0w = 1, b1w = 1, MH_var = 0.01, s2 = 1)
abs(start_time - Sys.time())


plot(apply(test_result$w, 1, sum), type = "l")
apply(test_result$w, 1, sum)[4990:5000]

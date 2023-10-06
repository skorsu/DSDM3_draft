rm(list = ls())

### Required Libraries: --------------------------------------------------------
library(foreach)
library(doParallel)
library(tidyverse)
library(LaplacesDemon)
library(ClusterZI)
library(cluster)
library(ecodist)

### Setting: -------------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/" ### path for saving the data and result
case_name <- "z_2_pi_90_J_100" ### the path that use for differentiating each case
z_case <- 2
pi_gm_case <- 0.9
J_case <- 100
data_rep <- 4
ncore <- 5

### User-defined functions: ----------------------------------------------------
### Function: Data Simulation
data_sim <- function(N, J, pi_gamma, z, aPhi = 1, bPhi = 9,
                     aLambda = 2, bLambda = 8, aNoise = 1, bNoise = 1){
  
  ### (a) Get the "marginal" probability
  pPhi <- rbeta(J/10, aPhi, bPhi)
  pLambda <- rbeta(J/10, aLambda, bLambda)
  pNoise <- rbeta(0.8 * J, aNoise, bNoise)
  
  ### (b) Shi's adjustment
  pPhi_A <- (1 - (z/5)) * pPhi
  pLambda_A <- ((sum(pLambda) + ((z/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  pPhi_B <- (1 + (z/5)) * pPhi
  pLambda_B <- ((sum(pLambda) - ((z/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  
  ### (c) At-risk indicator
  gamma_mat <- matrix(rbinom(N * J, 1, pi_gamma), nrow = N, ncol = J)
  
  ### (d) Calculate the probability matrix for each observations and normalize
  prob_mat <- rbind(t(matrix(rep(c(pPhi_A, pLambda_A, pNoise), N/2), ncol = N/2)),
                    t(matrix(rep(c(pPhi_B, pLambda_B, pNoise), N/2), ncol = N/2)))
  pg <- prob_mat * gamma_mat
  
  ### (f) Generate the Dirichlet random variables
  dat <- matrix(NA, ncol = J, nrow = N)
  for(i in 1:N){
    prob <- rdirichlet(1, 200 * pg[i, ]/sum(pg[i, ]))
    dat[i, ] <- rmultinom(1, 15000, prob)
  }
  
  dat
  
}

### Function: Summarize the simulated data
summarise_dat <- function(list_simDat){
  
  lapply(list_simDat, function(x){as.data.frame(x) |> 
      mutate(obs = paste0("OB", str_pad(1:50, 3, pad = "0"))) |>
      pivot_longer(cols = -obs) |>
      mutate(taxa_name = paste0("TX", str_pad(str_extract(name, "[:digit:]+$"), 3, pad = "0"))) |>
      ggplot(aes(taxa_name, obs, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low="white", high="palegreen3") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
      labs(x = "Variable (Taxa Count)", y = "Observation", fill = "Count")})
  
}

### Data Simulation: -----------------------------------------------------------
registerDoParallel(ncore)
dat <- foreach(t = 1:data_rep) %dopar% {
  set.seed(t)
  data_sim(N = 50, J = J_case, pi_gamma = pi_gm_case, z = z_case)
}
stopImplicitCluster()

summarise_dat(dat)

### Save the simulated data
saveRDS(dat, paste0(path, "data_", case_name, ".RData"))

### Clustering: ----------------------------------------------------------------

set.seed(1)
test <- ZIDM_ZIDM(iter = 25000, K_max = 10, z = dat[[1]], 
                  theta_vec = rep(1, 10), launch_iter = 5, MH_var = 1, mu = 0,
                  s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, print_iter = 2500)
apply(test$assign, 1, function(x){length(unique(x))}) |> plot(type = "l")
salso(test$assign[-(1:10000), ]) |> table(sort(rep(1:2, 25)))

set.seed(1)
test <- ZIDM_ZIDM(iter = 25000, K_max = 10, z = dat[[1]], 
                  theta_vec = rep(1, 10), launch_iter = 1, MH_var = 1, mu = 0,
                  s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, print_iter = 2500)
apply(test$assign, 1, function(x){length(unique(x))}) |> plot(type = "l")
salso(test$assign[-(1:10000), ]) |> table(sort(rep(1:2, 25)))

set.seed(1)
test <- ZIDM_ZIDM(iter = 25000, K_max = 10, z = dat[[1]], 
                  theta_vec = rep(1, 10), launch_iter = 10, MH_var = 1, mu = 0,
                  s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, print_iter = 2500)
apply(test$assign, 1, function(x){length(unique(x))}) |> plot(type = "l")
salso(test$assign[-(1:10000), ]) |> table(sort(rep(1:2, 25)))

set.seed(1)
test <- ZIDM_ZIDM(iter = 25000, K_max = 10, z = dat[[1]], 
                  theta_vec = rep(1, 10), launch_iter = 10, MH_var = 1, mu = 0,
                  s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 4, print_iter = 2500)
apply(test$assign, 1, function(x){length(unique(x))}) |> plot(type = "l")
salso(test$assign[-(1:10000), ]) |> table(sort(rep(1:2, 25)))

set.seed(1)
test <- ZIDM_ZIDM(iter = 25000, K_max = 10, z = dat[[2]], 
                  theta_vec = rep(1, 10), launch_iter = 1, MH_var = 1, mu = 0,
                  s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, print_iter = 2500)
apply(test$assign, 1, function(x){length(unique(x))}) |> plot(type = "l")
salso(test$assign[-(1:10000), ]) |> table(sort(rep(1:2, 25)))

### Finding: change in a and b does not change the result. 
### Launch Iteration has an effect on the result.

registerDoParallel(ncore)
assign <- foreach(t = 1:data_rep) %dopar% {
  
  set.seed(t)
  runtime <- rep(NA, 8)
  
  ### ZIDM-ZIDM
  start_time <- Sys.time()
  result_ZZ <- ZIDM_ZIDM(iter = 25000, K_max = 10, z = dat[[t]], 
                      theta_vec = rep(1, 10), launch_iter = 5, MH_var = 1, mu = 0,
                      s2 = 1, r0g = 1, r1g = 1, r0c = 1, r1c = 1, print_iter = 2500)
  runtime[1] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### DM-ZIDM
  start_time <- Sys.time()
  result_DZ <- DM_ZIDM(iter = 25000, K_max = 10, z = dat[[t]], 
                       theta_vec = rep(1, 10), launch_iter = 5, MH_var = 1, mu = 0,
                       s2 = 1, r0c = 1, r1c = 1, print_iter = 2500)
  runtime[2] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### DM-DM
  start_time <- Sys.time()
  result_DD <- DM_DM(iter = 25000, K_max = 2, z = dat[[t]], 
                     theta_vec = rep(1, 10), MH_var = 1, mu = 0,
                     s2 = 1, print_iter = 2500)
  runtime[3] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### DM-sDM
  start_time <- Sys.time()
  result_DsD <- DM_DM(iter = 25000, K_max = 10, z = dat[[t]], 
                      theta_vec = rep(0.001, 10), MH_var = 1, mu = 0,
                      s2 = 1, print_iter = 2500)
  runtime[4] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### pam (Bray-Curtis)
  start_time <- Sys.time()
  pam_bc <- lapply(2:10, 
                   function(x, dat_mat){pam(bcdist(dat_mat), x)$silinfo$avg.width}, 
                   dat_mat = dat[[t]]) |>
    unlist() |>
    which.max() + 1
  result_pam_BC <- pam(bcdist(dat[[t]]), pam_bc)$clustering
  runtime[5] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### pam (Aitchison)
  start_time <- Sys.time()
  pam_ait <- lapply(2:10, 
                    function(x, dat_mat){pam(adj_aitchison(dat_mat), x)$silinfo$avg.width}, 
                    dat_mat = dat[[t]]) |>
    unlist() |>
    which.max() + 1
  result_pam_AT <- pam(adj_aitchison(dat[[t]]), pam_ait)$clustering
  runtime[6] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### hclust (Bray-Curtis)
  start_time <- Sys.time()
  hclust_bc <- lapply(2:10, 
                      function(x, dat_mat){hclust(bcdist(dat_mat), method = "complete") |> 
                          cutree(x) |>
                          silhouette(bcdist(dat_mat)) |>
                          (\(.) mean(.[, 3]))()}, 
                      dat_mat = dat[[t]]) |>
    unlist() |>
    which.max() + 1
  result_hclust_BC <- hclust(bcdist(dat[[t]]), method = "complete") |> 
    cutree(hclust_bc)
  runtime[7] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  ### hclust (Aitchison)
  start_time <- Sys.time()
  hclust_ait <- lapply(2:10, 
                      function(x, dat_mat){hclust(adj_aitchison(dat_mat), method = "complete") |> 
                          cutree(x) |>
                          silhouette(adj_aitchison(dat_mat)) |>
                          (\(.) mean(.[, 3]))()}, 
                      dat_mat = dat[[t]]) |>
    unlist() |>
    which.max() + 1
  result_hclust_AT <- hclust(adj_aitchison(dat[[t]]), method = "complete") |> 
    cutree(hclust_ait)
  runtime[8] <- difftime(Sys.time(), start_time, units = "mins") |> as.numeric()
  
  list(result_ZZ = result_ZZ$assign, result_DZ = result_DZ$assign, 
       result_DD = result_DD, result_DsD = result_DsD, result_pam_BC = result_pam_BC,
       result_pam_AT = result_pam_AT, result_hclust_BC = result_hclust_BC, 
       result_hclust_AT = result_hclust_AT, runtime = runtime)
  
  }

stopImplicitCluster()

### Save the simulated data
saveRDS(assign, paste0(path, "result_", case_name, ".RData")) 

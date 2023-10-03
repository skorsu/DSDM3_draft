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
case_name <- "z_1_pi_90_J_100" ### the path that use for differentiating each case
z_case <- 1
pi_gm_case <- 0.9
J_case <- 100
data_rep <- 2
ncore <- 5

### User-defined functions: ----------------------------------------------------
### Function: Data Simulation
data_sim <- function(N, J, pi_gamma, z, aPhi = 1, bPhi = 9,
                     aLambda = 2, bLambda = 8){
  
  ### (a) Get the "marginal" probability
  pPhi <- rbeta(J/2, aPhi, bPhi)
  pLambda <- rbeta(J/2, aLambda, bLambda)
  
  ### (b) Shi's adjustment
  pPhi_A <- (1 - (z/5)) * pPhi
  pLambda_A <- ((sum(pLambda) + ((z/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  pPhi_B <- (1 + (z/5)) * pPhi
  pLambda_B <- ((sum(pLambda) - ((z/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  
  ### (c) At-risk indicator
  gamma_mat <- matrix(rbinom(N * J, 1, pi_gamma), nrow = N, ncol = J)
  
  ### (d) Calculate the probability matrix for each observations and normalize
  prob_mat <- rbind(t(matrix(rep(c(pPhi_A, pLambda_A), N/2), ncol = N/2)),
                    t(matrix(rep(c(pPhi_B, pLambda_B), N/2), ncol = N/2)))
  pg <- prob_mat * gamma_mat
  
  ### (f) Generate the Dirichlet random variables
  dat <- matrix(NA, ncol = J, nrow = N)
  for(i in 1:N){
    prob <- rdirichlet(1, 200 * pg[i, ]/sum(pg[i, ]))
    dat[i, ] <- rmultinom(1, 15000, prob)
  }
  
  dat
  
}

### Function: adjusted Aitchison distance
adj_aitchison <- function(simDat){
  simDat[simDat == 0] <- 1e-10 ### Replace 0 with the small number
  geomMean <- exp(rowMeans(log(simDat)))
  dist(log(simDat/geomMean), method = "euclidean")
}

### Data Simulation: -----------------------------------------------------------
registerDoParallel(ncore)
dat <- foreach(t = 1:data_rep) %dopar% {
  set.seed(t)
  data_sim(N = 50, J = J_case, pi_gamma = pi_gm_case, z = z_case)
}
stopImplicitCluster()

### Save the simulated data
saveRDS(dat, paste0(path, "data_", case_name, ".RData")) 

### Clustering: ----------------------------------------------------------------

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

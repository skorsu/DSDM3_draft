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
library(latex2exp)
library(sparseMbClust)
library(tidyverse)
library(pheatmap)
library(mixtools)
library(coda.base)

### The corresponding report on Overleaf is result_1222.tex

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
# sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

### Data Simulation followed Shi's paper ---------------------------------------
data_sim_shi <- function(N, J, pi_gamma, z_case, aPhi = 1, bPhi = 9,
                         aLambda = 2, bLambda = 8, zsum){
  
  ### (a) Get the "marginal" probability
  pPhi <- rbeta(J/2, aPhi, bPhi)
  pLambda <- rbeta(J/2, aLambda, bLambda)
  
  ### (b) Shi's adjustment
  pPhi_A <- (1 - (z_case/5)) * pPhi
  pLambda_A <- ((sum(pLambda) + ((z_case/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  pPhi_B <- (1 + (z_case/5)) * pPhi
  pLambda_B <- ((sum(pLambda) - ((z_case/5) * sum(pPhi)))/sum(pLambda)) * pLambda
  
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
    dat[i, ] <- rmultinom(1, zsum, prob)
  }
  
  list(dat = dat, gamma_mat = gamma_mat, prob_mat = prob_mat)
  
}

### Simulate the data based on the meeting -------------------------------------
datsim_new <- function(n, Jnoise, Jsignal, pi_gamma,
                       ZSumNoise, ### Noise
                       caseSignal, aPhi, bPhi, aLambda, bLambda, ZSumSignal ### Signal
                       ){
  
  ### Generate the noise
  gammaNoise <- matrix(rbinom(n * Jnoise, 1, pi_gamma), 
                       nrow = n, ncol = Jnoise)
  datNoise <- matrix(0, ncol = Jnoise, nrow = n)
  for(i in 1:n){
    probNoise <- rdirichlet(1, gammaNoise[i, ])
    datNoise[i, ] <- rmultinom(1, ZSumNoise, probNoise)
  }
  
  ### Generate the signal
  ListSignal <- data_sim_shi(n, Jsignal, pi_gamma, caseSignal, 
                             aPhi, bPhi, aLambda, bLambda, ZSumSignal)
  
  ### Create the matrix for storing the final result
  cbind(datNoise, ListSignal$dat)
  
}

### Make it to simulate more than 2 clusters -----------------------------------
### Simulate the data ----------------------------------------------------------
nData <- 20
K <- 5 ### The number that we need to run the datsim_new is ceiling(K/2).
nK <- 20 ### The number of the observations for each cluster
J <- 50 ### Number of the variables (0.8 will be noise and 0.2 will be the informative variables)
impMarker <- matrix(c(41:50, 41:50, 1:10, 1:10, 21:30), ncol = 10, byrow = TRUE) ### Indicate the location of the important variables

registerDoParallel(5)
datSimK <- foreach(t = 1:nData) %dopar% {
  
  ksim <- 1
  simDat_sep <- vector("list", ceiling(K/2))
  
  set.seed(t)
  while(ksim <= ceiling(K/2)){
    
    test_dat <- tryCatch({
      datsim_new(2 * nK, Jnoise = 0.8 * J, Jsignal = 0.2 * J, 
                 pi_gamma = 1, ZSumNoise = 20000, 
                 caseSignal = 3, aPhi = 1, bPhi = 1, 
                 aLambda = 1, bLambda = 1, ZSumSignal = 10000);
    }, error = function(e){"ERROR"})
    
    if(is.matrix(test_dat)){
      simDat_sep[[ksim]] <- test_dat
      ksim <- ksim + 1
    }
    
  }
  
  simDat <- do.call(rbind, simDat_sep)
  simDat <- simDat[1:(K * nK), ]
  actual_clus <- sort(rep(1:K, nK))
  
  simDat_new <- matrix(NA, nrow = K * nK, ncol = 50)
  
  for(i in 1:(K * nK)){
    impVar <- impMarker[actual_clus[i], ]
    simDat_new[i, impVar] <- simDat[i, 41:50] ### Replace the cluster indicating variables
    simDat_new[i, which(!(1:50 %in% impVar))] <- simDat[i, -(41:50)] ### Replace the other variables
  }
  
  ### Random order
  obs_index <- sample(1:(K * nK))
  actual_clus <- actual_clus[obs_index]
  simDat_new <- simDat_new[obs_index, ]
  
  list(dat = simDat_new, clus = actual_clus)
  
}
stopImplicitCluster()

### Save the simulated data
save_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/result_1222/"
# save_path <- "/Users/kevin-imac/Desktop/"
case_name <- "diffindex_3_K_5"
saveRDS(datSimK, paste0(save_path, case_name, "_simDat.RData"))

### Example of the data --------------------------------------------------------
rm(datSimK)
# save_path <- "/Users/kevin-imac/Desktop/result_cluster/"
# case_name <- "diffindex_15"
dat <- readRDS(paste0(save_path, case_name, "_simDat.RData")) ## Data
pheatmap(dat[[1]]$dat[sort(dat[[1]]$clus, index.return = TRUE)$ix, ], 
         display_numbers = F, color = colorRampPalette(c('white','springgreen4'))(100), 
         cluster_rows = F, cluster_cols = F)

#### Function: Calculate mean and SD -------------------------------------------
meanSD <- function(x, dplace = 5){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

### Run all models and save the result -----------------------------------------
### ZIDM-ZIDM
set.seed(3124, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
resultZZ <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- mod(iter = 100000, Kmax = 10, nbeta_split = 5, z = dat[[t]]$dat, 
                     atrisk_init = matrix(1, ncol = 50, nrow = 100), 
                     beta_init = matrix(0, ncol = 50, nrow = 10), 
                     ci_init = rep(0, 100), theta = 1, mu = 0, s2 = 1, s2_MH = 1, 
                     launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                     thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)
saveRDS(resultZZ, paste0(save_path, case_name, "_ZZ.RData"))

sapply(1:nData, function(x){resultZZ[[x]]$time}) %>% meanSD()

sapply(1:nData, 
       function(x){apply(resultZZ[[x]]$result$ci_result, 1, 
                         function(y){length(unique(y))})}) %>%
  matplot(type = "l", ylim = c(1, 10), main = "ZIDM-ZIDM",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")

sapply(1:nData, 
       function(x){apply(resultZZ[[x]]$result$ci_result, 1, 
                         function(y){length(unique(y))})}) %>%
  colMeans() %>% meanSD()

actual_clus <- sapply(1:nData, function(x){dat[[x]]$clus})

viZZ <- sapply(1:nData,
               function(x){as.numeric(salso(resultZZ[[x]]$result$ci_result[-(1:500), ]))})
apply(viZZ, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viZZ[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

bdZZ <- sapply(1:nData,
               function(x){as.numeric(salso(resultZZ[[x]]$result$ci_result[-(1:500), ], loss = binder()))})
apply(bdZZ, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdZZ[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##: ----------------------------------------------------------------------------








### DM-ZIDM
start_ova <- Sys.time()
registerDoParallel(5)
resultDZ <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- DMZIDM(iter = 100000, Kmax = 10, z = dat$dat[[t]], 
                        beta_mat = matrix(0, ncol = 50, nrow = 10), 
                        ci_init = rep(0, 50), 
                        theta = 10, launch_iter = 10, r0c = 1, r1c = 1, 
                        thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDZ, paste0(save_path, case_name, "_DZ.RData"))

# start_ova <- Sys.time()
# registerDoParallel(5)
# resultDZ <- foreach(t = 1:nData) %:%
#   foreach(k = 2:10) %dopar% {
#     start_time <- Sys.time()
#     clus_result <- DMZIDM(iter = 100000, Kmax = k, z = dat$dat[[t]], 
#                           beta_mat = matrix(0, ncol = 50, nrow = 10), 
#                           ci_init = rep(0, 50), 
#                           theta = 10, launch_iter = 10, r0c = 1, r1c = 1, 
#                           thin = 100)
#     tot_time <- difftime(Sys.time(), start_time, units = "secs")
#     list(time = tot_time, result = clus_result)
#   }
# stopImplicitCluster()
# difftime(Sys.time(), start_ova)
# 
# saveRDS(resultDZ, paste0(save_path, case_name, "_DZ.RData"))

### DM-DM
start_ova <- Sys.time()
registerDoParallel(5)
resultDD <- foreach(t = 1:nData) %:%
  foreach(k = 2:10) %dopar% {
  start_time <- Sys.time()
  clus_result <- DMDM(iter = 100000, Kmax = k, z = dat$dat[[t]], 
                      beta_mat = matrix(0, ncol = 50, nrow = 10), 
                      ci_init = rep(0, 50), 
                      theta = 10, thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDD, paste0(save_path, case_name, "_DD.RData"))

### DM-sDM
start_ova <- Sys.time()
registerDoParallel(5)
resultDsD <- foreach(t = 1:nData) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- DMDM(iter = 100000, Kmax = 10, z = dat$dat[[t]], 
                      beta_mat = matrix(0, ncol = 50, nrow = 10), 
                      ci_init = rep(0, 50), 
                      theta = 1e-10, thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultDsD, paste0(save_path, case_name, "_DsD.RData"))

### Begin the analysis ---------------------------------------------------------

#### (If needed) Load the result -----------------------------------------------
# save_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/result_1224/"
# case_name <- "easiest_case"

save_path <- "/Users/kevin-imac/Desktop/result_cluster/"
case_name <- "more_difficult_case"
nData <- 20

dat <- readRDS(paste0(save_path, case_name, "_simDat.RData")) ## Data

resultZZ <- readRDS(paste0(save_path, case_name, "_ZZ.RData")) ## ZIDM-ZIDM
resultDZ <- readRDS(paste0(save_path, case_name, "_DZ.RData")) ## DM-ZIDM
resultDD <- readRDS(paste0(save_path, case_name, "_DD.RData")) ## DM-DM
resultDsD <- readRDS(paste0(save_path, case_name, "_DsD.RData")) ## DM-sDM
resultDP <- readRDS(paste0(save_path, case_name, "_DP_new.RData")) ## DP
resultMFM <- readRDS(paste0(save_path, case_name, "_MFM_new.RData")) ## MFM


#### Analyze: Computational Time -----------------------------------------------

## resultDZ <- lapply(1:nData, function(x){resultDZ[[x]][[9]]})
test <- list(resultZZ, resultDZ, resultDsD, resultDP, resultMFM)
length(test)
sapply(1:length(test), 
       function(x){sapply(1:nData, function(y){test[[x]][[y]]$time})}) %>%
  `colnames<-` (c("ZIDM-ZIDM", "DM-ZIDM", "DM-sDM", "DP", "MFM")) %>%
  apply(2, meanSD)

### DM-DM
sapply(1:nData, 
       function(x){sapply(1:9, function(y){resultDD[[x]][[y]]$time})}) %>%
  colSums() %>%
  meanSD()

#### Number of active cluster via MCMC -----------------------------------------
##### ZIDM-ZIDM
sapply(1:nData, 
       function(x){apply(resultZZ[[x]]$result$ci_result, 1, function(y){length(unique(y))})}) %>%
  colMeans() %>%
  meanSD()

sapply(1:nData, 
       function(x){apply(resultZZ[[x]]$result$ci_result, 1, function(y){length(unique(y))})}) %>%
  matplot(type = "l", ylim = c(1, 10), main = "ZIDM-ZIDM: Difficult Level = 3",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")

##### DM-ZIDM
sapply(1:nData, 
       function(x){apply(resultDZ[[x]]$result$ci_result, 1, function(y){length(unique(y))})}) %>%
  colMeans() %>%
  meanSD()

sapply(1:nData, 
       function(x){apply(resultDZ[[x]]$result$ci_result, 1, function(y){length(unique(y))})}) %>%
  matplot(type = "l", ylim = c(1, 10), main = "DM-ZIDM: Difficult Level = 3",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")

##### DM-sDM
sapply(1:nData, 
       function(x){apply(resultDsD[[x]]$result, 1, function(y){length(unique(y))})}) %>%
  colMeans() %>%
  meanSD()

sapply(1:nData, 
       function(x){apply(resultDsD[[x]]$result, 1, function(y){length(unique(y))})}) %>%
  matplot(type = "l", ylim = c(1, 10), main = "DM-sDM: Difficult Level = 3",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")

##### DP
sapply(1:nData, 
       function(x){apply(resultDP[[x]]$result, 1, function(y){max(y) - min(y) + 1})}) %>%
  colMeans() %>%
  meanSD()

sapply(1:nData, 
       function(x){apply(resultDP[[x]]$result, 1, function(y){max(y) - min(y) + 1})}) %>%
  matplot(type = "l", ylim = c(1, 50), main = "DP: Difficult Level = 3 (with 100,000 iterations)",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")

##### MFM
sapply(1:nData, 
       function(x){apply(resultMFM[[x]]$result, 1, function(y){max(y) - min(y) + 1})}) %>%
  colMeans() %>%
  meanSD()

sapply(1:nData, 
       function(x){apply(resultMFM[[x]]$result, 1, function(y){max(y) - min(y) + 1})}) %>%
  matplot(type = "l", ylim = c(1, 50), main = "MFM: Difficult Level = 3 (with 100,000 iterations)",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")


#### Loss: VI ------------------------------------------------------------------
actual_clus <- sapply(1:nData, function(x){dat$clus[[x]]})

##### ZIDM-ZIDM
viZZ <- sapply(1:nData,
               function(x){as.numeric(salso(resultZZ[[x]]$result$ci_result[-(1:500), ]))})
apply(viZZ, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viZZ[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

table(apply(viZZ, 2, function(x){length(unique(x))}))
table(viZZ[, 2], actual_clus[, 2])
table(viZZ[, 3], actual_clus[, 3])
table(viZZ[, 10], actual_clus[, 10])
table(viZZ[, 12], actual_clus[, 12])
table(viZZ[, 15], actual_clus[, 15])
table(viZZ[, 20], actual_clus[, 20])
which(apply(viZZ, 2, function(x){length(unique(x))}) != 1)

##### DM-ZIDM
viDZ <- sapply(1:nData,
               function(x){as.numeric(salso(resultDZ[[x]]$result$ci_result[-(1:500), ]))})
apply(viDZ, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viDZ[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### DM-sDM
viDsD <- sapply(1:nData,
                function(x){as.numeric(salso(resultDsD[[x]]$result[-(1:500), ]))})
apply(viDsD, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viDsD[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### DP
viDP <- sapply(1:nData,
               function(x){as.numeric(salso(resultDP[[x]]$result[-(1:500), ]))})
apply(viDP, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viDP[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### MFM
viMFM <- sapply(1:nData,
                function(x){as.numeric(salso(resultMFM[[x]]$result[-(1:500), ]))})
apply(viMFM, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viMFM[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

#### Loss: binder --------------------------------------------------------------
##### ZIDM-ZIDM
bdZZ <- sapply(1:nData,
               function(x){as.numeric(salso(resultZZ[[x]]$result$ci_result[-(1:500), ], loss = binder()))})
apply(bdZZ, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdZZ[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

which(apply(bdZZ, 2, function(x){length(unique(x))}) == 1)

for(i in 1:20){
  if(! (i %in% c(4, 7))){
    print(paste0("i: ", i))
    print(table(bdZZ[, i], actual_clus[, i]))
  }
}

##### DM-ZIDM
bdDZ <- sapply(1:nData,
               function(x){as.numeric(salso(resultDZ[[x]]$result$ci_result[-(1:500), ], loss = binder()))})
apply(bdDZ, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdDZ[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### DM-sDM
bdDsD <- sapply(1:nData,
                function(x){as.numeric(salso(resultDsD[[x]]$result[-(1:500), ], loss = binder()))})
apply(bdDsD, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdDsD[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### DP
bdDP <- sapply(1:nData,
               function(x){as.numeric(salso(resultDP[[x]]$result[-(1:500), ], loss = binder()))})
apply(bdDP, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdDP[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### MFM
bdMFM <- sapply(1:nData,
                function(x){as.numeric(salso(resultMFM[[x]]$result[-(1:500), ], loss = binder()))})
apply(bdMFM, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdMFM[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

#### DM-DM ---------------------------------------------------------------------
##### Loss: VI
viDD <- sapply(1:nData,
               function(x){as.numeric(salso(resultDD[[x]][[1]]$result[-(1:500), ]))})
apply(viDD, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(viDD[, x], actual_clus[, x])[1, 2]}) %>% meanSD()

##### Loss: binder
bdDD <- sapply(1:nData,
               function(x){as.numeric(salso(resultDD[[x]][[1]]$result[-(1:500), ], loss = binder()))})
apply(bdDD, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(bdDD[, x], actual_clus[, x])[1, 2]}) %>% meanSD()


resultDD[[20]][[1]]
ci_result <- salso(resultDD[[1]][[1]]$result)
ci_result
for(i in 2:10){
  print(-2 * sum(logmar(dat$dat[[1]], matrix(1, ncol = 50, nrow = 50), matrix(0, ncol = 50, nrow = 5))[, 1]) + 
          ((i * 50) + i) * log(50))
}

### Non-parametric models ------------------------------------------------------
#### EM: Choose the number of cluster by using BIC

ciresult_EM <- matrix(0, ncol = nData, nrow = 50)

for(t in 1:nData){
  
  ### Perform EM
  seed_try <- 1
  EMresult <- tryCatch({set.seed(seed_try);
    nchoose <- multmixmodel.sel(dat$dat[[t]], 2:10); ### Find the optimal number of cluster
    kopt <- nchoose["BIC", "Winner"]; 
    set.seed(seed_try);
    multmixEM(y = dat$dat[[t]], k = kopt)
  }, error = function(e){"ERROR"})
  
  while(sum(EMresult == "ERROR") == 1){
    seed_try <- seed_try + 1
    EMresult <- tryCatch({set.seed(seed_try);
      nchoose <- multmixmodel.sel(dat$dat[[t]], 2:10); 
      kopt <- nchoose["BIC", "Winner"]; 
      set.seed(seed_try);
      multmixEM(y = dat$dat[[t]], k = kopt)
    }, error = function(e){"ERROR"})
  }
  
  ### Perform clustering assignment
  ciresult_EM[, t] <- apply(EMresult$posterior, 1, function(x){which(rmultinom_1(x) == 1)})
  
}

#### Number of cluster
apply(ciresult_EM, 2, function(x){length(unique(x))}) %>% meanSD()
sapply(1:nData, 
       function(x){mclustcomp(ciresult_EM[, x], dat$clus[[x]])[1, 2]}) %>% meanSD()

#### PAM: Bray-Curtiss
registerDoParallel(5)
resultBC <- foreach(t = 1:nData) %dopar% {
  
  silVec <- rep(NA, 9)
  bcDist <- bcdist(dat$dat[[t]])
  
  for(k in 2:10){
    silCal <- silhouette(pam(bcDist, k), dmatrix = bcDist)
    silVec[(k-1)] <- mean(silCal[, "sil_width"])
  }
  
  list(kopt = which.max(silVec) + 1, 
       result =  pam(bcDist, which.max(silVec) + 1)$clustering)
  
}
stopImplicitCluster()

sapply(1:nData, function(x){resultBC[[x]]$kopt}) %>% meanSD()
clusBC <- sapply(1:nData, function(x){resultBC[[x]]$result})
sapply(1:nData, 
       function(x){mclustcomp(clusBC[, x], dat$clus[[x]])[1, 2]}) %>% meanSD()

#### PAM: Aitchison
registerDoParallel(5)
resultAT <- foreach(t = 1:nData) %dopar% {
  
  silVec <- rep(NA, 9)
  
  datAdj <- dat$dat[[t]]
  datAdj[datAdj == 0] <- 1e-10
  atDist <- dist(datAdj, "aitchison")
  
  for(k in 2:10){
    silCal <- silhouette(pam(atDist, k), dmatrix = atDist)
    silVec[(k-1)] <- mean(silCal[, "sil_width"])
  }
  
  list(kopt = which.max(silVec) + 1, 
       result =  pam(atDist, which.max(silVec) + 1)$clustering)
  
}
stopImplicitCluster()
sapply(1:nData, function(x){resultAT[[x]]$kopt}) %>% meanSD()
clusAT <- sapply(1:nData, function(x){resultAT[[x]]$result})
sapply(1:nData, 
       function(x){mclustcomp(clusAT[, x], dat$clus[[x]])[1, 2]}) %>% meanSD()

#### PAM: Euclidean
registerDoParallel(5)
resultEU <- foreach(t = 1:nData) %dopar% {
  
  silVec <- rep(NA, 9)
  euDist <- dist(dat$dat[[t]], "euclidean")
  
  for(k in 2:10){
    silCal <- silhouette(pam(euDist, k), dmatrix = euDist)
    silVec[(k-1)] <- mean(silCal[, "sil_width"])
  }
  
  list(kopt = which.max(silVec) + 1, 
       result =  pam(euDist, which.max(silVec) + 1)$clustering)
  
}
stopImplicitCluster()

sapply(1:nData, function(x){resultEU[[x]]$kopt}) %>% meanSD()
clusEU <- sapply(1:nData, function(x){resultEU[[x]]$result})
sapply(1:nData, 
       function(x){mclustcomp(clusEU[, x], dat$clus[[x]])[1, 2]}) %>% meanSD()

### ----------------------------------------------------------------------------

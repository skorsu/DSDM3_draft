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

### Simulate the data ----------------------------------------------------------
nData <- 20
datsim <- vector(mode = 'list', length = nData)
clussim <- vector(mode = 'list', length = nData)
index <- 1
set.seed(31807)
while(index <= nData){
  
  ### Change the setting of the data here
  simDat <- tryCatch({
    datsim_new(n = 50, Jnoise = 40, Jsignal = 10, 
               pi_gamma = 1, ZSumNoise = 20000, 
               caseSignal = 3, aPhi = 1, bPhi = 1, 
               aLambda = 1, bLambda = 1, ZSumSignal = 10000);
  }, error = function(e){"ERROR"})
  
  if(is.matrix(simDat)){
    
    ### Add difficulty by shuffle the row
    order <- sample(1:50)
    datsim[[index]] <- simDat[order, ]
    clussim[[index]] <- sort(rep(1:2, 25))[order]
    index <- index + 1
  }
  
}

### Save the simulated data
datlist <- list(dat = datsim, clus = clussim)
save_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/sensitivity/"
# save_path <- "/Users/kevin-imac/Desktop/sensitivity_0119/"
case_name <- "diff_index_3"
saveRDS(datlist, paste0(save_path, case_name, "_simDat.RData"))

### Example of the data --------------------------------------------------------
rm(datsim, clussim, datlist)
dat <- readRDS(paste0(save_path, case_name, "_simDat.RData")) ## Data

pheatmap(dat$dat[[1]][sort(dat$clus[[1]], index.return = TRUE)$ix, ], 
         display_numbers = F, color = colorRampPalette(c('white','springgreen4'))(100), 
         cluster_rows = F, cluster_cols = F)

### Set of the hyperparameter --------------------------------------------------
hyperParam <- list(c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 20, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 100, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 0.1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 10, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 100, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 5, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 20, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 0.1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 10, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 0.1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 10, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 4, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 9, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 4),
                   c(Kmax = 10, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 9),
                   c(Kmax = 5, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 20, nbeta_split = 10, theta = 10, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1))

### Run all models and save the result -----------------------------------------
nHyperSet <- length(hyperParam)

start_ova <- Sys.time()
registerDoParallel(5)
result <- foreach(h = 1:nHyperSet) %:% ### Each set of the hyperparameter
  foreach(r = 1:nData) %dopar% { ### The number of replicated data
    start_time <- Sys.time()
    clus_result <- mod(iter = 100000, Kmax = hyperParam[[h]]["Kmax"], 
                       nbeta_split = hyperParam[[h]]["nbeta_split"], 
                       z = dat$dat[[r]], 
                       atrisk_init = matrix(1, ncol = 50, nrow = 50), 
                       beta_init = matrix(0, ncol = 50, nrow = hyperParam[[h]]["Kmax"]), 
                       ci_init = rep(0, 50), theta = hyperParam[[h]]["theta"], 
                       mu = 0, s2 = hyperParam[[h]]["s2"], 
                       s2_MH = hyperParam[[h]]["s2_MH"], 
                       launch_iter = hyperParam[[h]]["launch_iter"], 
                       r0g = hyperParam[[h]]["r0g"], 
                       r1g = hyperParam[[h]]["r1g"], 
                       r0c = hyperParam[[h]]["r0c"], 
                       r1c = hyperParam[[h]]["r1c"], 
                       thin = 100)
    tot_time <- difftime(Sys.time(), start_time, units = "secs")
    list(time = tot_time, result = clus_result$ci_result)
  }
stopImplicitCluster()
difftime(Sys.time(), start_ova)

### Save the result ------------------------------------------------------------
for(i in 1:nHyperSet){
  
  file_case <- paste(paste0(names(hyperParam[[i]]), "_", hyperParam[[i]]), 
                     collapse = "_")
  save_result <- list(hyper = hyperParam[[i]], 
                      result = result[[i]])
  saveRDS(save_result, 
          paste0(save_path, case_name, "_",file_case, "_MB.RData"))
  
}

### ----------------------------------------------------------------------------

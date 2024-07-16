library(tidyverse)
library(dirmult)
library(foreach)
library(doParallel)
library(salso)

### Directory path: ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/Data/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/Data/"
}

### Data Simulation followed Shi's paper: --------------------------------------
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

### Test: ----------------------------------------------------------------------
##### Case 1 -- n = 50 with 2 clusters
testDat <- datsim_new(n = 50, Jnoise = 200, Jsignal = 50, pi_gamma = 1,
                      ZSumNoise = 100000, caseSignal = 2.5, aPhi = 1, bPhi = 1, aLambda = 1, bLambda = 1, 
                      ZSumSignal = 25000)
randIndex <- sample(1:50)
testDat <- testDat[randIndex, ]
testDatPlot <- data.frame(testDat, Index = 1:50) %>%
  pivot_longer(!Index)

testDatPlot$name <- factor(testDatPlot$name, levels = paste0("X", 1:500))
ggplot(testDatPlot, aes(x = name, y = Index, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

mod <- mod_adaptive(iter = 2500, Kmax = 10, nbeta_split = 25, 
                    z = testDat, atrisk_init = matrix(1, nrow = 50, ncol = 250), 
                    beta_init = matrix(0, nrow = 10, ncol = 250), 
                    ci_init = rep(0, 50), theta = 1, mu = 0, s2 = 0.1, 
                    s2_MH = 1e-3, t_thres = 500, launch_iter = 30, 
                    r0g = 1, r1g = 1, r0c = 1, r1c = 4, thin = 1)

apply(mod$ci_result, 1, function(x){length(unique(x))}) %>% plot(type = "l")
as.numeric(salso(mod$ci_result)) %>% table(sort(rep(1:2, 25))[randIndex])
plot(mod$beta_result[6, 220, ], type = "l")

###### Simulated the data
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(1)
simu_data_case_I <- foreach(t = 1:20) %dopar% {
  testDat <- datsim_new(n = 50, Jnoise = 200, Jsignal = 50, pi_gamma = 1,
                        ZSumNoise = 100000, caseSignal = 2.5, aPhi = 1, bPhi = 1, aLambda = 1, bLambda = 1, 
                        ZSumSignal = 25000)
  randIndex <- sample(1:50)
  testDat <- testDat[randIndex, ]
  clusAssign <- sort(rep(1:2, 25))[randIndex]
  list(dat = testDat, clus = clusAssign)
}
stopImplicitCluster()
saveRDS(simu_data_case_I, file = paste0(path, "/Simulation Study/simu_data_case_I.rds"))

##### Case 2 -- n = 100 with 5 clusters
testDat_1 <- datsim_new(n = 50, Jnoise = 200, Jsignal = 50, pi_gamma = 1,
                      ZSumNoise = 100000, caseSignal = 3, aPhi = 1, bPhi = 1, aLambda = 1, bLambda = 1, 
                      ZSumSignal = 25000)
testDat_2 <- datsim_new(n = 50, Jnoise = 200, Jsignal = 50, pi_gamma = 1,
                        ZSumNoise = 100000, caseSignal = 3, aPhi = 1, bPhi = 1, aLambda = 1, bLambda = 1, 
                        ZSumSignal = 25000)

testDat <- rbind(cbind(testDat_1[1:25, 201:250], testDat_1[1:25, -(201:250)]),
                 cbind(testDat_1[26:50, 1:50], testDat_1[26:50, 201:250], testDat_1[26:50, 51:200]), 
                 cbind(testDat_2[1:25, 1:100], testDat_2[1:25, 201:250], testDat_2[1:25, 101:200]),
                 testDat_2[26:50, ])

testDatPlot <- data.frame(testDat, Index = 1:100) %>%
  pivot_longer(!Index)

testDatPlot$name <- factor(testDatPlot$name, levels = paste0("X", 1:500))
ggplot(testDatPlot, aes(x = name, y = Index, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")

mod <- mod_adaptive(iter = 2500, Kmax = 10, nbeta_split = 50, 
                    z = testDat, atrisk_init = matrix(1, nrow = 100, ncol = 250), 
                    beta_init = matrix(0, nrow = 10, ncol = 250), 
                    ci_init = rep(0, 100), theta = 1, mu = 0, s2 = 0.1, 
                    s2_MH = 1e-3, t_thres = 500, launch_iter = 30, 
                    r0g = 1, r1g = 1, r0c = 1, r1c = 9, thin = 1)

apply(mod$ci_result, 1, function(x){length(unique(x))}) %>% plot(type = "l")
mod$ci_result[500, ]
table(salso(mod$ci_result[-(1:500), ]), sort(rep(1:4, 25)))
plot(mod$beta_result[1, 220, ], type = "l")

###### Simulated the data
set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(1)
simu_data_case_II <- foreach(t = 1:20) %dopar% {
  testDat_1 <- datsim_new(n = 50, Jnoise = 200, Jsignal = 50, pi_gamma = 1,
                          ZSumNoise = 100000, caseSignal = 3, aPhi = 1, bPhi = 1, aLambda = 1, bLambda = 1, 
                          ZSumSignal = 25000)
  testDat_2 <- datsim_new(n = 50, Jnoise = 200, Jsignal = 50, pi_gamma = 1,
                          ZSumNoise = 100000, caseSignal = 3, aPhi = 1, bPhi = 1, aLambda = 1, bLambda = 1, 
                          ZSumSignal = 25000)
  
  testDat <- rbind(cbind(testDat_1[1:25, 201:250], testDat_1[1:25, -(201:250)]),
                   cbind(testDat_1[26:50, 1:50], testDat_1[26:50, 201:250], testDat_1[26:50, 51:200]), 
                   cbind(testDat_2[1:25, 1:100], testDat_2[1:25, 201:250], testDat_2[1:25, 101:200]),
                   testDat_2[26:50, ])
  randIndex <- sample(1:100)
  testDat <- testDat[randIndex, ]
  clusAssign <- sort(rep(1:4, 25))[randIndex]
  list(dat = testDat, clus = clusAssign)
}
stopImplicitCluster()
saveRDS(simu_data_case_II, file = paste0(path, "/Simulation Study/simu_data_case_II.rds"))


### Simulate the data: ---------------------------------------------------------
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

### Save the simulated data: ---------------------------------------------------
datlist <- list(dat = datsim, clus = clussim)
case_name <- "diffindex_3"
saveRDS(datlist, paste0(path, case_name, "_simDat.RData"))

### ----------------------------------------------------------------------------

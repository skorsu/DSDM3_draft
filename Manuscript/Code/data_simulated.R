library(dirmult)
library(foreach)
library(doParallel)

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

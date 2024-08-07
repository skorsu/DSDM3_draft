library(dirmult)

sim_Signal <- function(N, J, pi_gamma, z_case, aPhi, bPhi, aLambda, bLambda, zsum){
  
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
  
  N1 <- round(N/2)
  N2 <- N - N1
  
  prob_mat <- rbind(t(matrix(rep(c(pPhi_A, pLambda_A), N1), ncol = N1)),
                    t(matrix(rep(c(pPhi_B, pLambda_B), N2), ncol = N2)))
  pg <- prob_mat * gamma_mat
  
  ### (f) Generate the Dirichlet random variables
  dat <- matrix(NA, ncol = J, nrow = N)
  for(i in 1:N){
    prob <- rdirichlet(1, 200 * pg[i, ]/sum(pg[i, ]))
    dat[i, ] <- rmultinom(1, zsum, prob)
  }
  
  list(dat = dat, gamma_mat = gamma_mat, prob_mat = prob_mat)
  
}

sim_clusDat <- function(n, Jnoise, Jsignal, pZero, ZSumNoise, ZSumSignal, 
                        seed, shuffle = TRUE, caseSignal = 3, aPhi = 1, bPhi = 1, 
                        aLambda = 1, bLambda = 1){
  
  set.seed(seed)
  
  ### Check the required packages
  if(!require(dirmult)){
    stop("Missing Package: dirmult")
  }
  
  pi_gamma <- 1 - pZero
  
  ### Generate the noise
  gammaNoise <- matrix(rbinom(n * Jnoise, 1, pi_gamma), 
                       nrow = n, ncol = Jnoise)
  datNoise <- matrix(0, ncol = Jnoise, nrow = n)
  for(i in 1:n){
    probNoise <- rdirichlet(1, gammaNoise[i, ])
    datNoise[i, ] <- rmultinom(1, ZSumNoise, probNoise)
  }
  
  ### Generate the signal
  ListSignal <- sim_Signal(n, Jsignal, pi_gamma, caseSignal, 
                           aPhi, bPhi, aLambda, bLambda, ZSumSignal)
  
  ### Cluster Assignment
  simulated_ci <- c(rep(1, round(n/2)), rep(2, n - round(n/2)))
  simulated_dat <- cbind(datNoise, ListSignal$dat)
  
  ### Create the matrix for storing the final result
  if(shuffle == TRUE){
    index <- sample(1:n)
    simulated_ci <- simulated_ci[index]
    simulated_dat <- simulated_dat[index, ]
  }
  
  list(dat = simulated_dat, c = simulated_ci)
  
}

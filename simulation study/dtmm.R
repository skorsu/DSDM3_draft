## library('devtools')
## devtools::install_github('MaStatLab/DTMM')

library(DTMM)
library(phylo)
library(tidyverse)
library(ape)
library(salso)

### Import the data
# save_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/simulation study/result_1222/"
save_path <- "/Users/kevin-imac/Desktop/"
case_name <- "easiest_case"
nData <- 20
dat <- readRDS(paste0(save_path, case_name, "_simDat.RData"))

### Edit the R code for running DTMM
DTMM_edit <- function(Y, tree, tau_vec = 10 ^ seq(-1, 4, 0.5), nu_vec = 1, 
                      theta_vec = seq(0.01, 0.99, 0.08), init_c, beta = "default",
                      mcmc_iter = 2000){
  
  
  if(beta == "default"){
    temp = runif(1)
    alpha = temp/(1 - temp)
  }
  
  tree <- reorder(tree, order = "postorder")
  edge <- apply(tree$edge, 2, rev)
  p = dim(Y)[2]
  gamma_sample = rbinom(p - 1, 1, 0.5)
  
  HDTMcpp(Y, edge, tau_vec, nu_vec, theta_vec,
          init_c, gamma_sample, alpha, mcmc_iter, select = FALSE)
}

### Generate the tree structure for our model
tree_structure <- read.tree(text = paste0("(", paste(paste0("V", 1:50), collapse = ", "), ");"))
summary(tree_structure)

testMod <- DTMM_edit(Y = dat$dat[[2]], tree = tree_structure, 
                     init_c = rep(1, 50), mcmc_iter = 5000)

### Post-Processing
plot(apply(testMod$post_c, 2, function(x){length(unique(x))}), type = "l")
clusVI <- as.numeric(salso(t(testMod$post_c)))
table(clusVI, dat$clus[[2]])

clusBINDER <- as.numeric(salso(t(testMod$post_c)), loss = "binder")
table(clusBINDER, dat$clus[[2]])


??read.tree


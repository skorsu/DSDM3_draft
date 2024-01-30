### Required Library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(salso)
library(gridExtra)
library(stringr)
library(sparseMbClust)
library(cluster)
library(ecodist)
library(factoextra)
library(rbiom)
library(ggcorrplot)

### Import the data
path <- "/Users/kevinkvp/Desktop/"
# path <- "/Users/kevin-imac/Desktop/"
ni <- read.csv(paste0(path, "Nicaragua_12mo_Metadata_csv.csv"))
ml <- read.csv(paste0(path, "Mali_12mo_Metadata_csv.csv"))

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")

### Functions ------------------------------------------------------------------
### recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}

## Data Cleaning ---------------------------------------------------------------
### Create a new variable and join the two datasets
dat <- cbind("country" = c("NI"), ni[, intersect(colnames(ni), colnames(ml))]) %>%
  rbind(cbind("country" = c("ML"), ml[, intersect(colnames(ni), colnames(ml))]))
for(i in c(1, 2, 3, 5)){
  dat[, i] <- factor(dat[, i])
}

### Clean the data
dat <- dat %>% dplyr::select(-Age)
demo_dat <- dat[, 1:4]
taxa_dat <- dat[, -(1:4)]
taxa_dat <- taxa_dat[, colMeans(taxa_dat > 0) >= 0.1]

view(cbind(demo_dat, taxa_dat))

### Visualize the data ---------------------------------------------------------
pheatmap(taxa_dat, 
         display_numbers = FALSE, color = colorRampPalette(c('white','lightblue3'))(100), 
         cluster_rows = F, cluster_cols = F,
         show_rownames = FALSE, show_colnames = FALSE)

### Begin Sensitivity Analysis -------------------------------------------------
hyperParam <- list(c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 5, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 20, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 1, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 10, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 0.1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 10, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 0.1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 10, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 0.1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 10, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 5, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 20, r0g = 1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 0.1, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 10, r1g = 1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 0.1, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 10, r0c = 1, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 4, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 9, r1c = 1),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 4),
                   c(Kmax = 10, nbeta_split = 5, theta = 1, s2 = 1, s2_MH = 1, launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 9))


registerDoParallel(6)
result <- foreach(h = 1:length(hyperParam)) %dopar% {
  set.seed((2*h) + 1)
  start_time <- Sys.time()
  clus_result <- mod(iter = 100000, Kmax = hyperParam[[h]]["Kmax"], 
                     nbeta_split = hyperParam[[h]]["nbeta_split"], 
                     z = as.matrix(taxa_dat), 
                     atrisk_init = matrix(1, ncol = 46, nrow = 91), 
                     beta_init = matrix(0, ncol = 46, nrow = hyperParam[[h]]["Kmax"]), 
                     ci_init = rep(0, 91), theta = hyperParam[[h]]["theta"], 
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

# path <- "/Users/kevinkvp/Desktop/"
result_path <- "Github Repo/ClusterZI/simulation study/sensitivity/annika/"
saveRDS(result, paste0(path, result_path, "result.RData"))

### Post-Processing ------------------------------------------------------------
## plot(apply(clus_result$ci_result, 1, function(x){length(unique(x))}), type = "l")

### Pairwise salso
dat_mat <- matrix(runif(441, -1, 1), ncol = 21)
diag(dat_mat) <- 1
ggcorrplot(dat_mat, lab = TRUE, 
           colors = c("white", "white", "red"))
?ggcorrplot


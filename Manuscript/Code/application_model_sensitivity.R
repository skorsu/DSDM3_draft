### Load libraries -------------------------------------------------------------
library(Rcpp)
library(foreach)
library(doParallel)
library(tidyverse)

### Import the data ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

sourceCpp(paste0(path, "src/clusterZI.cpp"))

datpath <- "/Users/kevin-imac/Desktop/Annika/"
if(! file.exists(datpath)){
  datpath <- "/Users/kevinkvp/Desktop/Annika/Application/"
}

### Data: 6 and 8 Months
ni68 <- read.csv(paste0(datpath, "Data/Nicaragua_6mo_8mo_genus.csv"))
ml68 <- read.csv(paste0(datpath, "Data/Mali_6mo_8mo_genus.csv"))

### Data: 12 Months
ni12 <- read.csv(paste0(datpath, "Data/Nicaragua_12mo_Metadata_csv.csv"))
ml12 <- read.csv(paste0(datpath, "Data/Mali_12mo_Metadata_csv.csv"))

## Data Pre-processing ---------------------------------------------------------
### For each nationality, first split 6 and 8 from x68. Then, choose only the 
### common infants among three datasets.
ni06 <- ni68 %>% filter(Age..months. == 6)
ni08 <- ni68 %>% filter(Age..months. == 8)

niInfant <- intersect(intersect(str_extract(ni06$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"),
                                str_extract(ni08$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}")),
                      str_extract(ni12$ID, "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"))

ni06 <- ni06[which(str_extract(ni06$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %in% niInfant),]
ni06 <- ni06 %>% rename(ID = ID.)
ni08 <- ni08[which(str_extract(ni08$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %in% niInfant),]
ni08 <- ni08 %>% rename(ID = ID.)
ni12 <- ni12[which(str_extract(ni12$ID, "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %in% niInfant),]

ml06 <- ml68 %>% filter(Age..months. == 6)
ml08 <- ml68 %>% filter(Age..months. == 8)

mlInfant <- intersect(intersect(str_extract(ml06$X.SampleID, "^[:digit:]+\\."), 
                                str_extract(ml08$X.SampleID, "^[:digit:]+\\.")), 
                      str_extract(ml12$ID, "^[:digit:]+\\."))

ml06 <- ml06[which(str_extract(ml06$X.SampleID, "^[:digit:]+\\.") %in% mlInfant), ]
ml06 <- ml06 %>% rename(ID = X.SampleID)
ml08 <- ml08[which(str_extract(ml08$X.SampleID, "^[:digit:]+\\.") %in% mlInfant), ]
ml08 <- ml08 %>% rename(ID = X.SampleID)
ml12 <- ml12[which(str_extract(ml12$ID, "^[:digit:]+\\.") %in% mlInfant), ]

### For the columns, we will first join the Nicaraguan and Malian with the same timestamp
dat06 <- cbind("country" = c("NI"), ni06[, c(intersect(colnames(ni06), colnames(ml06)))]) %>%
  rbind(cbind("country" = c("ML"), ml06[, c(intersect(colnames(ni06), colnames(ml06)))]))
dat08 <- cbind("country" = c("NI"), ni08[, c(intersect(colnames(ni08), colnames(ml08)))]) %>%
  rbind(cbind("country" = c("ML"), ml08[, c(intersect(colnames(ni08), colnames(ml08)))]))
dat12 <- cbind("country" = c("NI"), ni12[, c(intersect(colnames(ni12), colnames(ml12)))]) %>%
  rbind(cbind("country" = c("ML"), ml12[, c(intersect(colnames(ni12), colnames(ml12)))]))

### For the taxa, select only the column with the non-negative count proportion is greater than or equal 0.1
dat06 <- dat06[, colMeans(dat06 > 0) >= 0.1]
dat08 <- dat08[, colMeans(dat08 > 0) >= 0.1]
dat12 <- dat12[, colMeans(dat12 > 0) >= 0.1]

commomTaxa <- intersect(intersect(colnames(dat06[, -(1:5)]), colnames(dat08[, -(1:5)])), 
                        colnames(dat12[, -(1:5)]))

dat06 <- dat06[, c(1:5, which(colnames(dat06[, -(1:5)]) %in% commomTaxa) + 5)]
dat08 <- dat08[, c(1:5, which(colnames(dat08[, -(1:5)]) %in% commomTaxa) + 5)]
dat12 <- dat12[, c(1:5, which(colnames(dat12[, -(1:5)]) %in% commomTaxa) + 5)]

identical(colnames(dat06), colnames(dat08))
identical(colnames(dat12[, -(1:5)]), colnames(dat08[, -(1:5)]))

### Run the models for sensitivity analysis ------------------------------------
#### 6-month
nRepli <- 20
nameResult <- "microbiome_result_at_risk_sensitivity_6mo.RData"

set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
resultZZ <- foreach(t = 1:nRepli) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- mod(iter = 100000, Kmax = 10, nbeta_split = 5, z = as.matrix(dat06[, -(1:5)]), 
                     atrisk_init = matrix(1, ncol = 38, nrow = 90), 
                     beta_init = matrix(0, ncol = 38, nrow = 10), 
                     ci_init = rep(0, 90), theta = 1, mu = 0, s2 = 1, s2_MH = 1, 
                     launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                     thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultZZ, paste0(path, "Manuscript/Result/", nameResult))

#### 8-month
nameResult <- "microbiome_result_at_risk_sensitivity_8mo.RData"

set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
resultZZ <- foreach(t = 1:nRepli) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- mod(iter = 100000, Kmax = 10, nbeta_split = 5, z = as.matrix(dat08[, -(1:5)]), 
                     atrisk_init = matrix(1, ncol = 38, nrow = 90), 
                     beta_init = matrix(0, ncol = 38, nrow = 10), 
                     ci_init = rep(0, 90), theta = 1, mu = 0, s2 = 1, s2_MH = 1, 
                     launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                     thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultZZ, paste0(path, "Manuscript/Result/", nameResult))

#### 12-month
nameResult <- "microbiome_result_at_risk_sensitivity_12mo.RData"

set.seed(1415, kind = "L'Ecuyer-CMRG")
start_ova <- Sys.time()
registerDoParallel(5)
resultZZ <- foreach(t = 1:nRepli) %dopar% {
  
  start_time <- Sys.time()
  clus_result <- mod(iter = 100000, Kmax = 10, nbeta_split = 5, z = as.matrix(dat12[, -(1:5)]), 
                     atrisk_init = matrix(1, ncol = 38, nrow = 90), 
                     beta_init = matrix(0, ncol = 38, nrow = 10), 
                     ci_init = rep(0, 90), theta = 1, mu = 0, s2 = 1, s2_MH = 1, 
                     launch_iter = 10, r0g = 1, r1g = 1, r0c = 1, r1c = 1, 
                     thin = 100)
  tot_time <- difftime(Sys.time(), start_time, units = "secs")
  list(time = tot_time, result = clus_result)
  
}
stopImplicitCluster()
difftime(Sys.time(), start_ova)

saveRDS(resultZZ, paste0(path, "Manuscript/Result/", nameResult))




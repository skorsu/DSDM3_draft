# Required Libraries
library(tidyverse)
library(foreach)
library(doParallel)
library(salso)
library(ggplot2)
library(mclustcomp)
library(latex2exp)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(ggpubr)

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

# Path
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}
datapath <- paste0(path, "Data/Application Data/")
resultpath <- paste0(path, "Result/zupancic/")
# file.exists(resultpath)

## Metadata
metadata <- read.delim(paste0(datapath, "zupancic/ob_zupancic.metadata.txt"))
View(metadata)

## Data
dat <- read.delim(paste0(datapath, "zupancic/ob_zupancic.otu_table.100.denovo.rdp_assigned"))
taxaName <- dat$X
dat <- dat %>% dplyr::select(-X)
dat <- dat[, which(colnames(dat) %in% metadata$X)]
dim(dat)
metadata <- metadata[which(metadata$X %in% colnames(dat)), ]
dim(metadata)

### Combined Duplicated Taxa
taxaSplit <- do.call(rbind.data.frame, 
                     lapply(taxaName, function(x){str_replace(str_split_fixed(x, "\\;", 8), "\\([:digit:]+\\)\\;*", "")}))
View(taxaSplit)
colnames(taxaSplit) <- c("k", "p", "c", "o", "f", "g", "s", "d")
taxaSplit <- taxaSplit[, 1:6]
distinctTaxa <- taxaSplit %>% distinct()

dim(dat)
dim(taxaSplit)
dim(distinctTaxa)

start_clean <- Sys.time()
cleanDat <- matrix(NA, nrow = 220, ncol = 169)
cleanDat_row <- rep(NA, 169)
j <- 1

for(i in 1:220){
  
  colIndex <- which(sapply(1:105998, function(x){sum(taxaSplit[x, ] == distinctTaxa[i, ])}) == 6)
  
  if(length(colIndex) != 0){
    
    if(length(colIndex) == 1){
      cleanDat[j, ] <- as.numeric(dat[colIndex, ])
    } else {
      cleanDat[j, ] <- colSums(dat[colIndex, ])
    }
    
    cleanDat_row[j] <- paste(distinctTaxa[i, ], collapse = ";")
    
    j <- j + 1
    
  }
  
  print(c(i, j - 1))
  
}
difftime(Sys.time(), start_clean)

rownames(cleanDat) <- cleanDat_row
colnames(cleanDat) <- colnames(dat)
cleanDat <- t(cleanDat)
View(cleanDat)
cleanDat <- cleanDat[, -which(colMeans(cleanDat > 0) < 0.1)]

saveRDS(cleanDat, paste0(datapath, "zupancic/clean_zupancic.rds"))

### Read the data
otuTab <- readRDS(paste0(datapath, "zupancic/clean_zupancic.rds"))
otuTab <- as.matrix(otuTab)
dim(otuTab)

# mod <- mod_adaptive(iter = 3000, Kmax = 20, nbeta_split = 10,
#                     z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 169, ncol = 80),
#                     beta_init = matrix(0, nrow = 20, ncol = 80),
#                     ci_init = rep(0, 169),
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-5,
#                     t_thres = 2000, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# 
# apply(mod$ci_result, 1, uniqueClus) %>% plot(type = "l")
# 
# mod$beta_result[, 1, ] %>% t() %>% as.data.frame() %>%
#   mutate(iter = 1:3000) %>%
#   pivot_longer(!iter) %>%
#   ggplot(aes(x = iter, y = value, color = name)) +
#   geom_line()
# 
# sumTab <- data.frame(ID = rownames(otuTab),
#            clus = as.numeric(salso(mod$ci_result[-(1:1500), ]))) %>%
#   inner_join(metadata, by = c("ID" = "X"))
# 
# table(sumTab$clus, sumTab$DiseaseState)

### RUN THE MODEL
set.seed(1)
ciInit <- matrix(0, nrow = 169, ncol = 12)
ciInit[, 4] <- sample(0:2, 169, replace = TRUE)
ciInit[, 5] <- sample(0:2, 169, replace = TRUE)
ciInit[, 6] <- sample(0:2, 169, replace = TRUE)
ciInit[, 7] <- sample(0:4, 169, replace = TRUE)
ciInit[, 8] <- sample(0:4, 169, replace = TRUE)
ciInit[, 9] <- sample(0:4, 169, replace = TRUE)
ciInit[, 10] <- sample(0:19, 169, replace = TRUE)
ciInit[, 11] <- sample(0:19, 169, replace = TRUE)
ciInit[, 12] <- sample(0:19, 169, replace = TRUE)

xiInitDum <- lapply(4:12, function(y){sapply(0:max(ciInit[, y]), function(x){
  pp <- otuTab[which(ciInit[, y] == x), ]
  if(is.matrix(pp)){
    p <- colSums(pp)/sum(pp)
  } else {
    p <- pp/sum(pp)
  }
  # p <- colSums(otuTab[which(ciInit[, y] == x), ])/sum(otuTab[which(ciInit[, y] == x), ])
  ifelse(is.infinite(log(p/(1-p))), -20, log(p/(1-p)))
}) %>% t()
})

is.matrix(otuTab[1, ])

xiInit <- vector("list", 12)
xiInit[[1]] <- matrix(0, nrow = 20, ncol = 80)
xiInit[[2]] <- matrix(0, nrow = 20, ncol = 80)
xiInit[[3]] <- matrix(0, nrow = 20, ncol = 80)

xiInit[[4]] <- rbind(xiInitDum[[1]], matrix(0, nrow = 17, ncol = 80))
xiInit[[5]] <- rbind(xiInitDum[[2]], matrix(0, nrow = 17, ncol = 80))
xiInit[[6]] <- rbind(xiInitDum[[3]], matrix(0, nrow = 17, ncol = 80))

xiInit[[7]] <- rbind(xiInitDum[[4]], matrix(0, nrow = 15, ncol = 303))
xiInit[[8]] <- rbind(xiInitDum[[5]], matrix(0, nrow = 15, ncol = 303))
xiInit[[9]] <- rbind(xiInitDum[[6]], matrix(0, nrow = 15, ncol = 303))

xiInit[[10]] <- xiInitDum[[7]]
xiInit[[11]] <- xiInitDum[[8]]
xiInit[[12]] <- xiInitDum[[9]]

resultName <- c(paste0("result_cleaned_zupancic_chain_", 1:3, "_init_oneClus_s2_1en1_s2MH_1en3.rds"),
                paste0("result_cleaned_zupancic_chain_", 1:3, "_init_3clus_s2_1en1_s2MH_1en3.rds"),
                paste0("result_cleaned_zupancic_chain_", 1:3, "_init_5clus_s2_1en1_s2MH_1en3.rds"),
                paste0("result_cleaned_zupancic_chain_", 1:3, "_init_20clus_s2_1en1_s2MH_1en3.rds"))

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:12) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 50000, Kmax = 20, nbeta_split = 10,
                      z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 169, ncol = 80),
                      beta_init = as.matrix(xiInit[[t]]),
                      ci_init = ciInit[, t],
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
                      t_thres = 10000, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

# mod <- mod_adaptive(iter = 3000, Kmax = 20, nbeta_split = 10,
#                     z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 169, ncol = 80),
#                     beta_init = matrix(0, nrow = 20, ncol = 80),
#                     ci_init = rep(0, 169),
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-5,
#                     t_thres = 2000, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)

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
resultpath <- paste0(path, "Result/baxter/")
# file.exists(resultpath)

## Metadata
metadata <- read.delim(paste0(datapath, "baxter/metadata.tsv"))
taxaList <- read.delim(paste0(datapath, "baxter/glne007.final.an.unique_list.0.03.cons.taxonomy"))
View(metadata)
table(metadata$Site, metadata$dx)
table(metadata$Site)

## Data
dat <- read.delim(paste0(datapath, "baxter/glne007.final.an.unique_list.0.03.subsample.0.03.filter.shared"))
dat <- dat %>% dplyr::select(-label, -X, -numOtus)
datID <- dat$Group
dat <- dat %>% dplyr::select(-Group)
rownames(dat) <- datID
# dat <- dat[, -which(colMeans(dat > 0) < 0.1)]

### Filter only one location
datUMich <- dat[which(rownames(dat) %>% as.numeric() %in% metadata$sample[which(metadata$Site == "U Michigan")]), ]
datUMich <- datUMich[, -which(colMeans(datUMich > 0) < 0.1)]
dim(datUMich)

### Combined Duplicated Taxa
taxaSplit <- do.call(rbind.data.frame, 
                     lapply(taxaList$Taxonomy, function(x){str_replace(str_split_fixed(x, "\\;", 6), "\\([:digit:]+\\)\\;*", "")}))
colnames(taxaSplit) <- c("k", "p", "c", "o", "f", "g")
rownames(taxaSplit) <- taxaList$OTU
distinctTaxa <- taxaSplit %>% distinct()

start_cleaned <- Sys.time()
cleanedUmich <- matrix(NA, nrow = 107, ncol = 309)
cleanedUmichCOL <- rep(NA, 309)
j <- 1

for(i in 1:309){
  
  colIndex <- (colnames(datUMich) %in% taxaList$OTU[which(sapply(1:9467, function(x){sum(taxaSplit[x, ] == distinctTaxa[i, ])}) == 6)]) %>%
    which()
  
  if(length(colIndex) != 0){
    
    if(length(colIndex) == 1){
      cleanedUmich[, j] <- datUMich[, colIndex]
    } else {
      cleanedUmich[, j] <- rowSums(datUMich[, colIndex])
    }
    
    cleanedUmichCOL[j] <- paste(distinctTaxa[i, ], collapse = ";")
    
    j <- j + 1
    
  }
  
  print(c(i, j))
  
}
difftime(Sys.time(), start_cleaned)

cleanedUmich <- cleanedUmich[, 1:81]
rownames(cleanedUmich) <- rownames(datUMich)
colnames(cleanedUmich) <- cleanedUmichCOL[1:81]

View(cleanedUmich)
saveRDS(cleanedUmich, paste0(datapath, "baxter/cleanUmich.rds"))

### Read the data
metadata <- read.delim(paste0(datapath, "baxter/metadata.tsv"))
otuTab <- readRDS(paste0(datapath, "baxter/cleanUmich.rds"))

data.frame(x = 1:81, p = sort(colSums(otuTab)/sum(otuTab), decreasing = TRUE)) %>%
  ggplot(aes(x = x, y = p)) +
  geom_bar(stat = "identity")

x <- rnorm(81, 0, sqrt(1))
data.frame(x = 1:81, p = sort(exp(x)/sum(exp(x)), decreasing = TRUE)) %>%
  ggplot(aes(x = x, y = p)) +
  geom_bar(stat = "identity")

mod <- mod_adaptive(iter = 3000, Kmax = 20, nbeta_split = 10,
                    z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 107, ncol = 81),
                    beta_init = matrix(0, nrow = 20, ncol = 81),
                    ci_init = rep(0, 107),
                    theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
                    t_thres = 2000, launch_iter = 30,
                    r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)

apply(mod$ci_result, 1, uniqueClus) %>% plot(type = "l")

mod$beta_result[, 1, ] %>% t() %>% as.data.frame() %>%
  mutate(iter = 1:3000) %>%
  pivot_longer(!iter) %>%
  ggplot(aes(x = iter, y = value, color = name)) +
  geom_line()

salsoSum <- data.frame(sample = as.numeric(rownames(otuTab)), clus = as.numeric(salso(mod$ci_result[-(1:1000), ]))) %>%
  inner_join(metadata)

table(salsoSum$clus)
table(salsoSum$clus, salsoSum$dx)
table(salsoSum$clus, salsoSum$Hx_Prev)
table(salsoSum$clus, salsoSum$Gender)
table(salsoSum$clus, salsoSum$Smoke)
table(salsoSum$clus, salsoSum$Diabetic)
table(salsoSum$clus, salsoSum$stage)

salsoSum %>%
  group_by(clus) %>%
  summarise(mean(Age), sd(Age), mean(BMI), sd(BMI))

### Run the models
## 1, nrow = 107, ncol = 81

set.seed(1)
ciInit <- matrix(0, nrow = 107, ncol = 12)
ciInit[, 4] <- sample(0:2, 107, replace = TRUE)
ciInit[, 5] <- sample(0:2, 107, replace = TRUE)
ciInit[, 6] <- sample(0:2, 107, replace = TRUE)
ciInit[, 7] <- sample(0:4, 107, replace = TRUE)
ciInit[, 8] <- sample(0:4, 107, replace = TRUE)
ciInit[, 9] <- sample(0:4, 107, replace = TRUE)
ciInit[, 10] <- sample(0:19, 107, replace = TRUE)
ciInit[, 11] <- sample(0:19, 107, replace = TRUE)
ciInit[, 12] <- sample(0:19, 107, replace = TRUE)

xiInitDum <- lapply(4:12, function(y){sapply(0:max(ciInit[, y]), function(x){
  p <- colSums(otuTab[which(ciInit[, y] == x), ])/sum(otuTab[which(ciInit[, y] == x), ])
  ifelse(is.infinite(log(p/(1-p))), -20, log(p/(1-p)))
}) %>% t()
})

xiInit <- vector("list", 12)
xiInit[[1]] <- matrix(0, nrow = 20, ncol = 303)
xiInit[[2]] <- matrix(0, nrow = 20, ncol = 303)
xiInit[[3]] <- matrix(0, nrow = 20, ncol = 303)

xiInit[[4]] <- rbind(xiInitDum[[1]], matrix(0, nrow = 17, ncol = 303))
xiInit[[5]] <- rbind(xiInitDum[[2]], matrix(0, nrow = 17, ncol = 303))
xiInit[[6]] <- rbind(xiInitDum[[3]], matrix(0, nrow = 17, ncol = 303))

xiInit[[7]] <- rbind(xiInitDum[[4]], matrix(0, nrow = 15, ncol = 303))
xiInit[[8]] <- rbind(xiInitDum[[5]], matrix(0, nrow = 15, ncol = 303))
xiInit[[9]] <- rbind(xiInitDum[[6]], matrix(0, nrow = 15, ncol = 303))

xiInit[[10]] <- xiInitDum[[7]]
xiInit[[11]] <- xiInitDum[[8]]
xiInit[[12]] <- xiInitDum[[9]]

resultName <- c(paste0("result_cleaned_UMich_chain_", 1:3, "_init_oneClus_s2_1en1_s2MH_1en3.rds"),
                paste0("result_cleaned_UMich_chain_", 1:3, "_init_3clus_s2_1en1_s2MH_1en3.rds"),
                paste0("result_cleaned_UMich_chain_", 1:3, "_init_5clus_s2_1en1_s2MH_1en3.rds"),
                paste0("result_cleaned_UMich_chain_", 1:3, "_init_20clus_s2_1en1_s2MH_1en3.rds"))

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
globalTime <- Sys.time()
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 25000, Kmax = 20, nbeta_split = 5,
                      z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 107, ncol = 81),
                      beta_init = as.matrix(xiInit[[t]]),
                      ci_init = ciInit[, t],
                      theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
                      t_thres = 5000, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
}
stopImplicitCluster()
difftime(Sys.time(), globalTime)

## DUMMY: ----------------------

data.frame(x = 1:303, y = colSums(datUMich)/sum(datUMich)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity")

xx <- rnorm(303, 0, sqrt(10))
data.frame(x = 1:303, y = sort(exp(xx)/sum(exp(xx)), decreasing = TRUE)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity")

start_time <- Sys.time()
mod <- mod_adaptive(iter = 3000, Kmax = 10, nbeta_split = 20,
                    z = as.matrix(datUMich), atrisk_init = matrix(1, nrow = 107, ncol = 303),
                    beta_init = matrix(0, nrow = 10, ncol = 303),
                    ci_init = rep(0, 107),
                    theta = 1, mu = 0, s2 = 0.1, s2_MH = 1e-5,
                    t_thres = 1000, launch_iter = 30,
                    r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
difftime(Sys.time(), start_time)

apply(mod$ci_result, 1, uniqueClus) %>% plot(type = "l")
plot(mod$beta_result[1, 1, ], type = "l")

combMeta <- data.frame(ID = as.numeric(rownames(datUMich)), 
                       clus = as.numeric(salso(mod$ci_result[-(1:500), ]))) %>%
  inner_join(metadata, by = join_by(ID == sample))

table(combMeta$clus, combMeta$Dx_Bin)
table(combMeta$clus, combMeta$Gender)
table(combMeta$clus, combMeta$Smoke)




# ### Metadata
# # metadat <- read.delim(paste0(datapath, "claesson/mapping-orig.txt"))
# metadat <- read.delim(paste0(datapath, "bacteremia/task.txt"))
# 
# ### OTU Table
# otuTab <- read.delim(paste0(datapath, "bacteremia/gg/taxatable.txt"))
# View(otuTab)
# taxaName <- otuTab[, 1]
# otuTab <- otuTab[, -1]
# otuTab <- t(otuTab)
# colnames(otuTab) <- taxaName
# dim(otuTab)
# otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
# # View(otuTab)
# dim(otuTab)
# View(otuTab)
# ### Filter only the observation code that appear only once
# 
# data.frame(j = 1:143, p = as.numeric(colSums(otuTab)/sum(otuTab))) %>%
#   ggplot(aes(x = j, y = p)) +
#   geom_bar(stat = "identity")
# 
# start_time <- Sys.time()
# mod <- mod_adaptive(iter = 10000, Kmax = 10, nbeta_split = 20,
#                     z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 28, ncol = 143),
#                     beta_init = matrix(0, nrow = 10, ncol = 143),
#                     ci_init = rep(0, 28),
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                     t_thres = 1000, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# difftime(Sys.time(), start_time)
# 
# apply(mod$ci_result, 1, uniqueClus) %>% plot(type = "l")
# salso(mod$ci_result[-(1:2000), ])
# 
# colSums(otuTab) %>% as.numeric() %>% which.max()
# apply(otuTab, 2, var) %>% as.numeric() %>% which.max()
# 
# t(mod$beta_result[, 93, ]) %>%
#   as.data.frame() %>%
#   mutate(Iter = 1:10000) %>%
#   pivot_longer(!Iter) %>%
#   ggplot(aes(x = Iter, y = value, color = name)) + 
#   geom_line()
# 
# rownames(otuTab)
# metadat$X.SampleID
# 
# ### OTU Table -- turnbaugh
# # otuTab <- read.delim(paste0(datapath, "turnbaugh/refseq/taxatable.txt"))
# # View(otuTab)
# # taxaName <- otuTab[, 1]
# # otuTab <- otuTab[, -1]
# # otuTab <- t(otuTab)
# # colnames(otuTab) <- taxaName
# # dim(otuTab)
# # otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
# # # View(otuTab)
# # dim(otuTab)
# 
# ### KARLSSON
# # otuTab <- read.delim(paste0(datapath, "karlsson/taxatable.txt"))
# # taxaName <- otuTab[, 1]
# # otuTab <- otuTab[, -1]
# # otuTab <- t(otuTab)
# # colnames(otuTab) <- taxaName
# # dim(otuTab)
# # otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
# # # View(otuTab)
# # dim(otuTab)
# 
# ### KARLSON DATASET
# # mod <- mod_adaptive(iter = 500, Kmax = 10, nbeta_split = 30, 
# #                     z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 144, ncol = 1460),
# #                     beta_init = matrix(0, nrow = 10, ncol = 1460),
# #                     ci_init = rep(0, 144),
# #                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-5,
# #                     t_thres = 100, launch_iter = 30,
# #                     r0g = 1, r1g = 1, r0c = 9, r1c = 1, thin = 1)
# 
# start_time <- Sys.time()
# mod <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 10,
#                     z = as.matrix(otuTab), atrisk_init = matrix(1, nrow = 281, ncol = 231),
#                     beta_init = matrix(0, nrow = 10, ncol = 231),
#                     ci_init = rep(0, 281),
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                     t_thres = 100, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# difftime(Sys.time(), start_time)
# saveRDS(mod, file = paste0(resultpath, "test_turnbaugh_1en3_5000iter_thres_100.rds"))
# 
# rm(mod)
# 
# ### Test: Result
# result <- readRDS(paste0(resultpath, "test_turnbaugh_1en10_5000iter_thres_100.rds"))
# clusVI <- as.numeric(salso(result$ci_result))
# clusBinder <- as.numeric(salso(result$ci_result, loss = "binder"))
# 
# table(clusVI)
# table(clusBinder)
# table(clusVI, clusBinder)
# 
# apply(result$ci_result, 1, uniqueClus) %>% plot(type = "l")
# 
# colSums(otuTab) %>% as.numeric() %>% which.max()
# apply(otuTab, 2, var) %>% as.numeric() %>% which.max()
# 
# t(result$beta_result[, 1, ]) %>%
#   as.data.frame() %>%
#   mutate(Iter = 1:5000) %>%
#   pivot_longer(!Iter) %>%
#   ggplot(aes(x = Iter, y = value, color = name)) + 
#   geom_line()
# 
# table(clusVI, metadat$OBESITYCAT)
# 
# boxplot(metadat$AGE)
# table(metadat$OBESITYCAT)
# mean(metadat$AGE)
# sd(metadat$AGE)
# 
# 
# 
# table(metadat$TWIN_MOTHER)
# table(metadat$ZYGOSITY)
# 
# data.frame(age = metadat$AGE, gr = paste0(metadat$TWIN_MOTHER, ifelse(is.na(metadat$ZYGOSITY), "", paste0(" - ", metadat$ZYGOSITY))))
# 
# paste0(metadat$TWIN_MOTHER, ifelse(is.na(metadat$ZYGOSITY), "", paste0(" - ", metadat$ZYGOSITY))) %>%
#   table(clusVI)
# 
# table(metadat$TWIN_MOTHER, clusSALSOtest[, 1])
# 
# data.frame(metadat$TWIN_MOTHER, metadat$ZYGOSITY) %>% View()
# 
# table(clusSALSOtest[, 1], metadat$OBESITYCAT)
# data.frame(age = metadat$AGE, clus = clusSALSOtest[, 1]) %>%
#   group_by(clus) %>%
#   summarise(n = n(), mean(age), sd(age))
# View(metadat)
# 
# 
# 
# ### Metadata
# # read.delim(paste0(datapath, "karlsson/task-impaired-diabetes.txt"))
# # read.delim(paste0(datapath, "karlsson/task-normal-diabetes.txt"))
# # View(metadat)
# # metadat$X.SampleID[which(metadat$Day == 4)]
# # 
# # dat <- read.delim(paste0(datapath, "david/refseq/taxatable.txt"))
# # taxaRef <- dat[, 1]
# # # which(colnames(dat) %in% metadat$X.SampleID[which(metadat$Day == 4)])
# # dat <- dat[, which(colnames(dat) %in% metadat$X.SampleID[which(metadat$Day == 4)])]
# # dat <- t(dat)
# # colnames(dat) <- taxaRef
# # dat <- dat[, -which(colMeans(dat > 0) < 0.1)]
# # 
# # data.frame(ID = 1:393, p = as.numeric(colSums(dat)/sum(dat))) %>%
# #   ggplot(aes(x = ID, y = p)) +
# #   geom_bar(stat = "identity")
# # 
# # set.seed(1)
# # x <- rnorm(400, 0, sqrt(10))
# # data.frame(ID = 1:400, p = exp(x)/sum(exp(x))) %>%
# #   ggplot(aes(x = ID, y = p)) +
# #   geom_bar(stat = "identity")
# # 
# # x[1] <- rnorm(1, 0, sqrt(10))
# # data.frame(ID = 1:400, p = exp(x)/sum(exp(x))) %>%
# #   ggplot(aes(x = ID, y = p)) +
# #   geom_bar(stat = "identity")
# 
# mod <- mod_adaptive(iter = 7500, Kmax = 10, nbeta_split = 30, 
#                     z = as.matrix(dat), atrisk_init = matrix(1, nrow = 17, ncol = 393),
#                     beta_init = matrix(0, nrow = 10, ncol = 393),
#                     ci_init = rep(0, 17),
#                     theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
#                     t_thres = 2700, launch_iter = 30,
#                     r0g = 1, r1g = 1, r0c = 9, r1c = 1, thin = 1)
# 
# apply(mod$ci_result, 1, uniqueClus) %>% plot(type = "l")
# 
# salso(mod$ci_result[-(1:200), ])
# plot(mod$beta_result[1, 1, ], type = "l")
# 
# dietInfo <- metadat[which(metadat$Day == 4), c("X.SampleID", "Diet")]
# colnames(dietInfo)[1] <- "ID"
# 
# # kclus <- data.frame(clus = kmeans(dat, 2)$cluster)
# # kclus %>%
# #   mutate(ID = rownames(kclus)) %>%
# #   inner_join(dietInfo) 
# 
# # dim(dat)
# # dat <- read.delim(paste0(datapath, "david/refseq/taxatable.txt"))
# # dim(dat)
# # 
# # dat <- dat[-which(rowMeans(dat[, -1] > 0) < 0.1), ]
# # taxaName <- dat[, 1]
# # dat <- t(dat)
# # dat <- dat[-1, ]
# # mode(dat) <- "numeric"
# # colnames(dat) <- taxaName
# 
# # # Demographic
# # demoBH <- read.delim(paste0(datapath, "ravel/task-black-hispanic.txt"))
# # demoWB <- read.delim(paste0(datapath, "ravel/task-white-black.txt"))
# # colnames(demoBH)[1] <- "ID"
# # colnames(demoWB)[1] <- "ID"
# # demoWBH <- full_join(demoBH, demoWB)
# # 
# # demoN <- read.delim(paste0(datapath, "ravel/task-nugent-category.txt"))
# # colnames(demoN)[1] <- "ID"
# # colnames(demoN)[2] <- "Nugent"
# # demoWBHN <- full_join(demoWBH, demoN)
# # colnames(demoWBHN) <- c("ID", "Ethnicity", "Nugent")
# # 
# # demoWBHN[is.na(demoWBHN)] <- "No Information"
# # demoWBHN$Ethnicity <- factor(demoWBHN$Ethnicity, levels = c("Black", "Hispanic", "White", "No Information"))
# # demoWBHN$Nugent <- factor(demoWBHN$Nugent, labels = c("High", "Low", "No Information"))
# # 
# # p1 <- data.frame(table(demoWBHN$Ethnicity)) %>%
# #   ggplot(aes(x = Var1, y = Freq)) +
# #   geom_bar(stat = "identity", width = 0.3) +
# #   geom_text(aes(label = Freq), position = position_dodge(width = 0.3), vjust = -0.5) +
# #   theme_bw() +
# #   theme(axis.title.x = element_blank()) +
# #   labs(title = "Ethnicity: ravel dataset from vangay",
# #        y = "Frequency")
# # 
# # p2 <- data.frame(table(demoWBHN$Nugent)) %>%
# #   ggplot(aes(x = Var1, y = Freq)) +
# #   geom_bar(stat = "identity", width = 0.3) +
# #   geom_text(aes(label = Freq), position = position_dodge(width = 0.3), vjust = -0.5) +
# #   theme_bw() +
# #   theme(axis.title.x = element_blank()) +
# #   labs(title = "Predicted Nugent group: ravel dataset from vangay",
# #        y = "Frequency")
# # 
# # grid.arrange(p1, p2, ncol = 2)
# # 
# # dat <- dat[which(rownames(dat) %in% intersect(rownames(dat), demoWBHN$ID)), ]
# # 
# # # Sequencing Dept
# # median(rowSums(dat))
# # min(rowSums(dat))
# # max(rowSums(dat))
# # 
# # which(rowSums(dat) == min(rowSums(dat)))
# # dat[375, ]
# # 
# # ggplot(data.frame(x = rowSums(dat)), aes(x = x)) +
# #   geom_histogram() +
# #   labs(x = "Total Read Count", title = "ravel dataset: Distribution of the Total Read Count") +
# #   theme_bw() +
# #   theme(axis.title.y = element_blank())
# # 
# # # Shannon Diversity and Simpson Index
# # shannon_d <- sapply(1:375, function(x){
# #   pi <- dat[x, dat[x, ] != 0]/sum(dat[x, dat[x, ] != 0])
# #   -sum(pi * log(pi))})
# # 
# # simpson_d <- sapply(1:375, function(x){
# #   pi <- dat[x, dat[x, ] != 0]/sum(dat[x, dat[x, ] != 0])
# #   1 - sum(pi^2)})
# # 
# # numericLonger <- data.frame(trc = rowSums(dat), shannon_d, simpson_d) %>%
# #   mutate(ID = 1:375, Ethnicity = demoWBHN$Ethnicity, Nugent = demoWBHN$Nugent) %>%
# #   pivot_longer(!c(ID, Ethnicity, Nugent))
# # 
# # numericLonger$name <- factor(numericLonger$name, labels = c("Shannon", "Simpson", "Total Read Count"))
# # 
# # ggplot(numericLonger, aes(x = Nugent, y = value)) +
# #   geom_boxplot() +
# #   facet_wrap( ~ name, scales = "free_y") +
# #   theme_bw() +
# #   labs(title = "Total Read Count and alpha-Diversity for each predicted Nugent group") +
# #   theme(axis.title.y = element_blank())
# # 
# # ggplot(numericLonger, aes(x = Ethnicity, y = value)) +
# #   geom_boxplot() +
# #   facet_wrap( ~ name, scales = "free_y") +
# #   theme_bw() +
# #   labs(title = "Total Read Count and alpha-Diversity for each Ethnicity group") +
# #   theme(axis.title.y = element_blank())
# # 
# # ovarallSum <- data.frame(trc = rowSums(dat), shannon_d, simpson_d) %>%
# #   summarise(n(), median(trc), min(trc), max(trc), median(shannon_d), min(shannon_d), 
# #             max(shannon_d), median(simpson_d), min(simpson_d), 
# #             max(simpson_d))
# # 
# # paste0(ovarallSum[, 2], " (", ovarallSum[, 3], ", ", ovarallSum[, 4], ")")
# # paste0(round(ovarallSum[, 5], 4), " (", round(ovarallSum[, 6], 4), ", ", round(ovarallSum[, 7], 4), ")")
# # paste0(round(ovarallSum[, 8], 4), " (", round(ovarallSum[, 9], 4), ", ", round(ovarallSum[, 10], 4), ")")
# # 
# # groupSum <- data.frame(trc = rowSums(dat), shannon_d, simpson_d) %>%
# #   mutate(Nugent = demoWBHN$Nugent) %>%
# #   group_by(Nugent) %>%
# #   summarise(n(), median(trc), min(trc), max(trc), median(shannon_d), min(shannon_d), 
# #             max(shannon_d), median(simpson_d), min(simpson_d), 
# #             max(simpson_d))
# # 
# # paste0(groupSum[1, 3], " (", groupSum[1, 4], ", ", groupSum[1, 5], ")")
# # paste0(groupSum[2, 3], " (", groupSum[2, 4], ", ", groupSum[2, 5], ")")
# # paste0(groupSum[3, 3], " (", groupSum[3, 4], ", ", groupSum[3, 5], ")")
# # 
# # paste0(round(groupSum[1, 6], 4), " (", round(groupSum[1, 7], 4), ", ", round(groupSum[1, 8], 4), ")")
# # paste0(round(groupSum[2, 6], 4), " (", round(groupSum[2, 7], 4), ", ", round(groupSum[2, 8], 4), ")")
# # paste0(round(groupSum[3, 6], 4), " (", round(groupSum[3, 7], 4), ", ", round(groupSum[3, 8], 4), ")")
# # 
# # paste0(round(groupSum[1, 9], 4), " (", round(groupSum[1, 10], 4), ", ", round(groupSum[1, 11], 4), ")")
# # paste0(round(groupSum[2, 9], 4), " (", round(groupSum[2, 10], 4), ", ", round(groupSum[2, 11], 4), ")")
# # paste0(round(groupSum[3, 9], 4), " (", round(groupSum[3, 10], 4), ", ", round(groupSum[3, 11], 4), ")")
# # 
# # which.min(shannon_d)
# # which.min(simpson_d)
# # which.min(rowSums(dat))
# # which.max(rowSums(dat))
# # which.max(shannon_d)
# # which.max(simpson_d)
# # 
# # # ### Default set of Hyperparameters
# # # set.seed(1, kind = "L'Ecuyer-CMRG")
# # # registerDoParallel(6)
# # # globalTime <- Sys.time()
# # # foreach(t = 1:6) %dopar% {
# # #   start_time <- Sys.time()
# # #   mod <- mod_adaptive(iter = 25000, Kmax = 10, nbeta_split = 5,
# # #                       z = as.matrix(dat), atrisk_init = matrix(1, nrow = 375, ncol = 56),
# # #                       beta_init = matrix(0, nrow = 10, ncol = 56),
# # #                       ci_init = rep(0, 375),
# # #                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
# # #                       t_thres = 2500, launch_iter = 30,
# # #                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# # #   comp_time <- difftime(Sys.time(), start_time, units = "secs")
# # #   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, "result_ravel_chain_", t, "_init_oneClus_defaultHyper.rds"))
# # # }
# # # stopImplicitCluster()
# # # difftime(Sys.time(), globalTime)
# # # 
# # # ### Different starting point: 1e-3 - PART 1
# # # set.seed(1)
# # # ciInit <- matrix(NA, nrow = 375, ncol = 6)
# # # ciInit[, 1] <- sample(0:4, 375, replace = TRUE)
# # # ciInit[, 2] <- sample(0:4, 375, replace = TRUE)
# # # ciInit[, 3] <- sample(0:4, 375, replace = TRUE)
# # # ciInit[, 4] <- sample(0:19, 375, replace = TRUE)
# # # ciInit[, 5] <- sample(0:19, 375, replace = TRUE)
# # # ciInit[, 6] <- sample(0:19, 375, replace = TRUE)
# # # 
# # # KmaxVec <- c(20, 20, 20, 50, 50, 50)
# # # 
# # # xiInit <- lapply(1:6, function(y){sapply(0:max(ciInit[, y]), function(x){
# # #   colSums(dat[which(ciInit[, y] == x), ])/sum(dat[which(ciInit[, y] == x), ])
# # # }) %>% t()
# # # })
# # # 
# # # xiInit[[1]] <- rbind(xiInit[[1]], matrix(0, nrow = 15, ncol = 56))
# # # xiInit[[2]] <- rbind(xiInit[[2]], matrix(0, nrow = 15, ncol = 56))
# # # xiInit[[3]] <- rbind(xiInit[[3]], matrix(0, nrow = 15, ncol = 56))
# # # xiInit[[4]] <- rbind(xiInit[[4]], matrix(0, nrow = 30, ncol = 56))
# # # xiInit[[5]] <- rbind(xiInit[[5]], matrix(0, nrow = 30, ncol = 56))
# # # xiInit[[6]] <- rbind(xiInit[[6]], matrix(0, nrow = 30, ncol = 56))
# # # 
# # # resultName <- c(paste0("result_ravel_chain_", 1:3, "_init_5clus_Kmax_20_defaultHyper.rds"),
# # #                 paste0("result_ravel_chain_", 1:3, "_init_20clus_Kmax_50_defaultHyper.rds"))
# # # 
# # # set.seed(1, kind = "L'Ecuyer-CMRG")
# # # registerDoParallel(6)
# # # globalTime <- Sys.time()
# # # foreach(t = 1:6) %dopar% {
# # #   start_time <- Sys.time()
# # #   mod <- mod_adaptive(iter = 25000, Kmax = KmaxVec[t], nbeta_split = 5,
# # #                       z = as.matrix(dat), atrisk_init = matrix(1, nrow = 375, ncol = 56),
# # #                       beta_init = as.matrix(xiInit[[t]]),
# # #                       ci_init = ciInit[, t],
# # #                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
# # #                       t_thres = 2500, launch_iter = 30,
# # #                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# # #   comp_time <- difftime(Sys.time(), start_time, units = "secs")
# # #   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
# # # }
# # # stopImplicitCluster()
# # # difftime(Sys.time(), globalTime)
# # # 
# # # ### Different starting point: 1e-3 - PART 2
# # # set.seed(1)
# # # ciInit <- matrix(NA, nrow = 375, ncol = 3)
# # # ciInit[, 1] <- sample(0:2, 375, replace = TRUE)
# # # ciInit[, 2] <- sample(0:2, 375, replace = TRUE)
# # # ciInit[, 3] <- sample(0:2, 375, replace = TRUE)
# # # 
# # # xiInit <- lapply(1:3, function(y){sapply(0:max(ciInit[, y]), function(x){
# # #   colSums(dat[which(ciInit[, y] == x), ])/sum(dat[which(ciInit[, y] == x), ])
# # # }) %>% t()
# # # })
# # # 
# # # xiInit[[1]] <- rbind(xiInit[[1]], matrix(0, nrow = 12, ncol = 56))
# # # xiInit[[2]] <- rbind(xiInit[[2]], matrix(0, nrow = 12, ncol = 56))
# # # xiInit[[3]] <- rbind(xiInit[[3]], matrix(0, nrow = 12, ncol = 56))
# # # 
# # # resultName <- paste0("result_ravel_chain_", 1:3, "_init_3clus_Kmax_15_defaultHyper.rds")
# # # 
# # # set.seed(1, kind = "L'Ecuyer-CMRG")
# # # registerDoParallel(3)
# # # globalTime <- Sys.time()
# # # foreach(t = 1:3) %dopar% {
# # #   start_time <- Sys.time()
# # #   mod <- mod_adaptive(iter = 25000, Kmax = 15, nbeta_split = 5,
# # #                       z = as.matrix(dat), atrisk_init = matrix(1, nrow = 375, ncol = 56),
# # #                       beta_init = as.matrix(xiInit[[t]]),
# # #                       ci_init = ciInit[, t],
# # #                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3,
# # #                       t_thres = 2500, launch_iter = 30,
# # #                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
# # #   comp_time <- difftime(Sys.time(), start_time, units = "secs")
# # #   saveRDS(list(time = comp_time, mod = mod), file = paste0(resultpath, resultName[t]))
# # # }
# # # stopImplicitCluster()
# # # difftime(Sys.time(), globalTime)
# # 
# # # ### Post Analysis: -------------------------------------------------------------
# # resultFilename <- c(paste0(resultpath, "result_ravel_chain_", c(1, 2, 5), "_init_oneClus_defaultHyper.rds"),
# #                     paste0(resultpath, "result_ravel_chain_", 1:3, "_init_3clus_Kmax_15_defaultHyper.rds"),
# #                     paste0(resultpath, "result_ravel_chain_", 1:3, "_init_5clus_Kmax_20_defaultHyper.rds"),
# #                     paste0(resultpath, "result_ravel_chain_", 1:3, "_init_20clus_Kmax_50_defaultHyper.rds"))
# # 
# # ### Computational time
# # registerDoParallel(3)
# # compTime <- foreach(t = 1:12, .combine = cbind) %dopar% {
# #   result <- readRDS(resultFilename[t])
# #   as.numeric(result$time)
# # }
# # stopImplicitCluster()
# # 
# # mean(compTime/3600)
# # sd(compTime/3600)
# # 
# # ### Active Cluster - Combine both two hyperparameters
# # registerDoParallel(4)
# # activeClusMat <- foreach(t = 1:12, .combine = cbind) %dopar% {
# #   result <- readRDS(resultFilename[t])
# #   apply(result$mod$ci_result, 1, uniqueClus)
# # }
# # stopImplicitCluster()
# # 
# # activeClusMatPlot <- activeClusMat %>%
# #   as.data.frame() %>%
# #   mutate(iter = 1:25000) %>%
# #   pivot_longer(!iter)
# # 
# # activeClusMatPlot$name <- factor(activeClusMatPlot$name, 
# #                                  levels = paste0("result.", 1:12), labels = paste0("Chain ", 1:12))
# # 
# # ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
# #   geom_line() +
# #   theme_bw() +
# #   theme(legend.position = "bottom") +
# #   labs(title = "Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusteres",
# #        color = "MCMC Chain") +
# #   guides(color = guide_legend(ncol = 12))
# # 
# # ### Check the convergence of xi
# # data.frame(colMeans(otuHIV), apply(otuHIV, 2, var))
# # 
# # registerDoParallel(2)
# # xiFirst <- foreach(t = 1:12) %dopar% {
# #   result <- readRDS(resultFilename[t])
# #   result$mod$beta_result[, 1, ] %>% t() %>%
# #     as.data.frame() %>%
# #     mutate(Iteration = 1:25000) %>%
# #     pivot_longer(!Iteration) %>%
# #     transmute(Iteration, xi = value,
# #               Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
# # }
# # stopImplicitCluster()
# # 
# # xiFirstLong <- lapply(1:12, function(x){data.frame(xiFirst[[x]], chain = paste0("Chain ", x))}) %>%
# #   bind_rows()
# # 
# # xiFirstLong$chain <- factor(xiFirstLong$chain, levels = paste0("Chain ", 1:12))
# # 
# # ggplot(xiFirstLong, aes(x = Iteration, y = xi, color = Cluster)) +
# #   geom_line() +
# #   facet_wrap(. ~ chain) +
# #   theme_bw() +
# #   theme(legend.position = "none") +
# #   labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k1}$")), x = "Iteration", y = TeX("\\xi"))
# # 
# # registerDoParallel(2)
# # xiThird <- foreach(t = 1:12) %dopar% {
# #   result <- readRDS(resultFilename[t])
# #   result$mod$beta_result[, 2, ] %>% t() %>%
# #     as.data.frame() %>%
# #     mutate(Iteration = 1:25000) %>%
# #     pivot_longer(!Iteration) %>%
# #     transmute(Iteration, xi = value,
# #               Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+")))
# # }
# # stopImplicitCluster()
# # 
# # xiThirdLong <- lapply(1:12, function(x){data.frame(xiThird[[x]], chain = paste0("Chain ", x))}) %>%
# #   bind_rows()
# # 
# # xiThirdLong$chain <- factor(xiThirdLong$chain, levels = paste0("Chain ", 1:12))
# # 
# # ggplot(xiThirdLong, aes(x = Iteration, y = xi, color = Cluster)) +
# #   geom_line() +
# #   facet_wrap(. ~ chain) +
# #   theme_bw() +
# #   theme(legend.position = "none") +
# #   labs(title = TeX(paste0("Trace plot: ", "$\\xi_{k2}$")), x = "Iteration", y = TeX("\\xi"))
# # 
# # ### Cluster Assignment - Individual
# # set.seed(1, kind = "L'Ecuyer-CMRG")
# # registerDoParallel(2)
# # clusSALSO <- foreach(t = 1:12, .combine = cbind) %dopar% {
# #   result <- readRDS(resultFilename[t])
# #   as.numeric(salso(result$mod$ci_result[-(1:5000), ]))
# # }
# # stopImplicitCluster()
# # 
# # ### Individual - Demographic
# # demoLong <- lapply(1:12, function(x){
# #   dumDat <- data.frame(ID = rownames(dat), clus = clusSALSO[, x]) %>%
# #     inner_join(demoWBHN)
# #   DemoPercent <- data.frame(table(dumDat$clus, dumDat$Nugent)) %>%
# #     group_by(Var1) %>%
# #     mutate(Percent = Freq/sum(Freq), Chain = paste0("Chain ", x)) %>%
# #     arrange(Var1)
# #   DemoPercent
# # }) %>%
# #   bind_rows()
# # 
# # demoLong$Chain <- factor(demoLong$Chain, levels = paste0("Chain ", 1:12))
# # demoLong$Var1 <- factor(demoLong$Var1, labels = paste0("Cluster ", 1:6))
# # 
# # demoLong %>%
# #   ggplot(aes(x = Var1, y = Percent, fill = Var2)) +
# #   geom_bar(stat = "identity") +
# #   geom_text(aes(label = paste0("n = ",Freq)), position = position_stack(vjust = 0.5), size = 2) +
# #   scale_y_continuous(labels = scales::percent) +
# #   theme_bw() +
# #   theme(legend.position = "bottom") +
# #   scale_fill_manual(values = c("springgreen3", "coral1", "gray90")) +
# #   facet_wrap(. ~ Chain, scales = "free_x") +
# #   labs(title = "The distribution of the predicted Nugent score for each cluster in every MCMC chain.",
# #        x = "Cluster", fill = "Predicted Nugent Score Group")
# # 
# # 
# # #### Individual - Relative Abundance
# # taxaName <- data.frame(c = str_remove(str_extract(colnames(dat), "c__[:alpha:]+"), "c__"),
# #                        o = str_remove(str_extract(colnames(dat), "o__[:alpha:]+"), "o__"),
# #                        f = str_remove(str_extract(colnames(dat), "f__[:alpha:]+"), "f__"),
# #                        g = str_remove(str_extract(colnames(dat), "g__[:alpha:]+"), "g__"),
# #                        s = str_remove(str_extract(colnames(dat), "s__[:alpha:]+[:space:]*[:alpha:]*"), "s__"))
# # 
# # highTaxa <- sapply(1:6, function(x){
# #   (dat[which(clusSALSO[, 1] == x), ]/rowSums(dat[which(clusSALSO[, 1] == x), ])) %>%
# #     colMeans() %>%
# #     sort(decreasing = TRUE, index.return = TRUE) %>%
# #     .$ix %>%
# #     .[1:8]})
# # 
# # otu_visual <- dat/rowSums(dat)
# # highTaxaIndex <- union(union(highTaxa[, 1], highTaxa[, 2]), highTaxa[, 3])
# # highTaxaIndexLabel <- ifelse(is.na(taxaName[highTaxaIndex, "s"]), taxaName[highTaxaIndex, "g"], taxaName[highTaxaIndex, "s"])
# # 
# # otuRela <- matrix(NA, ncol = length(highTaxaIndex) + 1, nrow = 375)
# # for(i in 1:length(highTaxaIndex)){
# #   otuRela[, i] <- otu_visual[, highTaxaIndex[i]]
# # }
# # otuRela[, length(highTaxaIndex) + 1] <- 1 - rowSums(otuRela[, 1:length(highTaxaIndex)])
# # colnames(otuRela) <- c(highTaxaIndexLabel, "Others")
# # 
# # otuRelaPlot <- otuRela %>%
# #   as.data.frame() %>% 
# #   mutate(ID = str_extract(rownames(dat), "[:digit:]+"), cluster = paste0("Cluster ", clusSALSO[, 1])) %>%
# #   pivot_longer(!c(ID, cluster))
# # 
# # otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(highTaxaIndexLabel, "Others"))
# # # otuRelaPlot$ID <- factor(otuRelaPlot$ID, levels = c(rownames(otuHIVvisual)[which(clusSALSO[, 1] == 1)][sort(otuHIVvisual[which(clusSALSO[, 1] == 1), 3], decreasing = TRUE, index.return = TRUE)$ix],
# # #                                                     rownames(otuHIVvisual)[which(clusSALSO[, 1] == 2)][sort(otuHIVvisual[which(clusSALSO[, 1] == 2), 6], decreasing = TRUE, index.return = TRUE)$ix],
# # #                                                     rownames(otuHIVvisual)[which(clusSALSO[, 1] == 3)][sort(otuHIVvisual[which(clusSALSO[, 1] == 3), 1], decreasing = TRUE, index.return = TRUE)$ix]))
# # # 
# # # rownames(otuHIVvisual)[which(clusSALSO[, 1] == 1)][sort(otuHIVvisual[which(clusSALSO[, 1] == 1), 3], decreasing = TRUE, index.return = TRUE)$ix]
# # # rownames(otuHIVvisual)[which(clusSALSO[, 1] == 2)][sort(otuHIVvisual[which(clusSALSO[, 1] == 2), 6], decreasing = TRUE, index.return = TRUE)$ix]
# # # rownames(otuHIVvisual)[which(clusSALSO[, 1] == 3)][sort(otuHIVvisual[which(clusSALSO[, 1] == 3), 1], decreasing = TRUE, index.return = TRUE)$ix]
# # 
# # ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
# #   geom_bar(position = "stack", stat = "identity") +
# #   facet_wrap(. ~ cluster, scales = "free_x") +
# #   theme_bw() +
# #   scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
# #   scale_y_continuous(labels = scales::percent) +
# #   theme(axis.text.x = element_text(angle = 90)) +
# #   theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
# #   guides(fill = guide_legend(nrow = 2)) +
# #   labs(fill = "Taxa", y = "Relative Abundance",
# #        title = paste0("Relative Abundance for each cluster"))
# # 
# # relaPlot <- lapply(1:12, function(x){
# #   otuRelaPlot <- otuRela %>%
# #     as.data.frame() %>% 
# #     mutate(ID = str_extract(rownames(dat), "[:digit:]+"), cluster = paste0("Cluster ", clusSALSO[, x])) %>%
# #     pivot_longer(!c(ID, cluster))
# #   otuRelaPlot$name <- factor(otuRelaPlot$name, levels = c(highTaxaIndexLabel, "Others"))
# #   ggplot(otuRelaPlot, aes(x = ID, y = value, fill = name)) +
# #     geom_bar(position = "stack", stat = "identity") +
# #     facet_grid(. ~ cluster, scales = "free_x") +
# #     theme_bw() +
# #     scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
# #     scale_y_continuous(labels = scales::percent) +
# #     theme(axis.text.x = element_text(angle = 90)) +
# #     theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
# #     guides(fill = guide_legend(nrow = 2)) +
# #     labs(fill = "Taxa", y = "Relative Abundance",
# #          title = paste0("Chain ", x, " - Relative Abundance for each cluster"))
# # })
# # 
# # ### Combine all chains
# # registerDoParallel(4)
# # MCMCCombine <- foreach(t = 1:12, .combine = rbind) %dopar% {
# #   result <- readRDS(resultFilename[t])
# #   result$mod$ci_result[15001:25000, ]
# # }
# # stopImplicitCluster()
# # 
# # set.seed(1)
# # clusComb <- as.numeric(salso(MCMCCombine))
# # 
# # #### Demographic - Combined Chain
# # combineDemo <- data.frame(ID = rownames(dat), clus = clusComb) %>%
# #   inner_join(demoWBHN)
# # 
# # table(combineDemo$clus, combineDemo$Nugent)
# # 
# # combineDemoPercent <- data.frame(table(combineDemo$clus, combineDemo$Nugent)) %>%
# #   group_by(Var1) %>%
# #   mutate(Percent = Freq/sum(Freq)) %>%
# #   arrange(Var1)
# # 
# # combineDemoPercent$Var1 <- factor(combineDemoPercent$Var1, labels = paste0("Cluster ", 1:5))
# # 
# # ggplot(combineDemoPercent, aes(x = Var1, y = Percent, fill = Var2)) +
# #   geom_bar(stat = "identity") +
# #   geom_text(aes(label = paste0("n = ",Freq)), position = position_stack(vjust = 0.5)) +
# #   scale_y_continuous(labels = scales::percent) +
# #   theme_bw() +
# #   theme(legend.position = "bottom") +
# #   scale_fill_manual(values = c("springgreen3", "coral1", "gray90")) +
# #   labs(title = "The distribution of the predicted Nugent score for each cluster when combining MCMC chain.",
# #        x = "Cluster", fill = "Predicted Nugent Score Group") 
# # 
# # #### Relative Abundance - Combined Chain
# # otuRelaCombPlot <- otuRela %>%
# #   as.data.frame() %>% 
# #   mutate(ID = str_extract(rownames(dat), "[:digit:]+"), cluster = paste0("Cluster ", clusComb)) %>%
# #   pivot_longer(!c(ID, cluster))
# # 
# # otuRelaCombPlot$name <- factor(otuRelaCombPlot$name, levels = c(highTaxaIndexLabel, "Others"))
# # 
# # ggplot(otuRelaCombPlot, aes(x = ID, y = value, fill = name)) +
# #   geom_bar(position = "stack", stat = "identity") +
# #   facet_wrap(. ~ cluster, scales = "free_x") +
# #   theme_bw() +
# #   scale_fill_manual(values = c(turbo(n = 15), "gray90")) + 
# #   scale_y_continuous(labels = scales::percent) +
# #   theme(axis.text.x = element_text(angle = 90)) +
# #   theme(legend.position = "bottom") +
# #   guides(fill = guide_legend(nrow = 2)) +
# #   labs(fill = "Taxa", y = "Relative Abundance",
# #        title = paste0("Relative Abundance for each cluster of the combined MCMC chain"))
# # 
# # data.frame(trc = rowSums(dat), shannon_d, simpson_d) %>%
# #   mutate(ID = 1:375, clusComb) %>%
# #   group_by(clusComb) %>%
# #   summarise(n(), median(trc), min(trc), max(trc), median(shannon_d), min(shannon_d), 
# #             max(shannon_d), median(simpson_d), min(simpson_d), 
# #             max(simpson_d)) %>% View()

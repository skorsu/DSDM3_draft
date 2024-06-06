# devtools::install_github("YushuShi/MicrobiomeCluster")

library(ape)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(salso)

### Global Objects: ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

### Papa 2012, IBD: ------------------------------------------------------------
# file.exists(paste0(path, "Manuscript/Data/Application Data/ibd_alm_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/ibd_alm_results.tar"))
info <- read.table(paste0(path, "ibd_alm_results/ibd_alm.metadata.txt"), sep = "\t", header = TRUE) %>%
  dplyr::select(-X)
info <- data.frame(ID = str_replace(str_extract(info$sample, pattern = "[:digit:]+\\-[:alpha:]"), "\\-", ""),
                   info %>% dplyr::select(-sample))
taxa <- read.table(paste0(path, "ibd_alm_results/RDP/ibd_alm.otu_table.100.denovo.rdp_assigned"))
taxa <- t(taxa)
taxa <- taxa[, -which(colMeans(taxa > 0) < 0.1)]

dim(taxa)

### K-means
plot(sapply(2:10, function(x){sqrt(sum(kmeans(dist(taxa), x)$withinss))}), type = "b")
kmean_clus <- kmeans(dist(taxa), 2)$cluster
kmean_info <- data.frame(ID = str_extract(names(kmean_clus), "[:digit:]+[:alpha:]"),
                         clus = kmean_clus) %>%
  inner_join(info)
table(kmean_info$clus, kmean_info$DiseaseState)

### Our model
set.seed(1)
start_time <- Sys.time()
ourMod_clus <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 100, 
                            z = as.matrix(taxa), 
                            atrisk_init = matrix(1, nrow = 91, ncol = 616), 
                            beta_init = matrix(0, nrow = 5, ncol = 616), 
                            ci_init = rep(0, 91), 
                            theta = 1, mu = 0, s2 = 1e-3, 
                            s2_MH = 1e-5, t_thres = 2500, launch_iter = 30, 
                            r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 10)
difftime(Sys.time(), start_time)

# salso(resultSmits$ci_result, loss = "binder") %>% table()

colSums(resultSmits$MH_accept == 1)/colSums(resultSmits$MH_accept != -1)




plot(apply(resultSmits$ci_result, 1, function(x){length(unique(x))}), type = "l")

table(summary(salso(resultSmits$ci_result, loss = "binder"))$estimate)

ourMod <- inner_join(data.frame(ID = str_extract(rownames(taxa), pattern = "[:digit:]+[:alpha:]"),
                                cluster = as.numeric(salso(resultSmits$ci_result, loss = "binder"))),
                     data.frame(ID = str_replace(str_extract(info$sample, pattern = "[:digit:]+\\-[:alpha:]"), "\\-", ""),
                                info %>% dplyr::select(-sample)))

table(ourMod$cluster, ourMod$DiseaseState)

str_replace(str_extract(info$sample, pattern = "[:digit:]+\\-[:alpha:]"), "\\-", "")

data.frame(ID = rownames(taxa), 
           cluster = as.numeric(salso(resultSmits$ci_result, loss = "binder")))

rownames(taxa)

rownames(info)

info$sample

table(info$gender)
table(info$DiseaseState)

table(resultSmits$ci_result[99, ], resultSmits$ci_result[100, ])

### Smits Dataset: -------------------------------------------------------------
load(paste0(path, "Manuscript/Data/Application Data/Smits/Smits.RData")) 
datSmits <- t(Smits$otutable)
dim(datSmits)

mean(colMeans(datSmits > 0) < 0.1) ### Non-zero proportion
datSmitsADJ <- datSmits[, -which(colMeans(datSmits > 0) < 0.9)] ### Filter: not less than 10% non-zero columns

dim(datSmitsADJ)

set.seed(1)
resultSmits <- mod_adaptive(iter = 1000, Kmax = 5, nbeta_split = 10, 
                            z = as.matrix(datSmitsADJ), 
                            atrisk_init = matrix(1, nrow = 259, ncol = 47), 
                            beta_init = matrix(0, nrow = 5, ncol = 47), 
                            ci_init = rep(0, 259), 
                            theta = 1, mu = 0, s2 = 1e-5, 
                            s2_MH = 1e-3, t_thres = 500, launch_iter = 30, 
                            r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)
apply(resultSmits$ci_result, 1, function(x){length(unique(x))})



table(Smits$sampleinfo, resultSmits$ci_result[12, ])

table(Smits$sampleinfo, salso(resultSmits$ci_result, loss = "binder"))

data.frame(LD = colMeans(datSmitsADJ[which(Smits$sampleinfo == "Late Dry"), ]),
           ED = colMeans(datSmitsADJ[-which(Smits$sampleinfo == "Late Dry"), ]))

names(Smits$sampleinfo) == rownames(datSmitsADJ)


table(Smits$sampleinfo)
datSmits <- datSmits[, -which(colMeans(datSmits == 0) == 1)]


mean(datSmits == 0) ### Proportion of zero
colMeans(datSmits == 0) %>%
  hist()

prop_zero <- seq(0, 1, 0.01)
zeroOVA <- rep(NA, length(prop_zero))
propCol <- rep(NA, length(prop_zero))

for(i in 1:length(prop_zero)){
  colInt <- colMeans(datSmits == 0) < prop_zero[i]
  propCol[i] <- mean(colInt)
  zeroOVA[i] <- mean(datSmits[, colInt] == 0)
}

data.frame(prop_zero, zeroOVA, propCol) %>%
  pivot_longer(!prop_zero) %>%
  dplyr::mutate(name = factor(name, levels = c("zeroOVA", "propCol"),
                              labels = c("Overall zero", "Columns included"))) %>%
  ggplot(aes(x = prop_zero, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "The maximum proportion of zero for each columns allowed into the analysis",
       y = "Proportion",
       title = "Smits: Taxa selection") +
  theme(legend.position = "bottom", legend.title = element_blank())

propCol[which(prop_zero == 0.9)] * dim(datSmits)

relaTab <- datSmits[, colMeans(datSmits == 0) < 0.75]/rowSums(datSmits[, colMeans(datSmits == 0) < 0.75])
data.frame(x = colnames(relaTab), y = relaTab[1, ], index = 1) %>%
  `rownames<-`(NULL) %>%
  ggplot(aes(x = index, y = y, fill = x)) +
  geom_bar(position = "fill", stat = "identity")

relaTab[, which(colMeans(relaTab) > 0.01)] %>%
  data.frame() %>%
  dplyr::mutate(Other = 1 - rowSums(relaTab[, which(colMeans(relaTab) > 0.01)]),
                Sample = rownames(relaTab[, which(colMeans(relaTab) > 0.01)])) %>%
  pivot_longer(!Sample) %>%
  ggplot(aes(x = factor(Sample), y = value, fill = factor(name))) +
  geom_bar(position = "fill", stat = "identity")

1 - rowSums(relaTab[, which(colMeans(relaTab) > 0.001)])

hist(colMeans(relaTab))

which()

propCol[prop_zero == 0.75] * 12000
table(kmeans(dist(datSmits[, colMeans(datSmits == 0) < 0.75]), 2)$cluster, Smits$sampleinfo)

data.frame(x = names(datSmits[1, colMeans(datSmits == 0) < 0.75]),
           y = datSmits[1, colMeans(datSmits == 0) < 0.75]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point()

colMeans(datSmits == 0) < 0.5
sum(colMeans(datSmits == 0) < 0.75)
datSmits[, as.numeric(which(colMeans(datSmits == 0) < 0.75))] %>% View()

# t1 <- read.csv(paste0(path, "Manuscript/Data/Application Data/Smits/aan4834_table_s4.csv"))
# t2 <- read.csv(paste0(path, "Manuscript/Data/Application Data/Smits/aan4834_table_s3.csv"))

tab <- read.csv(paste0(path, "Manuscript/Data/Application Data/Smits/aan4834_table_s1.csv"))
tab$X.SampleID

str_extract(names(Smits$sampleinfo), "^[:alpha:]+") %>%
  as.factor() %>%
  table()

names(Smits$sampleinfo) %in% t1$Sample

colnames(t1)
colnames(t2)

MicrobiomeCluster::convCounts2Abun()

# aan4834_table_s4.csv

Smits$sampleinfo

Smits$otutable


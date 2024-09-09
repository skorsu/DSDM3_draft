### Libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(salso)
library(selbal)
library(taxize)
library(colorspace)
# library(ClusterZI)

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

add <- function(x) Reduce("+", x)

### Dataset: -------------------------------------------------------------------
#### HIV dataset
##### Import the data from selbal packages
datHIV <- HIV
sexHIV <- factor(datHIV[, 61], labels = c("non-MSM", "MSM"))
statusHIV <- factor(datHIV[, 62], labels = c("Healthy", "HIV-patient"))
otuHIV <- datHIV[, 1:60]

### The column name are recorded as a genus. Obtain the family and phylum from
### NCBI database.

taxaName <- data.frame(c = str_remove(str_extract(colnames(otuHIV), "c_[:alpha:]+"), "c_"),
                       o = str_remove(str_extract(colnames(otuHIV), "o_[:alpha:]+"), "o_"),
                       f = str_remove(str_extract(colnames(otuHIV), "f_[:alpha:]+"), "f_"),
                       g = str_remove(str_extract(colnames(otuHIV), "g_[:alpha:]+"), "g_"))

searchTERM <- ifelse(! is.na(taxaName$f), taxaName$f, 
                     ifelse(taxaName$g == "unclassified", NA, taxaName$g))

groupDB <- data.frame(matrix(NA, ncol = 2, nrow = 60))
colnames(groupDB) <- c("Phylum", "Family")

set.seed(1)
for(i in 1:60){
  if(!is.na(searchTERM[i])){
    dbResult <- classification(searchTERM[i], db = "ncbi")
    if(is.data.frame(dbResult[[1]])){
      groupDB[i, ] <- c(subset(dbResult[[1]], rank == "phylum")$name,
                        subset(dbResult[[1]], rank == "family")$name)
    }
  }
}

# groupDB[which(is.na(searchTERM)), ]

### Basic Info
length(sexHIV)
table(sexHIV)
table(statusHIV)
dim(otuHIV)
mean(otuHIV == 0)

path <- "/Users/kevinkvp/Desktop/Github Repo/Manuscript/"
# path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/"

resultFilename <- c(paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1:2, "_init_oneClus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_3clus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_5clus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_20clus_JUL10_fixed.rds"))

##### Obtain the cluster assignment from combined chains.
registerDoParallel(5)
MCMCCombine <- foreach(t = 1:5, .combine = rbind) %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$ci_result[15001:25000, ]
}
stopImplicitCluster()

set.seed(1)
clusComb <- as.numeric(salso(MCMCCombine))

#### Obtain the estimated Relative Abundance
zDEPTH <- as.numeric(rowSums(otuHIV))
estSUMz_allCHAIN <- vector("list", 5)

for(z in 1:5){
  
  resultTEST <- readRDS(resultFilename[z])
  
  registerDoParallel(5)
  estSUMz_allCHAIN[[z]] <- foreach(i = 1:155, .combine = rbind) %dopar% {
    ci <- resultTEST$mod$ci_result[15001:25000, i] ### Cluster Assignment
    at_risk <- t(resultTEST$mod$atrisk_result[i, , 15001:25000]) ### At-risk indicator
    clus_conc <- t(sapply(1:10000, function(x){resultTEST$mod$beta_result[ci[x] + 1, , 15000 + x]})) ### Cluster Concentration
    
    set.seed(1)
    estUNNORMprob <- sapply(1:10000, function(t){
      sapply(1:60, function(j){ifelse(at_risk[t, j] == 0, 0, rgamma(1, exp(clus_conc[t, j]), 1))})
    }) %>% t()
    
    estNORMprob <- estUNNORMprob/rowSums(estUNNORMprob)
    
    set.seed(1)
    estZ <- sapply(1:10000, function(x){
      rmultinom(1, zDEPTH[i], estNORMprob[x, ])
    })
    
    rowSums(estZ)
    
  }
  stopImplicitCluster()
  
  rm(resultTEST)
  
}

### Get the estimated Relative Abundance
estSUMz <- add(estSUMz_allCHAIN)
estRela <- estSUMz/rowSums(estSUMz)

rownames(estRela) <- rownames(otuHIV)
groupDB[is.na(groupDB)] <- "Others"

### Obtain the estimated aggregate relative abundance
agg_estRela <- matrix(NA, nrow = 155, ncol = 32)
for(i in 1:32){
  
  jIndex <- which(sapply(1:60, function(x){sum(groupDB[x, ] == distinct(groupDB)[i, ]) == 2}))
  if(length(jIndex) == 1){
    agg_estRela[, i] <- estRela[, jIndex]
  } else {
    agg_estRela[, i] <- rowSums(estRela[, jIndex])
  }
  
}

colTab <- distinct(groupDB)
labNames <- c("Prevotellaceae", "Oscillospiraceae", "Bacteroidaceae", "Lachnospiraceae",
              "Succinivibrionaceae", "Tannerellaceae", "Erysipelotrichaceae", "Rikenellaceae",
              "Acidaminococcaceae", "Barnesiellaceae", "Coprobacillaceae", "Veillonellaceae",
              "Clostridiaceae", "Peptostreptococcaceae", "Odoribacteraceae", "Selenomonadaceae",
              "Victivallaceae", "Sutterellaceae", "Bifidobacteriaceae", "Enterobacteriaceae", 
              "Elusimicrobiaceae", "Streptococcaceae", "Desulfovibrionaceae", "Christensenellaceae", 
              "Bacillota", "Brachyspiraceae", "Thalassospiraceae", "Anaeroplasmataceae",
              "Defluviitaleaceae", "Porphyromonadaceae", "Coriobacteriaceae", "Others")

colnames(agg_estRela) <- colTab$Family

aggRELA_long <- data.frame(agg_estRela, ID = rownames(otuHIV)) %>%
  inner_join(data.frame(ID = rownames(otuHIV), Cluster = paste0("Cluster ", clusComb))) %>%
  # mutate(ID = singhEDD$SampleID, Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!c(ID, Cluster))

aggRELA_long$name <- factor(aggRELA_long$name, levels = labNames)
# finalAGG <- inner_join(aggRELA_long, colorDistinct)
# finalAGG$name <- factor(finalAGG$name, levels = colorDistinct[, 1])

testCOL <- sequential_hcl(29, "PuBuGn")
plotCOL <- c("red", testCOL[1], "yellow", testCOL[-1], "grey90")

aggRELA_long %>%
  ggplot(aes(x = ID, y = value, fill = name)) +
  geom_bar(stat = "identity", linewidth = 0.25) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 17.5),
        axis.text.y = element_text(size = 15), legend.title = element_text(size = 20)) +
  scale_fill_manual(values = plotCOL) +
  guides(color = guide_legend(order = 1), 
         fill = guide_legend(nrow = 4, order = 2)) +
  labs(y = "Relative Abundances", x = "Participants", fill = "Taxon", color = "Order") +
  facet_grid(. ~ Cluster, scales = "free")

# sapply(1:32, function(x){kruskal.test(agg_estRela[, x] ~ clusComb)$p.val})


# distinct(groupDB)


# resultTEST <- readRDS(resultFilename[1])
# 
# registerDoParallel(5)
# estSUMz_allCHAIN[[1]] <- foreach(i = 1:155, .combine = rbind) %dopar% {
#   ci <- resultTEST$mod$ci_result[15001:25000, i] ### Cluster Assignment
#   at_risk <- t(resultTEST$mod$atrisk_result[i, , 15001:25000]) ### At-risk indicator
#   clus_conc <- t(sapply(1:10000, function(x){resultTEST$mod$beta_result[ci[x] + 1, , 15000 + x]})) ### Cluster Concentration
#   
#   set.seed(1)
#   estUNNORMprob <- sapply(1:10000, function(t){
#     sapply(1:60, function(j){ifelse(at_risk[t, j] == 0, 0, rgamma(1, exp(clus_conc[t, j]), 1))})
#   }) %>% t()
#   
#   estNORMprob <- estUNNORMprob/rowSums(estUNNORMprob)
#   
#   set.seed(1)
#   estZ <- sapply(1:10000, function(x){
#     rmultinom(1, zDEPTH[i], estNORMprob[x, ])
#   })
#   
#   rowSums(estZ)
#   
# }
# stopImplicitCluster()
# 
# estSUMz/rowSums(estSUMz)
# 
# ci <- resultTEST$mod$ci_result[15001:25000, 1] ### Cluster Assignment
# at_risk <- t(resultTEST$mod$atrisk_result[1, , 15001:25000]) ### At-risk indicator
# clus_conc <- t(sapply(1:10000, function(x){resultTEST$mod$beta_result[ci[x] + 1, , 15000 + x]}))
# 
# set.seed(1)
# estUNNORMprob <- sapply(1:10000, function(i){
#   sapply(1:60, function(j){ifelse(at_risk[i, j] == 0, 0, rgamma(1, exp(clus_conc[i, j]), 1))})
# }) %>% t()
# 
# estNORMprob <- estUNNORMprob/rowSums(estUNNORMprob)
# 
# set.seed(1)
# estZ <- sapply(1:10000, function(x){
#   rmultinom(1, zDEPTH[1], estNORMprob[x, ])
# })
# 
# rowSums(estZ)
# 
# colSums(t(estZ))

#test <- rmultinom(1, zDEPTH[1], estNORMprob[1, ])
#data.frame(est = test/sum(test), act = as.numeric(otuHIV[1, ]/sum(otuHIV[1, ])))

### Basic Demo for Each Cluster
demoHIV <- data.frame(datHIV[, -(1:60)], Cluster = paste0("Cluster ", clusComb))
table(demoHIV$Cluster, paste0(demoHIV$HIV_Status, " ", demoHIV$MSM))
table(demoHIV$Cluster, demoHIV$MSM)

salsoCombDemo <- data.frame(table(clusComb, statusHIV, sexHIV)) %>%
  mutate(status_sex = paste0(statusHIV, ": ", sexHIV), 
         Cluster = paste0("Cluster ", clusComb))
salsoCombDemo$Cluster <- factor(salsoCombDemo$Cluster, levels = paste0("Cluster ", 1:3))

##### Obtain the relative abundance
acRela <- (otuHIV/rowSums(otuHIV))

### Richness and Shannon
rsHIV <- data.frame(Cluster = paste0("Cluster ", clusComb), 
           Richness = rowSums(otuHIV > 0), 
           Shannon = apply(acRela, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})) %>%
  pivot_longer(!Cluster)

rsHIV$name <- factor(rsHIV$name, labels = c("Richness", "Shannon Diversity"))

rsHIV %>%
  ggplot(aes(y = value, x = Cluster)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

##### Obtain the phylum/family
pf <- data.frame(f = str_extract(colnames(otuHIV), "f_[:alpha:]*"),
                 g = str_extract(colnames(otuHIV), "g_[:alpha:]*"))
pf[is.na(pf)] <- "No Info"

### Obtain the aggregate relative abundance (de novo level)
agg_acRela <- matrix(NA, nrow = 155, ncol = 56)
for(i in 1:56){
  
  jIndex <- which(sapply(1:60, function(x){sum(pf[x, ] == distinct(pf)[i, ]) == 2}))
  if(length(jIndex) == 1){
    agg_acRela[, i] <- acRela[, jIndex]
  } else {
    agg_acRela[, i] <- rowSums(acRela[, jIndex])
  }
  
}

### Choose only the highest bacteria
high_combTaxa <- union(sort(colMeans(agg_acRela[which(clusComb == 1), ]), decreasing = TRUE, index.return = TRUE)$ix[1:15],
                       sort(colMeans(agg_acRela[which(clusComb == 2), ]), decreasing = TRUE, index.return = TRUE)$ix[1:15]) %>%
  union(sort(colMeans(agg_acRela[which(clusComb == 3), ]), decreasing = TRUE, index.return = TRUE)$ix[1:15])

### Final Aggregation
colTab <- distinct(pf)[high_combTaxa, ]
aggRELA <- matrix(NA, ncol = length(high_combTaxa), nrow = 155)
for(i in 1:(length(high_combTaxa))){
  aggRELA[, i] <- agg_acRela[, high_combTaxa[i]]
}

colTab <- colTab %>%
  mutate(plotLab = ifelse(f == "No Info", ifelse(g == "g_unclassified", "<DROP>", str_remove(g, "g_")), str_remove(f, "f_")),
         groupLab = ifelse(f == "No Info", ifelse(g == "g_unclassified", "<DROP>", "Genus"), "Family"),
         index = 1:20)
### Drop 14, 17, 18

colTab <- colTab[-c(14, 17, 18), ]
colTab <- colTab %>% 
  add_row(f = "Other", g = "Other", plotLab = "Other", groupLab = "Other", index = 21) %>%
  mutate(index = 1:18)
aggRELA <- aggRELA[, -c(14, 17, 18)]

### 4 and 10 are the same => Combine to 4
aggRELA[, 4] <- aggRELA[, 4] + aggRELA[, 10]
aggRELA <- aggRELA[, -10]
colTab <- colTab[-10, ]

aggRELA <- aggRELA %>%
  as.data.frame() %>%
  mutate(XL = 1 - rowSums(aggRELA))
colnames(aggRELA) <- colTab$plotLab

colorDistinct <- cbind(c("Bacteroides", "Prevotella", "Faecalibacterium", "Lachnospiraceae", "Alistipes", 
                         "Parabacteroides", "Ruminococcaceae", "Blautia", "Ruminococcus", 
                         "Succinivibrio", "Barnesiella", "Lachnospira", "Sutterella", "Alloprevotella",
                         "Coprococcus", "Victivallis", "Other"),
                       c("Genus", "Genus", "Genus", "Family", "Genus", "Genus", "Family", "Genus", 
                         "Genus", "Genus", "Genus", "Genus", "Genus", "Genus", "Genus", 
                         "Genus", "Other"),
                       c("#C15252", "#39C161", "#FFD700", "#FFFACD", "#F9E79F", "#87CEEB", 
                         "#ADD8E6", "#ADD8E6", "#4682B4", "#00FFFF", "#6495ED", "#CCCCFF", "#8A2BE2", 
                         "#E6E6FA", "#9966CC", "#6A0DAD", "grey90"))

colnames(colorDistinct) <- c("name", "fg", "col")
colorDistinct <- data.frame(colorDistinct)
#colorDistinct$fg <- factor(colorDistinct$fg, levels = c("Family", "Genus", "Other"), labels = c("red", "black", "grey"))

aggRELA_long <- data.frame(aggRELA, ID = rownames(otuHIV)) %>%
  inner_join(data.frame(ID = rownames(otuHIV), Cluster = paste0("Cluster ", clusComb))) %>%
  # mutate(ID = singhEDD$SampleID, Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!c(ID, Cluster))

aggRELA_long$name <- factor(aggRELA_long$name, levels = colorDistinct[, 1])
finalAGG <- inner_join(aggRELA_long, colorDistinct)
finalAGG$name <- factor(finalAGG$name, levels = colorDistinct[, 1])

finalAGG %>%
  ggplot(aes(x = ID, y = value, fill = name, color = fg)) +
  geom_bar(stat = "identity", linewidth = 0.5, alpha = 0.75) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 17.5),
        axis.text.y = element_text(size = 15), legend.title = element_text(size = 20)) +
  scale_color_manual(values = c("red", "black", "grey")) +
  scale_fill_manual(values = colorDistinct[, 3]) +
  guides(color = guide_legend(order = 1), 
         fill = guide_legend(nrow = 2, order = 2)) +
  labs(y = "Relative Abundances", x = "Participants", fill = "Taxon", color = "Order") +
  facet_grid(. ~ Cluster, scales = "free")

### EDD Dataset - de novo level
## data("singhEDD")
otuTab <- readRDS("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/clean_singh_species.rds")
metaData <- read.delim("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/edd_singh.metadata.txt")

# otuTab <- readRDS("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/Data/Application Data/singh/clean_singh_species.rds")
# metaData <- read.delim("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/Data/Application Data/singh/edd_singh.metadata.txt")

### Import the cluster result
path <- "/Users/kevinkvp/Desktop/Github Repo/Manuscript/"
# path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/"

### Filename
filename <- paste0("Result/Result/result_cleaned_singh_species_chain_", 1:2, "_init_", 
                   c("oneClus", "oneClus", "3clus", "3clus", "5clus", "5clus", "20clus", "20clus"), 
                   "_s2_1_s2MH_1en5_FIXED.rds")

#### Obtain the cluster assignment via salso packages
lastMCMC <- vector("list", 8)

for(i in c(1, 2, 3, 5, 7)){
  
  if(file.exists(paste0(path, filename[i]))){
    set.seed(1)
    result <- readRDS(paste0(path, filename[i]))
    lastMCMC[[i]] <- result$mod$ci_result[40001:50000, ]
    rm(result)
  }
  
}

combClus <- lapply(1:8, function(x){
  
  if(is.matrix(lastMCMC[[x]])){
    data.frame(lastMCMC[[x]])
  } else {
    NULL
  }
  
}) %>%
  bind_rows() %>%
  as.matrix() %>%
  salso() %>%
  as.numeric()

data.frame(SampleID = rownames(otuTab), combClus) %>%
  inner_join(metaData) %>%
  group_by(combClus, Status) %>%
  summarise(n = n())

sumTab <- data.frame(SampleID = rownames(otuTab), clus = combClus) %>%
  inner_join(metaData)
table(sumTab$clus, sumTab$Status)

### Obtain the relative abundance (de novo level)
acRela <- (otuTab/rowSums(otuTab))

sumTab$SampleID == rownames(acRela)

rsDAT <- data.frame(Cluster = paste0("Cluster ", sumTab$clus), 
                    Richness = rowSums(otuTab > 0), 
                    Shannon = apply(acRela, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})) %>%
  pivot_longer(!Cluster)

rsDAT$name <- factor(rsDAT$name, labels = c("Richness", "Shannon Diversity"))

rsDAT %>%
  ggplot(aes(y = value, x = Cluster)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

### Obtain the phylum/family
pf <- data.frame(p = str_extract(colnames(otuTab), "p__[:alpha:]*"),
                 f = str_extract(colnames(otuTab), "f__[:alpha:]*"))

### Obtain the aggregate relative abundance (de novo level)
agg_acRela <- matrix(NA, nrow = 303, ncol = 26)
for(i in 1:26){
  
  jIndex <- which(sapply(1:213, function(x){sum(pf[x, ] == distinct(pf)[i, ]) == 2}))
  if(length(jIndex) == 1){
    agg_acRela[, i] <- acRela[, jIndex]
  } else {
    agg_acRela[, i] <- rowSums(acRela[, jIndex])
  }
  
}

### Choose only the highest bacteria
high_combTaxa <- union(sort(colMeans(agg_acRela[which(combClus == 1), ]), decreasing = TRUE, index.return = TRUE)$ix[1:13],
                       sort(colMeans(agg_acRela[which(combClus == 2), ]), decreasing = TRUE, index.return = TRUE)$ix[1:13]) %>%
  union(sort(colMeans(agg_acRela[which(combClus == 3), ]), decreasing = TRUE, index.return = TRUE)$ix[1:13])

### Final Aggregation
colTab <- distinct(pf)[high_combTaxa, ]
aggRELA <- matrix(NA, ncol = length(high_combTaxa) + 1, nrow = 303)
colTab <- colTab %>%
  mutate(famName = str_extract(colTab$f, "[:upper:][:alpha:]*")) %>%
  add_row(p = "Other", f = "Other", famName = "Other")

for(i in 1:(length(high_combTaxa))){
  aggRELA[, i] <- agg_acRela[, high_combTaxa[i]]
}
aggRELA[, (length(high_combTaxa) + 1)] <- 1 - rowSums(aggRELA[, -(length(high_combTaxa) + 1)])
colnames(aggRELA) <- colTab$famName

colorDistinct <- rbind(c("Lachnospiraceae", "#c99930", "Firmicutes"),
                       c("Ruminococcaceae", "#f9bb00", "Firmicutes"),
                       c("Veillonellaceae", "#ffe184", "Firmicutes"),
                       c("Streptococcaceae", "#bf9000", "Firmicutes"),
                       c("Acidaminococcaceae", "#f1c232", "Firmicutes"),
                       c("Erysipelotrichaceae", "#ffe599", "Firmicutes"),
                       c("Clostridiaceae", "#fff2cc", "Firmicutes"),
                       c("Enterococcaceae", "#ffd966", "Firmicutes"),
                       c("Bacteroidaceae", "#39C161", "Bacteroidetes"),
                       c("Porphyromonadaceae", "#2B9F4D", "Bacteroidetes"),
                       c("Rikenellaceae", "#75D792", "Bacteroidetes"),
                       c("Prevotellaceae", "#54D97B", "Bacteroidetes"),
                       c("Enterobacteriaceae", "#E23232", "Proteobacteria"),
                       c("Sutterellaceae", "#C15252", "Proteobacteria"),
                       c("Pasteurellaceae", "#A41D00", "Proteobacteria"),
                       c("Coriobacteriaceae", "#565DFF", "Actinobacteria"),
                       c("Other","grey90" , "Other"))
                    
aggRELA_long <- aggRELA %>% 
  as.data.frame() %>%
  mutate(ID = rownames(otuTab)) %>% 
  inner_join(data.frame(ID = rownames(otuTab), Cluster = paste0("Cluster ", combClus))) %>%
  # mutate(ID = singhEDD$SampleID, Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!c(ID, Cluster))

aggRELA_long$name <- factor(aggRELA_long$name, levels = colorDistinct[, 1])
aggRELA_long %>%
  ggplot(aes(x = ID, y = value, fill = name)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 17.5),
        axis.text.y = element_text(size = 15), legend.title = element_text(size = 20)) +
  scale_fill_manual(values = colorDistinct[, 2]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(y = "Relative Abundances", x = "Participants", fill = "Family") +
  facet_grid(. ~ Cluster, scales = "free")

### EDD Dataset - grouped genus level
# otuTab <- readRDS("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/clean_singh.rds")
# metaData <- read.delim("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/edd_singh.metadata.txt")
# path <- "/Users/kevinkvp/Desktop/Github Repo/Manuscript/"

otuTab <- readRDS("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/Data/Application Data/singh/clean_singh.rds")
metaData <- read.delim("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/Data/Application Data/singh/edd_singh.metadata.txt")
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI - Manuscript Code - Draft/Manuscript/"

resultFilename <- c(paste0(path, "Result/singh/result_cleaned_singh_chain_", 1:2, "_init_oneClus_s2_1en1_s2MH_1en3.rds"),
                    paste0(path, "Result/singh/result_cleaned_singh_chain_", 1, "_init_3clus_s2_1en1_s2MH_1en3.rds"),
                    paste0(path, "Result/singh/result_cleaned_singh_chain_", 1, "_init_5clus_s2_1en1_s2MH_1en3.rds"),
                    paste0(path, "Result/singh/result_cleaned_singh_chain_", 1, "_init_20clus_s2_1en1_s2MH_1en3.rds"))

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
combineClus <- foreach(t = 1:5, .combine = "rbind") %dopar% {
  result <- readRDS(resultFilename[t])
  result$mod$ci_result[15001:25000, ]
}
stopImplicitCluster()

set.seed(1)
combSalso <- as.numeric(salso(combineClus))
table(combSalso)
sumTab <- data.frame(SampleID = rownames(otuTab), clus = combSalso) %>%
  inner_join(metaData)
table(sumTab$clus, sumTab$Status)

### Obtain the relative abundance (de novo level)
acRela <- (otuTab/rowSums(otuTab))

### Richness and Shannon
# sumTab$SampleID == rownames(acRela)

rsDAT <- data.frame(Cluster = paste0("Cluster ", sumTab$clus), 
                    Richness = rowSums(otuTab > 0), 
                    Shannon = apply(acRela, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})) %>%
  pivot_longer(!Cluster)

rsDAT$name <- factor(rsDAT$name, labels = c("Richness", "Shannon Diversity"))

rsDAT %>%
  ggplot(aes(y = value, x = Cluster)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

### Obtain the phylum/family
pf <- data.frame(p = str_extract(colnames(otuTab), "p__[:alpha:]*"),
                 f = str_extract(colnames(otuTab), "f__[:alpha:]*"))

### Obtain the aggregate relative abundance (de novo level)
agg_acRela <- matrix(NA, nrow = 303, ncol = 32)
for(i in 1:32){
  
  jIndex <- which(sapply(1:213, function(x){sum(pf[x, ] == distinct(pf)[i, ]) == 2}))
  if(length(jIndex) == 1){
    agg_acRela[, i] <- acRela[, jIndex]
  } else {
    agg_acRela[, i] <- rowSums(acRela[, jIndex])
  }
  
}

### Choose only the highest bacteria
high_combTaxa <- union(sort(colMeans(agg_acRela[which(combSalso == 1), ]), decreasing = TRUE, index.return = TRUE)$ix[1:10],
                       sort(colMeans(agg_acRela[which(combSalso == 2), ]), decreasing = TRUE, index.return = TRUE)$ix[1:10])

### Final Aggregation
colTab <- distinct(pf)[high_combTaxa, ]
aggRELA <- matrix(NA, ncol = length(high_combTaxa) + 1, nrow = 303)
colTab <- colTab %>%
  mutate(famName = str_extract(colTab$f, "[:upper:][:alpha:]*")) %>%
  add_row(p = "Other", f = "Other", famName = "Other")

for(i in 1:(length(high_combTaxa))){
  aggRELA[, i] <- agg_acRela[, high_combTaxa[i]]
}
aggRELA[, (length(high_combTaxa) + 1)] <- 1 - rowSums(aggRELA[, -(length(high_combTaxa) + 1)])
colnames(aggRELA) <- colTab$famName

colorDistinct <- rbind(c("Lachnospiraceae", "#c99930", "Firmicutes"),
                       c("Ruminococcaceae", "#f9bb00", "Firmicutes"),
                       c("Veillonellaceae", "#ffe184", "Firmicutes"),
                       c("Streptococcaceae", "#bf9000", "Firmicutes"),
                       c("Acidaminococcaceae", "#f1c232", "Firmicutes"),
                       c("Bacteroidaceae", "#39C161", "Bacteroidetes"),
                       c("Porphyromonadaceae", "#2B9F4D", "Bacteroidetes"),
                       c("Rikenellaceae", "#75D792", "Bacteroidetes"),
                       c("Prevotellaceae", "#54D97B", "Bacteroidetes"),
                       c("Enterobacteriaceae", "#E23232", "Proteobacteria"),
                       c("Sutterellaceae", "#C15252", "Proteobacteria"),
                       c("Pseudomonadaceae", "#A41D00", "Proteobacteria"),
                       c("Fusobacteriaceae", "#751197", "Fusobacteria"),
                       c("Other","grey90" , "Other"))

set.seed(1)
combSalso <- as.numeric(salso(combineClus))
table(combSalso)
sumTab <- data.frame(SampleID = rownames(otuTab), clus = combSalso) %>%
  inner_join(metaData)
table(sumTab$clus, sumTab$Status)

aggRELA_long <- aggRELA %>% 
  as.data.frame() %>%
  mutate(ID = rownames(otuTab)) %>% 
  inner_join(data.frame(ID = rownames(otuTab), Cluster = paste0("Cluster ", combSalso))) %>%
  # mutate(ID = singhEDD$SampleID, Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!c(ID, Cluster))

aggRELA_long$name <- factor(aggRELA_long$name, levels = colorDistinct[, 1])
aggRELA_long %>%
  ggplot(aes(x = ID, y = value, fill = name)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 20),
        axis.text.y = element_text(size = 15), legend.title = element_text(size = 20)) +
  scale_fill_manual(values = colorDistinct[, 2]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(y = "Relative Abundances", x = "Participants", fill = "Family") +
  facet_grid(. ~ Cluster, scales = "free")


# registerDoParallel(5)
# activeClusMat <- foreach(t = 1:5, .combine = cbind) %dopar% {
#   result <- readRDS(resultFilename[t])
#   apply(result$mod$ci_result, 1, uniqueClus)
# }
# stopImplicitCluster()
# 
# activeClusMatPlot <- activeClusMat %>%
#   as.data.frame() %>%
#   mutate(iter = 1:25000) %>%
#   pivot_longer(!iter)
# 
# activeClusMatPlot$name <- factor(activeClusMatPlot$name,
#                                  levels = paste0("result.", 1:12), labels = paste0("Chain ", 1:12))
# 
# ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
#   geom_line() +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   labs(title = "Singh (Genus Level) - Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusters",
#        color = "MCMC Chain") +
#   guides(color = guide_legend(ncol = 5)) +
#   scale_y_continuous(breaks = seq(0, 12, 1))





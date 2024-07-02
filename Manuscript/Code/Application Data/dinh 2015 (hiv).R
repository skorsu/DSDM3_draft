# devtools::install_github("YushuShi/MicrobiomeCluster")

library(ape)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(salso)
library(ecodist)
library(foreach)
library(doParallel)
library(cluster)
library(ggplot2)
library(ecodist)
library(coda.base)
library(mclustcomp)
library(vcd)
library(gridExtra)

### User-defined functions: ----------------------------------------------------
uniqueClus <- function(x){
  length(unique(x))
}

### Global Objects: ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

### Dinh 2015, HIV: ------------------------------------------------------------
### Import Data and Cleaning OTU: ----------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/hiv_dinh_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/hiv_dinh_results.tar"))

#### Metadata
metData <- read.table(paste0(path, "hiv_dinh_results/hiv_dinh.metadata.txt"), sep = "\t", header = TRUE)
metData[metData == "<not provided>"] <- NA
dim(metData)

#### OTU Table
otuTab <- read.table(paste0(path, "hiv_dinh_results/RDP/hiv_dinh.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
dim(otuTab)

##### Check the samples: All samples in OTU has the metadata.
str_extract(rownames(otuTab), "[^X]+") %in% metData$Sample_Name_s
metData$Sample_Name_s %in% str_extract(rownames(otuTab), "[^X]+") 
rownames(otuTab) <- str_extract(rownames(otuTab), "[^X]+")

##### Filter OTU table
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)

### EDA: -----------------------------------------------------------------------
# View(metData)

metNew <- metData %>%
  transmute(ID = Sample_Name_s,
            Disease = disease_s,
            Gender = ifelse(sex_s == "Female", "Female", ifelse(is.na(is_gay_s), "Male: Hetorosexual", "Male: Homosexual")),
            Education = ifelse(Education_s == "college", "College", ifelse(Education_s == "high_school", "Up to High School", "No High School")),
            Drinker = ifelse(Alcohol_s == "Non_drinker", "Non-Drinker", ifelse(Bingeing_s == "binge_drinker", "Drinker: Binge", "Drinker: Non-Binge")),
            Smoking = ifelse(Current_smoking_s == "Not_currently_smoking", "No", "Yes"), 
            BMI = BMI_s, Age = age_s, Cholesterol = as.numeric(lvCholesterol_s), 
            IFNgamma = as.numeric(IFNgamma_s), IL1Beta = as.numeric(IL1Beta_s), 
            IL6 = as.numeric(IL6_s), Kcal = as.numeric(Kcal_s), TNFalpha = as.numeric(TNFalpha_s))

metNew$Education <- factor(metNew$Education, levels = c("College", "Up to High School", "No High School"))

# View(metNew)

metNew %>%
  group_by(Disease) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Disease, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Disease Status")

metNew %>%
  group_by(Gender) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Gender, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Gender")

metNew %>%
  group_by(Education) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Education, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Education Level")

metNew %>%
  group_by(Drinker) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Drinker, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Alcohol Drinker")

metNew %>%
  group_by(Smoking) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Smoking, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Smoking Behavior")

metNew %>%
  dplyr::select(ID, BMI, Age, Cholesterol, IFNgamma, IL1Beta, IL6, Kcal, TNFalpha) %>%
  pivot_longer(!ID) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(. ~ name, scales = "free_y")

### Taxa for visualization:-----------------------------------------------------
colTaxa <- data.frame(Kingdom = str_replace(str_extract(colnames(otuTab), "k\\_\\_[:alpha:]+"), "k\\_\\_", ""),
                      Phylum = str_replace(str_extract(colnames(otuTab), "p\\_\\_[:alpha:]+"), "p\\_\\_", ""),
                      Class = str_replace(str_extract(colnames(otuTab), "c\\_\\_[:alpha:]+"), "c\\_\\_", ""),
                      Order = str_replace(str_extract(colnames(otuTab), "o\\_\\_[:alpha:]+"), "o\\_\\_", ""),
                      Family = str_replace(str_extract(colnames(otuTab), "f\\_\\_[:alpha:]+"), "f\\_\\_", ""),
                      Genus = str_replace(str_extract(colnames(otuTab), "g\\_\\_[:alpha:]+"), "g\\_\\_", ""),
                      Species = str_replace(str_extract(colnames(otuTab), "s\\_\\_[:alpha:]+"), "s\\_\\_", ""))

colTaxa[which(is.na(colTaxa$Class)), "Class"] <- "Other Firmicutes"
# table(colTaxa$Phylum)
# table(colTaxa$Class)
# table(colTaxa$Family)
# table(colTaxa$Phylum, colTaxa$Family) %>% sum()

ClassName <- names(table(colTaxa$Class))
ClassOTU <- sapply(1:9, function(x){rowSums(otuTab[, which(colTaxa$Class == ClassName[x])])}) %>%
  `colnames<-`(ClassName)
  
ClassOTU <- ClassOTU/rowSums(ClassOTU)
# View(ClassOTU)

### Visualization
ClassOTUlonger <- ClassOTU %>%
  as.data.frame() %>%
  mutate(ID = rownames(ClassOTU)) %>%
  pivot_longer(!ID, names_to = "Class")

ClassRank <- c("Bacteroidia", "Bacilli", "Clostridia", "Erysipelotrichia", "Negativicutes", "Other Firmicutes",
               "Fusobacteriia", "Betaproteobacteria", "Gammaproteobacteria")
LabelRank <- c("Bacteroidetes - Bacteroidia", "Firmicutes - Bacilli", 
               "Firmicutes - Clostridia", "Firmicutes - Erysipelotrichia", 
               "Firmicutes - Negativicutes", "Firmicutes - Other Firmicutes",
               "Fusobacteria - Fusobacteriia", "Proteobacteria - Betaproteobacteria", 
               "Proteobacteria - Gammaproteobacteria")
ColorPhyla <- c("#FF9999", "#99FF99", "#66FF66", "#33FF33", "#00FF00", "#00CC00",
                "#81EFFF", "#CC99FF", "#9966FF")

ClassOTUlonger$Class <- factor(ClassOTUlonger$Class, levels = ClassRank, labels = LabelRank)

ggplot(ClassOTUlonger, aes(x = ID, y = value, fill = Class)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = ColorPhyla) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  labs(y = "Relative Abundance")

### PAM: -----------------------------------------------------------------------
#### Choosing the number of clusters for each distance metric
data.frame(k = 2:10,
           Euclidean = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(otuTab), y))[, 3])}),
           BrayCurtis = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(bcdist(otuTab), y))[, 3])}),
           Aitchison = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(otuTab + 1e-20, method = "aitchison"), y))[, 3])})) %>%
  pivot_longer(cols = !k) %>%
  ggplot(aes(x = k, y = value)) +
  geom_line() +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_bw() +
  labs(y = "Average Silhouette", x = "Number of clusters")

##### PAM - Cluster Assignment with the metadata
set.seed(1)
metPAM <- data.frame(ID = names(pam(bcdist(otuTab), 7)$clustering),
                     clusEu = pam(dist(otuTab), 7)$clustering, 
                     clusBC = pam(bcdist(otuTab), 7)$clustering,
                     clusAT = pam(dist(otuTab + 1e-20, method = "aitchison"), 2)$clustering) %>%
  `rownames<-`(NULL) %>%
  inner_join(metNew)

transmute(metPAM, ID, Cluster = paste0("Cluster ", clusEu)) %>%
  inner_join(ClassOTUlonger) %>%
  ggplot(aes(x = ID, y = value, fill = Class)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = ColorPhyla) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  labs(y = "Relative Abundance", title = "Relative Abundance - Cluster by PAM with Euclidean distance") +
  facet_grid(. ~ Cluster, scales = "free_x")

transmute(metPAM, ID, Cluster = paste0("Cluster ", clusBC)) %>%
  inner_join(ClassOTUlonger) %>%
  ggplot(aes(x = ID, y = value, fill = Class)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = ColorPhyla) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  labs(y = "Relative Abundance", title = "Relative Abundance - Cluster by PAM with Bray-Curtis distance") +
  facet_grid(. ~ Cluster, scales = "free_x")

transmute(metPAM, ID, Cluster = paste0("Cluster ", clusAT)) %>%
  inner_join(ClassOTUlonger) %>%
  ggplot(aes(x = ID, y = value, fill = Class)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = ColorPhyla) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  labs(y = "Relative Abundance", title = "Relative Abundance - Cluster by PAM with Aitchison distance") +
  facet_grid(. ~ Cluster, scales = "free_x")

##### Metadata - Euclidean -----------------------------------------------------
###### Conclusion - Only Cluster 2 and 3 can be interpret as an HIV groups. Provide 4 singletons
structable(Disease ~ clusEu, metPAM)
structable(Gender ~ clusEu, metPAM)
structable(Education ~ clusEu, metPAM)
structable(Drinker ~ clusEu, metPAM)
structable(Smoking ~ clusEu, metPAM)

#### Cluster 2 - Male, Education HS or higher, Smoking
#### Cluster 3 - Male, Education HS or higher, Drinker

##### Metadata - Bray-Curtis ---------------------------------------------------
###### Conclusion - Cluster 4 and 5 can be interpret as an HIV groups. Provide 2 singletons
structable(Disease ~ clusBC, metPAM)
structable(Gender ~ clusBC, metPAM)
structable(Education ~ clusBC, metPAM)
structable(Drinker ~ clusBC, metPAM)
structable(Smoking ~ clusBC, metPAM)

#### Cluster 4 - Education HS or higher
#### Cluster 5 - Male, Education HS or higher, Drinker

### ZIDM-ZIDM: -----------------------------------------------------------------
#### Import the result form our model
result <- readRDS(paste0(path, "Manuscript/Result/Application Data/Dinh/hiv_dinh_result_split_10_theta_1_rn_ADAP500.rds"))

#### Run for 5000 iterations, thinning every 5th iterations
#### AMH - Apply the adpative proposal after the 500th iteration

#### Computational Time
sapply(1:6, function(x){as.numeric(result[[x]]$time)})/3600

#### Convergence Checking: -----------------------------------------------------
##### Active Cluster in MCMC iterations
actClus <- sapply(1:6, function(x){apply(result[[x]]$mod$ci_result, 1, uniqueClus)}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:6)) %>%
  mutate(Iteration = 1:1000) %>%
  pivot_longer(!Iteration)

actClus$name <- factor(actClus$name, labels = c(paste0("One Cluster, Kmax = 10 - Chain ", 1:2),
                                                paste0("Three Clusters, Kmax = 20 - Chain ", 1:2),
                                                paste0("Singletons, Kmax = 36 - Chain ", 1:2)))

ggplot(actClus, aes(x = Iteration, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Thinned Iteration", y = "# Active Clusters",
       title = "Dinh: Number of active clusters via ZIDM-ZIDM",
       color = "Initialization") +
  theme(legend.position = "bottom")

##### Acceptance Rate - Split-Merge Step
sapply(1:6, function(x){mean(result[[x]]$mod$sm_accept)})

##### Acceptance Rate - Cluster Concentration
Kmax <- c(10, 10, 20, 20, 36, 36)
# lapply(1:6, function(y){sapply(1:Kmax[y], function(x){c(sum(result[[y]]$mod$MH_accept[, x] != -1),
#                                                    sum(result[[y]]$mod$MH_accept[, x] == 1)/sum(result[[y]]$mod$MH_accept[, x] != -1),
#                                                    sum(result[[y]]$mod$MH_accept[(1:500), x] == 1)/sum(result[[y]]$mod$MH_accept[(1:500), x] != -1),
#                                                    sum(result[[y]]$mod$MH_accept[-(1:500), x] == 1)/sum(result[[y]]$mod$MH_accept[-(1:500), x] != -1))}) %>%
#     t()})

#### The average acceptance rate for the cluster concentration for every chain
lapply(1:6, function(y){sapply(1:Kmax[y], function(x){sum(result[[y]]$mod$MH_accept[, x] == 1)/sum(result[[y]]$mod$MH_accept[, x] != -1)})}) %>%
  lapply(mean)

plot(result[[1]]$mod$beta_result[1, 1, ], type = "l")


#### Cluster Assignment: -------------------------------------------------------
set.seed(1)
clusZZ <- lapply(1:6, function(x){as.data.frame(result[[x]]$mod$ci_result[seq(250, 1000, 10), ])}) %>%
  bind_rows() %>%
  as.matrix() %>%
  salso(loss = "binder")

##### Estimated Relative Abundance
estProb <- lapply(1:6, function(y){
  
  sapply(1:36, function(x){
    ci <- result[[y]]$mod$ci_result[seq(250, 1000, 10), x]
    dumProb <- matrix(NA, nrow = 76, ncol = 1024)
    for(i in 1:76){
      estClusConc <- exp(result[[y]]$mod$beta_result[ci[i] + 1, , 1250 + ((i - 1) * 50)])
      atrisk <- result[[y]]$mod$atrisk_result[x, , 250 + ((i - 1) * 10)]
      dumProb[i, ] <- (estClusConc * atrisk)/sum(estClusConc * atrisk)
    }
    
    colSums(dumProb)/sum(dumProb)
    
  }) %>% t() %>% as.data.frame()
  
  
})

estProb <- (estProb[[1]] + estProb[[2]] + estProb[[3]] + estProb[[4]] + estProb[[5]] + estProb[[6]])/sum(estProb[[1]] + estProb[[2]] + estProb[[3]] + estProb[[4]] + estProb[[5]] + estProb[[6]])
colnames(estProb) <- colnames(otuTab)

ClassEstOTU <- sapply(1:9, function(x){rowSums(estProb[, which(colTaxa$Class == ClassName[x])])}) %>%
  `colnames<-`(ClassName)

ClassEstOTU <- ClassEstOTU/rowSums(ClassEstOTU)

ClassEstlonger <- ClassEstOTU %>%
  as.data.frame() %>%
  mutate(ID = rownames(ClassOTU)) %>%
  pivot_longer(!ID, names_to = "Class")

ClassEstlonger$Class <- factor(ClassEstlonger$Class, levels = ClassRank, labels = LabelRank)

data.frame(ID = rownames(otuTab), clusterZZ = paste0("Cluster ", clusZZ)) %>%
  inner_join(ClassOTUlonger) %>%
  ggplot(aes(x = ID, y = value, fill = Class)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = ColorPhyla) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  labs(y = "Relative Abundance", title = "Relative Abundance - Cluster by PAM with ZIDM-ZIDM") +
  facet_grid(. ~ clusterZZ, scales = "free_x")

data.frame(ID = rownames(otuTab), clusterZZ = paste0("Cluster ", clusZZ)) %>%
  inner_join(ClassEstlonger) %>%
  ggplot(aes(x = ID, y = value, fill = Class)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = ColorPhyla) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  labs(y = "Estimated Relative Abundance", title = "Estimated Relative Abundance - Cluster by PAM with ZIDM-ZIDM") +
  facet_grid(. ~ clusterZZ, scales = "free_x")

##### ZIDM-ZIDM: Profile Analysis ----------------------------------------------
metZZ <- data.frame(ID = rownames(otuTab), clusZZ) %>%
  inner_join(metNew)

###### Conclusion - 9 Clusters in total (No singleton)
###### Cluster 1, 8 - non-HIV while Cluster 5, 6, 9 are HIV.
structable(Disease ~ clusZZ, metZZ)
structable(Gender ~ clusZZ, metZZ)
structable(Education ~ clusZZ, metZZ)
structable(Drinker ~ clusZZ, metZZ)
structable(Smoking ~ clusZZ, metZZ)

##### Cluster 1 - Not Smoking
##### Cluster 5 - Male, Education HS or higher, smoking
##### Cluster 6 - Male, College, Drinker
##### Cluster 8 - Education HS or higher, smoking
##### Cluster 9 - Drinker, smoking

















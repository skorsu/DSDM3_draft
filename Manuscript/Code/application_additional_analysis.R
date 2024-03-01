### Required Library
library(readxl)
library(stringr)
library(tidyverse)

library(salso)
library(foreach)
library(doParallel)
library(xtable)

library(mclustcomp)
library(ggplot2)
library(gridExtra)

library(reshape2)
library(salso)
library(gridExtra)
library(sparseMbClust)
library(cluster)
library(ecodist)
library(factoextra)
library(rbiom)
library(ggcorrplot)

library(pheatmap)
library(gridExtra)

### Import the data and result -------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Additional Data
addml <- read_excel(paste0(path, "Data/Application/Mali_RB Metadata.xlsx"))
addni <- read_excel(paste0(path, "Data/Application/Nicaragua_RB Metadata.xlsx"))

### Result
annikaZZ <- readRDS(paste0(path, "Result/microbiome_result.RData"))

### Data: 6 and 8 Months
ni68 <- read.csv(paste0(path, "Data/Application/Nicaragua_6mo_8mo_genus.csv"))
ml68 <- read.csv(paste0(path, "Data/Application/Mali_6mo_8mo_genus.csv"))

### Data: 12 Months
ni12 <- read.csv(paste0(path, "Data/Application/Nicaragua_12mo_Metadata_csv.csv"))
ml12 <- read.csv(paste0(path, "Data/Application/Mali_12mo_Metadata_csv.csv"))

## Retrived the infant ID in the analysis --------------------------------------
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

obsID <- c(str_extract(dat06[1:45, 2], "[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"), 
           str_extract(dat06[-(1:45), 2], "^[:digit:]+"))

## Cleaning the additional data ------------------------------------------------
#### (!) Breast Feeding is all yes.
addml %>%
  filter(Child_ID %in% obsID) %>%
  arrange(Child_ID) %>%
  filter(`Monthly_Visit#` %in% c("6MO", "8MO", "12MO")) %>%
  dplyr::select(Child_ID, `Monthly_Visit#`)

addml$Vegetables

addni %>%
  filter(Child_ID %in% gsub("\\.", "-", obsID)) %>%
  arrange(Child_ID) %>%
  filter(`Monthly_Visit#` %in% c("6 MO", "8 MO", "12 MO")) %>%
  dplyr::select(Child_ID, `Monthly_Visit#`, Cow_milk_5.2) %>%
  mutate(cowMilk = ifelse(Cow_milk_5.2 == "NEVER", "NO", "YES")) %>%
  dplyr::select(-Cow_milk_5.2) %>%
  pivot_wider(names_from = `Monthly_Visit#`, values_from = cowMilk)



### Analyze --------------------------------------------------------------------
#### Computational Time 
c(annikaZZ[[1]]$time, annikaZZ[[2]]$time, annikaZZ[[3]]$time)/60

### Cluster Assignment
clusBinder <- sapply(1:3, 
                     function(x){as.numeric(salso(annikaZZ[[x]]$result$ci_result[-(1:500), ],
                                                  loss = "binder"))})

obsID <- c(str_extract(dat06[1:45, 2], "[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"), 
           str_extract(dat06[-(1:45), 2], "^[:digit:]+"))

#### First, define the bacteria that will show in the plot, while the rest will represents as Other ...
datList <- list(dat06, dat08, dat12)

registerDoParallel(5)
impTaxa <- foreach(t = 1:3) %dopar% {
  
  prop <- as.numeric(colMeans(datList[[t]][, -(1:5)]/rowSums(datList[[t]][, -(1:5)])))
  which(prop %in% boxplot.stats(prop)$out)
  
}
stopImplicitCluster()
impTaxa <- union(impTaxa[[1]], union(impTaxa[[2]], impTaxa[[3]]))

#### str_match_all(colnames(dat06[-(1:5)]), "[:alpha:]{1}\\_{2}\\.?[:alpha:]+")
#### Go for the phylum and genus
taxaName <- cbind(str_match(colnames(dat06[-(1:5)]), "p\\_{2}\\.?[:alpha:]+") %>%
                    str_match("[:upper:]{1}[:alpha:]+"), 
                  str_match(colnames(dat06[-(1:5)]), "g\\_{2}\\.?[:alpha:]+") %>%
                    str_match("[:upper:]{1}[:alpha:]+"))

impTaxaInd <- rep(FALSE, 38)
impTaxaInd[impTaxa] <- TRUE

#### Get the column name, group the other taxa
showName <- data.frame(taxaName, impTaxaInd) %>%
  replace_na(list(X1 = "Bacteria", X2 = "Other")) %>%
  mutate(showName = ifelse(impTaxaInd, X2, 
                           ifelse(X1 == "Bacteria", "Others", paste0("Other ", X1)))) %>%
  .$showName

### Preprocess: Reindex --------------------------------------------------------
#### Calculate the relative Taxa count
##### 6-Month Data
registerDoParallel(5)
dat06rela <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(showName, ID = rep(obsID[t], 38), 
             clus = rep(paste0("Cluster ", clusBinder[t, 1]), 38),
             taxa = as.numeric(dat06[t, -(1:5)])) %>%
    group_by(ID, clus, showName) %>%
    summarise(taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(relaTaxa = taxa/sum(taxa)) %>%
    dplyr::select(-taxa)
  
}
stopImplicitCluster()
dat06rela$showName <- factor(dat06rela$showName,
                             levels = c("Bifidobacterium", "Other Actinobacteria",
                                        "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                        "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                        "Streptococcus", "Veillonella", "Other Firmicutes", 
                                        "Others"))

colours <- c("#A4C639", "#BAC394", "#DDA661", "#BD8B46", "#FFD197",
             "#BC6D62", "#DE6A5D", "#DE998B", "#FF6558", "#FFC6B8", 
             "#FFE4DA", "#E5E5E5")

dat06rela %>%
  as.data.frame() %>%
  ggplot(aes(x = ID, y = relaTaxa, fill = factor(showName))) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +  # Vertical x-axis tick labels
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance", title = "6-Month Infants", x = "Sample") +
  facet_grid(~clus, scales = "free_x")

##### 8-Month Data
registerDoParallel(5)
dat08rela <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(showName, ID = rep(obsID[t], 38), 
             clus = rep(paste0("Cluster ", clusBinder[t, 2]), 38),
             taxa = as.numeric(dat08[t, -(1:5)])) %>%
    group_by(ID, clus, showName) %>%
    summarise(taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(relaTaxa = taxa/sum(taxa)) %>%
    dplyr::select(-taxa)
  
}
stopImplicitCluster()
dat08rela$showName <- factor(dat08rela$showName,
                             levels = c("Bifidobacterium", "Other Actinobacteria",
                                        "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                        "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                        "Streptococcus", "Veillonella", "Other Firmicutes", 
                                        "Others"))

dat08rela %>%
  as.data.frame() %>%
  ggplot(aes(x = ID, y = relaTaxa, fill = factor(showName))) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +  # Vertical x-axis tick labels
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +  # Vertical x-axis tick labels
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance", title = "8-Month Infants", x = "Sample") +
  facet_grid(~clus, scales = "free_x")

##### 12-Month Data
registerDoParallel(5)
dat12rela <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(showName, ID = rep(obsID[t], 38), 
             clus = rep(paste0("Cluster ", clusBinder[t, 3]), 38),
             taxa = as.numeric(dat12[t, -(1:5)])) %>%
    group_by(ID, clus, showName) %>%
    summarise(taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(relaTaxa = taxa/sum(taxa)) %>%
    dplyr::select(-taxa)
  
}
stopImplicitCluster()
dat12rela$showName <- factor(dat12rela$showName,
                             levels = c("Bifidobacterium", "Other Actinobacteria",
                                        "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                        "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                        "Streptococcus", "Veillonella", "Other Firmicutes", 
                                        "Others"))

dat12rela %>%
  as.data.frame() %>%
  ggplot(aes(x = ID, y = relaTaxa, fill = factor(showName))) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +  # Vertical x-axis tick labels
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance", title = "12-Month Infants", x = "Sample") +
  facet_grid(~clus, scales = "free_x")

dat06rela %>%
  dplyr::select(-clus) %>%
  pivot_wider(names_from = showName, values_from = relaTaxa) %>%
  mutate(clus = clusBinder[, 1]) %>%
  group_by(clus) %>%
  summarise(meanBifidobacterium = mean(Bifidobacterium),
            meanBacteroides = mean(Bacteroides),
            meanPrevotella = mean(Prevotella))

dat08rela %>%
  dplyr::select(-clus) %>%
  pivot_wider(names_from = showName, values_from = relaTaxa) %>%
  mutate(clus = clusBinder[, 2]) %>%
  group_by(clus) %>%
  summarise(meanBifidobacterium = mean(Bifidobacterium),
            meanBacteroides = mean(Bacteroides),
            meanPrevotella = mean(Prevotella))

dat12rela %>%
  dplyr::select(-clus) %>%
  pivot_wider(names_from = showName, values_from = relaTaxa) %>%
  mutate(clus = clusBinder[, 3]) %>%
  group_by(clus) %>%
  summarise(meanBifidobacterium = mean(Bifidobacterium),
            meanBacteroides = mean(Bacteroides),
            meanPrevotella = mean(Prevotella))

#### Reindex by using the mean Bifidobacterium
clusBinderN <- matrix(NA, nrow = 90, ncol = 3)
##### 6-Month
clusBinderN[which(clusBinder[, 1] == 4), 1] <- 1
clusBinderN[which(clusBinder[, 1] == 3), 1] <- 2
clusBinderN[which(clusBinder[, 1] == 5), 1] <- 3
clusBinderN[which(clusBinder[, 1] == 1), 1] <- 4
clusBinderN[which(clusBinder[, 1] == 2), 1] <- 5
##### 8-Month
clusBinderN[which(clusBinder[, 2] == 2), 2] <- 1
clusBinderN[which(clusBinder[, 2] == 4), 2] <- 2
clusBinderN[which(clusBinder[, 2] == 3), 2] <- 3
clusBinderN[which(clusBinder[, 2] == 1), 2] <- 4
##### 12-Month
clusBinderN[which(clusBinder[, 3] == 1), 3] <- 1
clusBinderN[which(clusBinder[, 3] == 2), 3] <- 2
clusBinderN[which(clusBinder[, 3] == 4), 3] <- 3
clusBinderN[which(clusBinder[, 3] == 3), 3] <- 4

### Start analysis based on the new index cluster ------------------------------
dat06[, 5] <- factor(dat06[, 5], labels = c("Control", "Rice Bran", "Rice Bran"))
dat08[, 5] <- factor(dat08[, 5], labels = c("Control", "Rice Bran", "Rice Bran"))

table(clusBinderN[, 3], dat12[, 5])
table(dat06[, 5], dat06[, 3], clusBinderN[, 1])
table(dat08[, 5], dat08[, 3], clusBinderN[, 2])
table(dat12[, 5], dat12[, 3], clusBinderN[, 3])

#### Calculate the relative Taxa count with adjusted cluster index
##### 6-Month Data
registerDoParallel(5)
dat06Nrela <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(showName, ID = rep(obsID[t], 38), 
             clus = rep(paste0("Cluster ", clusBinderN[t, 1]), 38),
             taxa = as.numeric(dat06[t, -(1:5)])) %>%
    group_by(ID, clus, showName) %>%
    summarise(taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(relaTaxa = taxa/sum(taxa)) %>%
    dplyr::select(-taxa)
  
}
stopImplicitCluster()
dat06Nrela$showName <- factor(dat06Nrela$showName,
                              levels = c("Bifidobacterium", "Other Actinobacteria",
                                         "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                         "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                         "Streptococcus", "Veillonella", "Other Firmicutes", 
                                         "Others"))

colours <- c("#A4C639", "#BAC394", "#DDA661", "#BD8B46", "#FFD197",
             "#BC6D62", "#DE6A5D", "#DE998B", "#FF6558", "#FFC6B8", 
             "#FFE4DA", "#E5E5E5")

dat06Nrela %>%
  as.data.frame() %>%
  ggplot(aes(x = ID, y = relaTaxa, fill = factor(showName))) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +  # Vertical x-axis tick labels
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance", title = "6-Month Infants", x = "Sample") +
  facet_grid(~clus, scales = "free_x")

##### 8-Month Data
registerDoParallel(5)
dat08Nrela <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(showName, ID = rep(obsID[t], 38), 
             clus = rep(paste0("Cluster ", clusBinderN[t, 2]), 38),
             taxa = as.numeric(dat08[t, -(1:5)])) %>%
    group_by(ID, clus, showName) %>%
    summarise(taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(relaTaxa = taxa/sum(taxa)) %>%
    dplyr::select(-taxa)
  
}
stopImplicitCluster()
dat08Nrela$showName <- factor(dat08Nrela$showName,
                              levels = c("Bifidobacterium", "Other Actinobacteria",
                                         "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                         "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                         "Streptococcus", "Veillonella", "Other Firmicutes", 
                                         "Others"))

dat08Nrela %>%
  as.data.frame() %>%
  ggplot(aes(x = ID, y = relaTaxa, fill = factor(showName))) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +  # Vertical x-axis tick labels
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +  # Vertical x-axis tick labels
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance", title = "8-Month Infants", x = "Sample") +
  facet_grid(~clus, scales = "free_x")

##### 12-Month Data
registerDoParallel(5)
dat12Nrela <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(showName, ID = rep(obsID[t], 38), 
             clus = rep(paste0("Cluster ", clusBinderN[t, 3]), 38),
             taxa = as.numeric(dat12[t, -(1:5)])) %>%
    group_by(ID, clus, showName) %>%
    summarise(taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(relaTaxa = taxa/sum(taxa)) %>%
    dplyr::select(-taxa)
  
}
stopImplicitCluster()
dat12Nrela$showName <- factor(dat12Nrela$showName,
                              levels = c("Bifidobacterium", "Other Actinobacteria",
                                         "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                         "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                         "Streptococcus", "Veillonella", "Other Firmicutes", 
                                         "Others"))

dat12Nrela %>%
  as.data.frame() %>%
  ggplot(aes(x = ID, y = relaTaxa, fill = factor(showName))) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "bottom") +  # Vertical x-axis tick labels
  guides(fill = guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative abundance", title = "12-Month Infants", x = "Sample") +
  facet_grid(~clus, scales = "free_x")

### Trajectory Plot ------------------------------------------------------------
c06 <- data.frame(obs = 1:90, clus6 = factor(clusBinderN[, 1])) %>%
  arrange(clus6) %>% 
  mutate(index6 = 90:1)
c08 <- data.frame(obs = 1:90, clus8 = factor(clusBinderN[, 2])) %>%
  arrange(clus8) %>% 
  mutate(index8 = 90:1)
c12 <- data.frame(obs = 1:90, clus12 = factor(clusBinderN[, 3])) %>%
  arrange(clus12) %>% 
  mutate(index12 = 90:1)

obsChange <- NULL
for(i in 1:90){
  
  if(c06[which(c06$obs == i), "clus6"] != c08[which(c08$obs == i), "clus8"]){
    obsChange <- rbind(obsChange,
                       c(c08[which(c08$obs == i), "clus8"], c06[which(c06$obs == i), "index6"], c08[which(c08$obs == i), "index8"]))
  }
  
}

obsChange2 <- NULL
for(i in 1:90){
  
  if(c08[which(c08$obs == i), "clus8"] != c12[which(c12$obs == i), "clus12"]){
    obsChange2 <- rbind(obsChange2,
                        c(c12[which(c12$obs == i), "clus12"], c08[which(c08$obs == i), "index8"], c12[which(c12$obs == i), "index12"]))
  }
  
}

dim(obsChange2)

ggplot() +
  geom_text(data = c06, aes(x = rep(1, 90), y = index6, label = obs, color = clus6), size = 2.5) +
  geom_text(data = c08, aes(x = rep(2, 90), y = index8, label = obs, color = clus8), size = 2.5) +
  geom_text(data = c12, aes(x = rep(3, 90), y = index12, label = obs, color = clus12), size = 2.5) +
  geom_segment(aes(x = rep(1.01, 47), y = obsChange[, 2], 
                   xend = rep(1.99, 47), yend = obsChange[, 3],
                   color = factor(obsChange[, 1])), alpha = 0.4, linetype = "dashed") +
  geom_segment(aes(x = rep(2.01, 62), y = obsChange2[, 2], 
                   xend = rep(2.99, 62), yend = obsChange2[, 3],
                   color = factor(obsChange2[, 1])), alpha = 0.4, linetype = "dashed") +
  theme_minimal() + 
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("6 months", "8 months", "12 months")) +
  labs(x = "Timestamps", y = "", title = "The change of cluster behavior for each infants across three timestamps")

### ----------------------------------------------------------------------------



registerDoParallel(5)
dat06ggPlot <- foreach(t = 1:5) %dopar% {
  
  main_cap <- paste0("Cluster ", t, ": 6-Month Infants")
  dat06Plot %>%
    dplyr::filter(ID %in% obsID[which(clusBinder[, 1] == t)]) %>%
    ggplot(aes(x = ID, y = prop_taxa, fill = factor(g))) +
    geom_bar(position = "stack", stat="identity") +
    scale_fill_manual("", values=colours) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5)) +  # Vertical x-axis tick labels
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Relative abundance", title = main_cap)
  
}
stopImplicitCluster()

grid.arrange(grobs = dat06ggPlot)

registerDoParallel(5)
dat08Plot <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(taxon, ID = rep(obsID[t], 38), 
             taxa = as.numeric(dat08[t, -(1:5)])) %>%
    group_by(ID, p, g) %>%
    summarise(sum_taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(prop_taxa = sum_taxa/sum(sum_taxa)) %>%
    dplyr::select(-sum_taxa)
  
}
stopImplicitCluster()

dat08Plot$g <- factor(dat08Plot$g,
                      levels = c("Bifidobacterium", "Other Bifidobacterium",
                                 "Bacteroides", "Other Bacteroidetes", "Prevotella",
                                 "Faecalibacterium", "Megasphaera", "Other Firmicutes",
                                 "Ruminococcus", "Streptococcus", "Veillonella", 
                                 "Others"))

registerDoParallel(5)
dat08ggPlot <- foreach(t = 1:4) %dopar% {
  
  main_cap <- paste0("Cluster ", t, ": 8-Month Infants")
  dat08Plot %>%
    dplyr::filter(ID %in% obsID[which(clusBinder[, 2] == t)]) %>%
    ggplot(aes(x = ID, y = prop_taxa, fill = factor(g))) +
    geom_bar(position = "stack", stat="identity") +
    scale_fill_manual("", values=colours) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5)) +  # Vertical x-axis tick labels
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Relative abundance", title = main_cap)
  
}
stopImplicitCluster()

grid.arrange(grobs = dat08ggPlot)

registerDoParallel(5)
dat12Plot <- foreach(t = 1:90, .combine = "rbind") %dopar% {
  
  data.frame(taxon, ID = rep(obsID[t], 38), 
             taxa = as.numeric(dat12[t, -(1:5)])) %>%
    group_by(ID, p, g) %>%
    summarise(sum_taxa = sum(taxa)) %>%
    ungroup() %>%
    mutate(prop_taxa = sum_taxa/sum(sum_taxa)) %>%
    dplyr::select(-sum_taxa)
  
}
stopImplicitCluster()

dat12Plot$g <- factor(dat12Plot$g,
                      levels = c("Bifidobacterium", "Other Bifidobacterium",
                                 "Bacteroides", "Other Bacteroidetes", "Prevotella",
                                 "Faecalibacterium", "Megasphaera", "Other Firmicutes",
                                 "Ruminococcus", "Streptococcus", "Veillonella", 
                                 "Others"))

registerDoParallel(5)
dat12ggPlot <- foreach(t = 1:4) %dopar% {
  
  main_cap <- paste0("Cluster ", t, ": 12-Month Infants")
  dat12Plot %>%
    dplyr::filter(ID %in% obsID[which(clusBinder[, 3] == t)]) %>%
    ggplot(aes(x = ID, y = prop_taxa, fill = factor(g))) +
    geom_bar(position = "stack", stat="identity") +
    scale_fill_manual("", values=colours) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5)) +  # Vertical x-axis tick labels
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Relative abundance", title = main_cap)
  
}
stopImplicitCluster()

grid.arrange(grobs = dat12ggPlot)



sapply(1:3, 
       function(x){apply(annikaZZ[[x]]$result$ci_result, 1, 
                         function(y){length(unique(y))})}) %>%
  matplot(type = "l", ylim = c(1, 10), main = " ",
          xlab = "Iteration (Thinning)", ylab = "Active Clusters")

clusVI <- sapply(1:3, 
                 function(x){as.numeric(salso(annikaZZ[[x]]$result$ci_result[-(1:500), ]))})

mclustcomp(clusVI[, 1], clusBinder[, 1]); table(VI = clusVI[, 1], Binder = clusBinder[, 1])
mclustcomp(clusVI[, 2], clusBinder[, 2]); table(VI = clusVI[, 2], Binder = clusBinder[, 2])
mclustcomp(clusVI[, 3], clusBinder[, 3]); table(VI = clusVI[, 3], Binder = clusBinder[, 3])

table(clusBinder[, 1])
table(clusBinder[, 2])
table(clusBinder[, 3])



### Visualzation ---------------------------------------------------------------
#### Line plot: Change of the cluster
c6 <- data.frame(obs = 1:90, clus6 = factor(clusBinder[, 1])) %>%
  arrange(clus6) %>% 
  mutate(index6 = 90:1)
c8 <- data.frame(obs = 1:90, clus8 = factor(clusBinder[, 2])) %>%
  arrange(clus8) %>% 
  mutate(index8 = 90:1)
c12 <- data.frame(obs = 1:90, clus12 = factor(clusBinder[, 3])) %>%
  arrange(clus12) %>% 
  mutate(index12 = 90:1)

obsChange <- NULL
for(i in 1:90){
  
  if(c6[which(c6$obs == i), "clus6"] != c8[which(c8$obs == i), "clus8"]){
    obsChange <- rbind(obsChange,
                       c(c8[which(c8$obs == i), "clus8"], c6[which(c6$obs == i), "index6"], c8[which(c8$obs == i), "index8"]))
  }
  
}

obsChange2 <- NULL
for(i in 1:90){
  
  if(c8[which(c8$obs == i), "clus8"] != c12[which(c12$obs == i), "clus12"]){
    obsChange2 <- rbind(obsChange2,
                        c(c12[which(c12$obs == i), "clus12"], c8[which(c8$obs == i), "index8"], c12[which(c12$obs == i), "index12"]))
  }
  
}

ggplot() +
  geom_text(data = c6, aes(x = rep(1, 90), y = index6, label = obs, color = clus6), size = 2.5) +
  geom_text(data = c8, aes(x = rep(2, 90), y = index8, label = obs, color = clus8), size = 2.5) +
  geom_text(data = c12, aes(x = rep(3, 90), y = index12, label = obs, color = clus12), size = 2.5) +
  geom_segment(aes(x = rep(1.01, 77), y = obsChange[, 2], 
                   xend = rep(1.99, 77), yend = obsChange[, 3],
                   color = factor(obsChange[, 1])), alpha = 0.4, linetype = "dashed") +
  geom_segment(aes(x = rep(2.01, 62), y = obsChange2[, 2], 
                   xend = rep(2.99, 62), yend = obsChange2[, 3],
                   color = factor(obsChange2[, 1])), alpha = 0.4, linetype = "dashed") +
  theme_minimal() + 
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("6 months", "8 months", "12 months")) +
  labs(x = "Timestamps", y = "", title = "The change of cluster behavior for each infants across three timestamps")


### Venn Diagram
VennData <- vector("list", 13)
#### 6-month
VennData[[1]] <- which(clusBinder[, 1] == 1) 
VennData[[2]] <- which(clusBinder[, 1] == 2) 
VennData[[3]] <- which(clusBinder[, 1] == 3) 
VennData[[4]] <- which(clusBinder[, 1] == 4) 
VennData[[5]] <- which(clusBinder[, 1] == 5)
#### 8-month
VennData[[6]] <- which(clusBinder[, 2] == 1) 
VennData[[7]] <- which(clusBinder[, 2] == 2) 
VennData[[8]] <- which(clusBinder[, 2] == 3) 
VennData[[9]] <- which(clusBinder[, 2] == 4)
#### 12-month
VennData[[10]] <- which(clusBinder[, 3] == 1) 
VennData[[11]] <- which(clusBinder[, 3] == 2) 
VennData[[12]] <- which(clusBinder[, 3] == 3) 
VennData[[13]] <- which(clusBinder[, 3] == 4)

names(VennData) <- c(paste0("06M: Cluster ", 1:5), 
                     paste0("08M: Cluster ", 1:4),
                     paste0("12M: Cluster ", 1:4))

ggVennDiagram(VennData, order.set.by = "name", relative_height = 0.35)

c6V1 <- ggVennDiagram(VennData[c(1, 6:9)], label_alpha = 0,
                      label = "count", edge_lty = 1, edge_size = 0.5,
                      category.names = c("Original","Cluster 1", 
                                         "Cluster 2","Cluster 3", "Cluster 4")) +
  scale_fill_gradientn(colours = c('grey75', "white", 'white'),
                       values = c(0, .Machine$double.eps, 1)) +
  theme(legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  labs(title = "Originally in Cluster 1 (6m-infant)")

c6V2 <- ggVennDiagram(VennData[c(2, 6:9)], label_alpha = 0,
                      label = "count", edge_lty = 1, edge_size = 0.5,
                      category.names = c("6M: Cluster 2","8M: Cluster 1", 
                                         "8M: Cluster 2","8M: Cluster 3", "8M: Cluster 4")) +
  scale_fill_gradientn(colours = c('grey75', "white", 'white'),
                       values = c(0, .Machine$double.eps, 1)) +
  theme(legend.position = "none")

grid.arrange(c6V1, c6V1, c6V1, c6V1, c6V1)

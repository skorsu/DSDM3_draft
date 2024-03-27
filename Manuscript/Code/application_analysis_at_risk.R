### Required Library
library(stringr)
library(tidyverse)
library(salso)
library(foreach)
library(doParallel)
library(xtable)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(rbiom)
library(pheatmap)
library(readxl)
library(MCMCprecision)

### User-defined Functions -----------------------------------------------------
meanSD <- function(x, dplace = 3){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

### Import the data and result -------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

datpath <- "/Users/kevin-imac/Desktop/Annika/"
if(! file.exists(datpath)){
  datpath <- "/Users/kevinkvp/Desktop/Annika/Application/"
}

### Additional Data
addml <- read_excel(paste0(datpath, "Data/Mali_RB Metadata.xlsx"))
addni <- read_excel(paste0(datpath, "Data/Nicaragua_RB Metadata.xlsx"))

### Data: 6 and 8 Months
ni68 <- read.csv(paste0(datpath, "Data/Nicaragua_6mo_8mo_genus.csv"))
ml68 <- read.csv(paste0(datpath, "Data/Mali_6mo_8mo_genus.csv"))

### Data: 12 Months
ni12 <- read.csv(paste0(datpath, "Data/Nicaragua_12mo_Metadata_csv.csv"))
ml12 <- read.csv(paste0(datpath, "Data/Mali_12mo_Metadata_csv.csv"))

### Result with at-risk indicator
annikaZZ <- readRDS(paste0(path, "Result/microbiome_result_at_risk.RData"))

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

## Cluster from salso: ---------------------------------------------------------
salso_clus <- sapply(1:3, 
                     function(x){as.numeric(salso(annikaZZ[[x]]$result$ci_result[-(1:500), ]))})

sapply(1:3, function(x){apply(annikaZZ[[x]]$result$ci_result, 1, function(y){length(unique(y))})}) %>%
  as.data.frame() %>%
  `colnames<-`(c("6MO", "8MO", "12MO")) %>%
  mutate(iter = 1:1000) %>%
  pivot_longer(!(iter), names_to = "Month", values_to = "Cluster") %>%
  ggplot(aes(x = iter, y = Cluster, 
             color = factor(Month, levels = c("6MO", "8MO", "12MO"),
                            labels = c("6 Month", "8 Month", "12 Month")))) +
  geom_line() +
  theme_minimal() +
  labs(x = "Thinned Iteration", y = "Number of Active Clusters",
       color = "Timestamp", title = "Active clusters over MCMC iterations for each timestamp") +
  scale_y_continuous(limits = c(1, 10), breaks = seq(1, 10, 1)) +
  theme(legend.position = "bottom", legend.box = "horizontal")

### Need to add: gender, nationality, rice bran

table(salso_clus[, 1], dat06[, 1])

### Group Taxa: ----------------------------------------------------------------

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

showName <- factor(showName,
                   levels = c("Bifidobacterium", "Other Actinobacteria",
                              "Bacteroides", "Prevotella", "Other Bacteroidetes",
                              "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                              "Streptococcus", "Veillonella", "Other Firmicutes", 
                              "Others"),
                   labels = c("Bifidobacterium", "Other_Actinobacteria",
                              "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                              "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                              "Streptococcus", "Veillonella", "Other_Firmicutes", 
                              "Others"))

### Estimates: -----------------------------------------------------------------
#### Check the convergence



#### Calculate the estimate relative abundance for each observations
datList <- list(dat06, dat08, dat12)
obsID <- c(str_replace_all(str_extract(dat08$ID[1:45], "^([:alpha:]{2}\\.){2}[:digit:]{2}"), "\\.", "\\-"),
           str_extract(dat08$ID[-(1:45)], "^[:digit:]{7}"))
bactList <- c("Bifidobacterium", "Other_Actinobacteria", "Bacteroides", 
              "Prevotella", "Other_Bacteroidetes", "Faecalibacterium", 
              "Megasphaera", "Ruminococcus", "Streptococcus", "Veillonella", 
              "Other_Firmicutes", "Others")
colours <- c("#A4C639", "#BAC394", "#DDA661", "#BD8B46", "#FFD197",
             "#BC6D62", "#DE6A5D", "#DE998B", "#FF6558", "#FFC6B8", 
             "#FFE4DA", "#E5E5E5")

relaPlot <- function(timestamp_index, actual_month){
  
  ### Get the estimated relative abundance
  set.seed(1415, kind = "L'Ecuyer-CMRG")
  registerDoParallel(5)
  estiRela <- foreach(i = 1:90, .combine = "rbind") %dopar% {
    
    beteVec <- matrix(NA, nrow = 1000, ncol = 38)
    
    for(t in 1:1000){
      ### Beta vector
      ci <- annikaZZ[[timestamp_index]]$result$ci_result[t, i] + 1
      beteVec[t, ] <- annikaZZ[[timestamp_index]]$result$beta_result[ci, ,t]
    }
    
    atriskBetaE <- exp(beteVec) * t(annikaZZ[[timestamp_index]]$result$atrisk_result[i, , ])
    colMeans(rdirichlet(1000, colMeans(atriskBetaE[-(1:500), ])))
  }
  stopImplicitCluster()
  
  ### Group the bacteria
  actualDat <- datList[[timestamp_index]]
  estiRelaGroup <- matrix(NA, nrow = 90, ncol = 12)
  actualRela <- actualDat[, -(1:5)]/rowSums(actualDat[, -(1:5)])
  actualRelaGroup <- matrix(NA, nrow = 90, ncol = 12)
  
  for(i in 1:length(bactList)){
    index <- which(showName == bactList[i])
    if(length(index) == 1){
      estiRelaGroup[, i] <- estiRela[, which(showName == bactList[i])]
      actualRelaGroup[, i] <- actualRela[, which(showName == bactList[i])]
    } else {
      estiRelaGroup[, i] <- rowSums(estiRela[, which(showName == bactList[i])])
      actualRelaGroup[, i] <- rowSums(actualRela[, which(showName == bactList[i])])
    }
  }
  
  colnames(estiRelaGroup) <- bactList
  colnames(actualRelaGroup) <- bactList
  
  ### Actual
  actual_title <- paste0(actual_month, " Infants: Actual Relative Abundance")
  actPlot <- data.frame(actualRelaGroup) %>%
    mutate(id = obsID, Cluster = paste0("Cluster ", salso_clus[, timestamp_index])) %>%
    pivot_longer(!c(id, Cluster), names_to = "Bacteria", values_to = "Prop") %>%
    ggplot(aes(fill = factor(Bacteria, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                  "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                  "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                  "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                  "Others"),
                             labels = c("Bifidobacterium", "Other Actinobacteria",
                                        "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                        "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                        "Streptococcus", "Veillonella", "Other Firmicutes", 
                                        "Others")), y = Prop, x = id)) + 
    geom_bar(position="fill", stat = "identity") +
    scale_fill_manual("", values = colours) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Relative Abundance", title = actual_title, x = "Sample") +
    facet_grid(~Cluster, scales = "free_x")
  
  ### Estimate
  estimate_title <- paste0(actual_month, " Infants: Estimated Relative Abundance")
  estPlot <- data.frame(estiRelaGroup) %>%
    mutate(id = obsID, Cluster = paste0("Cluster ", salso_clus[, timestamp_index])) %>%
    pivot_longer(!c(id, Cluster), names_to = "Bacteria", values_to = "Prop") %>%
    ggplot(aes(fill = factor(Bacteria, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                  "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                  "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                  "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                  "Others"),
                             labels = c("Bifidobacterium", "Other Actinobacteria",
                                        "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                        "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                        "Streptococcus", "Veillonella", "Other Firmicutes", 
                                        "Others")), y = Prop, x = id)) + 
    geom_bar(position="fill", stat = "identity") +
    scale_fill_manual("", values = colours) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Estimated Relative Abundance", title = estimate_title, x = "Sample") +
    facet_grid(~Cluster, scales = "free_x")
  
  list(estiRelaGroup = estiRelaGroup, 
       actualRelaGroup = actualRelaGroup, plot = grid.arrange(actPlot, estPlot))
  
}

mo6plot <- relaPlot(timestamp_index = 1, actual_month = "6-Month")
mo8plot <- relaPlot(timestamp_index = 2, actual_month = "8-Month")
mo12plot <- relaPlot(timestamp_index = 3, actual_month = "12-Month")

###

data.frame(dat06[, c(1, 3, 5)], cluster = salso_clus[, 1]) %>%
  group_by(cluster, country) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = factor(cluster, labels = paste0("Cluster", 1:3)), y = prop, fill = country)) +
  geom_bar(position="fill", stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = " ", y = "Proportion", 
       title = "6-Month: Proportion of the infant by their nationality")

data.frame(dat06[, c(1, 3, 5)], cluster = salso_clus[, 1]) %>%
  group_by(cluster, country) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = factor(cluster, labels = paste0("Cluster", 1:3)), y = prop, fill = country)) +
  geom_bar(position="fill", stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = " ", y = "Proportion", 
       title = "6-Month: Proportion of the infant by their nationality")


### Draft: ---------------------------------------------------------------------
###### Plot

relaPlot(3, "12-Month")

set.seed(1415, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
estiRela <- foreach(i = 1:90, .combine = "rbind") %dopar% {
  
  beteVec <- matrix(NA, nrow = 1000, ncol = 38)
  
  for(t in 1:1000){
    ### Beta vector
    ci <- annikaZZ[[1]]$result$ci_result[t, i] + 1
    beteVec[t, ] <- annikaZZ[[1]]$result$beta_result[ci, ,t]
  }
  
  atriskBetaE <- exp(beteVec) * t(annikaZZ[[1]]$result$atrisk_result[i, , ])
  colMeans(rdirichlet(1000, colMeans(atriskBetaE[-(1:500), ])))
}
stopImplicitCluster()

### Collapse to the important taxa
obsID <- c(str_replace_all(str_extract(dat08$ID[1:45], "^([:alpha:]{2}\\.){2}[:digit:]{2}"), "\\.", "\\-"),
           str_extract(dat08$ID[-(1:45)], "^[:digit:]{7}"))
bactList <- c("Bifidobacterium", "Other_Actinobacteria", "Bacteroides", 
              "Prevotella", "Other_Bacteroidetes", "Faecalibacterium", 
              "Megasphaera", "Ruminococcus", "Streptococcus", "Veillonella", 
              "Other_Firmicutes", "Others")

estiRelaGroup <- matrix(NA, nrow = 90, ncol = 12)
actualRela <- dat06[, -(1:5)]/rowSums(dat06[, -(1:5)])
actualRelaGroup <- matrix(NA, nrow = 90, ncol = 12)

for(i in 1:length(bactList)){
  index <- which(showName == bactList[i])
  if(length(index) == 1){
    estiRelaGroup[, i] <- estiRela[, which(showName == bactList[i])]
    actualRelaGroup[, i] <- actualRela[, which(showName == bactList[i])]
  } else {
    estiRelaGroup[, i] <- rowSums(estiRela[, which(showName == bactList[i])])
    actualRelaGroup[, i] <- rowSums(actualRela[, which(showName == bactList[i])])
  }
}

colnames(estiRelaGroup) <- bactList
colnames(actualRelaGroup) <- bactList

colours <- c("#A4C639", "#BAC394", "#DDA661", "#BD8B46", "#FFD197",
             "#BC6D62", "#DE6A5D", "#DE998B", "#FF6558", "#FFC6B8", 
             "#FFE4DA", "#E5E5E5")

### Actual
actPlot <- data.frame(actualRelaGroup) %>%
  mutate(id = obsID, Cluster = paste0("Cluster ", salso_clus[, 1])) %>%
  pivot_longer(!c(id, Cluster), names_to = "Bacteria", values_to = "Prop") %>%
  ggplot(aes(fill = factor(Bacteria, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                "Others"),
                           labels = c("Bifidobacterium", "Other Actinobacteria",
                                      "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                      "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                      "Streptococcus", "Veillonella", "Other Firmicutes", 
                                      "Others")), y = Prop, x = id)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Relative Abundance", title = "6-Month Infants: Actual Relative Abundance", x = "Sample") +
  facet_grid(~Cluster, scales = "free_x")


### Estimate
estPlot <- data.frame(estiRelaGroup) %>%
  mutate(id = obsID, Cluster = paste0("Cluster ", salso_clus[, 1])) %>%
  pivot_longer(!c(id, Cluster), names_to = "Bacteria", values_to = "Prop") %>%
  ggplot(aes(fill = factor(Bacteria, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                "Others"),
                           labels = c("Bifidobacterium", "Other Actinobacteria",
                                      "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                      "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                      "Streptococcus", "Veillonella", "Other Firmicutes", 
                                      "Others")), y = Prop, x = id)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual("", values = colours) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Estimated Relative Abundance", title = "6-Month Infants: Estimated Relative Abundance", x = "Sample") +
  facet_grid(~Cluster, scales = "free_x")

grid.arrange(actPlot, estPlot)

######

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

data.frame(est_prob_group, iter = 1:1000) %>%
  pivot_longer(!iter, names_to = "Bacteria", values_to = "Prop") %>%
  ggplot(aes(x = iter, y = Prop, 
             color = factor(Bacteria, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                 "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                 "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                 "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                 "Others"),
                            labels = c("Bifidobacterium", "Other Actinobacteria",
                                       "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                       "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                       "Streptococcus", "Veillonella", "Other Firmicutes", 
                                       "Others")))) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "vertical", legend.margin = margin(),
        legend.title = element_blank()) + 
  scale_color_manual(values = c("#A4C639", "#BAC394", "#DDA661", "#BD8B46", "#FFD197",
                                "#BC6D62", "#DE6A5D", "#DE998B", "#FF6558", "#FFC6B8", 
                                "#FFE4DA", "#E5E5E5")) +
  labs(x = "Thinned Iteration", y = "Relative Abundance")





































# ------------------------------------------------------------------------------
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

lineAdd_12 <- data.frame(x = rep(1, 20), y = sort(rep(1:5, 4)), 
                      xend = rep(2, 20), yend = rep(1:4, 5), 
                      f = (table(clusBinderN[, 1], clusBinderN[, 2])/rowSums(table(clusBinderN[, 1], clusBinderN[, 2]))) %>%
                        t() %>%
                        as.numeric()) %>%
  filter(f != 0)
lineAdd_23 <- data.frame(x = rep(2, 16), y = sort(rep(1:4, 4)), 
                         xend = rep(3, 16), yend = rep(1:4, 4), 
                         f = (table(clusBinderN[, 2], clusBinderN[, 3])/rowSums(table(clusBinderN[, 2], clusBinderN[, 3]))) %>%
                           t() %>%
                           as.numeric()) %>%
  filter(f != 0)

rbind(data.frame(ci = 1:5, m = rep("6 Month", 5), 
                 meanB = dat06Nrela %>%
                   dplyr::select(-clus) %>%
                   pivot_wider(names_from = showName, values_from = relaTaxa) %>%
                   mutate(clus = clusBinderN[, 1]) %>%
                   group_by(clus) %>%
                   summarise(meanBifidobacterium = mean(Bifidobacterium),
                             meanBacteroides = mean(Bacteroides),
                             meanPrevotella = mean(Prevotella)) %>%
                   .$meanBifidobacterium),
      data.frame(ci = 1:4, m = rep("8 Month", 4),
                 meanB = dat08Nrela %>%
                   dplyr::select(-clus) %>%
                   pivot_wider(names_from = showName, values_from = relaTaxa) %>%
                   mutate(clus = clusBinderN[, 2]) %>%
                   group_by(clus) %>%
                   summarise(meanBifidobacterium = mean(Bifidobacterium),
                             meanBacteroides = mean(Bacteroides),
                             meanPrevotella = mean(Prevotella)) %>%
                   .$meanBifidobacterium),
      data.frame(ci = 1:4, m = rep("12 Month", 4),
                 meanB = dat12Nrela %>%
                   dplyr::select(-clus) %>%
                   pivot_wider(names_from = showName, values_from = relaTaxa) %>%
                   mutate(clus = clusBinderN[, 3]) %>%
                   group_by(clus) %>%
                   summarise(meanBifidobacterium = mean(Bifidobacterium),
                             meanBacteroides = mean(Bacteroides),
                             meanPrevotella = mean(Prevotella)) %>%
                   .$meanBifidobacterium)) %>%
  mutate(sizeN = unlist(lapply(1:3, function(x){as.numeric(table(clusBinderN[, x]))}))) %>%
  ggplot(aes(x = factor(m, levels = c("6 Month", "8 Month", "12 Month")), 
             y = ci, label = sizeN)) +
  geom_point(aes(color = meanB), size = 15) +
  geom_text() +
  geom_segment(data = lineAdd_12, 
               aes(x = x + 0.075, y = y, xend = xend - 0.075, yend = yend, linewidth = f), 
               alpha = 0.5, inherit.aes = FALSE) +
  geom_segment(data = lineAdd_23, 
               aes(x = x + 0.075, y = y, xend = xend - 0.075, yend = yend, linewidth = f), 
               alpha = 0.5, inherit.aes = FALSE) +
  scale_y_reverse() +
  # geom_text(data = lineAdd_12, aes(x = x + 0.15, y = y + yend/4, label = f), inherit.aes = FALSE) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  guides(linewidth = "none") +
  labs(y = "Cluster", x = "Timestamp", 
       color = "The mean of the Bifidobacterium proportion", 
       title = "The change of cluster assignment across three timestamps") +
  scale_color_gradient(low = "grey90", high = "springgreen3")

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

bar1 <- data.frame(x = rep(1, 5), y = 90 - as.numeric(cumsum(table(clusBinderN[, 1]))))
bar2 <- data.frame(x = rep(2, 4), y = 90 - as.numeric(cumsum(table(clusBinderN[, 2]))))
bar3 <- data.frame(x = rep(3, 4), y = 90 - as.numeric(cumsum(table(clusBinderN[, 3]))))

ggplot() +
  geom_text(data = c06, aes(x = rep(1, 90), y = index6, label = obs, color = clus6), size = 4) +
  geom_text(data = c08, aes(x = rep(2, 90), y = index8, label = obs, color = clus8), size = 4) +
  geom_text(data = c12, aes(x = rep(3, 90), y = index12, label = obs, color = clus12), size = 4) +
  geom_segment(aes(x = rep(1.01, 55), y = obsChange[, 2], 
                   xend = rep(1.99, 55), yend = obsChange[, 3]),
               color = "grey50", alpha = 0.6, linewidth = 0.30) +
  geom_segment(aes(x = rep(2.01, 52), y = obsChange2[, 2], 
                   xend = rep(2.99, 52), yend = obsChange2[, 3]), 
               color = "grey50", alpha = 0.6, linewidth = 0.30) +
  geom_segment(data = bar1, 
               aes(x = x - 0.05, y = y, xend = x + 0.05, yend = y), 
               alpha = 0.5, inherit.aes = FALSE) +
  geom_segment(data = bar2, 
               aes(x = x - 0.05, y = y, xend = x + 0.05, yend = y), 
               alpha = 0.5, inherit.aes = FALSE) +
  geom_segment(data = bar3, 
               aes(x = x - 0.05, y = y, xend = x + 0.05, yend = y), 
               alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() + 
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("6 months", "8 months", "12 months")) +
  labs(x = "Timestamps", y = "", title = "The change of cluster behavior for each infants across three timestamps")



### ----------------------------------------------------------------------------

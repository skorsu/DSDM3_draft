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
activeClus <- function(x){
  length(unique(x))
}

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

### Result for sensitivity analysis
annika6M <- readRDS(paste0(path, "Result/microbiome_result_at_risk_sensitivity_6mo.RData"))
annika8M <- readRDS(paste0(path, "Result/microbiome_result_at_risk_sensitivity_8mo.RData"))
annika12M <- readRDS(paste0(path, "Result/microbiome_result_at_risk_sensitivity_12mo.RData"))

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

## Additional Dataset: ---------------------------------------------------------
addDat <- rbind(addml %>%
                  mutate(across(where(is.character), tolower)) %>%
                  mutate(Diarrhea = ifelse(Diarrhea_Episode == "yes", 
                                           ifelse(difftime(Today_Date, Diarrhea_Episode_Date_End) <= 14, "yes", "no"), 
                                           Diarrhea_Episode)) %>%
        rowwise %>% 
        mutate(null_message_GR = as.integer(all(c(Rice, Maize, Cowpea, Porridge) %in% NA)) > 0) %>%
        mutate(allNo_GR = all(c(Rice, Maize, Cowpea, Porridge) %in% "no")) %>%
        mutate(GR = ifelse(null_message_GR == TRUE, NA, ifelse(allNo_GR == TRUE, "no", "yes"))) %>%
        mutate(null_message_PR = as.integer(all(c(`Chicken / fish / meat`, Eggs) %in% NA)) > 0) %>%
        mutate(allNo_PR = all(c(`Chicken / fish / meat`, Eggs) %in% "no")) %>%
        mutate(Protein = ifelse(null_message_PR == TRUE, NA, ifelse(allNo_PR == TRUE, "no", "yes"))) %>%
        dplyr::select(Child_ID, Month =`Monthly_Visit#`, Breastfed = Breastfed_Yesterday, Cow_milk, 
                      Fruits_Natural_Juices, Vegetables, Soups, Antibiotics = Receive_Antibiotics, 
                      Weight, Height, Diarrhea, GR, Protein) %>%
        mutate(Month = toupper(Month)) %>%
        mutate(Nationality = "ML"),
      
      addni %>%
        mutate(across(where(is.character), tolower)) %>%
        rowwise %>% 
        mutate(null_message_GR= as.integer(all(c(Rice_cereal_5.3, Gallo_Pinto_5.7) %in% NA)) > 0) %>%
        mutate(allNo_GR = all(c(Rice_cereal_5.3, Gallo_Pinto_5.7) %in% c("used to consume, but not now", "never"))) %>%
        mutate(GR = ifelse(null_message_GR == TRUE, NA, ifelse(allNo_GR == TRUE, "no", "yes"))) %>%
        mutate(null_message_PR = as.integer(all(c(Chicken_5.6, Cheese_5.8, Yogurt_5.10) %in% NA)) > 0) %>%
        mutate(allNo_PR = all(c(Chicken_5.6, Cheese_5.8, Yogurt_5.10) %in% c("used to consume, but not now", "never"))) %>%
        mutate(Protein = ifelse(null_message_PR == TRUE, NA, ifelse(allNo_PR == TRUE, "no", "yes"))) %>%
        dplyr::select(Child_ID, Month =`Monthly_Visit#`, Breastfed = Breastfed_Yesterday_4.1, Cow_milk = Cow_milk_5.2,
                      Fruits_Natural_Juices = Fruits_Natural_Juices_5.4, Vegetables = Vegetables_5.5,
                      Soups = Soups_5.9, Antibiotics = Receive_Antibiotics_6, 
                      Weight = Weight_11, Height = Height_12,
                      Diarrhea = Diarrhea_Episode_Past_14_Days_3, GR, Protein) %>% 
        mutate(Cow_milk = ifelse(is.na(Cow_milk), NA, ifelse(Cow_milk %in% c("used to consume, but not now", "never"), "no", "yes")),
               Fruits_Natural_Juices = ifelse(is.na(Fruits_Natural_Juices), NA, ifelse(Fruits_Natural_Juices %in% c("used to consume, but not now", "never"), "no", "yes")),
               Vegetables = ifelse(is.na(Vegetables), NA, ifelse(Vegetables %in% c("used to consume, but not now", "never"), "no", "yes")),
               Soups = ifelse(is.na(Soups), NA, ifelse(Soups %in% c("used to consume, but not now", "never"), "no", "yes"))) %>%
        mutate(Child_ID = toupper(Child_ID), Month = str_replace_all(toupper(Month), " ", "")) %>%
        mutate(Nationality = "NI"))

### Number of active clusters --------------------------------------------------
resultList <- list(annika6M, annika8M, annika12M)
activeList <- lapply(1:3, function(x){sapply(1:20, function(y){apply(resultList[[x]][[y]]$result$ci_result, 1, activeClus)})})

rbind(activeList[[1]] %>%
        `colnames<-`(paste0("Chain ", 1:20)) %>%
        as.data.frame() %>%
        mutate(Iteration = 1:1000, Timestamp = "6-Month") %>%
        pivot_longer(!c(Iteration, Timestamp), names_to = "Chain", values_to = "aClus"),
      activeList[[2]] %>%
        `colnames<-`(paste0("Chain ", 1:20)) %>%
        as.data.frame() %>%
        mutate(Iteration = 1:1000, Timestamp = "8-Month") %>%
        pivot_longer(!c(Iteration, Timestamp), names_to = "Chain", values_to = "aClus"),
      activeList[[3]] %>%
        `colnames<-`(paste0("Chain ", 1:20)) %>%
        as.data.frame() %>%
        mutate(Iteration = 1:1000, Timestamp = "12-Month") %>%
        pivot_longer(!c(Iteration, Timestamp), names_to = "Chain", values_to = "aClus")) %>%
  ggplot(aes(x = Iteration, y = aClus, color = Chain)) +
  geom_line() +
  theme_bw() +
  guides(color = "none") +
  facet_grid(. ~ factor(Timestamp, levels = c("6-Month", "8-Month", "12-Month"))) +
  scale_y_continuous(name = "Active Clusters", limits = c(1, 10), breaks = seq(1, 10)) +
  labs(title = "Active Clusters via MCMC for each timestamp")

### Result from salso function: ------------------------------------------------
salsoList <- lapply(1:3, function(x){sapply(1:20, function(y){salso(resultList[[x]][[y]]$result$ci_result[-(1:500), ])})})
timeStampIndex <- c("6-Month", "8-Month", "12-Month")
lapply(1:3, function(x){apply(salsoList[[x]], 2, activeClus) %>%
    table() %>%
    as.data.frame() %>%
    `colnames<-`(c("nClus", "Freq")) %>%
    mutate(Timestamp = timeStampIndex[x])})

sapply(1:20, function(x){salso(annika6M[[x]]$result$ci_result[-(1:500), ])}) %>%
  apply(2, activeClus) %>%
  table()

### Estimated Relative Abundance: ----------------------------------------------


### Cluster from salso: --------------------------------------------------------
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
       actualRelaGroup = actualRelaGroup, plot = grid.arrange(actPlot, estPlot),
       estiRelaUngroup = estiRela, actualRelaUngroup = actualRela)
  
}

mo6plot <- relaPlot(timestamp_index = 1, actual_month = "6-Month")
mo8plot <- relaPlot(timestamp_index = 2, actual_month = "8-Month")
mo12plot <- relaPlot(timestamp_index = 3, actual_month = "12-Month")

### Table: Relative Taxa: ------------------------------------------------------
#### Bacteria of interest: c("Bifidobacterium", "Bacteroides", "Prevotella", "Streptococcus")
list(mo6plot$actualRelaGroup, mo8plot$actualRelaGroup, mo12plot$actualRelaGroup)

reindexFn <- function(listDat, matClus, bactSort){
  
  clusMat <- matrix(NA, ncol = 3, nrow = 90)
  
  for(t in 1:3){
    reindex_cluster <- data.frame(listDat[[t]], Cluster = matClus[, t]) %>%
      group_by(Cluster) %>%
      summarise(meanInt = mean({{bactSort}})) %>%
      arrange(-meanInt) %>%
      .$Cluster
    
    k <- length(reindex_cluster)
    for(i in 1:k){
      clusMat[which(matClus[, t] == reindex_cluster[i]), t] <- i
    }
    
  }
  
  clusMat
    
}

newClus <- reindexFn(list(mo6plot$actualRelaGroup, mo8plot$actualRelaGroup, mo12plot$actualRelaGroup),
                     salso_clus, Bifidobacterium)

### Cluster Size: --------------------------------------------------------------
rbind(c(as.numeric(table(newClus[, 1])), rep(NA, 2)),
      as.numeric(table(newClus[, 2])),
      as.numeric(table(newClus[, 3]))) %>%
  `rownames<-`(c("6-Month", "8-Month", "12-Month")) %>%
  `colnames<-`(paste0("Cluster ", 1:5)) %>%
  xtable()

### Average Relative Abundance (Visualization): --------------------------------
heatPlotAVG <- function(relaPlotObj, cluster_assign, monthLab){
  
  k <- length(unique(cluster_assign))
  meanEst <- matrix(NA, nrow = k, ncol = 12)
  meanSDEst <- matrix(NA, nrow = k, ncol = 12)
  meanAct <- matrix(NA, nrow = k, ncol = 12)
  meanSDAct <- matrix(NA, nrow = k, ncol = 12)
  bactName <- colnames(relaPlotObj$actualRelaGroup)
  
  for(i in 1:k){
    
    actPop <- data.frame(Cluster = cluster_assign, relaPlotObj$actualRelaGroup) %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::select(-Cluster)
    meanAct[i, ] <- colMeans(actPop)
    meanSDAct[i, ] <- apply(actPop, 2, meanSD)
    
    EstPop <- data.frame(Cluster = cluster_assign, relaPlotObj$estiRelaGroup) %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::select(-Cluster)
    meanEst[i, ] <- colMeans(EstPop)
    meanSDEst[i, ] <- apply(EstPop, 2, meanSD)
    
  }
  
  actTitle <- paste0(monthLab, ": Actual Average relative abundance")
  estTitle <- paste0(monthLab, ": Estimated Average relative abundance")
  
  actHeat <- inner_join(t(meanAct) %>%
    `colnames<-`(paste0("Cluster_", 1:k)) %>%
    data.frame(Bact = bactName) %>%
    pivot_longer(!Bact, values_to = "prop_col"),
    t(meanSDAct) %>%
      `colnames<-`(paste0("Cluster_", 1:k)) %>%
      data.frame(Bact = bactName) %>%
      pivot_longer(!Bact, values_to = "prop_print")) %>%
    ggplot(aes(x = factor(name, levels = paste0("Cluster_", 1:k), 
                          labels = paste0("Cluster ", 1:k)), 
               y = forcats::fct_rev(factor(Bact, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                            "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                            "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                            "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                            "Others"),
                                           labels = c("Bifidobacterium", "Other Actinobacteria",
                                                      "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                                      "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                      "Streptococcus", "Veillonella", "Other Firmicutes", 
                                                      "Others"))), fill = prop_col, label = prop_print)) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = " ", y = "Bacteria", title = actTitle, 
         fill = "Average Relative Abundance")
  
  estHeat <- inner_join(t(meanEst) %>%
                          `colnames<-`(paste0("Cluster_", 1:k)) %>%
                          data.frame(Bact = bactName) %>%
                          pivot_longer(!Bact, values_to = "prop_col"),
                        t(meanSDEst) %>%
                          `colnames<-`(paste0("Cluster_", 1:k)) %>%
                          data.frame(Bact = bactName) %>%
                          pivot_longer(!Bact, values_to = "prop_print")) %>%
    ggplot(aes(x = factor(name, levels = paste0("Cluster_", 1:k), 
                          labels = paste0("Cluster ", 1:k)), 
               y = forcats::fct_rev(factor(Bact, levels = c("Bifidobacterium", "Other_Actinobacteria",
                                                            "Bacteroides", "Prevotella", "Other_Bacteroidetes",
                                                            "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                            "Streptococcus", "Veillonella", "Other_Firmicutes", 
                                                            "Others"),
                                           labels = c("Bifidobacterium", "Other Actinobacteria",
                                                      "Bacteroides", "Prevotella", "Other Bacteroidetes",
                                                      "Faecalibacterium", "Megasphaera", "Ruminococcus", 
                                                      "Streptococcus", "Veillonella", "Other Firmicutes", 
                                                      "Others"))), fill = prop_col, label = prop_print)) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = " ", y = "Bacteria", title = estTitle, 
         fill = "Average Relative Abundance")
  
  list(actHeat, estHeat)
  
}

grid.arrange(grobs = heatPlotAVG(mo6plot, newClus[, 1], "6-Month"))
grid.arrange(grobs = heatPlotAVG(mo8plot, newClus[, 2], "8-Month"))
grid.arrange(grobs = heatPlotAVG(mo12plot, newClus[, 3], "12-Month"))

### Descriptive Statistic Plot: ------------------------------------------------
descPlot <- function(descDat_type, intVar, capt, monthLab, g1level = "yes", g2level = "no", 
                     g1label = "Yes", g2label = "No", g1Color = "olivedrab3", g2Color = "rosybrown1"){
  
  k <- length(unique(descDat_type$Cluster))
  ptitle <- paste0(monthLab, ": ", capt)
  
  descDat_type %>%
    group_by(Cluster, {{intVar}}) %>%
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = factor(Cluster, labels = paste0("Cluster ", 1:k)), 
               y = freq, fill = factor({{intVar}}, levels = c(g1level, g2level, "N/A"),
                                       label = c(g1label, g2label, "Not Available")))) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(capt, values = c(g1Color, g2Color, "grey90")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = " ", y = "Percentage", title = ptitle)
  
}

plotData <- function(monthabbv, monthLabel, monthDat, monthCluster){
  
  descDat <- data.frame(Child_ID = obsID, monthDat[, c(1, 3, 5)], Cluster = monthCluster) %>%
    left_join(addDat[which(addDat$Month == monthabbv), ]) %>%
    replace(is.na(.), "N/A") %>%
    mutate(Group = str_to_title(Group))
  
  plotList <- list(descPlot(descDat, Sex, "Sex", monthLabel, "M", "F", "Male", "Female", "lightblue", "lightpink"),
                   descPlot(descDat, country, "Nationality", monthLabel, "NI", "ML", "Nicaraguan", "Malian", "darkblue", "lightgreen"),
                   descPlot(descDat, Group, "Group", monthLabel, "Rice Bran", "Control", "Rice Bran", "Control"),
                   descPlot(descDat, Breastfed, "Breastfed", monthLabel),
                   descPlot(descDat, Cow_milk, "Cow milk", monthLabel),
                   descPlot(descDat, Fruits_Natural_Juices, "Fruits Natural Juices", monthLabel),
                   descPlot(descDat, Vegetables, "Vegetables", monthLabel),
                   descPlot(descDat, Soups, "Soups", monthLabel),
                   descPlot(descDat, Antibiotics, "Antibiotics", monthLabel),
                   descPlot(descDat, Diarrhea, "Diarrhea", monthLabel),
                   descPlot(descDat, GR, "Grains and Legumes", monthLabel),
                   descPlot(descDat, Protein, "Protein", monthLabel))
  
  plotList
  
}

plot6mo <- plotData("6MO", "6-Month", dat06, newClus[, 1])
plot8mo <- plotData("8MO", "8-Month", dat08, newClus[, 2])
plot12mo <- plotData("12MO", "12-Month", dat12, newClus[, 3])

i <- 2
grid.arrange(grobs = list(plot6mo[[i]], plot8mo[[i]], plot12mo[[i]]))

grid.arrange(grobs = list(plot6mo[[1]], plot8mo[[1]], plot12mo[[1]],
                          plot6mo[[2]], plot8mo[[2]], plot12mo[[2]],
                          plot6mo[[3]], plot8mo[[3]], plot12mo[[3]]))

### Alpha-Diversity: -----------------------------------------------------------
alphaSumm <- function(relaPlotObj, cluster_assign, moLabel){
  
  ### Calculate the Alpha Diversity
  #### Richness
  actualRich <- as.numeric(rowSums(relaPlotObj$actualRelaUngroup != 0))
  estiRich <- as.numeric(rowSums(relaPlotObj$estiRelaUngroup != 0))
  
  #### Simpson
  actualSS <- as.numeric(1 - rowSums(relaPlotObj$actualRelaUngroup^2))
  estiSS <- as.numeric(1 - rowSums(relaPlotObj$estiRelaUngroup^2))
  
  #### Inverse Simpson
  actualISS <- as.numeric(1/rowSums(relaPlotObj$actualRelaUngroup^2))
  estiISS <- as.numeric(1/rowSums(relaPlotObj$estiRelaUngroup^2))
  
  #### Shannon
  actualSN <- as.matrix(relaPlotObj$actualRelaUngroup)
  actualSN[actualSN == 0] <- 1e-200
  actualSN <- as.numeric(-rowSums(actualSN * log(actualSN)))
  
  estiSN <- as.matrix(relaPlotObj$estiRelaUngroup)
  estiSN[estiSN == 0] <- 1e-200
  estiSN <- as.numeric(-rowSums(estiSN * log(estiSN)))
  
  ### Plot
  k <- length(unique(cluster_assign))
  plotTitle = paste0(moLabel, ": Alpha Diversity")
  
  pplot <- rbind(data.frame(Actual = actualRich, Estimate = estiRich, cluster = cluster_assign) %>%
                   pivot_longer(c(Actual, Estimate)) %>%
                   mutate(Quantity = "Richness"),
                 data.frame(Actual = actualSS, Estimate = estiSS, cluster = cluster_assign) %>%
                   pivot_longer(c(Actual, Estimate)) %>%
                   mutate(Quantity = "Simpson"),
                 data.frame(Actual = actualISS, Estimate = estiISS, cluster = cluster_assign) %>%
                   pivot_longer(c(Actual, Estimate)) %>%
                   mutate(Quantity = "Inverse Simpson"),
                 data.frame(Actual = actualSN, Estimate = estiSN, cluster = cluster_assign) %>%
                   pivot_longer(c(Actual, Estimate)) %>%
                   mutate(Quantity = "Shannon")) %>%
    ggplot(aes(x = name, y = value, fill = factor(cluster, labels = paste0("Cluster ", 1:k)))) +
    geom_boxplot() +
    facet_grid(factor(Quantity, levels = c("Richness", "Shannon", "Simpson", "Inverse Simpson")) ~ ., scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = " ", y = " ", fill = "Cluster", title = plotTitle)
  
  ### xtable
  actualQ <- list(actualRich, actualSN, actualSS, actualISS)
  estiQ <- list(estiRich, estiSN, estiSS, estiISS)
  report_table_actual <- data.frame(matrix(NA, ncol = 4, nrow = k))
  report_table_estimate <- data.frame(matrix(NA, ncol = 4, nrow = k))
  
  for(i in 1:k){
    report_table_actual[i, ] <- sapply(1:4, function(y){meanSD(lapply(1:4, function(x){actualQ[[x]][which(cluster_assign == i)]})[[y]])})
    report_table_estimate[i, ] <- sapply(1:4, function(y){meanSD(lapply(1:4, function(x){estiQ[[x]][which(cluster_assign == i)]})[[y]])})
  }
  
  colnames(report_table_actual) <- c("Richness", "Shannon", "Simpson", "Inverse Simpson")
  colnames(report_table_estimate) <- c("Richness", "Shannon", "Simpson", "Inverse Simpson")
  rownames(report_table_actual) <- paste0("Cluster ", 1:k)
  rownames(report_table_estimate) <- paste0("Cluster ", 1:k)
  
  printList <- list(report_table_actual, report_table_estimate)
  attr(printList, "subheadings") <- c("Actual", "Estimate")
  
  xtableObj <- xtableList(printList)
  
  list(plot = pplot, xtable = xtableObj, actualQ = actualQ, estiQ = estiQ)

}

alpha6mo <- alphaSumm(mo6plot, newClus[, 1], "6-Month")
alpha8mo <- alphaSumm(mo8plot, newClus[, 2], "8-Month")
alpha12mo <- alphaSumm(mo12plot, newClus[, 3], "12-Month")
grid.arrange(grobs = list(alpha6mo$plot, alpha8mo$plot, alpha12mo$plot))
print.xtableList(alpha6mo$xtable, booktabs = TRUE)

### Observations movement: -----------------------------------------------------
obsMovement <- function(clusAssign){
  
  ### Get the bar separating each cluster
  bar1 <- data.frame(x = rep(1, 3), y = 90 - as.numeric(cumsum(table(clusAssign[, 1]))))
  bar2 <- data.frame(x = rep(2, 5), y = 90 - as.numeric(cumsum(table(clusAssign[, 2]))))
  bar3 <- data.frame(x = rep(3, 5), y = 90 - as.numeric(cumsum(table(clusAssign[, 3]))))
  
  ### Get the observation index
  c06 <- data.frame(obs = 1:90, clus6 = factor(clusAssign[, 1])) %>%
    arrange(clus6) %>% 
    mutate(index6 = 90:1)
  c08 <- data.frame(obs = 1:90, clus8 = factor(clusAssign[, 2])) %>%
    arrange(clus8) %>% 
    mutate(index8 = 90:1)
  c12 <- data.frame(obs = 1:90, clus12 = factor(clusAssign[, 3])) %>%
    arrange(clus12) %>% 
    mutate(index12 = 90:1)
  
  ### Get the linking line
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
  
  d1 <- dim(obsChange)[1]
  d2 <- dim(obsChange2)[1]
  
  ggplot() +
    geom_text(data = c06, aes(x = rep(1, 90), y = index6, label = obs, color = clus6), size = 4) +
    geom_text(data = c08, aes(x = rep(2, 90), y = index8, label = obs, color = clus8), size = 4) +
    geom_text(data = c12, aes(x = rep(3, 90), y = index12, label = obs, color = clus12), size = 4) +
    geom_segment(aes(x = rep(1.01, d1), y = obsChange[, 2], 
                     xend = rep(1.99, d1), yend = obsChange[, 3]),
                 color = "grey50", alpha = 0.6, linewidth = 0.30) +
    geom_segment(aes(x = rep(2.01, d2), y = obsChange2[, 2], 
                     xend = rep(2.99, d2), yend = obsChange2[, 3]), 
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
  
  
}

#### Reindex based on the alpha diversity
reindexFnA <- function(alphaList, matClus){
  
  clusMat <- matrix(NA, ncol = 3, nrow = 90)
  
  for(t in 1:3){
    reindex_cluster <- data.frame(alp = alphaList[[t]], Cluster = matClus[, t]) %>%
      group_by(Cluster) %>%
      summarise(meanInt = mean(alp)) %>%
      arrange(-meanInt) %>%
      .$Cluster
    
    k <- length(reindex_cluster)
    for(i in 1:k){
      clusMat[which(matClus[, t] == reindex_cluster[i]), t] <- i
    }
    
  }
  
  clusMat
  
}

#### Bacteria of interest: c("Bifidobacterium", "Bacteroides", "Prevotella", "Streptococcus")
reindexFn(list(mo6plot$actualRelaGroup, mo8plot$actualRelaGroup, mo12plot$actualRelaGroup),
          salso_clus, Streptococcus) %>%
  obsMovement()

j <- 2
reindexFnA(list(alpha6mo$actualQ[[j]], alpha8mo$actualQ[[j]], alpha12mo$actualQ[[j]]), newClus) %>%
  obsMovement()

### ----------------------------------------------------------------------------

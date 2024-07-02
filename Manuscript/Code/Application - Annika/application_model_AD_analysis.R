### Load libraries: ------------------------------------------------------------
library(foreach)
library(doParallel)
library(stringr)
library(tidyverse)
library(salso)
library(ggplot2)
library(readxl)
library(latex2exp)
library(MCMCprecision)
library(gridExtra)
library(ggpubr)
library(ggh4x)

### User-defined Functions: ----------------------------------------------------
meanSD <- function(x, dplace = 3){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

uniqueClus <- function(x){
  length(unique(x))
}

### Import the data: -----------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

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

### Additional Data
addml <- read_excel(paste0(datpath, "Data/Mali_RB Metadata.xlsx"))
addni <- read_excel(paste0(datpath, "Data/Nicaragua_RB Metadata.xlsx"))

listData <- list(m6 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_6m_combined.RData")), 
                 m8 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_8m_combined.RData")), 
                 m12 = readRDS(paste0(path, "Manuscript/Result/microbiome_result_12m_combined.RData")))

### Data Pre-processing: -------------------------------------------------------
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

### Additional Dataset: --------------------------------------------------------
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

### Run the salso: -------------------------------------------------------------
salsoClus <- sapply(1:3, function(x){salso(listData[[x]]$ci_result)})

### Actual relative abundance for each observation: ----------------------------
actRela <- lapply(list(dat06[, -(1:5)], dat08[, -(1:5)], dat12[, -(1:5)]),
                  function(x){as.matrix(x/rowSums(x))})

### Estimate the relative abundance for each observation: ----------------------
set.seed(1)
estRela <- lapply(1:3, function(z){sapply(1:90, function(y){ci <- listData[[z]]$ci_result[, y] + 1
betaVec <- t(sapply(1:2000, function(x){listData[[z]]$beta_result[[x]][ci[x], ]}))
expBetaAR <- exp(betaVec) * t(listData[[z]]$atrisk_result[y, , ])
simRelaAbun <- t(sapply(1:2000, function(x){as.numeric(rdirichlet(1, expBetaAR[x, ]))}))
colMeans(simRelaAbun)}) %>%
    t()})

### Sort by some bacteria: -----------------------------------------------------
### Group Taxa
datList <- list(dat06, dat08, dat12)

registerDoParallel(5)
impTaxa <- foreach(t = 1:3) %dopar% {
  
  prop <- as.numeric(colMeans(datList[[t]][, -(1:5)]/rowSums(datList[[t]][, -(1:5)])))
  which(prop %in% boxplot.stats(prop)$out)
  
}
stopImplicitCluster()
impTaxa <- union(impTaxa[[1]], union(impTaxa[[2]], impTaxa[[3]]))

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

obsID <- c(str_replace_all(str_extract(dat08$ID[1:45], "^([:alpha:]{2}\\.){2}[:digit:]{2}"), "\\.", "\\-"),
           str_extract(dat08$ID[-(1:45)], "^[:digit:]{7}"))
bactList <- c("Bifidobacterium", "Other_Actinobacteria", "Bacteroides", 
              "Prevotella", "Other_Bacteroidetes", "Faecalibacterium", 
              "Megasphaera", "Ruminococcus", "Streptococcus", "Veillonella", 
              "Other_Firmicutes", "Others")
colours <- c("#A4C639", "#BAC394", "#DDA661", "#BD8B46", "#FFD197",
             "#BC6D62", "#DE6A5D", "#DE998B", "#FF6558", "#FFC6B8", 
             "#FFE4DA", "#E5E5E5")

### Group the bacteria
groupRelaTaxa <- lapply(1:3, function(x){
  
  estiRelaGroup <- matrix(NA, nrow = 90, ncol = 12)
  actualRelaGroup <- matrix(NA, nrow = 90, ncol = 12)
  
  for(i in 1:length(bactList)){
    index <- which(showName == bactList[i])
    if(length(index) == 1){
      estiRelaGroup[, i] <- estRela[[x]][, which(showName == bactList[i])]
      actualRelaGroup[, i] <- actRela[[x]][, which(showName == bactList[i])]
    } else {
      estiRelaGroup[, i] <- rowSums(estRela[[x]][, which(showName == bactList[i])])
      actualRelaGroup[, i] <- rowSums(actRela[[x]][, which(showName == bactList[i])])
    }
  }
  
  colnames(estiRelaGroup) <- bactList
  colnames(actualRelaGroup) <- bactList
  
  list(actualRelaGroup = actualRelaGroup, estiRelaGroup = estiRelaGroup)
  
})

### Rearrange the cluster assignment: ------------------------------------------
# colnames(groupRelaTaxa[[1]]$actualRelaGroup)
# [1] "Bifidobacterium"     
# [2] "Other_Actinobacteria"
# [3] "Bacteroides"         
# [4] "Prevotella"          
# [5] "Other_Bacteroidetes" 
# [6] "Faecalibacterium"    
# [7] "Megasphaera"         
# [8] "Ruminococcus"        
# [9] "Streptococcus"       
# [10] "Veillonella"         
# [11] "Other_Firmicutes"    
# [12] "Others"  

orderedClus <- sapply(1:3, function(x){
  newIndex <- cbind(y = groupRelaTaxa[[x]]$estiRelaGroup[, "Bifidobacterium"], 
                    clus = salsoClus[, x]) %>%
    as.data.frame() %>%
    group_by(clus) %>%
    summarise(mean = mean(y)) %>%
    .$mean %>%
    sort(decreasing = TRUE, index.return=TRUE) %>%
    .$ix
  newCLUS <- rep(NA, 90)
  for(i in 1:2){
    newCLUS[which(salsoClus[, x] == newIndex[i])] <- i
  }
  
  newCLUS
  
})

table(orderedClus[, 1])
table(orderedClus[, 2])
table(orderedClus[, 3])

### (Updated): Profiling Analysis: ---------------------------------------------
# datProf <- data.frame(Child_ID = obsID, clus = orderedClus[, 3], dat12[, c(1, 3, 5)]) %>%
#   left_join(addDat[which(addDat$Month == "12MO"), ]) %>%
#   mutate_all(~replace(., is.na(.), "not available"))
# 
# colnames(datProf)
# # [1] "Child_ID"              "clus"                 
# # [3] "country"               "Sex"                  
# # [5] "Group"                 "Month"                
# # [7] "Breastfed"             "Cow_milk"             
# # [9] "Fruits_Natural_Juices" "Vegetables"           
# # [11] "Soups"                 "Antibiotics"          
# # [13] "Weight"                "Height"               
# # [15] "Diarrhea"              "GR"                   
# # [17] "Protein"               "Nationality"   
# 
# profAnaPlot <- function()
# 
# ### Breast Milk
# bm <- datProf %>%
#   group_by(clus, Breastfed) %>%
#   summarise(n = n()) %>%
#   mutate(prop = n/sum(n)) %>%
#   ggplot(aes(x = factor(clus, labels = paste0("Cluster ", 1:2)), y = prop, fill = factor(Breastfed, levels = c("yes", "no", "not available"),
#                                                labels = c("Yes", "No", "Not Available")))) +
#   geom_bar(position = "stack", stat = "identity") +
#   scale_fill_manual(values = c("darkseagreen3", "coral", "grey90")) +
#   theme_minimal() +
#   theme(legend.position = "bottom", legend.title = element_blank()) +
#   scale_y_continuous(labels = scales::percent_format()) +
#   labs(x = "Cluster", y = "Percentage", title = "Breast Milk")

### Relative Plots: ------------------------------------------------------------
relaPlot <- function(group_relaAbunList, clus_assign, actual_month){
  
  ### Actual
  actual_title <- paste0(actual_month, " Infants: Actual Relative Abundance")
  actPlot <- data.frame(group_relaAbunList$actualRelaGroup) %>%
    mutate(id = obsID, Cluster = paste0("Cluster ", clus_assign)) %>%
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5.5), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Relative Abundance", title = actual_title, x = "Sample") +
    facet_grid(~ Cluster, scales = "free_x")
  
  ### Estimate
  estimate_title <- paste0(actual_month, " Infants: Estimated Relative Abundance")
  estPlot <- data.frame(group_relaAbunList$estiRelaGroup) %>%
    mutate(id = obsID, Cluster = paste0("Cluster ", clus_assign)) %>%
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5.5), legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "Estimated Relative Abundance", title = estimate_title, x = "Sample") +
    facet_grid(~ Cluster, scales = "free_x")
  
  list(actPlot = actPlot, estPlot = estPlot)
  
}

abr6m <- relaPlot(groupRelaTaxa[[1]], orderedClus[, 1], "6-Month")
abr8m <- relaPlot(groupRelaTaxa[[2]], orderedClus[, 2], "8-Month")
abr12m <- relaPlot(groupRelaTaxa[[3]], orderedClus[, 3], "12-Month")

ggarrange(abr6m$actPlot, abr8m$actPlot, abr12m$actPlot,
          abr6m$estPlot, abr8m$estPlot, abr12m$estPlot,
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")


### Descriptive of the relative bacteria: --------------------------------------
heatPlotAVG <- function(group_relaAbunList, cluster_assign, monthLab){
  
  k <- length(unique(cluster_assign))
  meanEst <- matrix(NA, nrow = k, ncol = 12)
  meanSDEst <- matrix(NA, nrow = k, ncol = 12)
  meanAct <- matrix(NA, nrow = k, ncol = 12)
  meanSDAct <- matrix(NA, nrow = k, ncol = 12)
  bactName <- colnames(group_relaAbunList$actualRelaGroup)
  
  for(i in 1:k){
    
    actPop <- data.frame(Cluster = cluster_assign, group_relaAbunList$actualRelaGroup) %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::select(-Cluster)
    meanAct[i, ] <- colMeans(actPop)
    meanSDAct[i, ] <- apply(actPop, 2, meanSD)
    
    EstPop <- data.frame(Cluster = cluster_assign, group_relaAbunList$estiRelaGroup) %>%
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
  
  list(actHeat = actHeat, estHeat = estHeat)
  
}

h6m <- heatPlotAVG(groupRelaTaxa[[1]], orderedClus[, 1], "6-Month")
h8m <- heatPlotAVG(groupRelaTaxa[[2]], orderedClus[, 2], "8-Month")
h12m <- heatPlotAVG(groupRelaTaxa[[3]], orderedClus[, 3], "12-Month")

ggarrange(h6m$actHeat, h8m$actHeat, h12m$actHeat, h6m$estHeat, h8m$estHeat, h12m$estHeat,
          ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")

### Trajectory Plot: -----------------------------------------------------------
obsMovement <- function(clusAssign){
  
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
    geom_text(data = c06, aes(x = rep(1, 90), y = index6, label = obs, color = clus6), size = 3.5) +
    geom_text(data = c08, aes(x = rep(2, 90), y = index8, label = obs, color = clus8), size = 3.5) +
    geom_text(data = c12, aes(x = rep(3, 90), y = index12, label = obs, color = clus12), size = 3.5) +
    geom_segment(aes(x = rep(1.01, d1), y = obsChange[, 2], 
                     xend = rep(1.99, d1), yend = obsChange[, 3]),
                 color = "grey50", alpha = 0.6, linewidth = 0.30) +
    geom_segment(aes(x = rep(2.01, d2), y = obsChange2[, 2], 
                     xend = rep(2.99, d2), yend = obsChange2[, 3]), 
                 color = "grey50", alpha = 0.6, linewidth = 0.30) +
    theme_minimal() + 
    theme(legend.position = "none", axis.ticks = element_blank(), axis.text.y = element_blank()) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("6 months", "8 months", "12 months")) +
    labs(x = "Timestamps", y = "", title = "The change of cluster behavior for each infants across three timestamps")
  
}

obsMovement(orderedClus)

### Descriptive Statistic Plot: ------------------------------------------------
descPlot <- function(descDat_type, intVar, capt, monthLab, g1level = "yes", g2level = "no", 
                     g1label = "Yes", g2label = "No", g1Color = "darkseagreen3", g2Color = "coral"){
  
   
  
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
    scale_fill_manual(capt, values = c("Yes" = g1Color, "No" = g2Color, "Not Available" = "grey90")) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
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

plot6mo <- plotData("6MO", "6-Month", dat06, orderedClus[, 1])
plot8mo <- plotData("8MO", "8-Month", dat08, orderedClus[, 2])
plot12mo <- plotData("12MO", "12-Month", dat12, orderedClus[, 3])

grid.arrange(grobs = plot6mo)
grid.arrange(grobs = plot8mo)
grid.arrange(grobs = plot12mo)

plot12mo[[10]]

ggarrange(plot6mo[[2]], plot8mo[[2]], plot12mo[[2]],
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

ggarrange(plot6mo[[6]], plot6mo[[7]], plot6mo[[8]], plot6mo[[12]],
          ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")

ggarrange(plot8mo[[6]], plot8mo[[7]], plot8mo[[8]], plot8mo[[12]],
          ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")

ggarrange(plot12mo[[6]], plot12mo[[7]], plot12mo[[8]], plot12mo[[12]],
          ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")

ggarrange(plot6mo[[1]], plot8mo[[1]], plot12mo[[1]],
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(plot6mo[[3]], plot8mo[[3]], plot12mo[[3]],
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(plot6mo[[4]], plot6mo[[5]], plot6mo[[9]], plot6mo[[10]],
          ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(plot8mo[[4]], plot8mo[[5]], plot8mo[[9]], plot8mo[[10]],
          ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(plot12mo[[4]], plot12mo[[5]], plot12mo[[9]], plot12mo[[10]],
          ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")

### Alpha-Diversity: -----------------------------------------------------------
ts_label <- c("6-Month", "8-Month", "12-Month")
alpMat <- lapply(1:3, function(x){rbind(data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Richness", stat = "Actual", val = as.numeric(rowSums(actRela[[x]] != 0))),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Richness", stat = "Estimated", val = rowSums(estRela[[x]] != 0)),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Shannon", stat = "Actual", val = as.numeric(-rowSums(ifelse(actRela[[x]] == 0, 1e-200, actRela[[x]]) * log(ifelse(actRela[[x]] == 0, 1e-200, actRela[[x]]))))),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Shannon", stat = "Estimated", val = as.numeric(-rowSums(ifelse(estRela[[x]] == 0, 1e-200, estRela[[x]]) * log(ifelse(estRela[[x]] == 0, 1e-200, estRela[[x]]))))),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Simpson", stat = "Actual", val = 1 - as.numeric(rowSums(actRela[[x]]^2))),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Simpson", stat = "Estimated", val = 1 - rowSums(estRela[[x]]^2)),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Inverse Simpson", stat = "Actual", val = 1/as.numeric(rowSums(actRela[[x]]^2))),
                                        data.frame(time = ts_label[x], clus = paste0("Cluster ", orderedClus[, x]), alp = "Inverse Simpson", stat = "Estimated", val = 1/rowSums(estRela[[x]]^2)))})
alpMat %>%
  bind_rows(.id = NULL) %>%
  ggplot(aes(x = stat, y = val, fill = clus)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "bottom", axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ggh4x::facet_grid2(factor(time, levels = paste0(c(6, 8, 12), "-Month")) ~ factor(alp, levels = c("Richness", "Shannon", "Simpson", "Inverse Simpson")), 
                     scales = "free_y", independent = "y") +
  labs(fill = "Cluster")

###: ---------------------------------------------------------------------------

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

### Average Relative Abundance: ------------------------------------------------
heatPlotAVG <- function(relaPlotObj, cluster_assign, monthLab){
  
  k <- length(unique(cluster_assign))
  meanEst <- matrix(NA, nrow = k, ncol = 12)
  meanAct <- matrix(NA, nrow = k, ncol = 12)
  bactName <- colnames(relaPlotObj$actualRelaGroup)
  
  for(i in 1:k){
    
    meanAct[i, ] <- data.frame(Cluster = cluster_assign, relaPlotObj$actualRelaGroup) %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::select(-Cluster) %>%
      colMeans()
    
    meanEst[i, ] <- data.frame(Cluster = cluster_assign, relaPlotObj$estiRelaGroup) %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::select(-Cluster) %>%
      colMeans()
    
  }
  
  actTitle <- paste0(monthLab, ": Actual Average relative abundance")
  estTitle <- paste0(monthLab, ": Estimated Average relative abundance")
  
  actHeat <- t(meanAct) %>%
    `colnames<-`(paste0("Cluster_", 1:k)) %>%
    data.frame(Bact = bactName) %>%
    pivot_longer(!Bact) %>%
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
                                                      "Others"))), fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 3)), color = "black") +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = " ", y = "Bacteria", title = actTitle, 
         fill = "Average Relative Abundance")
  
  estHeat <- t(meanEst) %>%
    `colnames<-`(paste0("Cluster_", 1:k)) %>%
    data.frame(Bact = bactName) %>%
    pivot_longer(!Bact) %>%
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
                                                      "Others"))), fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 3)), color = "black") +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = " ", y = "Bacteria", title = estTitle, 
         fill = "Average Relative Abundance")
  
  list(actHeat, estHeat)
  
}

grid.arrange(grobs = heatPlotAVG(mo6plot, salso_clus[, 1], "6-Month"))
grid.arrange(grobs = heatPlotAVG(mo8plot, salso_clus[, 2], "8-Month"))
grid.arrange(grobs = heatPlotAVG(mo12plot, salso_clus[, 3], "12-Month"))

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
                   descPlot(descDat, country, "Nationality", monthLabel, "NI", "ML", "Nicaraguan", "Malian", "darkblue", "green"),
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
  
  grid.arrange(grobs = plotList)
  
}

plotData("6MO", "6-Month", dat06, salso_clus[, 1])
plotData("8MO", "8-Month", dat08, salso_clus[, 2])
plotData("12MO", "12-Month", dat12, salso_clus[, 3])



### ----------------------------------------------------------------------------

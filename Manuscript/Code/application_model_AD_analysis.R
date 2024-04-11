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

### Estimate the relative abundance for each observation: ----------------------
set.seed(1)
ci <- listData[[1]]$ci_result[, 1] + 1
betaVec <- t(sapply(1:2000, function(x){listData[[1]]$beta_result[[x]][ci[x], ]}))
expBetaAR <- exp(betaVec) * t(listData[[1]]$atrisk_result[1, , ])
simRelaAbun <- t(sapply(1:2000, function(x){as.numeric(rdirichlet(1, expBetaAR[x, ]))}))
cbind(colMeans(simRelaAbun), as.numeric(dat06[1, -(1:5)])/sum(as.numeric(dat06[1, -(1:5)])))

### Trajectory Plot: -----------------------------------------------------------

obsMovement <- function(clusAssign){
  
  ### Get the bar separating each cluster
  bar1 <- data.frame(x = rep(1, 2), y = 90 - as.numeric(cumsum(table(clusAssign[, 1]))))
  bar2 <- data.frame(x = rep(2, 2), y = 90 - as.numeric(cumsum(table(clusAssign[, 2]))))
  bar3 <- data.frame(x = rep(3, 2), y = 90 - as.numeric(cumsum(table(clusAssign[, 3]))))
  
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

obsMovement(salsoClus)





###: ---------------------------------------------------------------------------

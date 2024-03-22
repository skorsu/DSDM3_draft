### Required Library
library(readxl)
library(stringr)
library(tidyverse)
library(lubridate)
library(VennDiagram)

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

datpath <- "/Users/kevin-imac/Desktop/Annika/"
if(! file.exists(path)){
  datpath <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}

### Additional Data
addml <- read_excel(paste0(datpath, "Data/Mali_RB Metadata.xlsx"))
addni <- read_excel(paste0(datpath, "Data/Nicaragua_RB Metadata.xlsx"))

### Result
annikaZZ <- readRDS(paste0(path, "Result/microbiome_result.RData"))

### Data: 6 and 8 Months
ni68 <- read.csv(paste0(datpath, "Data/Nicaragua_6mo_8mo_genus.csv"))
ml68 <- read.csv(paste0(datpath, "Data/Mali_6mo_8mo_genus.csv"))

### Data: 12 Months
ni12 <- read.csv(paste0(datpath, "Data/Nicaragua_12mo_Metadata_csv.csv"))
ml12 <- read.csv(paste0(datpath, "Data/Mali_12mo_Metadata_csv.csv"))

## Retrieve the infant ID in the analysis --------------------------------------
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

## Cluster Assignment ----------------------------------------------------------
clusBinder <- sapply(1:3, 
                     function(x){as.numeric(salso(annikaZZ[[x]]$result$ci_result[-(1:500), ],
                                                  loss = "binder"))})
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

clusObs <- data.frame(ChildID = gsub("\\.", "-", obsID), clusBinderN)

### Interested Aspect ----------------------------------------------------------
mlVar <- str_replace_all(colnames(addml), "\\_", " ")[-c(1:6, 31, 32)]
ngVar <- str_replace_all(str_extract(colnames(addni), "\\S*(?=\\_[:digit:])"), "\\_", " ")[-(1:3)]

vennCol <- venn.diagram(list(Malian = mlVar, Nicaraguan = ngVar), filename = NULL)

vennCol[[5]]$label  <- paste(setdiff(mlVar, ngVar), collapse="\n")  
vennCol[[6]]$label <- paste(setdiff(ngVar, mlVar)  , collapse="\n")  
vennCol[[7]]$label <- paste(intersect(mlVar, ngVar), collapse="\n")

grid.newpage()
grid.draw(vennCol)

### Transform the variables ----------------------------------------------------
#### Diarrhea
addml %>%
  mutate(across(where(is.character), tolower)) %>%
  dplyr::select(Today_Date, Diarrhea_Episode_Date_Start, Diarrhea_Episode_Date_End, Diarrhea_Episode) %>%
  mutate(Diarrhea = ifelse(Diarrhea_Episode == "yes", 
                           ifelse(difftime(Today_Date, Diarrhea_Episode_Date_End) <= 14, "yes", "no"), 
                           Diarrhea_Episode)) %>%
  view()


addni %>%
  dplyr::select(Diarrhea_Episode_Past_14_Days_3)

#### GR
##### ML: rice, maize, cowpea, porridge
##### NI: rice cereal, Gallo Pinto, Brown Rice

colnames(addml)
addml %>%
  mutate(across(where(is.character), tolower)) %>%
  dplyr::select(Rice, Maize, Cowpea, Porridge) %>%
  rowwise %>% 
  mutate(null_message_GR = as.integer(all(c(Rice, Maize, Cowpea, Porridge) %in% NA)) > 0) %>%
  mutate(allNo_GR = all(c(Rice, Maize, Cowpea, Porridge) %in% "no")) %>%
  mutate(GR = ifelse(null_message_GR == TRUE, NA, ifelse(allNo_GR == TRUE, "no", "yes"))) %>%
  view()

colnames(addni)

table(factor(addni$Rice_cereal_5.3))
table(factor(addni$Gallo_Pinto_5.7))
table(factor(addni$Brown_Rice_RB_5.11))

addni %>%
  mutate(across(where(is.character), tolower)) %>% 
  dplyr::select(Rice_cereal_5.3, Gallo_Pinto_5.7) %>%
  rowwise %>% 
  mutate(null_message_GR= as.integer(all(c(Rice_cereal_5.3, Gallo_Pinto_5.7) %in% NA)) > 0) %>%
  mutate(allNo_GR = all(c(Rice_cereal_5.3, Gallo_Pinto_5.7) %in% c("used to consume, but not now", "never"))) %>%
  mutate(GR = ifelse(null_message_GR == TRUE, NA, ifelse(allNo_GR == TRUE, "no", "yes"))) %>%
  view()

#### Protein
##### ML: Chicken/fish/meat, Eggs
##### NI: Chicken, Cheese, and Yogurt

colnames(addml)
addml %>%
  mutate(across(where(is.character), tolower)) %>%
  dplyr::select(`Chicken / fish / meat`, Eggs) %>%
  rowwise %>% 
  mutate(null_message_PR = as.integer(all(c(`Chicken / fish / meat`, Eggs) %in% NA)) > 0) %>%
  mutate(allNo_PR = all(c(`Chicken / fish / meat`, Eggs) %in% "no")) %>%
  mutate(Protein = ifelse(null_message_PR == TRUE, NA, ifelse(allNo_PR == TRUE, "no", "yes"))) %>%
  view()

colnames(addni)

table(factor(addni$Chicken_5.6))
table(factor(addni$Cheese_5.8))
table(factor(addni$Yogurt_5.10))

addni %>%
  mutate(across(where(is.character), tolower)) %>% 
  dplyr::select(Chicken_5.6, Cheese_5.8, Yogurt_5.10) %>%
  rowwise %>% 
  mutate(null_message_PR = as.integer(all(c(Chicken_5.6, Cheese_5.8, Yogurt_5.10) %in% NA)) > 0) %>%
  mutate(allNo_PR = all(c(Chicken_5.6, Cheese_5.8, Yogurt_5.10) %in% c("used to consume, but not now", "never"))) %>%
  mutate(Protein = ifelse(null_message_PR == TRUE, NA, ifelse(allNo_PR == TRUE, "no", "yes"))) %>%
  view()

#### Combine the dataset -------------------------------------------------------
rbind(addml %>%
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
        mutate(Nationality = "NI")) -> addInfo

view(addInfo)

### Check the observations -----------------------------------------------------
#### Note: ML ID is number xxxxxxx
nat <- c(rep("Nicaraguan", 45), rep("Malian", 45))
check6M <- str_replace_all(obsID, "\\.", "\\-") %in% addInfo$Child_ID[addInfo$Month == "6MO"]
check8M <- str_replace_all(obsID, "\\.", "\\-") %in% addInfo$Child_ID[addInfo$Month == "8MO"]
check12M <- str_replace_all(obsID, "\\.", "\\-") %in% addInfo$Child_ID[addInfo$Month == "12MO"]

data.frame(Nationality = nat, ID = str_replace_all(obsID, "\\.", "\\-"),
           check6M, check8M, check12M) %>%
  pivot_longer(cols = c(check6M, check8M, check12M)) %>%
  mutate(Month = factor(name, levels = c("check12M", "check8M", "check6M"),
                        labels = c("12 Month", "8 Month", "6 Month")),
         Exists = ifelse(value, 1, 0)) %>%
  ggplot(aes(x = ID, y = Month)) + 
  geom_tile(aes(alpha = Exists, fill = Nationality), colour = "white") + 
  theme_minimal() + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Join the additional info with cluster information --------------------------
allClus <- rbind(left_join(clusObs %>%
                             dplyr::select(Child_ID = ChildID, cluster = X1) %>%
                             mutate(Month = "6MO"),
                           addInfo[addInfo$Month == "6MO", ] %>%
                             dplyr::select(-Month)),
                 left_join(clusObs %>%
                             dplyr::select(Child_ID = ChildID, cluster = X2) %>%
                             mutate(Month = "8MO"),
                           addInfo[addInfo$Month == "8MO", ] %>%
                             dplyr::select(-Month)),
                 left_join(clusObs %>%
                             dplyr::select(Child_ID = ChildID, cluster = X3) %>%
                             mutate(Month = "12MO"),
                           addInfo[addInfo$Month == "12MO", ] %>%
                             dplyr::select(-Month)))
view(allClus)

### One additional variable at a time ------------------------------------------
oneAddPlot <- function(intVar, capt, 
                       yesColor = "olivedrab3", noColor = "rosybrown1"){
  
  allClus %>%
    replace(is.na(.), "N/A") %>%
    group_by(Month, cluster, {{intVar}}) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = cluster, 
               y = freq, fill = factor({{intVar}}, levels = c("yes", "no", "N/A"),
                                       label = c("Yes", "No", "Not Available")))) +
    geom_bar(position = "stack", stat = "identity") +
    facet_wrap(. ~ factor(Month, levels = c("6MO", "8MO", "12MO"), 
                          label = c("6 Month", "8 Month", "12 Month")),
               scales = "free_x") +
    scale_fill_manual(capt, values = c(yesColor, noColor, "grey90")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    labs(x = "Cluster", y = "Proportion", title = paste0(capt, " by Cluster"))

}

colnames(allClus)

oneAddPlot(Breastfed, "Breast Milk")
oneAddPlot(Cow_milk, "Cow Milk")
oneAddPlot(Fruits_Natural_Juices, "Natural Juice")
oneAddPlot(Vegetables, "Vegetables")
oneAddPlot(Breastfed, "Breast Milk")
oneAddPlot(Soups, "Soups")
oneAddPlot(Antibiotics, "Antibiotics")
oneAddPlot(Diarrhea, "Diarrhea")
oneAddPlot(GR, "Grains and Legumes")
oneAddPlot(Protein, "Protein")

ggplot(allClus, aes(x = Weight, y = Height, color = factor(cluster, 
                                                           labels = paste0("Cluster ", 1:5)))) +
  geom_point() +
  facet_wrap(. ~ factor(Month, levels = c("6MO", "8MO", "12MO"), 
                        label = c("6 Month", "8 Month", "12 Month")),
             scales = "free_x") +
  theme_bw() +
  labs(title = "Weight and Height of the infants in each cluster") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values=c("purple", "#E69F00", "#56B4E9", "green", "pink"))
  
### Alpha Diversity ------------------------------------------------------------
alpCalc <- function(dat0, clus0, mo, treh = 1e-323){
  
  
  result <- as.data.frame(matrix(NA, ncol = 7, nrow = 90))
  result[, 1] <- str_replace_all(obsID, "\\.", "\\-")
  result[, 2] <- clus0
  
  #### Richness
  result[, 3] <- apply(dat0[, -(1:5)], 1, function(x){sum(x != 0)})
  
  #### Transform to the relative abundance
  rTaxa <- dat0[, -(1:5)]/rowSums(dat0[, -(1:5)])
  
  #### Simpson Index
  result[, 4] <- 1 - rowSums(rTaxa^2)
  
  #### Inverse Simpson Index
  result[, 5] <- 1/rowSums(rTaxa^2)
  
  #### Shannon
  rTaxa[rTaxa == 0] <- treh ### The lowest possible
  result[, 6] <- -rowSums(rTaxa * log(rTaxa))
  result[, 7] <- mo
  
  colnames(result) <- c("Child_ID", "Cluster", "Richness", "Simpson_Index",
                        "Inverse_Simpson_Index", "Shannon", "Month")
  
  result %>%
    mutate(Cluster = factor(paste0("Cluster ", Cluster)))

}

alpPlot <- function(clusIndex, threshold, alQuan, alDesc){
  
  rbind(alpCalc(dat06, clusIndex[, 1], "6 Month", treh = threshold),
        alpCalc(dat08, clusIndex[, 2], "8 Month", treh = threshold),
        alpCalc(dat12, clusIndex[, 3], "12 Month", treh = threshold)) %>%
  ggplot(aes(x = Cluster, y = {{alQuan}}, group = Cluster)) +
  geom_boxplot() +
  facet_wrap(. ~ factor(Month, levels = c("6 Month", "8 Month", "12 Month"), 
                        label = c("6 Month", "8 Month", "12 Month")),
             scales = "free_x") +
  theme_bw() +
  labs(x = "", y = "", title = alDesc)
}

alpPlot(clusBinderN, 1e-323, Richness, "Richness")
alpPlot(clusBinderN, 1e-323, Simpson_Index, "Simpson Index")
alpPlot(clusBinderN, 1e-323, Inverse_Simpson_Index, "Inverse Simpson Index")
alpPlot(clusBinderN, 1e-323, Shannon, "Shannon")


rbind(alpCalc(dat06, clusBinderN[, 1], "6 Month"),
      alpCalc(dat08, clusBinderN[, 2], "8 Month"),
      alpCalc(dat12, clusBinderN[, 3], "12 Month")) %>%
  dplyr::select(-Child_ID) %>%
  group_by(Month, Cluster) %>%
  summarise(across(everything(), mean)) %>%
  arrange(Month)

### ----------------------------------------------------------------------------

### Visualization --------------------------------------------------------------
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

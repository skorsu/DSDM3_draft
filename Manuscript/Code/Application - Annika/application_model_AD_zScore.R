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
library(MVTests)

### User-defined Functions: ----------------------------------------------------
meanSD <- function(x, dplace = 3){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

uniqueClus <- function(x){
  length(unique(x))
}

multiplesheets <- function(fname) { 
  
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  
  # assigning names to data frames 
  names(data_frame) <- sheets 
  
  # print data frame 
  print(data_frame) 
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

### Growth Velocity: -----------------------------------------------------------
ni_GV <- read.csv(paste0(datpath, "Data/Nicaragua_GV.csv"))
ni_GV[ni_GV == "-"] <- NA

niGV_8m <- ni_GV %>%
  filter(Study.ID. %in% (str_extract(dat06$ID[1:45], "[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %>%
  str_replace_all("\\.", "\\-"))) %>%
  filter(Age..month. == 8) %>%
  dplyr::select(-Length.velocity.z.score..from.8.months., -Weight.velocity.z.score..from.8.months.,
                -Diet.Group, -Age..month., -Sex) %>%
  rename(Child_ID = Study.ID., wz_6m_8m = Weight.velocity.z.score..from.6.months.,
         lz_6m_8m = Length.velocity.z.score..from.6.months.)

niGV_12m <- ni_GV %>%
  filter(Study.ID. %in% (str_extract(dat06$ID[1:45], "[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %>%
                           str_replace_all("\\.", "\\-"))) %>%
  filter(Age..month. == 12) %>%
  dplyr::select(-Diet.Group, -Age..month., -Sex) %>%
  rename(Child_ID = Study.ID., wz_6m_12m = Weight.velocity.z.score..from.6.months.,
         lz_6m_12m = Length.velocity.z.score..from.6.months., 
         wz_8m_12m = Weight.velocity.z.score..from.8.months.,
         lz_8m_12m = Length.velocity.z.score..from.8.months.)

inner_join(niGV_8m, niGV_12m)


ni_ML <- multiplesheets(paste0(datpath, "Data/Mali_GV.xlsx"))
niID <- intersect(ni_ML$`6 Month`$Child_ID, ni_ML$`8 Month`$Child_ID) %>%
  intersect(ni_ML$`12 Month`$Child_ID) %>%
  intersect(as.numeric(str_extract(dat06$ID[-(1:45)], "[:digit:]{7}")))
colnames(ni_ML$`6 Month`)[2:3] <- c("l6", "w6")
colnames(ni_ML$`8 Month`)[2:3] <- c("l8", "w8")
colnames(ni_ML$`12 Month`)[2:3] <- c("l12", "w12")

datGV <- rbind(inner_join(niGV_8m, niGV_12m),
               inner_join(ni_ML$`6 Month`, ni_ML$`8 Month`) %>%
                 inner_join(ni_ML$`12 Month`) %>%
                 as.data.frame() %>%
                 filter(Child_ID %in% niID) %>%
                 transmute(Child_ID, wz_6m_8m = w8 - w6, lz_6m_8m = l8 - l6,
                           wz_6m_12m = w12 - w6, lz_6m_12m = l12 - l6,
                           wz_8m_12m = w12 - w8, lz_8m_12m = l12 - l8))

datGV <- as.data.frame(datGV)

### Run the salso: -------------------------------------------------------------
salsoClus <- sapply(1:3, function(x){salso(listData[[x]]$ci_result)})
salsoClus <- abs(salsoClus - 3)

datGV_clus <- data.frame(Child_ID = c(str_extract(dat06$ID[1:45], "[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %>%
                                        str_replace_all("\\.", "\\-"),
                                      str_extract(dat06$ID[-(1:45)], "[:digit:]{7}")), salsoClus) %>%
  left_join(datGV)

datGV_clus <- datGV_clus %>%
  mutate(Cluster_M8 = paste0("Cluster ", X2), Cluster_M12 = paste0("Cluster ", X3))

### 8M to 6M: ------------------------------------------------------------------
dat8m <- datGV_clus %>%
  dplyr::select(Cluster_M8, lz_6m_8m, wz_6m_8m) %>%
  pivot_longer(!Cluster_M8)

dat8m$name <- factor(dat8m$name, labels = c("Length from 6 month", "Weight from 6 month"))

g1 <- dat8m %>%
  filter(Cluster_M8 == "Cluster 1" & name == "Length from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in length of each infant between 6 months and 8 months. (Cluster 1)")

g2 <- dat8m %>%
  filter(Cluster_M8 == "Cluster 2" & name == "Length from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in length of each infant between 6 months and 8 months. (Cluster 2)")

g3 <- dat8m %>%
  filter(Cluster_M8 == "Cluster 1" & name == "Weight from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in weight of each infant between 6 months and 8 months. (Cluster 1)")

g4 <- dat8m %>%
  filter(Cluster_M8 == "Cluster 2" & name == "Weight from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in weight of each infant between 6 months and 8 months. (Cluster 2)")

grid.arrange(g1, g2, g3, g4)

### Length
m8_l6 <- dat8m %>%
  filter(name == "Length from 6 month") %>%
  dplyr::select(-name) %>%
  transmute(group = factor(Cluster_M8), value = as.numeric(value))
  
var.test(value ~ group, m8_l6) ### Equal Variance Assumed
t.test(value ~ group, m8_l6, var.equal = TRUE) ### p-value is not less than 0.05 (Same)
wilcox.test(value ~ group, m8_l6, exact = FALSE, conf.int = TRUE) ### FTR - 

m8_l6 %>%
  group_by(group) %>%
  summarise(mn = mean(value, na.rm = TRUE), md = median(value, na.rm = TRUE))

### Weight
m8_w6 <- dat8m %>%
  filter(name == "Weight from 6 month") %>%
  dplyr::select(-name) %>%
  transmute(group = factor(Cluster_M8), value = as.numeric(value))

var.test(value ~ group, m8_w6) ### Equal Variance NOT Assumed
t.test(value ~ group, m8_w6, var.equal = FALSE) ### p-value is less than 0.05 (Difference)
wilcox.test(value ~ group, m8_w6, exact = FALSE, conf.int = TRUE)

m8_w6 %>%
  group_by(group) %>%
  summarise(mn = mean(value, na.rm = TRUE), md = median(value, na.rm = TRUE))

### Length and Weight
d8 <- datGV_clus %>%
  transmute(Cluster_M8, lz_6m_8m = as.numeric(lz_6m_8m), wz_6m_8m = as.numeric(wz_6m_8m)) %>%
  drop_na()

summary(BoxM(d8[, 2:3], d8[, 1])) ### Reject - Covariance matrix are not HOMO.
TwoSamplesHT2(d8[, 2:3], ifelse(d8[, 1] == "Cluster 2", 2, 1), alpha = 0.05, Homogenity = FALSE)

### Plots
b1 <- ggplot(dat8m, aes(x = Cluster_M8, y = as.numeric(value))) +
  geom_boxplot() +
  facet_grid(. ~ name) +
  theme_bw() + 
  labs(x = "Cluster (Timestamp: 8-Month)", y = "Difference in z-score",
       title = "The boxplot of the difference in z-scores for length and weight of individual infants, comparing 8-month to 6-month measurements.")

s1 <- datGV_clus %>%
  dplyr::transmute(Cluster = Cluster_M8, Length = as.numeric(lz_6m_8m), 
                   Weight = as.numeric(wz_6m_8m)) %>%
  ggplot(aes(x = Weight, y = Length, color = Cluster)) +
  geom_point() +
  theme_bw() +
  labs(title = "The scatterplot of the difference in z-scores for length and weight of individual infants, comparing 8-month to 6-month measurements.")

### 12M to 6M: ------------------------------------------------------------------
dat12m <- datGV_clus %>%
  dplyr::select(Cluster_M12, lz_6m_12m, wz_6m_12m) %>%
  pivot_longer(!Cluster_M12)

dat12m$name <- factor(dat12m$name, labels = c("Length from 6 month", "Weight from 6 month"))

g1 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 1" & name == "Length from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in length of each infant between 6 months and 12 months. (Cluster 1)")

g2 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 2" & name == "Length from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in length of each infant between 6 months and 12 months. (Cluster 2)")

g3 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 1" & name == "Weight from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in weight of each infant between 6 months and 12 months. (Cluster 1)")

g4 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 2" & name == "Weight from 6 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in weight of each infant between 6 months and 12 months. (Cluster 2)")

grid.arrange(g1, g2, g3, g4)

### Length
m12_l6 <- dat12m %>%
  filter(name == "Length from 6 month") %>%
  dplyr::select(-name) %>%
  transmute(group = factor(Cluster_M12), value = as.numeric(value))

var.test(value ~ group, m12_l6) ### Equal Variance NOT Assumed
t.test(value ~ group, m12_l6, var.equal = FALSE) ### p-value is not less than 0.05 (Same)
wilcox.test(value ~ group, m12_l6, exact = FALSE, conf.int = TRUE)

m12_l6 %>%
  group_by(group) %>%
  summarise(mn = mean(value, na.rm = TRUE), md = median(value, na.rm = TRUE))

### Weight
m12_w6 <- dat12m %>%
  filter(name == "Weight from 6 month") %>%
  dplyr::select(-name) %>%
  transmute(group = factor(Cluster_M12), value = as.numeric(value))

var.test(value ~ group, m12_w6) ### Equal Variance Assumed
t.test(value ~ group, m12_w6, var.equal = TRUE) ### p-value is less than 0.05 (Difference)
wilcox.test(value ~ group, m12_w6, exact = FALSE, conf.int = TRUE)

m12_w6 %>%
  group_by(group) %>%
  summarise(mn = mean(value, na.rm = TRUE), md = median(value, na.rm = TRUE))

### Length and Weight
d12m6 <- datGV_clus %>%
  transmute(Cluster_M12, lz_6m_12m = as.numeric(lz_6m_12m), wz_6m_12m = as.numeric(wz_6m_12m)) %>%
  drop_na()

summary(BoxM(d12m6[, 2:3], d12m6[, 1])) ### Reject - Covariance matrix are not HOMO.
TwoSamplesHT2(d12m6[, 2:3], ifelse(d12m6[, 1] == "Cluster 2", 2, 1), alpha = 0.05, Homogenity = FALSE)

### Plots
b2 <- ggplot(dat12m, aes(x = Cluster_M12, y = as.numeric(value))) +
  geom_boxplot() +
  facet_grid(. ~ name) +
  theme_bw() + 
  labs(x = "Cluster (Timestamp: 12-Month)", y = "Difference in z-score",
       title = "The boxplot of the difference in z-scores for length and weight of individual infants, comparing 12-month to 6-month measurements.")

s2 <- datGV_clus %>%
  dplyr::transmute(Cluster = Cluster_M12, Length = as.numeric(lz_6m_12m), 
                   Weight = as.numeric(wz_6m_12m)) %>%
  ggplot(aes(x = Weight, y = Length, color = Cluster)) +
  geom_point() +
  theme_bw() +
  labs(title = "The scatterplot of the difference in z-scores for length and weight of individual infants, comparing 12-month to 6-month measurements.")

### 12M to 8M: ------------------------------------------------------------------
dat12m <- datGV_clus %>%
  dplyr::select(Cluster_M12, lz_8m_12m, wz_8m_12m) %>%
  pivot_longer(!Cluster_M12)

dat12m$name <- factor(dat12m$name, labels = c("Length from 8 month", "Weight from 8 month"))

g1 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 1" & name == "Length from 8 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in length of each infant between 8 months and 12 months. (Cluster 1)")

g2 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 2" & name == "Length from 8 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in length of each infant between 8 months and 12 months. (Cluster 2)")

g3 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 1" & name == "Weight from 8 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in weight of each infant between 8 months and 12 months. (Cluster 1)")

g4 <- dat12m %>%
  filter(Cluster_M12 == "Cluster 2" & name == "Weight from 8 month") %>%
  .$value %>%
  as.numeric() %>%
  ggqqplot(conf.int = FALSE) +
  theme_bw() +
  labs(title = "QQ-Plot: The difference in weight of each infant between 8 months and 12 months. (Cluster 2)")

grid.arrange(g1, g2, g3, g4)

### Length
m12_l8 <- dat12m %>%
  filter(name == "Length from 8 month") %>%
  dplyr::select(-name) %>%
  transmute(group = factor(Cluster_M12), value = as.numeric(value))

var.test(value ~ group, m12_l8) ### Equal Variance Assumed
t.test(value ~ group, m12_l8, var.equal = TRUE) ### p-value is less than 0.05 (different)
wilcox.test(value ~ group, m12_l8, exact = FALSE, conf.int = TRUE)

m12_l8 %>%
  group_by(group) %>%
  summarise(mn = mean(value, na.rm = TRUE), md = median(value, na.rm = TRUE))

### Weight
m12_w8 <- dat12m %>%
  filter(name == "Weight from 8 month") %>%
  dplyr::select(-name) %>%
  transmute(group = factor(Cluster_M12), value = as.numeric(value))

var.test(value ~ group, m12_w8) ### Equal Variance Assumed
t.test(value ~ group, m12_w8, var.equal = TRUE) ### p-value is less than 0.05 (Difference)
wilcox.test(value ~ group, m12_w8, exact = FALSE, conf.int = TRUE)

m12_w8 %>%
  group_by(group) %>%
  summarise(mn = mean(value, na.rm = TRUE), md = median(value, na.rm = TRUE))

### Length and Weight
d12m8 <- datGV_clus %>%
  transmute(Cluster_M12, lz_8m_12m = as.numeric(lz_8m_12m), wz_8m_12m = as.numeric(wz_8m_12m)) %>%
  drop_na()

summary(BoxM(d12m8[, 2:3], d12m8[, 1])) ### FTR - Covariance matrix are HOMO.
TwoSamplesHT2(d12m8[, 2:3], ifelse(d12m8[, 1] == "Cluster 2", 2, 1), alpha = 0.05, Homogenity = TRUE)

### Plots
b3 <- ggplot(dat12m, aes(x = Cluster_M12, y = as.numeric(value))) +
  geom_boxplot() +
  facet_grid(. ~ name) +
  theme_bw() + 
  labs(x = "Cluster (Timestamp: 12-Month)", y = "Difference in z-score",
       title = "The boxplot of the difference in z-scores for length and weight of individual infants, comparing 12-month to 8-month measurements.")

s3 <- datGV_clus %>%
  dplyr::transmute(Cluster = Cluster_M12, Length = as.numeric(lz_8m_12m), 
                   Weight = as.numeric(wz_8m_12m)) %>%
  ggplot(aes(x = Weight, y = Length, color = Cluster)) +
  geom_point() +
  theme_bw() +
  labs(title = "The scatterplot of the difference in z-scores for length and weight of individual infants, comparing 12-month to 6-month measurements.")
  

grid.arrange(b1, b2, b3)
grid.arrange(s1, s2, s3)



###: ---------------------------------------------------------------------------

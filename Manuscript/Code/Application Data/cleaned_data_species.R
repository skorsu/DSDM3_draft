# devtools::install_github("YushuShi/MicrobiomeCluster")

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(foreach)
library(doParallel)

# Settings and Path: -----------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

# Clean the data: --------------------------------------------------------------
## Vincent 2013, CDI: ----------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/cdi_vincent_v3v5_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/cdi_vincent_v3v5_results.tar"))

### Check the metadata
metData <- read.table(paste0(path, "cdi_vincent_v3v5_results/cdi_vincent_v3v5.metadata.txt"), sep = "\t", header = TRUE) 

### OTU Table
otuTab <- read.table(paste0(path, "cdi_vincent_v3v5_results/RDP/cdi_vincent_v3v5.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]

### Check the sample
str_extract(metData$X, "[:alpha:]+[:digit:]+\\.[:alpha:]+[:digit:]+") %in% str_extract(rownames(otuTab), "[:alpha:]+[:digit:]+\\.[:alpha:]+[:digit:]+")
str_extract(rownames(otuTab), "[:alpha:]+[:digit:]+\\.[:alpha:]+[:digit:]+") %in% str_extract(metData$X, "[:alpha:]+[:digit:]+\\.[:alpha:]+[:digit:]+")

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/cdi_vincent_v3v5_cleaned_species.rds"))

## Ross 2015, Obesity: ---------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/ob_ross_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/ob_ross_results.tar"))

### Check the metadata
metData <- read.table(paste0(path, "ob_ross_results/ob_ross.metadata.txt"), sep = "\t", header = TRUE) 
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "ob_ross_results/RDP/ob_ross.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)

### Check the sample
sum(rownames(otuTab) %in% metData$sampleID) 
sum(metData$sampleID %in% rownames(otuTab)) 

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/ob_ross_cleaned_species.rds"))

## Chen 2012, CRC: -------------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/crc_xiang_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/crc_xiang_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "crc_xiang_results/crc_xiang.metadata.txt")) 
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "crc_xiang_results/RDP/crc_xiang.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

### Check the sample
sum(metData$X.SampleID %in% rownames(otuTab))
sum(rownames(otuTab) %in% metData$X.SampleID)

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/crc_xiang_cleaned_species.rds"))

### Wang 2012, CRC: --------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/crc_zhao_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/crc_zhao_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "crc_zhao_results/crc_zhao.metadata.txt")) 
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "crc_zhao_results/RDP/crc_zhao.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

### Check the sample
sum(metData$X.SampleID %in% rownames(otuTab))
sum(rownames(otuTab) %in% metData$X.SampleID)

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/crc_zhao_cleaned_species.rds"))

## Singh 2015, EDD: ------------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/edd_singh_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/edd_singh_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "edd_singh_results/edd_singh.metadata.txt")) 
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "edd_singh_results/RDP/edd_singh.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

### Check the sample
metData <- metData[which(metData$SampleID %in% intersect(metData$SampleID, rownames(otuTab))), ]
otuTab <- otuTab[which(rownames(otuTab) %in% intersect(metData$SampleID, rownames(otuTab))), ]

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/edd_singh_cleaned_species.rds"))

## Zhang 2013, LIV: ------------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/mhe_zhang_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/mhe_zhang_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "mhe_zhang_results/mhe_zhang.metadata.txt")) 
dim(metData)
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "mhe_zhang_results/RDP/mhe_zhang.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

### Check the sample
sum(metData$Sample_Name_s %in% rownames(otuTab))
sum(rownames(otuTab) %in% metData$Sample_Name_s)

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/mhe_zhang_cleaned_species.rds"))

## Schubert 2014, CDI: ---------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/cdi_schubert_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/cdi_schubert_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "cdi_schubert_results/cdi_schubert.metadata.txt")) 
dim(metData)
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "cdi_schubert_results/RDP/cdi_schubert.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

### Check the sample
otuTab <- otuTab[which(rownames(otuTab) %in% intersect(rownames(otuTab), metData$sample_id)), ]
dim(otuTab)

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/cdi_schubert_cleaned_species.rds"))

## Papa 2012, IBD: -------------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/ibd_alm_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/ibd_alm_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "ibd_alm_results/ibd_alm.metadata.txt")) 
dim(metData)
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "ibd_alm_results/RDP/ibd_alm.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

sum(str_extract(rownames(otuTab), "[:digit:]+[:alpha:]") %in% metData$X)
sum(metData$X %in% str_extract(rownames(otuTab), "[:digit:]+[:alpha:]"))

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/ibd_alm_cleaned_species.rds"))

## Morgan 2012, IBD: -----------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/ibd_huttenhower_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/ibd_huttenhower_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "ibd_huttenhower_results/ibd_huttenhower.metadata.txt")) 
dim(metData)
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "ibd_huttenhower_results/RDP/ibd_huttenhower.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

### Select sample
otuTab <- otuTab[which(str_extract(rownames(otuTab), "[:digit:]+") %in% intersect(str_extract(rownames(otuTab), "[:digit:]+"), metData$X.SampleID)), ]
dim(otuTab)

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/ibd_huttenhower_cleaned_species.rds"))

## Turnbaugh 2009, OB: ---------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/ob_gordon_2008_v2_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/ob_gordon_2008_v2_results.tar"))

### Check the metadata
metData <- read.delim(paste0(path, "ob_gordon_2008_v2_results/ob_gordon_2008_v2.metadata.txt")) 
dim(metData)
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "ob_gordon_2008_v2_results/RDP/ob_gordon_2008_v2.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

sum(metData$SampleID %in% rownames(otuTab))
sum(rownames(otuTab) %in% metData$SampleID)

saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/ob_gordon_2008_v2_cleaned_species.rds"))

## Zeller 2014, CRC: -----------------------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/crc_zeller_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/crc_zeller_results.tar"))

### Check the metadata
metData <- read.table(paste0(path, "crc_zeller_results/crc_zeller.metadata.txt"), sep = "\t") 
colnames(metData) <- metData[1, ]
metData <- metData[-1, ]
dim(metData)
View(metData)

### OTU Table
otuTab <- read.table(paste0(path, "crc_zeller_results/RDP/crc_zeller.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)
View(otuTab)

otuTab <- otuTab[which(rownames(otuTab) %in% intersect(str_replace(metData[, 1], "\\-", "\\."), rownames(otuTab))), ]
dim(otuTab)
saveRDS(otuTab, file = paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/crc_zeller_cleaned_species.rds"))
#: -----------------------------------------------------------------------------




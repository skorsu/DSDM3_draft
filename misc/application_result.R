### Libraries
library(foreach)
library(doParallel)
library(tidyverse)
library(salso)
library(ClusterZI)

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

### HIV Dataset
path <- "/Users/kevinkvp/Desktop/Github Repo/Manuscript/"
resultFilename <- c(paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1:2, "_init_oneClus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_3clus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_5clus_JUL10_fixed.rds"),
                    paste0(path, "Result/selbal_hiv/result_selbal_HIV_chain_", 1, "_init_20clus_JUL10_fixed.rds"))

registerDoParallel(5)
activeClusMat <- foreach(t = 1:5, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  uniqueCLUS(result$mod)
}
stopImplicitCluster()



result <- readRDS(resultFilename[1])
result$mod$ci_result

dim(result$mod$beta_result)


### EDD Dataset - de novo level
data("singhEDD")
otuTab <- readRDS("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/clean_singh_species.rds")
metaData <- read.delim("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/edd_singh.metadata.txt")

### Import the cluster result
path <- "/Users/kevinkvp/Desktop/Github Repo/Manuscript/"

### Filename
filename <- paste0("Result/Result/result_cleaned_singh_species_chain_", 1:2, "_init_", 
                   c("oneClus", "oneClus", "3clus", "3clus", "5clus", "5clus", "20clus", "20clus"), 
                   "_s2_1_s2MH_1en5_FIXED.rds")

#### Obtain the cluster assignment via salso packages
lastMCMC <- vector("list", 8)

for(i in c(1, 2, 3, 5, 7)){
  
  if(file.exists(paste0(path, filename[i]))){
    set.seed(1)
    result <- readRDS(paste0(path, filename[i]))
    lastMCMC[[i]] <- result$mod$ci_result[40001:50000, ]
    rm(result)
  }
  
}

combClus <- lapply(1:8, function(x){
  
  if(is.matrix(lastMCMC[[x]])){
    data.frame(lastMCMC[[x]])
  } else {
    NULL
  }
  
}) %>%
  bind_rows() %>%
  as.matrix() %>%
  salso() %>%
  as.numeric()

data.frame(singhEDD$SampleID, combClus)

data.frame(SampleID = rownames(otuTab), combClus) %>%
  inner_join(metaData) %>%
  group_by(combClus, Status) %>%
  summarise(n = n())

### Obtain the relative abundance (de novo level)
acRela <- (singhEDD[, -(1:2)]/rowSums(singhEDD[, -(1:2)]))

### Obtain the phylum/family
pf <- data.frame(p = str_extract(colnames(singhEDD[, -(1:2)]), "p__[:alpha:]*"),
                 f = str_extract(colnames(singhEDD[, -(1:2)]), "f__[:alpha:]*"))

### Obtain the aggregate relative abundance (de novo level)
agg_acRela <- matrix(NA, nrow = 303, ncol = 26)
for(i in 1:26){
  
  jIndex <- which(sapply(1:213, function(x){sum(pf[x, ] == distinct(pf)[i, ]) == 2}))
  if(length(jIndex) == 1){
    agg_acRela[, i] <- acRela[, jIndex]
  } else {
    agg_acRela[, i] <- rowSums(acRela[, jIndex])
  }
  
}

### Choose only the highest bacteria
high_combTaxa <- union(sort(colMeans(agg_acRela[which(combClus == 1), ]), decreasing = TRUE, index.return = TRUE)$ix[1:13],
                       sort(colMeans(agg_acRela[which(combClus == 2), ]), decreasing = TRUE, index.return = TRUE)$ix[1:13]) %>%
  union(sort(colMeans(agg_acRela[which(combClus == 3), ]), decreasing = TRUE, index.return = TRUE)$ix[1:13])

### Final Aggregation
colTab <- distinct(pf)[high_combTaxa, ]
aggRELA <- matrix(NA, ncol = length(high_combTaxa) + 1, nrow = 303)
colTab <- colTab %>%
  mutate(famName = str_extract(colTab$f, "[:upper:][:alpha:]*")) %>%
  add_row(p = "Other", f = "Other", famName = "Other")

for(i in 1:(length(high_combTaxa))){
  aggRELA[, i] <- agg_acRela[, high_combTaxa[i]]
}
aggRELA[, (length(high_combTaxa) + 1)] <- 1 - rowSums(aggRELA[, -(length(high_combTaxa) + 1)])
colnames(aggRELA) <- colTab$famName

colorDistinct <- rbind(c("Lachnospiraceae", "#c99930", "Firmicutes"),
                       c("Ruminococcaceae", "#f9bb00", "Firmicutes"),
                       c("Veillonellaceae", "#ffe184", "Firmicutes"),
                       c("Streptococcaceae", "#bf9000", "Firmicutes"),
                       c("Acidaminococcaceae", "#f1c232", "Firmicutes"),
                       c("Erysipelotrichaceae", "#ffe599", "Firmicutes"),
                       c("Clostridiaceae", "#fff2cc", "Firmicutes"),
                       c("Enterococcaceae", "#ffd966", "Firmicutes"),
                       c("Bacteroidaceae", "#39C161", "Bacteroidetes"),
                       c("Porphyromonadaceae", "#2B9F4D", "Bacteroidetes"),
                       c("Rikenellaceae", "#75D792", "Bacteroidetes"),
                       c("Prevotellaceae", "#54D97B", "Bacteroidetes"),
                       c("Enterobacteriaceae", "#E23232", "Proteobacteria"),
                       c("Sutterellaceae", "#C15252", "Proteobacteria"),
                       c("Coriobacteriaceae", "#565DFF", "Actinobacteria"),
                       c("Other","grey90" , "Other"))
                    
aggRELA_long <- aggRELA %>% 
  as.data.frame() %>%
  mutate(ID = singhEDD$SampleID) %>% 
  inner_join(data.frame(ID = rownames(otuTab), Cluster = paste0("Cluster ", combClus))) %>%
  # mutate(ID = singhEDD$SampleID, Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!c(ID, Cluster))

aggRELA_long$name <- factor(aggRELA_long$name, levels = colorDistinct[, 1])
aggRELA_long %>%
  ggplot(aes(x = ID, y = value, fill = name)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", fill = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 15)) +
  scale_fill_manual(values = colorDistinct[, 2]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(y = "Relative Abundances", x = "Participants") +
  facet_grid(. ~ Cluster, scales = "free")

### EDD Dataset - grouped genus level
otuTab <- readRDS("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/clean_singh.rds")
metaData <- read.delim("/Users/kevinkvp/Desktop/Github Repo/Manuscript/Data/Application Data/singh/edd_singh.metadata.txt")
path <- "/Users/kevinkvp/Desktop/Github Repo/Manuscript/"

resultFilename <- c(paste0(path, "Result/singh/result_cleaned_singh_chain_", 1:3, "_init_oneClus_s2_1en1_s2MH_1en3.rds"),
                    paste0(path, "Result/singh/result_cleaned_singh_chain_", 1:3, "_init_3clus_s2_1en1_s2MH_1en3.rds"),
                    paste0(path, "Result/singh/result_cleaned_singh_chain_", 1:3, "_init_5clus_s2_1en1_s2MH_1en3.rds"))

dim(otuTab)

registerDoParallel(5)
activeClusMat <- foreach(t = 1:9, .combine = cbind) %dopar% {
  result <- readRDS(resultFilename[t])
  apply(result$mod$ci_result, 1, uniqueClus)
}
stopImplicitCluster()

activeClusMatPlot <- activeClusMat %>%
  as.data.frame() %>%
  mutate(iter = 1:25000) %>%
  pivot_longer(!iter)

activeClusMatPlot$name <- factor(activeClusMatPlot$name,
                                 levels = paste0("result.", 1:12), labels = paste0("Chain ", 1:12))

ggplot(activeClusMatPlot, aes(x = iter, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Singh (Genus Level) - Active Clusters via MCMC Iterations", x = "Iteration", y = "Number of the active clusters",
       color = "MCMC Chain") +
  guides(color = guide_legend(ncol = 5)) +
  scale_y_continuous(breaks = seq(0, 12, 1))





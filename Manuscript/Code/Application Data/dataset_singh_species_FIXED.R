library(tidyverse)
library(salso)
library(latex2exp)

### Function
uniqueClus <- function(x){
  length(unique(x))
}

### Import the data
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/Manuscript/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/Manuscript/"
}
datapath <- paste0(path, "Data/Application Data/")
otuTab <- readRDS(paste0(datapath, "singh/clean_singh_species.rds"))
metaData <- read.delim(paste0(datapath, "singh/edd_singh.metadata.txt"))

mean(otuTab == 0)

### Which OTU for represetn taxap traceplot
colSums(otuTab) %>% sort(decreasing = TRUE, index.return = TRUE) %>% .$ix %>% .[1:5]
apply(otuTab, 2, var) %>% sort(decreasing = TRUE, index.return = TRUE) %>% .$ix %>% .[1:5]

### Filename
filename <- paste0("Result/Result/result_cleaned_singh_species_chain_", 1:2, "_init_", 
                   c("oneClus", "oneClus", "3clus", "3clus", "5clus", "5clus", "20clus", "20clus"), 
                   "_s2_1_s2MH_1en5_FIXED.rds")

for(i in c(1, 2, 3, 5, 7)){
  file.exists(paste0(path, filename[i])) %>% print()
}

### Obtain the result
activeClus <- matrix(NA, ncol = 8, nrow = 50000)
salsoClus <- matrix(NA, ncol = 8, nrow = 303)
xiImp <- vector("list", 8)
lastMCMC <- vector("list", 8)

start_time <- Sys.time()
for(i in c(1, 2, 3, 5, 7)){
  
  if(file.exists(paste0(path, filename[i]))){
    set.seed(1)
    result <- readRDS(paste0(path, filename[i]))
    activeClus[, i] <- apply(result$mod$ci_result, 1, uniqueClus)
    salsoClus[, i] <- as.numeric(salso(result$mod$ci_result[-(1:10000), ]))
    xiImp[[i]] <- lapply(c(73, 71, 3, 72, 74), function(x){result$mod$beta_result[, x, ] %>% t()})
    lastMCMC[[i]] <- result$mod$ci_result[40001:50000, ]
    rm(result)
  }
  
}
difftime(Sys.time(), start_time)

lapply(1:8, function(x){
  
  data.frame(SampleID = rownames(otuTab), clus = salsoClus[, x]) %>%
    inner_join(metaData) %>%
    group_by(clus, Status) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = Status, values_from = n)
  
})

### Active Cluster for each chain
activeLong <- as.data.frame(activeClus) %>%
  mutate(Iteration = 1:50000) %>%
  pivot_longer(!Iteration)

# chainName <- paste0(sort(rep(c(1, 3, 5, 20), 2)), c(rep(" Cluster", 2), rep(" Clusters", 6)),
#                     ": Chain ", 1:2)

chainName <- paste0("Chain ", c(1, 2, 3, 6, 4, 7, 5, 8))

activeLong$name <- factor(activeLong$name, labels = chainName)

activeLong %>% 
  filter(!is.na(value)) %>%
  ggplot(aes(x = Iteration, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1)) +
  labs(y = "Active Clusters", color = "")

### Cluster Concentration c(73, 71, 3, 72, 74)
xiTaxaIndex <- c(73, 71, 3, 72, 74)
xiPlotList <- lapply(1:5, function(y){
  
  lapply(1:8, function(x){
    
    if(is.matrix(xiImp[[x]][[y]])){
      data.frame(xiImp[[x]][[y]], Iteration = 1:50000, Chain = chainName[x]) %>%
        pivot_longer(!c(Iteration, Chain))
    } else {
      NULL
    }
    
  }) %>%
    bind_rows() %>%
    ggplot(aes(x = Iteration, y = value, color = name)) +
    geom_line() +
    geom_vline(xintercept = 40000, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none") + 
    facet_wrap(. ~ Chain) +
    labs(y = TeX(paste0("$\\xi_{", xiTaxaIndex[y], "}$")))
  
})

xiPlotList[[2]]

### Combined MCMC for the final cluster assignment
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

table(combClus)

data.frame(SampleID = rownames(otuTab), combClus) %>%
  inner_join(metaData) %>%
  group_by(combClus, Status) %>%
  summarise(n = n())

data.frame(SampleID = rownames(otuTab), combClus) %>%
  inner_join(metaData) %>%
  group_by(combClus, Status) %>%
  summarise(n = n()) %>%
  mutate(Cluster = paste0("Cluster ", combClus), Percent = n/sum(n)) %>%
  ggplot(aes(x = Cluster, y = Percent, fill = Status)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab2", "tomato1")) +
  labs(x = "", fill = "") +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1))

### Richness and Shannon
data.frame(richness = rowSums(otuTab > 0), combClus) %>%
  group_by(combClus) %>%
  summarise(mean(richness), sd(richness))

relaTab <- otuTab/rowSums(otuTab)
data.frame(shannon = apply(relaTab, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))}), combClus) %>%
  group_by(combClus) %>%
  summarise(mean(shannon), sd(shannon))

data.frame(Cluster = paste0("Cluster ", combClus), 
           Richness = rowSums(otuTab > 0), 
           Shannon = apply(relaTab, 1, function(x){-sum(ifelse(x == 0, 0, x * log(x)))})) %>%
  pivot_longer(!Cluster) %>%
  ggplot(aes(y = value, x = Cluster)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(. ~ name, scales = "free_y") +
  labs(y = "") +
  theme_bw()

### Phylum difference
taxoLevel <- data.frame(p = str_remove(str_extract(colnames(otuTab), "p__[^;]*"), "p__"),
                        c = str_remove(str_extract(colnames(otuTab), "c__[^;]*"), "c__"),
                        o = str_remove(str_extract(colnames(otuTab), "o__[^;]*"), "o__"),
                        f = str_remove(str_extract(colnames(otuTab), "f__[^;]*"), "f__"),
                        g = str_remove(str_extract(colnames(otuTab), "g__[^;]*"), "g__"))

pDistinct <- unique(taxoLevel$p)[1:4]

sapply(1:4, function(x){rowSums(relaTab[, which(taxoLevel$p == pDistinct[x])])}) %>%
  as.data.frame() %>%
  `colnames<-`(pDistinct) %>%
  mutate(Verrucomicrobia = relaTab[, which(taxoLevel$p == "Verrucomicrobia")],
         Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!Cluster) %>%
  ggplot(aes(x = Cluster, y = value, fill = name)) +
  geom_boxplot(width = 0.4) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  # scale_fill_manual(values = c("olivedrab2", "tomato1")) +
  labs(x = "", fill = "Phyla", y = "Relative Abundance") +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1))

sapply(1:4, function(x){rowSums(relaTab[, which(taxoLevel$p == pDistinct[x])])}) %>%
  as.data.frame() %>%
  `colnames<-`(pDistinct) %>%
  mutate(Verrucomicrobia = relaTab[, which(taxoLevel$p == "Verrucomicrobia")],
         Cluster = paste0("Cluster ", combClus)) %>%
  mutate(bf_ratio = Proteobacteria/Bacteroidetes) %>%
  group_by(Cluster) %>% summarise(mean(bf_ratio), sd(bf_ratio))

testTab <- sapply(1:4, function(x){rowSums(relaTab[, which(taxoLevel$p == pDistinct[x])])}) %>%
  as.data.frame() %>%
  `colnames<-`(pDistinct) %>%
  mutate(Verrucomicrobia = relaTab[, which(taxoLevel$p == "Verrucomicrobia")],
         Cluster = combClus) %>%
  filter(Cluster != 1)

wilcox.test(testTab$Bacteroidetes/testTab$Firmicutes, testTab$Cluster)
wilcox.test(testTab$Proteobacteria/testTab$Firmicutes, testTab$Cluster)
wilcox.test(testTab$Proteobacteria/testTab$Bacteroidetes, testTab$Cluster)

### Taxa Plot
DistincttaxoLevel <- distinct(taxoLevel)
dim(DistincttaxoLevel)
taxoLevelMatch <- sapply(1:60, function(y){sapply(1:213, function(x){sum(taxoLevel[x, ] == DistincttaxoLevel[y, ])})})
comb_relaTab <- sapply(1:60, function(x){
  
  if(length(which(taxoLevelMatch[, x] == 5)) == 1){
    relaTab[, which(taxoLevelMatch[, x] == 5)]
  } else {
    rowSums(relaTab[, which(taxoLevelMatch[, x] == 5)])
  }
  
})

colnames(comb_relaTab) <- ifelse(DistincttaxoLevel$g != "", paste0("g_", DistincttaxoLevel$g), 
                                 ifelse(DistincttaxoLevel$f != "", paste0("f_", DistincttaxoLevel$f), 
                                        ifelse(DistincttaxoLevel$o != "", paste0("o_", DistincttaxoLevel$o), paste0("p_", DistincttaxoLevel$p))))

high_combTaxa <- union(sort(colMeans(comb_relaTab[which(combClus == 1), ]), decreasing = TRUE, index.return = TRUE)$ix[1:10],
                       sort(colMeans(comb_relaTab[which(combClus == 2), ]), decreasing = TRUE, index.return = TRUE)$ix[1:10]) %>%
  union(sort(colMeans(comb_relaTab[which(combClus == 3), ]), decreasing = TRUE, index.return = TRUE)$ix[1:10])

comb_relaTab_VIZ <- data.frame(comb_relaTab[, high_combTaxa], OTHER = rowSums(comb_relaTab[, -high_combTaxa]),
                               ID = rownames(comb_relaTab), Cluster = paste0("Cluster ", combClus)) %>%
  pivot_longer(!c(ID, Cluster))

colorDistinct <- rbind(c("f_Lachnospiraceae", "#c99930", "Firmicutes"),
                       c("f_Ruminococcaceae", "#f9bb00", "Firmicutes"),
                       c("g_Clostridium_sensu_stricto", "#ffe184", "Firmicutes"),
                       c("g_Clostridium_XlVa", "#bf9000", "Firmicutes"),
                       c("g_Faecalibacterium", "#f1c232", "Firmicutes"),
                       c("g_Lachnospiracea_incertae_sedis", "#ffe599", "Firmicutes"),
                       c("g_Streptococcus", "#fff2cc", "Firmicutes"),
                       c("g_Veillonella", "#ffd966", "Firmicutes"),
                       c("g_Alistipes", "#39C161", "Bacteroidetes"),
                       c("g_Bacteroides", "#2B9F4D", "Bacteroidetes"),
                       c("g_Parabacteroides", "#75D792", "Bacteroidetes"),
                       c("g_Prevotella", "#54D97B", "Bacteroidetes"),
                       c("g_Cronobacter", "#E23232", "Proteobacteria"),
                       c("g_Salmonella", "#C15252", "Proteobacteria"),
                       c("f_Enterobacteriaceae", "#C10000", "Proteobacteria"),
                       c("OTHER","grey90" , "OTHER"))

comb_relaTab_VIZ$name <- factor(comb_relaTab_VIZ$name, levels = colorDistinct[, 1])

ggplot(comb_relaTab_VIZ, aes(x = ID, y = value, fill = name)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  # scale_fill_manual(values = c("olivedrab2", "tomato1")) +
  labs(x = "", fill = "") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 5, angle = 90)) +
  facet_wrap(. ~ Cluster, scale = "free_x", ncol = 2) +
  scale_fill_manual(values = colorDistinct[, 2]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(y = "Relative Abundance", x = "ID")






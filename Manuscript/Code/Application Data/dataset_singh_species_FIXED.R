library(tidyverse)
library(salso)

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

### Which OTU for represetn taxap traceplot
colSums(otuTab) %>% sort(decreasing = TRUE, index.return = TRUE) %>% .$ix %>% .[1:5]
apply(otuTab, 2, var) %>% sort(decreasing = TRUE, index.return = TRUE) %>% .$ix %>% .[1:5]

### Filename
filename <- paste0("Result/Result/result_cleaned_singh_species_chain_", 1:2, "_init_", 
                   c("oneClus", "oneClus", "3clus", "3clus", "5clus", "5clus", "20clus", "20clus"), 
                   "_s2_1_s2MH_1en5_FIXED.rds")

for(i in 1:8){
  file.exists(paste0(path, filename[i])) %>% print()
}

### Obtain the result
activeClus <- matrix(NA, ncol = 8, nrow = 50000)
salsoClus <- matrix(NA, ncol = 8, nrow = 303)
xiImp <- vector("list", 8)
lastMCMC <- vector("list", 8)

start_time <- Sys.time()
for(i in 1:8){
  
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

activeLong <- as.data.frame(activeClus) %>%
  mutate(Iteration = 1:50000) %>%
  pivot_longer(!Iteration)

activeLong$name <- factor(activeLong$name, 
                          labels = paste0(sort(rep(c(1, 3, 5, 20), 2)), c(rep(" Cluster", 2), rep(" Clusters", 6)),
                                          ": Chain ", 1:2))

activeLong %>% 
  filter(!is.na(value)) %>%
  ggplot(aes(x = Iteration, y = value, color = name)) +
  geom_line()

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

data.frame(SampleID = rownames(otuTab), combClus) %>%
  inner_join(metaData) %>%
  group_by(combClus, Status) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Status, values_from = n)

str_remove(str_extract(colnames(otuTab), "p__[^;]*"), "p__")
str_remove(str_extract(colnames(otuTab), "c__[^;]*"), "c__")
str_remove(str_extract(colnames(otuTab), "o__[^;]*"), "o__")
str_remove(str_extract(colnames(otuTab), "f__[^;]*"), "f__")
str_remove(str_extract(colnames(otuTab), "g__[^;]*"), "g__")

dim(otuTab)


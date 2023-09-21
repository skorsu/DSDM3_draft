### Required Library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(salso)
library(gridExtra)
library(sparseMbClust)

### Import the data
path <- "/Users/kevin-imac/Desktop/"
ni <- read.csv(paste0(path, "Nicaragua_12mo_Metadata_csv.csv"))
ml <- read.csv(paste0(path, "Mali_12mo_Metadata_csv.csv"))

### Create a new variable and join the two datasets
dat <- cbind("country" = c("NI"), ni[, intersect(colnames(ni), colnames(ml))]) %>%
  rbind(cbind("country" = c("ML"), ml[, intersect(colnames(ni), colnames(ml))]))
for(i in c(1, 2, 3, 5)){
  dat[, i] <- factor(dat[, i])
}

### Clean the data
dat <- dat %>% dplyr::select(-Age)
demo_dat <- dat[, 1:4]
taxa_dat <- dat[, -(1:4)]
taxa_dat <- taxa_dat[, colMeans(taxa_dat > 0) >= 0.1]

dim(taxa_dat)

rowSums(taxa_dat)

### Apply the ZIDM-ZIDM
set.seed(1)
start_time <- Sys.time()
result <- ZIDM_ZIDM(iter = 20000, K_max = 10, z = as.matrix(taxa_dat),
                    theta_vec = rep(1, 10), launch_iter = 5,
                    MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                    r0c = 1, r1c = 1, print_iter = 2000)
difftime(Sys.time(), start_time)

### Overall
nclus <- apply(result$assign, 1, function(x){length(unique(x))})
data.frame(iter = 1:20000, nclus) %>%
  ggplot(aes(x = iter, y = nclus)) +
  geom_line() +
  theme_bw() +
  labs(title = "Number of the active cluster for each MCMC iteration",
       x = "Iteration", y = "Number of active cluster")

### Analysis for each cluster

t1 <- salso(result$assign[-(1:10000), ], maxNClusters = 10, loss = binder()) %>%
  as.numeric()
t2 <- salso(result$assign[-(1:10000), ], maxNClusters = 10, loss = VI()) %>%
  as.numeric()
table(t1, t2)

#### Line plot for each cluster
index_name <- data.frame(colnames(taxa_dat), 
                         index = paste0("TX", str_pad(1:46, 2, pad = "0")))

colnames(taxa_dat) <- paste0("TX", str_pad(1:46, 2, pad = "0"))

list_plot <- vector(mode = "list", length =  4)

for(i in 1:2){
  title_text <- paste0("Cluster ", i)
  list_plot[[i]] <- rowid_to_column(taxa_dat[t1 == i, ]/rowSums(taxa_dat[t1 == i, ])) %>% 
    melt(id.vars = "rowid") %>%
    ggplot(aes(x = variable, y = value)) +
    geom_line(aes(color = rowid, group = rowid)) +
    ylim(0, 1) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title_text, x = "Taxa Count", y = "Proportion")
}

grid.arrange(grobs = list_plot)

### Demographic (1D)

summary(demo_dat[t1 == 1, c(1, 3, 4)])
summary(demo_dat[t1 == 2, c(1, 3, 4)])
summary(demo_dat[t1 == 3, c(1, 3, 4)])
summary(demo_dat[t1 == 4, c(1, 3, 4)])

### Try Shi's
start_time <- Sys.time()
shi_result <- PerformClustering(t(taxa_dat), "MFM", w = 1, thin = 1)
difftime(Sys.time(), start_time)

saveRDS(list(zz_result = t1, shi_result = shi_result),
        paste0(path, "result_mcb.RData"))
?saveRDS

salso(shi_result$crec, loss = binder(), maxNClusters = 10)

































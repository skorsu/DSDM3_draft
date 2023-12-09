### Required Library
library(tidyverse)
library(ggplot2)
library(reshape2)
library(salso)
library(gridExtra)
library(stringr)
library(sparseMbClust)
library(cluster)
library(ecodist)
library(factoextra)
library(rbiom)
## library(Bio)
## library(phyloseq)

### Import the data
# path <- "/Users/kevinkvp/Desktop/"
path <- "/Users/kevin-imac/Desktop/"
ni <- read.csv(paste0(path, "Nicaragua_12mo_Metadata_csv.csv"))
ml <- read.csv(paste0(path, "Mali_12mo_Metadata_csv.csv"))

### Create a new variable and join the two datasets
dat <- cbind("country" = c("NI"), ni[, intersect(colnames(ni), colnames(ml))]) %>%
  rbind(cbind("country" = c("ML"), ml[, intersect(colnames(ni), colnames(ml))]))
for(i in c(1, 2, 3, 5)){
  dat[, i] <- factor(dat[, i])
}

dim(dat)

dat[, 1:5]

### Function
## recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}

### Clean the data
dat <- dat %>% dplyr::select(-Age)
demo_dat <- dat[, 1:4]
taxa_dat <- dat[, -(1:4)]
taxa_dat <- taxa_dat[, colMeans(taxa_dat > 0) >= 0.1]

view(cbind(demo_dat, taxa_dat))

dim(dat)

### Apply: ZIDM-ZIDM
set.seed(1)
start_time <- Sys.time()
result_ZZ <- ZIDM_ZIDM(iter = 20000, K_max = 10, z = as.matrix(taxa_dat),
                    theta_vec = rep(1, 10), launch_iter = 5,
                    MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                    r0c = 1, r1c = 1, print_iter = 2000)
time_ZZ <- difftime(Sys.time(), start_time, units = "mins")

### Apply: DM-ZIDM
set.seed(1)
start_time <- Sys.time()
result_DZ <- DM_ZIDM(iter = 20000, K_max = 10, z = as.matrix(taxa_dat),
                     theta_vec = rep(1, 10), launch_iter = 5,
                     MH_var = 1, mu = 0, s2 = 1, r0c = 1, r1c = 1, 
                     print_iter = 2000)
time_DZ <- difftime(Sys.time(), start_time, units = "mins")

### Apply: DM-DM
set.seed(1)
start_time <- Sys.time()
result_DD <- DM_DM(iter = 20000, K_max = 5, z = as.matrix(taxa_dat),
                   theta_vec = rep(1, 10), MH_var = 1, mu = 0, s2 = 1, 
                   print_iter = 2000)
time_DD <- difftime(Sys.time(), start_time, units = "mins")

### Apply: DM-sDM
set.seed(1)
start_time <- Sys.time()
result_DsD <- DM_DM(iter = 20000, K_max = 10, z = as.matrix(taxa_dat),
                    theta_vec = rep(0.001, 10), MH_var = 1, mu = 0, s2 = 1, 
                    print_iter = 2000)
time_DsD <- difftime(Sys.time(), start_time, units = "mins")

### Apply: Shi (DP)
set.seed(1)
start_time <- Sys.time()
result_DP <- PerformClustering(t(taxa_dat), "DP", w = 1, thin = 1)
time_DP <- difftime(Sys.time(), start_time, units = "mins")

### Apply: Shi (MFM)
set.seed(1)
start_time <- Sys.time()
result_MFM <- PerformClustering(t(taxa_dat), "MFM", w = 1, thin = 1)
time_MFM <- difftime(Sys.time(), start_time, units = "mins")

### Save the result
saveRDS(list(result_ZZ = result_ZZ, result_DZ = result_DZ, 
             result_DD = result_DD, result_DsD = result_DsD, 
             result_DP = result_DP, result_MFM = result_MFM), 
        paste0(path, "result_mcb.RData"))

### Import the result
result <- readRDS(paste0(path, "result_mcb.RData"))

### Number of the Active cluster (MCMC)
n_unique <- function(x){
  length(unique(x))
}

data.frame(iter = 10001:20000, 
           zz = apply(result_ZZ$assign[-(1:10000), ], 1, n_unique),
           dz = apply(result_DZ$assign[-(1:10000), ], 1, n_unique),
           dp = apply(result_DP$crec, 1, n_unique),
           mfm = apply(result_MFM$crec, 1, n_unique)) %>%
  melt(id = "iter") %>%
  ggplot(aes(x = iter, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_color_discrete(name = "Method", 
                      labels = c("ZIDM-ZIDM", "DM-ZIDM", 
                                 "sparseMBClust (DP)", "sparseMBClust (MFM)")) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 22)) +
  labs(x = "Iteration (After Burn-in)", y = "Number of active clusters",
       title = "Number of the active cluster for each semiparametric methods.")

### ZIDM-ZIDM: Analysis for each cluster
t1 <- salso(result_ZZ$assign[-(1:10000), ], maxNClusters = 10, loss = binder()) %>%
  as.numeric()
t2 <- salso(result_ZZ$assign[-(1:10000), ], maxNClusters = 10, loss = VI()) %>%
  as.numeric()
table(t1, t2)

### ZIDM-ZIDM: Line plot for each cluster
taxa_index <- paste0("TX", str_pad(1:46, ceiling(log10(46)), pad = "0"))

tt <- str_extract(colnames(taxa_dat), pattern = regex("f__[:punct:]?[:alpha:]*[e]"))
tt <- ifelse(substr(tt, 4, 4) == ".", substr(tt, 5, 100), substr(tt, 4, 100))
taxa_name <- data.frame(k = substr(str_extract(colnames(taxa_dat), pattern = regex("k__[:alpha:]*")), 4, 1000),
                        p = substr(str_extract(colnames(taxa_dat), pattern = regex("p__[:alpha:]*")), 4, 1000),
                        c = substr(str_extract(colnames(taxa_dat), pattern = regex("c__[:alpha:]*")), 4, 1000),
                        o = substr(str_extract(colnames(taxa_dat), pattern = regex("o__[:alpha:]*")), 4, 1000),
                        f = tt,
                        g = substr(str_extract(colnames(taxa_dat), pattern = regex("g__[:alpha:]*")), 4, 1000))
taxa_name[taxa_name == ""] <- NA

plot_ticks <- ifelse(!is.na(taxa_name[, "g"]), taxa_name[, "g"], 
                     ifelse(!is.na(taxa_name[, "f"]), taxa_name[, "f"], taxa_index))

plot_tick <- c(taxa_index[1], 
               (ifelse(plot_ticks == lag(plot_ticks), taxa_index, plot_ticks))[-1])

list_plot <- vector(mode = "list", length =  4)
for(i in 1:4){
  title_text <- paste0("Cluster ", i)
  list_plot[[i]] <- rowid_to_column(taxa_dat[t1 == i, ]/rowSums(taxa_dat[t1 == i, ])) %>% 
    melt(id.vars = "rowid") %>%
    ggplot(aes(x = variable, y = value)) +
    geom_line(aes(color = rowid, group = rowid)) +
    ylim(0, 1) +
    scale_x_discrete(guide = guide_axis(angle = 90), 
                     labels = plot_tick) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title_text, x = "Taxa Count", y = "Proportion")
}

grid.arrange(grobs = list_plot)

### Demographic with Cluster
demo_clus <- cbind(t2, demo_dat)
table(t2, demo_clus$country)
table(t2, demo_clus$Sex)
table(t2, demo_clus$Group, demo_clus$country)
table(t2, demo_clus$Sex, demo_clus$country)
table(t2, demo_clus$Group, demo_clus$Sex)

### Shi (DP): salso
t1_dp <- salso(result_DP$crec, maxNClusters = 4, loss = VI()) %>% as.numeric()
t2_dp <- salso(result_DP$crec, maxNClusters = 4, loss = binder()) %>% as.numeric()

table(t1_dp, t2_dp)
table(t1_dp)
table(t2_dp)

t1dp_demo <- cbind(t1_dp, demo_dat)
table(t1_dp, demo_dat$country)
table(t1_dp, demo_dat$Group, demo_dat$country)
table(t2_dp, demo_dat$country)
table(t2_dp, demo_dat$Group, demo_dat$country)

### Shi (MFM)
t1_mfm <- salso(result_MFM$crec, maxNClusters = 4, loss = VI()) %>% as.numeric()
t2_mfm <- salso(result_MFM$crec, maxNClusters = 4, loss = binder()) %>% as.numeric()
table(t1_mfm, t2_mfm)
table(t1_mfm)
table(t2_mfm)
table(t1_mfm, demo_dat$country)
table(t1_mfm, demo_dat$Group, demo_dat$country)
table(t2_mfm, demo_dat$country)
table(t2_mfm, demo_dat$Group, demo_dat$country)

### Bray-Curtis
bc_mat <- bcdist(as.matrix(taxa_dat))
bc_hcut <- fviz_nbclust(taxa_dat, FUNcluster = hcut, method = "silhouette", 
             diss = bc_mat) +
  theme_bw() +
  labs(title = "Elbow Method: hcut (Bray-Curtis)")
bc_pam <- fviz_nbclust(taxa_dat, FUNcluster = pam, method = "silhouette", 
             diss = bc_mat) +
  theme_bw() +
  labs(title = "Elbow Method: pam (Bray-Curtis)")
grid.arrange(bc_hcut, bc_pam)

hclust(bc_mat)
bc_hcut <- hcut(bc_mat, k = 2, isdiss = TRUE, hc_method = "complete")
bc_hcut_2 <- bc_hcut$cluster
bc_pam <- pam(bc_mat, 2)$clustering

table(bc_hcut_2, demo_dat$country)
table(bc_hcut_2, demo_dat$Group)
table(bc_hcut_2, demo_dat$Group, demo_dat$country)

table(bc_pam, demo_dat$country)
table(bc_pam, demo_dat$Group)
table(bc_pam, demo_dat$Group, demo_dat$country)



### UniFrac
PhyDat


rm(list = ls())
df2newick(as.data.frame(cbind(tt_tab[2:19, 1:3],  x = colnames(taxa_dat)[2:19])))

str(taxa_name[2:19, -1])
str(df)

bact_name <- colnames(taxa_dat)
bact_name <- gsub("\\.\\.", ".", bact_name)
str_extract(bact_name, pattern = regex("p__[:alpha:]*"))
str_extract(bact_name, pattern = regex("c__[:alpha:]*"))
str_extract(bact_name, pattern = regex("o__[:alpha:]*"))

hclust(dist(as.matrix(t(taxa_dat)))) %>% plot()


read.dendrogram(text = "(1,(2, 3));") %>% plot()
colnames(taxa_dat)

test_mat <- t(taxa_dat)
rownames(test_mat) <- 1:46
tree <- ape::as.phylo(hclust(dist(test_mat)))
rownames()
plot(tree)

unifrac(test_mat, weighted = FALSE, tree = tree)

UniFrac(pp_dat)
phy_tree()
esophagus@otu_table
esophagus@phy_tree %>% plot()






























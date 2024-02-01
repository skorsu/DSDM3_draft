### Required Library
library(tidyverse)
library(stringr)


library(ggplot2)
library(reshape2)
library(salso)
library(gridExtra)
library(sparseMbClust)
library(cluster)
library(ecodist)
library(factoextra)
library(rbiom)
library(ggcorrplot)
library(mclustcomp)
library(pheatmap)
library(gridExtra)

sourceCpp("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/src/clusterZI.cpp")
# sourceCpp("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/src/clusterZI.cpp")

## Import the data -------------------------------------------------------------
path <- "/Users/kevinkvp/Desktop/Annika/"
# path <- "/Users/kevin-imac/Desktop/Annika/"

### 6 and 8 Months
ni68 <- read.csv(paste0(path, "Nicaragua_6mo_8mo_genus.csv"))
ml68 <- read.csv(paste0(path, "Mali_6mo_8mo_genus.csv"))

### 12 Months
ni12 <- read.csv(paste0(path, "Nicaragua_12mo_Metadata_csv.csv"))
ml12 <- read.csv(paste0(path, "Mali_12mo_Metadata_csv.csv"))

## Data Pre-processing ---------------------------------------------------------
### For each nationality, first split 6 and 8 from x68. Then, choose only the 
### common infants among three datasets.
ni06 <- ni68 %>% filter(Age..months. == 6)
ni08 <- ni68 %>% filter(Age..months. == 8)

identical(str_extract(ni06$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"),
          str_extract(ni08$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"))
### For Nicaraguan, 6 an 8 months dataset contain the same infants.

ni06 <- ni06[which(str_extract(ni06$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %in% str_extract(ni12$ID, "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}")),]
ni08 <- ni08[which(str_extract(ni08$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}") %in% str_extract(ni12$ID, "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}")),]

identical(str_extract(ni06$ID., "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"),
          str_extract(ni12$ID, "^[:alpha:]{2}\\.[:alpha:]{2}\\.[:digit:]{2}"))

ml06 <- ml68 %>% filter(Age..months. == 6)
ml08 <- ml68 %>% filter(Age..months. == 8)

str_extract(ml06$X.SampleID, "^[:digit:]+\\.") %in% str_extract(ml08$X.SampleID, "^[:digit:]+\\.")
str_extract(ml08$X.SampleID, "^[:digit:]+\\.") %in% str_extract(ml06$X.SampleID, "^[:digit:]+\\.")

ml06 <- ml06[which(str_extract(ml06$X.SampleID, "^[:digit:]+\\.") %in% str_extract(ml08$X.SampleID, "^[:digit:]+\\.")), ]
identical(str_extract(ml06$X.SampleID, "^[:digit:]+\\."), str_extract(ml08$X.SampleID, "^[:digit:]+\\."))

str_extract(ml12$ID, "^[:digit:]+\\.") %in% str_extract(ml08$X.SampleID, "^[:digit:]+\\.")
str_extract(ml08$X.SampleID, "^[:digit:]+\\.") %in% str_extract(ml12$ID, "^[:digit:]+\\.")
dim(ml06)
dim(ml08)
dim(ml12)



### Functions ------------------------------------------------------------------
### recursion function
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

meanSD <- function(x, dplace = 3){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (", ss, ")")
}

## Data Cleaning ---------------------------------------------------------------
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

view(cbind(demo_dat, taxa_dat))

### Visualize the data ---------------------------------------------------------
pheatmap(taxa_dat, 
         display_numbers = FALSE, color = colorRampPalette(c('white','lightblue3'))(100), 
         cluster_rows = F, cluster_cols = F,
         show_rownames = FALSE, show_colnames = FALSE)

result_path <- "Github Repo/ClusterZI/simulation study/sensitivity/annika/"
result_path <- NULL

### Post-Processing ------------------------------------------------------------
result <- readRDS(paste0(path, result_path, "result.RData"))

#### Check that VI and binder give the same result or not?
clusVI <- as.numeric(salso(result[[1]]$result[-(1:500), ]))
clusBinder <- as.numeric(salso(result[[1]]$result[-(1:500), ]), loss = "binder")
mclustcomp(clusVI, clusBinder) ### According to the ARI, these two loss functions give the same result.

table(clusVI, clusBinder)

### Visualization of the bacteria
relaTaxa <- taxa_dat/rowSums(taxa_dat)

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
colnames(relaTaxa) <- plot_tick

pheatmap(relaTaxa[which(clusVI == 1), ], 
         display_numbers = FALSE, color = colorRampPalette(c('white','maroon3'))(100), 
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Cluster 1: Proportion of the Taxa count")

pheatmap(relaTaxa[which(clusVI == 2), ], 
         display_numbers = FALSE, color = colorRampPalette(c('white','maroon3'))(100), 
         cluster_rows = F, cluster_cols = F,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Cluster 2: Proportion of the Taxa count")

pheatmap(relaTaxa[which(clusVI == 3), ], 
         display_numbers = FALSE, color = colorRampPalette(c('white','maroon3'))(100), 
         cluster_rows = F, cluster_cols = F,
         show_rownames = FALSE, show_colnames = TRUE,
         main = "Cluster 3: Proportion of the Taxa count")

### Demographic: Nationality (/)
table(demo_dat[which(clusVI == 1), "country"])
table(demo_dat[which(clusVI == 2), "country"])
table(demo_dat[which(clusVI == 3), "country"])

### Demographic: Sex (No)
table(demo_dat[which(clusVI == 1), "Sex"])
table(demo_dat[which(clusVI == 2), "Sex"])
table(demo_dat[which(clusVI == 3), "Sex"])

### Demographic: Group (No)
table(demo_dat[which(clusVI == 1), "Group"])
table(demo_dat[which(clusVI == 2), "Group"])
table(demo_dat[which(clusVI == 3), "Group"])

### 2 Demographics: Nationality x Sex
table(demo_dat[which(clusVI == 1), c("country", "Sex")])
table(demo_dat[which(clusVI == 2), c("country", "Sex")])
table(demo_dat[which(clusVI == 3), c("country", "Sex")])

### 2 Demographics: Nationality x Group
table(demo_dat[which(clusVI == 1), c("Group", "country")])
table(demo_dat[which(clusVI == 2), c("Group", "country")])
table(demo_dat[which(clusVI == 3), c("Group", "country")])

### 2 Demographics: Sex x Group
table(demo_dat[which(clusVI == 1), c("Sex", "Group")])
table(demo_dat[which(clusVI == 2), c("Sex", "Group")])
table(demo_dat[which(clusVI == 3), c("Sex", "Group")])

### 3 Demographics
table(demo_dat[which(clusVI == 1), c("Group", "Sex", "country")])
table(demo_dat[which(clusVI == 2), c("Group", "Sex", "country")])
table(demo_dat[which(clusVI == 3), c("Group", "Sex", "country")])

sapply(1:length(hyperParam),  function(x){result[[x]]$time})
sapply(1:length(hyperParam),
       function(x){apply(result[[x]]$result, 1, function(y){length(unique(y))})}) %>%
  matplot(type = "l")

sapply(1:length(hyperParam),
       function(x){apply(result[[x]]$result, 1, function(y){length(unique(y))})}) %>%
  apply(2, meanSD)

clusVI <- sapply(1:length(hyperParam), function(x){as.numeric(salso(result[[x]]$result[-(1:500), ]))})
clusBinder <- sapply(1:length(hyperParam), function(x){as.numeric(salso(result[[x]]$result[-(1:500), ], loss = "binder"))})

sapply(1:length(hyperParam), function(x){length(unique(clusVI[, x]))})
sapply(1:length(hyperParam), function(x){length(unique(clusBinder[, x]))})

ariVI <- diag(21)
ariBinder <- diag(21)
for(i in 1:length(hyperParam)){
  
  for(j in 1:length(hyperParam)){
    ariVI[i, j] <- mclustcomp(clusVI[, i], clusVI[, j], types = "adjrand")[, 2]
    ariBinder[i, j] <- mclustcomp(clusBinder[, i], clusBinder[, j], types = "adjrand")[, 2]
  }
  
}

ggcorrplot(ariVI, lab = TRUE, type = "upper", show.diag = TRUE, 
           title = "Post-processing: using VI as a loss function.",
           show.legend = FALSE)

ggcorrplot(ariBinder, lab = TRUE, type = "upper", show.diag = TRUE, 
           title = "Post-processing: using binder as a loss function.",
           show.legend = FALSE)



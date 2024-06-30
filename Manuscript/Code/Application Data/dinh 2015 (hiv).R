# devtools::install_github("YushuShi/MicrobiomeCluster")

library(ape)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(salso)
library(ecodist)
library(foreach)
library(doParallel)
library(cluster)
library(ggplot2)
library(ecodist)
library(coda.base)
library(mclustcomp)
library(vcd)

### User-defined functions: ----------------------------------------------------
uniqueClus <- function(x){
  length(unique(x))
}

### Global Objects: ------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

### Dinh 2015, HIV: ------------------------------------------------------------
### Import Data and Cleaning OTU: ----------------------------------------------
file.exists(paste0(path, "Manuscript/Data/Application Data/hiv_dinh_results.tar"))
untar(paste0(path, "Manuscript/Data/Application Data/hiv_dinh_results.tar"))

#### Metadata
metData <- read.table(paste0(path, "hiv_dinh_results/hiv_dinh.metadata.txt"), sep = "\t", header = TRUE)
metData[metData == "<not provided>"] <- NA
dim(metData)

#### OTU Table
otuTab <- read.table(paste0(path, "hiv_dinh_results/RDP/hiv_dinh.otu_table.100.denovo.rdp_assigned"))
otuTab <- t(otuTab)
dim(otuTab)

##### Check the samples: All samples in OTU has the metadata.
str_extract(rownames(otuTab), "[^X]+") %in% metData$Sample_Name_s
metData$Sample_Name_s %in% str_extract(rownames(otuTab), "[^X]+") 
rownames(otuTab) <- str_extract(rownames(otuTab), "[^X]+")

##### Filter OTU table
otuTab <- otuTab[, -which(colMeans(otuTab > 0) < 0.1)]
dim(otuTab)

### EDA: -----------------------------------------------------------------------
View(metData)
table(metData$Current_smoking_s)

metNew <- metData %>%
  transmute(ID = Sample_Name_s,
            Disease = disease_s,
            Gender = ifelse(sex_s == "Female", "Female", ifelse(is.na(is_gay_s), "Male: Hetorosexual", "Male: Homosexual")),
            Education = ifelse(Education_s == "college", "College", ifelse(Education_s == "high_school", "Up to High School", "No High School")),
            Drinker = ifelse(Alcohol_s == "Non_drinker", "Non-Drinker", ifelse(Bingeing_s == "binge_drinker", "Drinker: Binge", "Drinker: Non-Binge")),
            Smoking = ifelse(Current_smoking_s == "Not_currently_smoking", "No", "Yes"), 
            BMI = BMI_s, Age = age_s, Cholesterol = as.numeric(lvCholesterol_s), 
            IFNgamma = as.numeric(IFNgamma_s), IL1Beta = as.numeric(IL1Beta_s), 
            IL6 = as.numeric(IL6_s), Kcal = as.numeric(Kcal_s), TNFalpha = as.numeric(TNFalpha_s))

metNew$Education <- factor(metNew$Education, levels = c("College", "Up to High School", "No High School"))

metNew %>%
  group_by(Disease) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Disease, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Disease Status")

metNew %>%
  group_by(Gender) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Gender, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Gender")

metNew %>%
  group_by(Education) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Education, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Education Level")

metNew %>%
  group_by(Drinker) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Drinker, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Alcohol Drinker")

metNew %>%
  group_by(Smoking) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = Smoking, y = n)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() +
  labs(y = "Frequency", title = "Smoking Behavior")

metNew %>%
  dplyr::select(ID, BMI, Age, Cholesterol, IFNgamma, IL1Beta, IL6, Kcal, TNFalpha) %>%
  pivot_longer(!ID) %>%
  ggplot(aes(y = value)) +
  geom_boxplot() +
  facet_wrap(. ~ name, scales = "free_y")

### PAM: -----------------------------------------------------------------------
#### CHoosing the number of clusters for each distance metric
data.frame(k = 2:10,
           Euclidean = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(otuTab), y))[, 3])}),
           BrayCurtis = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(bcdist(otuTab), y))[, 3])}),
           Aitchison = sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(otuTab + 1e-20, method = "aitchison"), y))[, 3])})) %>%
  pivot_longer(cols = !k) %>%
  ggplot(aes(x = k, y = value)) +
  geom_line() +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_bw() +
  labs(y = "Average Silhouette", x = "Number of clusters")

##### PAM - Cluster Assignment with the metadata
metPAM <- data.frame(ID = names(pam(bcdist(otuTab), 7)$clustering),
                     clusEu = pam(dist(otuTab), 7)$clustering, 
                     clusBC = pam(bcdist(otuTab), 7)$clustering,
                     clusAT = pam(dist(otuTab + 1e-20, method = "aitchison"), 2)$clustering) %>%
  `rownames<-`(NULL) %>%
  inner_join(metNew)

View(metPAM)

structable(Disease ~ clusEu, metPAM)
structable(Gender ~ clusEu, metPAM)
structable( ~ clusEu, metPAM)



table(metPAM$clusEu, metPAM$Disease)

table(pam(dist(otuTab), 7)$clustering, metData$DiseaseState)
table(pam(dist(otuTab), 7)$clustering, metData$Alcohol_s)
table(pam(dist(otuTab), 7)$clustering, metData$Education_s)
table(pam(dist(otuTab), 7)$clustering, metData$sex_s)

##### Cluster: PAM with Bray-Curtis
table(pam(bcdist(otuTab), 7)$clustering, metData$DiseaseState)
table(pam(bcdist(otuTab), 7)$clustering, metData$Alcohol_s)
table(pam(bcdist(otuTab), 7)$clustering, metData$Education_s)
table(pam(bcdist(otuTab), 7)$clustering, metData$sex_s)

### ZIDM-ZIDM: -----------------------------------------------------------------
#### Import the result form our model
result <- readRDS(paste0(path, "Manuscript/Result/Application Data/Dinh/hiv_dinh_result_split_10_theta_1_rn.rds"))

#### Convergence Checking: -----------------------------------------------------
##### Active Cluster in MCMC iterations
actClus <- sapply(1:6, function(x){apply(result[[x]]$mod$ci_result, 1, uniqueClus)}) %>%
  as.data.frame() %>%
  `colnames<-`(paste0("Chain ", 1:6)) %>%
  mutate(Iteration = 1:1000) %>%
  pivot_longer(!Iteration)

actClus$name <- factor(actClus$name, labels = c(paste0("One Cluster, Kmax = 10 - Chain ", 1:2),
                                                paste0("Three Clusters, Kmax = 20 - Chain ", 1:2),
                                                paste0("Singletons, Kmax = 36 - Chain ", 1:2)))

ggplot(actClus, aes(x = Iteration, y = value, color = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Thinned Iteration", y = "# Active Clusters",
       title = "Dinh: Number of active clusters via ZIDM-ZIDM",
       color = "Initialization") +
  theme(legend.position = "bottom")

##### Acceptance Rate
sapply(1:6, function(x){mean(result[[x]]$mod$sm_accept)})

table(result[[1]]$mod$MH_accept[, 1])

table(result[[1]]$mod$MH_accept[-(1:2500), 1])

334/(1919 + 334)
320/(880 + 320)

sapply(1:6, function(x){table(result[[x]]$mod$sm_status, result[[x]]$mod$sm_accept)})
result[[6]]$time/(3600)

table(result[[1]]$mod$sm_status, result[[1]]$mod$sm_accept)

salso(result[[1]]$mod$ci_result[-(1:100), ], loss = "binder")

clusIndex <- sapply(1:6, function(x){set.seed(1); salso(result[[x]]$mod$ci_result[-(1:100), ], loss = "binder")})

registerDoParallel(5)
clusAssign <- foreach(t = 1:10, .combine = "cbind") %dopar% {
 as.numeric(salso(result[[1]]$mod$ci_result[-(1:100), ], loss = "binder"))
}
stopImplicitCluster()

as.matrix(expand.grid(1:10, 1:10)) %>%
  apply(1, function(x){mclustcomp(clusAssign[, x[1]], clusAssign[, x[2]])[1, 2]
}) %>%
  mean()


table(clusIndex[, 6], metData$DiseaseState)

plot(apply(result[[1]]$mod$ci_result, 1, uniqueClus), type = "l")

plot(result[[1]]$mod$beta_result[10, 834, ], type = "l")

which.max(colSums(otuTab))

set.seed(1)
clusAgg <- lapply(1:6, function(x){as.data.frame(result[[x]]$mod$ci_result[seq(100, 1000, 1), ])}) %>%
  bind_rows() %>%
  as.matrix() %>%
  salso(loss = "binder")
table(clusAgg, metData$DiseaseState)

plot(apply(result[[1]]$mod$ci_result, 1, uniqueClus), type = "l")

set.seed(1)
salso(result[[1]]$mod$ci_result[-(1:200), ], loss = "binder")

clusIndex <- as.numeric(salso(result[[1]]$mod$ci_result[-(1:200), ], loss = "binder"))
table(clusIndex, metData$DiseaseState)
table(clusIndex, metData$Alcohol_s)
table(clusIndex, metData$sex_s)

table(clusIndex, paste0(metData$sex_s, ": ", metData$is_gay_s))

data.frame(clus = paste0("Cluster ", clusIndex), val = metData$BMI_s) %>%
  ggplot(aes(x = clus, y = val, group = clus)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "BMI")

data.frame(clus = paste0("Cluster ", clusIndex), val = metData$age_s) %>%
  ggplot(aes(x = clus, y = val, group = clus)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Age")

data.frame(clus = paste0("Cluster ", clusIndex), val = metData$Fat_s) %>%
  ggplot(aes(x = clus, y = val, group = clus)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Fat")


data.frame(clusIndex, val = metData$Fat_s) %>%
  group_by(clusIndex) %>%
  summarise(mean(val), sd(val))

table(pam(dist(otuTab), 7)$clustering, metData$Alcohol_s)
table(clusIndex, metData$Education_s)

### Try: Random Starting Point (Theta = 1)
### Note: Dihn 6/27
set.seed(1)
clusInit <- matrix(0, ncol = 6, nrow = 36)
clusInit[, 3] <- sample(0:2, size = 36, replace = TRUE)
clusInit[, 4] <- sample(0:2, size = 36, replace = TRUE)
clusInit[, 5] <- 0:35
clusInit[, 6] <- 0:35

KmaxVec <- c(10, 10, 20, 20, 36, 36)

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
start_time_GLOBAL <- Sys.time()

result_split_10_theta_1_rn <- foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 5000, Kmax = KmaxVec[t], nbeta_split = 10,
                      z = as.matrix(otuTab),
                      atrisk_init = matrix(1, nrow = 36, ncol = 1024),
                      beta_init = matrix(0, nrow = KmaxVec[t], ncol = 1024),
                      ci_init = clusInit[, t],
                      theta = 1, mu = 0, s2 = 1e-3,
                      s2_MH = 1e-3, t_thres = 2500, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 5)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  return(list(time = comp_time, mod = mod))
}
stopImplicitCluster()
difftime(Sys.time(), start_time_GLOBAL)


# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(5)
# start_time_GLOBAL <- Sys.time()
# s2Mat <- data.frame(expand.grid(c(1e-3, 1e-4, 1e-5), c(1e-3, 1e-4, 1e-5)))
# 
# result_split_10_r0g_4 <- foreach(t = 1:9) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 3000, Kmax = 10, nbeta_split = 10, 
#                       z = as.matrix(otuTab), 
#                       atrisk_init = matrix(1, nrow = 36, ncol = 1024), 
#                       beta_init = matrix(0, nrow = 10, ncol = 1024), 
#                       ci_init = rep(0, 36), 
#                       theta = 1, mu = 0, s2 = s2Mat[t, 1], 
#                       s2_MH = s2Mat[t, 2], t_thres = 1500, launch_iter = 30, 
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 3)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   return(list(time = comp_time, mod = mod))
# }
# stopImplicitCluster()
# difftime(Sys.time(), start_time_GLOBAL)
# 
# saveRDS(result_split_10_r0g_4, file = paste0(path, "Manuscript/Data/Application Data/hiv_dinh_result_split_10_r0g_4.rds"))
# 
# # set.seed(1)
# # startTime <- Sys.time()
# # modResult <- mod_adaptive(iter = 1000, Kmax = 10, nbeta_split = 10, 
# #                           z = as.matrix(otuTab), 
# #                           atrisk_init = matrix(1, nrow = 36, ncol = 1024), 
# #                           beta_init = matrix(0, nrow = 10, ncol = 1024), 
# #                           ci_init = rep(0, 36), theta = 1, mu = 0, s2 = 1e-3, 
# #                           s2_MH = 1e-5, t_thres = 500, launch_iter = 30, 
# #                           r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 5)
# # totTime <- difftime(Sys.time(), startTime)
# 
# 
# 
# 
# plot(apply(modResult$ci_result, 1, activeClus), type = "l")
# 
# 
# #clus_assign <- data.frame(ID = str_extract(rownames(otuTab), "[^X]+"),
# #                          cluster = as.numeric(salso(modResult$ci_result, loss = "binder")))
# 
# ## This is the result when using 1e-3 for both s2.
# table(metData$sex_s, clus_assign$cluster)

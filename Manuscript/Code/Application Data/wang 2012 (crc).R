# Libraries: -------------------------------------------------------------------
library(salso)
library(tidyverse)
library(foreach)
library(doParallel)
library(cluster)

# User-defined functions: ------------------------------------------------------
uniqueClus <- function(x){
  length(unique(x))
}


# Path: ------------------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/"
}

data_path <- paste0(path, "Manuscript/Data/Application Data/Cleaned Data/Species/")
result_path <- paste0(path, "Manuscript/Result/")

# Import the data: -------------------------------------------------------------
dat <- readRDS(paste0(data_path, "crc_zhao_cleaned_species.rds")) 
dim(dat)

# Run the models: --------------------------------------------------------------
## The case description: Note July 7/2
hyperParam <- rbind(c(4, 10, 1e-3), c(9, 10, 1e-3), c(1, 5, 1e-3), 
                    c(1, 10, 1e-5), c(4, 10, 1e-5), c(9, 10, 1e-3))

hyperParam_2 <- rbind(c(1e-3, 1e-5, 1, 1, 30), c(1e-5, 1e-3, 1, 1, 30), 
                      c(1e-3, 1e-3, 1, 99, 30), c(1e-3, 1e-3, 99, 1, 30), 
                      c(1e-3, 1e-3, 1, 1, 1), c(1e-5, 1e-5, 1, 1, 1))

hyperParam_3 <- rbind(c(1e-3, 1e-5, 30, 1), c(1e-3, 1e-5, 30, 5), 
                      c(1e-3, 1e-5, 30, 25), c(1e-3, 1e-4, 30, 10), 
                      c(1e-3, 1e-10, 30, 10), c(1e-3, 1e-5, 50, 10))

hyperParam_4 <- rbind(c(5, 1, 1), c(5, 1e-5, 1), c(10, 1, 1), 
                      c(10, 1e-5, 1), c(10, 1, 9), c(10, 1, 9))

# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 5000, Kmax = hyperParam[t, 2], nbeta_split = 10,
#                       z = as.matrix(dat),
#                       atrisk_init = matrix(1, nrow = 102, ncol = 123),
#                       beta_init = matrix(0, nrow = hyperParam[t, 2], ncol = 123),
#                       ci_init = rep(0, 102),
#                       theta = 1e-5, mu = 0, s2 = 1e-3,
#                       s2_MH = hyperParam[t, 3], t_thres = 500, launch_iter = 30,
#                       r0g = 1, r1g = 1, r0c = 1, r1c = hyperParam[t, 1], thin = 2)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, result = mod), file = paste0(result_path, "result_crc_zhao_JULY02_case_", t, ".rds"))
# }
# stopImplicitCluster()

# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 5000, Kmax = 5, nbeta_split = 10,
#                       z = as.matrix(dat),
#                       atrisk_init = matrix(1, nrow = 102, ncol = 123),
#                       beta_init = matrix(0, nrow = 5, ncol = 123),
#                       ci_init = rep(0, 102),
#                       theta = 1e-10, mu = 0, s2 = hyperParam_2[t, 1],
#                       s2_MH = hyperParam_2[t, 2], t_thres = 500, launch_iter = hyperParam_2[t, 5],
#                       r0g = 1, r1g = 1, r0c = hyperParam_2[t, 3], r1c = hyperParam_2[t, 4], thin = 2)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, result = mod), file = paste0(result_path, "result_crc_zhao_JULY02_case_", 6+t, ".rds"))
# }
# stopImplicitCluster()

# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 5000, Kmax = 5, nbeta_split = hyperParam_3[t, 4],
#                       z = as.matrix(dat),
#                       atrisk_init = matrix(1, nrow = 102, ncol = 123),
#                       beta_init = matrix(0, nrow = 5, ncol = 123),
#                       ci_init = rep(0, 102),
#                       theta = 1e-10, mu = 0, s2 = hyperParam_3[t, 1],
#                       s2_MH = hyperParam_3[t, 2], t_thres = 500, launch_iter = hyperParam_3[t, 3],
#                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 2)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, result = mod), file = paste0(result_path, "result_crc_zhao_JULY02_case_", 12+t, ".rds"))
# }
# stopImplicitCluster()

# set.seed(1, kind = "L'Ecuyer-CMRG")
# registerDoParallel(6)
# foreach(t = 1:6) %dopar% {
#   start_time <- Sys.time()
#   mod <- mod_adaptive(iter = 5000, Kmax = hyperParam_4[t, 1], nbeta_split = 1,
#                       z = as.matrix(dat),
#                       atrisk_init = matrix(1, nrow = 102, ncol = 123),
#                       beta_init = matrix(0, nrow = hyperParam_4[t, 1], ncol = 123),
#                       ci_init = rep(0, 102),
#                       theta = 1e-10, mu = 0, s2 = hyperParam_4[t, 2],
#                       s2_MH = 1e-10, t_thres = 500, launch_iter = 30,
#                       r0g = 1, r1g = 1, r0c = 1, r1c = hyperParam_4[t, 3], thin = 2)
#   comp_time <- difftime(Sys.time(), start_time, units = "secs")
#   saveRDS(list(time = comp_time, result = mod), file = paste0(result_path, "result_crc_zhao_JULY02_case_", 18+t, ".rds"))
# }
# stopImplicitCluster()

set.seed(1, kind = "L'Ecuyer-CMRG")
registerDoParallel(6)
foreach(t = 1:6) %dopar% {
  start_time <- Sys.time()
  mod <- mod_adaptive(iter = 10000, Kmax = 5, nbeta_split = 1,
                      z = as.matrix(dat),
                      atrisk_init = matrix(1, nrow = 102, ncol = 123),
                      beta_init = matrix(0, nrow = 5, ncol = 123),
                      ci_init = rep(0, 102),
                      theta = 1e-10, mu = 0, s2 = 1e-3,
                      s2_MH = 1e-5, t_thres = 1000, launch_iter = 30,
                      r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 2)
  comp_time <- difftime(Sys.time(), start_time, units = "secs")
  saveRDS(list(time = comp_time, result = mod), file = paste0(result_path, "result_crc_zhao_JULY02_case_13_chain", t, ".rds"))
}
stopImplicitCluster()

# Run the result: --------------------------------------------------------------
result <- readRDS(paste0(result_path, "result_crc_zhao_JULY02_case_24.rds"))
apply(result$result$ci_result, 1, uniqueClus) %>% plot(type = "l")
mean(result$result$sm_accept)
table(factor(result$result$sm_status, labels = c("Merge", "Split")),
      factor(result$result$sm_accept, labels = c("Reject", "Accept")))
result$time/60

set.seed(1)
clusZZ <- as.numeric(salso(result$result$ci_result, loss = "binder"))
table(diseaseIND, clusZZ)

sapply(1:5, function(x){sum(result$result$MH_accept[, x] == 1)/sum(result$result$MH_accept[, x] != -1)})
result$result$beta_result[2, 1, ] %>% plot(type = "l")

# PAM: -------------------------------------------------------------------------
sapply(2:10, function(y){mean(cluster::silhouette(x = pam(dist(dat), y))[, 3])})
pam_EU <- pam(dist(dat), 2)$cluster
diseaseIND <- str_extract(names(pam_EU), "^[:alpha:]")   

table(diseaseIND, clusZZ)

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


# result$mod$beta_result[1, 1, ] %>%
#   plot(type = "l")
# 
# # result_crc_zhao_xiSplit_25_theta_1e3_ADP_1000_iter_10000
# result$time/(3600)
# result <- readRDS(paste0(result_path, "result_crc_zhao_xiSplit_10_theta_1e3_ADP_1000_iter_10000", ".rds"))
# apply(result$mod$ci_result, 1, uniqueClus) %>% plot(type = "l")
# salso(result$mod$ci_result, loss = "binder")



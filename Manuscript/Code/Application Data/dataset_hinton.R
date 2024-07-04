# Required Libraries
library(selbal) # Dataset
library(tidyverse)
library(salso)
library(ggplot2)

# User-defined function
uniqueClus <- function(x){
  length(unique(x))
}

# HIV dataset
datHIV <- HIV
statusHIV <- datHIV[, 62]
otuHIV <- datHIV[, 1:60]
## otuHIV[, -which(colMeans(otuHIV > 0) < 0.1)]
# log(colMeans(otuHIV)) %>% var()
View(otuHIV)
dim(otuHIV)

(otuHIV/rowSums(otuHIV)) %>%
  colMeans() %>%
  as.numeric() %>%
  var()

(otuHIV/rowSums(otuHIV)) %>%
  colMeans() %>%
  as.numeric() %>%
  sort(decreasing = TRUE)

otuHIV %>%
  colMeans() %>%
  as.numeric() %>%
  var()

set.seed(1)
result <- mod_adaptive(iter = 5000, Kmax = 10, nbeta_split = 5, 
                       z = as.matrix(otuHIV), atrisk_init = matrix(1, nrow = 155, ncol = 60), 
                       beta_init = matrix(0, nrow = 10, ncol = 60), 
                       ci_init = rep(0, 155), 
                       theta = 1, mu = 0, s2 = 1, s2_MH = 1e-3, 
                       t_thres = 500, launch_iter = 30, 
                       r0g = 1, r1g = 1, r0c = 1, r1c = 1, thin = 1)

# x <- rnorm(60, 0, (10)^2)
# exp(x)/sum(exp(x))

### Check the convergence of xi
result$beta_result[, 1, ] %>% t() %>%
  as.data.frame() %>%
  mutate(Iteration = 1:5000) %>%
  pivot_longer(!Iteration) %>%
  transmute(Iteration, xi = value,
            Cluster = paste0("Cluster ", str_extract(name, "[:digit:]+"))) %>%
  ggplot(aes(x = Iteration, y = xi, color = Cluster)) +
  geom_line() +
  theme_bw()

as.numeric(salso(result$ci_result[-(1:250), ])) %>% table(statusHIV)
apply(result$ci_result, 1, uniqueClus) %>% plot(type = "l")
sapply(1:10, function(x){sum(result$MH_accept[, x] == 1)/sum(result$MH_accept[, x] != -1)})


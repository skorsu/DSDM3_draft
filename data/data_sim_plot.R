rm(list = ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(gridExtra)

### Simulate the data
# source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")
source("/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/data_sim.R")
set.seed(31807)
dat_test <- data_sim(n = 100, K = 2, J_imp = 4, 
                     pi_gm_mat = matrix(c(0.75, 0.5), nrow = 2, ncol = 10), 
                     xi_scale = 5, sum_zi_imp = 80, sum_zi_unimp = 20)
dat_test$z

### Heat map
dat_heat <- dat_test$z
rownames(dat_heat) <- str_pad(1:dim(dat_heat)[1], 3, pad = "0")
colnames(dat_heat) <- str_pad(1:dim(dat_heat)[2], 3, pad = "0")

dat_heat %>% 
  as.data.frame() %>%
  rownames_to_column("observation") %>%
  pivot_longer(-c(observation), names_to = "variable", values_to = "counts") %>%
  ggplot(aes(x = variable, y = observation, fill = counts)) + 
  geom_raster() +
  scale_fill_gradient2(low = "white", mid = "grey", high = "black", midpoint = 50) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 9), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "bottom", plot.title = element_text(size = 30),
        legend.title = element_text(size = 20), legend.text = element_text(size = 15)) +
  labs(x = "Variable", y = "Observation", title = "Example of the Simulated Data",
       fill = TeX("$z_{ijk}$"))

### Line plot
dat <- cbind(index = 1:100, dat_test$z)
colnames(dat)[2:11] <- str_pad(1:10, 3, pad = "0")
t2 <- data.frame(dat[which(dat_test$ci == 0), ])
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)
p0 <- ggplot(tt2, aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "", title = "Simulated Data: First Cluster")
t2 <- data.frame(dat[which(dat_test$ci == 1), ])
tt2 <- melt(t2, id.vars = "index")
tt2$variable <- substr(tt2$variable, 2, 4)
p1 <- ggplot(tt2, aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "", title = "Simulated Data: Second Cluster")

grid.arrange(p0, p1)


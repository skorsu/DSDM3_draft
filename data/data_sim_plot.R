rm(list = ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(gridExtra)

### Simulate the data
set.seed(21)
source("/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/data_sim.R")
dat_test <- data_sim(N = 500, K = 4, J = 200, J_imp = 20, pi_gk = rep(0.25, 5), 
                     pi_g_ova = 0.75, xi_conc = 5, U_imp = 100, U_unimp = 100,
                     shuffle = FALSE)
dat_test$z

### Heat map
dat_heat <- dat_test$z
rownames(dat_heat) <- paste0("obs ", seq(1, dim(dat_heat)[1]))
colnames(dat_heat) <- str_pad(1:dim(dat_heat)[2], 3, pad = "0")

dat_heat %>% 
  as.data.frame() %>%
  rownames_to_column("observation") %>%
  pivot_longer(-c(observation), names_to = "variable", values_to = "counts") %>%
  ggplot(aes(x = variable, y = observation, fill = counts)) + 
  geom_raster() +
  scale_fill_gradient2(low = "white", mid = "grey", high = "black", midpoint = 50) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "bottom", plot.title = element_text(size = 30),
        legend.title = element_text(size = 20), legend.text = element_text(size = 15)) +
  labs(x = "Variable", y = "Observation", title = "Simulated Data",
       fill = TeX("$z_{ijk}$"))

### Line plot
dat <- cbind(index = 1:500, dat_test$z)
t2 <- data.frame(dat[which(dat_test$ci == 0), ])
p0 <- ggplot(melt(t2, id.vars = "index"), aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "Observation", title = "Simulated Data: First Cluster")
t2 <- data.frame(dat[which(dat_test$ci == 1), ])
melt(t2, id.vars = "index")
p1 <- ggplot(melt(t2, id.vars = "index"), aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "Observation", title = "Simulated Data: Second Cluster")
t2 <- data.frame(dat[which(dat_test$ci == 2), ])
melt(t2, id.vars = "index")
p2 <- ggplot(melt(t2, id.vars = "index"), aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "Observation", title = "Simulated Data: Third Cluster")
t2 <- data.frame(dat[which(dat_test$ci == 3), ])
melt(t2, id.vars = "index")
p3 <- ggplot(melt(t2, id.vars = "index"), aes(x = variable, y = value)) + 
  geom_line(aes(color = index, group = index)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.position = "hide", plot.title = element_text(size = 30)) +
  labs(x = "Variable", y = "Observation", title = "Simulated Data: Forth Cluster")

grid.arrange(p0, p1, p2, p3)


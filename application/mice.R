library(tidyverse)
library(salso)

### Import the data
load("/Users/kevinkvp/Desktop/mice_janet/genus.Rdata")
dat <- as.matrix(t(otu_table(phy.genus)))

### Clean the data
dim(dat)
m <- 10
prop_zero <- rep(m + 1)
for(i in 0:m){
  prop_zero[i + 1] <- sum(colMeans(dat > 0) >= i/m)
}

plot(x = (0:m)/m, y = prop_zero, type = "b",
     xlab = "Threshold", ylab = "Number of the remaining variables")

clean_dat <- dat[, colMeans(dat > 0) >= 0.5]

mod <- ZIDM_ZIDM(iter = 20000, K_max = 10, z = clean_dat, theta_vec = rep(1, 10),
                 launch_iter = 5,
                 MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                 r0c = 1, r1c = 1, print_iter = 2000)

zz_assign <- as.numeric(salso(mod$assign[-(1:10000), ], maxNClusters = 10))

### Demo Analysis
table(zz_assign, phy.genus@sam_data$Treatment)


otu_table(phy.genus) %>% dim()
tax_table(phy.genus)

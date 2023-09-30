### Library
library(tidyverse)
library(rjson)

### Read the data
json_dat <- fromJSON(file = "/Users/kevinkvp/Desktop/animal_parsa/cods.json")
unlist(json_dat[[1]]) %>% sum()


dat <- read.csv("/Users/kevinkvp/Desktop/animal_parsa/species_normed.csv")
dat[1, -1] %>% as.numeric()
apply(dat[-1, ], 1, mean)

?colMeans

view(dat[1:2, ])

###
dat[1, ] %>% as.numeric()
view(dat)
colnames(dat)

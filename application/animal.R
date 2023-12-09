### Library
library(tidyverse)
library(ggplot2)
library(rjson)

### Read the data
dat <- read.csv("/Users/kevin-imac/Desktop/animal/species_normed.csv")
dat <- t(dat)
colnames(dat) <- dat[1, ]
dat <- dat[-1, ]
dat <- as.data.frame(dat)
dat <- mutate_all(dat, function(x) as.numeric(as.character(x)))

animal_dat <- read.csv("/Users/kevin-imac/Desktop/animal/sample_metadata.csv")
dim(animal_dat)

dat <- inner_join(animal_dat, 
                  dat %>% mutate(animal_ID = rownames(dat)),
                  by = c("Sample" = "animal_ID"))
dim(dat)

### There are some NA
dat <- dat[rowSums(is.na(dat)) == 0, ]
which(as.numeric(colMeans(dat[, -(1:2)] > 0)) >= 0.10)

dat <- dat %>% dplyr::select(-X)
dim(dat)
table(str_to_title(dat$Host)) %>%
  as.data.frame() %>%
  rename(Count = Freq) %>%
  mutate(Animal = fct_reorder(Var1, -Count)) %>%
  arrange(-Count) %>%
  ggplot(aes(x = Animal, y = Count)) +
  geom_bar(stat = "identity", fill = "#6f6ec8") +
  geom_text(aes(label = Count), vjust = -0.15, color = "black",
            position = position_dodge(0.5), size = 5) +
  theme_bw() + 
  labs(title = "Frequencies of the animal in the dataset") +
  theme(axis.text.x = element_text(size = 20, angle = 0),
        axis.text.y = element_text(size = 15, angle = 0),
        plot.title = element_text(size = 20))
  
### plot for type of animals

dim(dat)

table(dat$Host)
view(dat[1:5, ])

dat[1, -1] %>% as.numeric()
apply(dat[-1, ], 1, mean)

?colMeans

view(dat[1:2, ])

###
dat[1, ] %>% as.numeric()
view(dat)
colnames(dat)

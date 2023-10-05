rm(list = ls())
options(scipen = 10)

### Required Libraries: --------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(salso)
library(mclustcomp)

### Setting: -------------------------------------------------------------------
path <- "/Users/kevin-imac/Desktop/simu_study/" ### path for saving the data and result
## path <- "/Users/kevinkvp/Desktop/simulation study/"
case_name <- "z_1_pi_90_J_100" ### the path that use for differentiating each case

### Read the RData file: -------------------------------------------------------
dat <- readRDS(paste0(path, "data_", case_name, ".RData")) 
result <- readRDS(paste0(path, "result_", case_name, ".RData")) 

### User-defined functions: ----------------------------------------------------
### Function: Summarize the simulated data
summarise_dat <- function(list_simDat){
  
  lapply(list_simDat, function(x){as.data.frame(x) |> 
      mutate(obs = paste0("OB", str_pad(1:50, 3, pad = "0"))) |>
      pivot_longer(cols = -obs) |>
      mutate(taxa_name = paste0("TX", str_pad(str_extract(name, "[:digit:]+$"), 3, pad = "0"))) |>
      ggplot(aes(taxa_name, obs, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low="white", high="palegreen3") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
      labs(x = "Variable (Taxa Count)", y = "Observation", fill = "Count")})
  
}

### Function: Create the "mean (sd)"
ms_format <- function(val_vec, r = 4){
  paste0(round(mean(val_vec), r), " (", round(sd(val_vec), r), ")")
}

### Function: apply salso for retriving the final cluster assignment
clus_list <- function(dat_list, burn_in, loss_type){
  ### Apply for ZIDM-ZIDM, DM-ZIDM, DM-DM, and DM-sDM
  dat_list[c('result_ZZ', 'result_DZ', 'result_DD', 'result_DsD')] |>
    lapply( `[`, -c(1:burn_in), ) |>
    lapply(salso, loss = loss_type) |>
    unlist() |>
    matrix(ncol = 50, byrow = TRUE) |>
    t() |>
    cbind(dat_list[c('result_pam_AT', 'result_pam_BC', 
                     'result_hclust_AT', 'result_hclust_BC')] |>
            unlist() |>
            matrix(ncol = 50, byrow = TRUE) |>
            t())
}

### Function: Calculate the mclustcomp
mclustcomp_across <- function(x){
  apply(x, 2, function(x){mclustcomp(x, sort(rep(1:2, 25)))[c(1, 5, 22), 2]}) |>
    t()
}

### Function: Calculate the unique clusters
nclus <- function(dat_list, burn_in){
  apply(dat_list$result_ZZ[-c(1:burn_in), ], 1, function(x){length(unique(x))}) |>
    cbind(apply(dat_list$result_DZ[-c(1:burn_in), ], 1, function(x){length(unique(x))}))
}

### Function: Proportion of Zero
pzero <- function(){}

### Plot the data: -------------------------------------------------------------
plot_list <- summarise_dat(dat)
plot_list[[5]]

### Time Analysis: -------------------------------------------------------------
#### Create the matrix for storing computational time
lapply(result, function(x){x$runtime}) |>
  unlist() |>
  matrix(nrow = 8) |>
  t() |>
  apply(2, ms_format, r = 4)

### Cluster Analysis (using VI() loss): ----------------------------------------
set.seed(1)
clus_assign_VI <- lapply(result, clus_list, burn_in = 10000, loss_type = "VI")
mclust_VI <- lapply(clus_assign_VI, mclustcomp_across)
lapply(1:8, function(x){simplify2array(mclust_VI)[x, , ] |> t()}) |>
  lapply(function(x){apply(x, 2, ms_format)})

lapply(clus_assign_VI, `[`, , 1:2) |>
  lapply(function(x){apply(x, 2, function(y){length(unique(y))})}) |>
  unlist() |>
  matrix(ncol = 30) |>
  t() |>
  apply(2, ms_format)

### Cluster Analysis (using binder() loss): ----------------------------------------
set.seed(1)
clus_assign_BD <- lapply(result, clus_list, burn_in = 10000, loss_type = "binder")
mclust_BD <- lapply(clus_assign_BD, mclustcomp_across)
lapply(1:8, function(x){simplify2array(mclust_BD)[x, , ] |> t()}) |>
  lapply(function(x){apply(x, 2, ms_format)})

lapply(clus_assign_BD, `[`, , 1:2) |>
  lapply(function(x){apply(x, 2, function(y){length(unique(y))})}) |>
  unlist() |>
  matrix(ncol = 30) |>
  t() |>
  apply(2, ms_format)

### Number of clusters for the non-parametric: ---------------------------------
lapply(clus_assign_BD, "[", , 5:8) |>
  lapply(function(x){apply(x, 2, function(y){length(unique(y))})}) |>
  unlist() |>
  matrix(ncol = 30) |>
  t() |>
  apply(2, ms_format)





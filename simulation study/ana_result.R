rm(list = ls())

### Required Libraries: --------------------------------------------------------
library(tidyverse)
library(ggplot2)

### Function: Summarize the simulated data
summarise_dat <- function(list_simDat){
  dat_plot <- vector("list", length(list_simDat))
  
  plot_list <- lapply(list_simDat, function(x){as.data.frame(x) |> 
      mutate(obs = paste0("OB", str_pad(1:50, 3, pad = "0"))) |>
      pivot_longer(cols = -obs) |>
      mutate(taxa_name = paste0("TX", str_pad(str_extract(name, "[:digit:]+$"), 3, pad = "0"))) |>
      ggplot(aes(taxa_name, obs, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low="white", high="palegreen3") +
      theme_minimal() +
      ## theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
      labs(x = "Variable (Taxa Count)", y = "Observation", fill = "Count")})
  
  prop_zero <- unlist(lapply(list_simDat, function(x){mean(x == 0)}))
  
  list(prop_zero = prop_zero, plot_list = plot_list)
  
}
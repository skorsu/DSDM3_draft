### Import function
path_mb <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data"
path_pc <- NULL

tryCatch({
  message("Try MB")
  source(paste0(path_mb, "/new_data_sim.R"))
  }, error=function(cond){
    message("Try PC")
    source(paste0(path_pc, "/new_data_sim.R"))
  })

### Create a plot
set.seed(12)
sim_plot(n = 100, J = 10, K = 5, scenario = 1, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_plot(n = 100, J = 10, K = 5, scenario = 2, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_plot(n = 100, J = 10, K = 5, scenario = 3, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_plot(n = 100, J = 10, K = 5, scenario = 4, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)

### Try with our model
set.seed(12)
sim_list <- dat_sim(n = 100, J = 10, K = 5, scenario = 1, xi_conc = 0.1, pi_gm = 0.75, sum_z = 2500)
sim_list$ci

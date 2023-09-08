library(salso)
library(mcclust)
library(gridExtra)

### Import the function
im_path <- "/Users/kevin-imac/Desktop/Github - Repo/ClusterZI/data/"
mb_path <- "/Users/kevinkvp/Desktop/Github Repo/ClusterZI/data/"
path <- NULL
if(dir.exists(mb_path)){
  path <- mb_path
} else {
  path <- im_path
}

source(paste0(path, "data_sim_DM.R"))

### Data Simulation
pat_mat1 <- diag(20)[c(2, 6, 10, 1, 3), ]
set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 10, pi_gm = c(0.75), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(754)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 1, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time
table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:10000), ], maxNClusters = 10)))
plot(apply(tt$assign, 1, function(x){length(unique(x))}), type = "l")
rm(pat_mat1, sim_dat, summa, tt)

### Data Simulation
pat_mat1 <- diag(20)[c(2, 4, 6, 10, 15), ]
set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 100, pi_gm = c(0.85), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(61)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 10, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time



table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:5000), ], maxNClusters = 10)))

rm(pat_mat1, sim_dat, summa, tt)

### Data Simulation
pat_mat1 <- matrix(0, nrow = 3, ncol = 20)
pat_mat1[1, 1:7] <- 1
pat_mat1[2, 8:14] <- 1
pat_mat1[3, 15:20] <- 1
set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 10, pi_gm = c(0.75), 
                 pi_c = c(1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(2500)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 10, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time

table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:5000), ], maxNClusters = 10)))

rm(pat_mat1, sim_dat, summa, tt)

### Data Simulation
pat_mat1 <- matrix(0, nrow = 3, ncol = 20)
pat_mat1[1, 1:7] <- 1
pat_mat1[2, 8:14] <- 1
pat_mat1[3, 15:20] <- 1
set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 100, pi_gm = c(0.85), 
                 pi_c = c(1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(2500)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 10, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time

table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:5000), ], maxNClusters = 10)))

rm(pat_mat1, sim_dat, summa, tt)

### Data Simulation
pat_mat1 <- diag(20)[1:5, ]
pat_mat1[3, 3] <- 0
pat_mat1[4, 4] <- 0
pat_mat1[5, 5] <- 0
pat_mat1[3, c(10, 12)] <- 1
pat_mat1[4, c(11, 12)] <- 1
pat_mat1[5, c(10, 11)] <- 1

set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 10, pi_gm = c(0.75), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(2500)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 10, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time

table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:5000), ], maxNClusters = 10)))

rm(pat_mat1, sim_dat, summa, tt)

### Data Simulation
pat_mat1 <- diag(20)[1:5, ]
pat_mat1[3, 3] <- 0
pat_mat1[4, 4] <- 0
pat_mat1[5, 5] <- 0
pat_mat1[3, c(10, 12)] <- 1
pat_mat1[4, c(11, 12)] <- 1
pat_mat1[5, c(10, 11)] <- 1

set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 100, pi_gm = c(0.85), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 500, z_sum_U = 1000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(2500)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 10, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time

table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:5000), ], maxNClusters = 10)))

rm(pat_mat1, sim_dat, summa, tt)

### Data Simulation
pat_mat1 <- diag(50)[1:5, ]

set.seed(72)
sim_dat <- simDM(n = 200, pattern = pat_mat1, xi_conc = 100, pi_gm = c(0.85), 
                 pi_c = c(1, 1, 1, 1, 1), z_sum_L = 3000, z_sum_U = 5000, 
                 theta = 0.01)
summa <- simDM_sum(sim_dat)
grid.arrange(grobs = summa$plot)

### Run the model
set.seed(2500)
start_time <- Sys.time()
tt <- ZIDM_ZIDM(iter = 15000, K_max = 10, z = sim_dat$z, theta_vec = rep(1, 10), 
                launch_iter = 5, MH_var = 10, mu = 0, s2 = 1, r0g = 1, r1g = 1, 
                r0c = 1, r1c = 1, print_iter = 1000)
Sys.time() - start_time

table(sim_dat$ci, as.numeric(salso(tt$assign[-c(1:5000), ], maxNClusters = 10)))

plot(apply(tt$assign, 1, function(x){length(unique(x))}), type = "l")


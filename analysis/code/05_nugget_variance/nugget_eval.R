library(argofda)
library(Matrix)
library(fda)
library(sparseinv)
library(dplyr)
library(fields)
library(parallel)
setwd('/home/dyarger/argofda')
load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

latvals <- seq(-79.5, 79.5, by = 1)
lonvals <- seq(20.5, 379.5, by = 1)

longvals <- ifelse(lonvals > 180, lonvals - 360, lonvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$order <- 1:nrow(grid_to_compute)

knots <- seq(0, 2000, length.out = 100)
penalty_matrix <- fda::bsplinepen(basisobj = create.bspline.basis(knots))

knots_pc <- seq(0, 2000, length.out = 100)
h_space <- 550

profYearAggr <- as.factor(profYearAggr)

cv <- function(index, h_space, K1, K2) {
  a <- proc.time()

  lat <- grid_to_compute[index,2]
  long <- grid_to_compute[index,1]

  df <- get_profile_data_subset(lat = lat, long = long, day = 45.25, h_time = 45.25, RG_def = FALSE,
                                h_space = h_space, mode = 'all',  min_prof=  5,
                                exclusion= TRUE)

  if (!(lat %in% grid_pc_temp$lat[grid_pc_temp$long == long] &
        lat %in% grid_pc_psal$lat[grid_pc_psal$long == long]  ) |
      is.null(df) | is.null(nrow(df))) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }

  pressure_min <- df %>% group_by(profile) %>% summarise(min(pressure), .groups = 'drop')
  pressure_max <- df %>% group_by(profile) %>% summarise(max(pressure), .groups = 'drop')
  short_profiles <- (1:nrow(pressure_max))[which(pressure_max[,2] - pressure_min[,2] < 800)]
  prof_lengths <- df %>% group_by(profile) %>% summarise(n(), .groups = 'drop')
  small_profiles <-(1:nrow(prof_lengths))[which(prof_lengths[,2] < 15)]
  df <- df[!(df$profile %in% short_profiles) & !(df$profile %in% small_profiles),]
  df$profile <- as.numeric(as.factor(df$profile))

  if (nrow(df) == 0 ) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }

  temp_pc <- coefs_pc_temp[[which(grid_pc_temp$long == long &
                                    grid_pc_temp$lat == lat)]][[1]][,1:K1]
  psal_pc <- coefs_pc_psal[[which(grid_pc_psal$long == long &
                                    grid_pc_psal$lat == lat)]][[1]][,1:K2]

  # get first few pc at this location
  # Calculate scores for each profile
  t_t <- get_scores(df, temp_pc, K1, type = 'temperature', knots = knots_pc)
  df$resids_temp <- unlist(t_t[[2]])
  rm(t_t)
  max_prof <- max(df$profile)

  s_s <- get_scores(df, psal_pc, K2, type = 'salinity', knots = knots_pc)
  df$resids_psal <- unlist(s_s[[2]])
  rm(s_s)
  basis <- create.bspline.basis(knots)
  # Fit spline to the variance of the residuals
  nr <- nrow(df)
  gc()
  if (nr %%2 ==0) {
    phi1 <- as(fda::eval.basis(basis, evalarg = df$pressure[1:(nr/2)]), 'sparseMatrix')
    phi2 <- as(fda::eval.basis(basis, evalarg = df$pressure[(nr/2 + 1):nr]), 'sparseMatrix')
  } else {
    phi1 <- as(fda::eval.basis(basis, evalarg = df$pressure[1:((nr-1)/2)]), 'sparseMatrix')
    phi2 <- as(fda::eval.basis(basis, evalarg = df$pressure[((nr-1)/2 + 1):nr]), 'sparseMatrix')
  }
  phi <- rbind(phi1, phi2)
  rm(phi1, phi2)
  gc()
 # phi <- as(fda::eval.basis(basis = create.bspline.basis(knots), evalarg = df$pressure), 'sparseMatrix')
  prof_lengths <- aggregate(df$pressure, list(df$profile), length)
  W <- argofda::construct_working_cov(p = df$pressure,
                                   weights = df$wt/rep(prof_lengths$x, prof_lengths$x),
                                   rho = .9,
                                   profile = df$profile, specify_cor = T)
  model_set <- model_set_up(covariates = matrix(nrow = nrow(df), ncol = 1, 1),
                            lambdas = 10^3, Phi = phi, Psi = penalty_matrix, W = W,
                            Y = log(df$resids_temp^2))
  var_temp <- complete_optim(U = model_set$U, V = model_set$V,W = W, Omega = model_set$Omega,
                             Y = log(df$resids_temp^2), Phi_all =model_set$Phi_all, lower_optim = 10^-5, upper_optim = 10^5,
                             selection_vector = model_set$reorder, knots =knots,
                             n_basis = length(knots)+2 )

  delayed_indexes <- df$mode == 'D'
  if (sum(delayed_indexes) == 0) {
    b <- proc.time()
    final_list <- list(var_temp, NULL, NULL,
                       c(index, (b-a)[3],
                         max_prof))
    print(final_list[[4]])
    return(final_list)
  }
  phi_psal <- phi[delayed_indexes,]
  model_set <- model_set_up(covariates = matrix(nrow = sum(delayed_indexes), ncol = 1, 1),
                            lambdas = 10^3, Phi = phi_psal, Psi = penalty_matrix, W = W[delayed_indexes, delayed_indexes],
                            Y = log(df$resids_psal[delayed_indexes]^2))
  var_psal <- complete_optim(U = model_set$U, V = model_set$V,W = W[delayed_indexes, delayed_indexes], Omega = model_set$Omega,
                             Y =log(df$resids_psal[delayed_indexes]^2), Phi_all =model_set$Phi_all, lower_optim = 10^-5, upper_optim = 10^5,
                             selection_vector = model_set$reorder, knots =knots,
                             n_basis = length(knots)+2 )

  b <- proc.time()
  final_list <- list(var_temp, var_psal, NULL,
                     c(index, (b-a)[3],
                       max_prof))
  print(final_list[[4]])
  return(final_list)
}
sequence <- as.list(seq(1, nrow(grid_to_compute), by = 1))
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
increments <- 500

profile_start <- 1 + (array_id - 1) * increments
#indexes <- seq(1, nrow(grid_to_compute), by = 1)
indexes <- profile_start:(profile_start+increments-1)

grid_to_compute_only_ours <- grid_to_compute[grid_to_compute$order %in% indexes,]

load('analysis/results/psal_cov_pca.RData')
coefs_pc_psal <- coefs_pc_psal[grid_pc_psal$long %in% grid_to_compute_only_ours$long &
                                 grid_pc_psal$lat %in% grid_to_compute_only_ours$lat]
grid_pc_psal <- grid_pc_psal[grid_pc_psal$long %in% grid_to_compute_only_ours$long &
                               grid_pc_psal$lat %in% grid_to_compute_only_ours$lat,]
load('analysis/results/temp_cov_pca.RData')
coefs_pc_temp <- coefs_pc_temp[grid_pc_temp$long %in% grid_to_compute_only_ours$long &
                            grid_pc_temp$lat %in% grid_to_compute_only_ours$lat]
grid_pc_temp <- grid_pc_temp[grid_pc_temp$long %in% grid_to_compute_only_ours$long &
                               grid_pc_temp$lat %in% grid_to_compute_only_ours$lat,]

if (array_id == ceiling(nrow(grid_to_compute)/increments)) {
  indexes <- indexes[sapply(sequence[indexes], function(x) {!is.null(x)})]
}
value <- indexes[length(indexes)]
# Run for a large number of grid points
object_name <- 'nugget_var'
command <- paste(object_name, ' <- lapply(indexes, cv,',
                 'K1 = 10, K2 = 10,',
               #  'mc.preschedule = F, mc.cores = 4,',
                 'h_space = h_space',
                 ')', sep="")
a <- proc.time()
eval(parse(text=command))
b <- proc.time();b-a

if (nchar(value) == 2) {
  value <- paste0('000',value)
} else if (nchar(value) == 3) {
  value <- paste0('00',value)
} else if(nchar(value) == 4) {
  value <- paste0('0',value)
}
file.to.save <- paste0('analysis/results/nugget_variance/nugget_var_', value, '.RData')
eval(parse(text=paste('save(', object_name, ', file = file.to.save)', sep = '')))

q()




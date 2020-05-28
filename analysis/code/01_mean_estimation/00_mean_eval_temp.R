
# implement the mean estimation on a 1 by 1 grid
library(argofda)
library(Matrix)
library(fda)
library(sparseinv)

load('analysis/data/jan_march_data.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

# create grid
load('analysis/data/RG_Defined.RData')
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$order <- 1:nrow(grid_to_compute)
grid_to_compute_merge <- merge(grid_to_compute, RG_defined_long, all = T, sort = FALSE)
grid_to_compute_merge$RG_defined <- ifelse(is.na(grid_to_compute_merge$value),
                                           FALSE,
                                           TRUE)
grid_to_compute_merge <- grid_to_compute_merge[order(grid_to_compute_merge$order), c(1,2,5)]

# build function for mclapply
# index determines which grid point we use
knots <- seq(0, 2000, length.out = 200)
penalty_mat <- as(fda::bsplinepen(create.bspline.basis(knots)), 'sparseMatrix')
h_space <- 900
mean_est <- function(index, knots, penalty_mat, h_space) {
  a <- proc.time()
  lat <- grid_to_compute[index,2]
  long <- grid_to_compute[index,1]
  # Get data
  RG_def <- grid_to_compute_merge[index,3]
  df <- get_profile_data_subset(lat = lat, long = long, day = 45.25, RG_def = RG_def,
                                h_space = h_space, h_time = 45.25, mode = 'all', #months = 1:3,
                                min_prof= 1, exclusion= TRUE)
  if (is.null(df)) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }
  n_temp <- nrow(df)
  n_prof_temp <- max(df$profile)

  df$pressure <- round(as.numeric(df$pressure), digits = 2)
  h_space_evaluated <- df$h_space[1]
  basis_eval <- as(fda::eval.basis(create.bspline.basis(knots), evalarg = df$pressure), 'sparseMatrix') #penalty matrix

  t_s <- df[, 'fracofyear'] - 45.25
  df$year <- factor(df$year, levels = 2007:2016)
  yearly_term <- sparse.model.matrix(model.frame(~df$year-1, df ))
  covariates <- cbind(yearly_term, df$EW_dist_km, df$NS_dist_km, df$EW_dist_km^2, df$NS_dist_km^2,
                      df$NS_dist_km*df$EW_dist_km,
                      t_s, t_s^2)
  W <- construct_working_cov(p = df$pressure, weights = df$wt, profile = df$profile,
                             rho = .999,
                             specify_cor = T)
  mod <- model_set_up(covariates = covariates, lambdas = 10^c(rep(2, 10), 10,10,15,15,15,7,11)* 10^-4,
                      Phi =basis_eval,Psi = penalty_mat, W = W,Y = df$temperature)
  spline_results <- complete_optim(U = mod[['U']], V = mod[['V']], W = W,
                                   Omega = mod[['Omega']], Y = df$temperature,
                                   Phi_all = mod[['Phi_all']],
                                   lower_optim = 10^-5, upper_optim = 10^5, tol = 10^-2,
                                   selection_vector = mod[['reorder']],
                                   n_basis = ncol(basis_eval), knots = knots)
  df_prof <- df[!duplicated(df$profile),]
  n_years <- table(df_prof$year)
  b <- proc.time()
  final_list <- list(spline_results,
                     n_years,
                     c(index, (b-a)[3], n_temp,
                       n_prof_temp, h_space_evaluated))
  print(final_list[[3]])
  return(final_list)
}

# set up array job issues and how things are proportioned
library(parallel)
sequence <- as.list(seq(1, nrow(grid_to_compute), by = 1))
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# array id from 1-116, each with 500 grid points (57590 overall)
increments <- 500

profile_start <- 1 + (array_id - 1) * increments
indexes <- profile_start:(profile_start+increments-1)

if (array_id == ceiling(nrow(grid_to_compute)/500)) {
  indexes <- indexes[sapply(sequence[indexes], function(x) {!is.null(x)})]
}

if (array_id == ceiling(nrow(grid_to_compute)/500)) {
  value <- sequence[[indexes[length(indexes)]]]-1
} else {
  value <- sequence[[indexes[length(indexes)]+1]]-1
}
object_name <-  paste('Grid_Pred', value, sep = '')

command <- paste(object_name, ' <- mclapply(indexes, mean_est, penalty_mat = penalty_mat, knots = knots, h_space = 900,',
                 'mc.preschedule = FALSE, mc.cores = 4',
                 ')', sep="")
a <- proc.time()
eval(parse(text=command)) # calculations are run here
b <- proc.time();b-a

# save the results

if (nchar(value) == 3) {
  value <- paste0('00',value)
} else if(nchar(value) == 4) {
  value <- paste0('0',value)
}
file.to.save <- paste0('analysis/results/mean_one_stage/mean_one_stage_temp_', value, '.Rdata')
eval(parse(text=paste0('save(', object_name, ', file = file.to.save)')))


library(gsw)
library(fda)
library(dplyr)
library(argofda)
load('analysis/results/joint_TS_20/preds.RData') # conditional predictions
source('analysis/code/08_mld/mld_source.R') # load mld funs

# set seed for each job array
rseed <- 348
set.seed(rseed)
all_seeds <- sample(1:(10^6), 116)
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(all_seeds[array_id])


df_preds <- df_preds[df_preds$index %in% (500*(array_id-1)+ 1):(500*array_id),]
gc()
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

grid_to_compute <- expand.grid('long' = longvals,'lat' = latvals)
grid_to_compute$index <- 1:nrow(grid_to_compute)
grid_to_compute <- grid_to_compute[(500*(array_id-1)+ 1):(500*array_id),]
### MLD ###

# load pcs and only get the ones we need
load('analysis/results/psal_one_stage_cov_pca.RData')
grid_merge <- left_join(grid_pc_psal, grid_to_compute)
coefs_pc_psal <- coefs_pc_psal[!is.na(grid_merge$index)]
grid_pc_psal <- grid_pc_psal[!is.na(grid_merge$index),]
load('analysis/results/temp_one_stage_cov_pca.RData')
grid_merge <- left_join(grid_pc_temp, grid_to_compute)
coefs_pc_temp <- coefs_pc_temp[!is.na(grid_merge$index)]
grid_pc_temp <- grid_pc_temp[!is.na(grid_merge$index),]

### get mean files
mean_files_temp <- list.files('analysis/results/mean_estimation/', pattern = 'temp')
mean_files_psal <- list.files('analysis/results/mean_estimation/', pattern = 'psal')
final_list <- list()

delta <- 2
pressure_vec <- seq(10, 2000, by = delta)
# for each mean file, do conditional simulations
pb_mld <- function(j, mean_files_temp, mean_files_psal, pressure_vec, K1, K2) {
  # load mean functions
  load(paste0('analysis/results/mean_estimation/', mean_files_temp[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  temp_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[1]])
    }
  })
  eval(parse(text=paste0('rm(', file, ')')))
  load(paste0('analysis/results/mean_estimation/', mean_files_psal[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  psal_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[1]])
    }
  })
  eval(parse(text=paste0('rm(', file, ')')))

  # set up relevant calculations for each grid point and year
  have_results <- which(sapply(temp_funs, length) != 0  & sapply(psal_funs, length) != 0)
  indexes_have_mean <- (1:500 + (j-1)*500)[have_results]
  indexes_have_both <- unique(df_preds$index[df_preds$index %in% indexes_have_mean])
  have_results_both <- which((1:500 + (j-1)*500) %in% indexes_have_both)
  df_indexes <- df_preds[df_preds$index %in% indexes_have_both,]
  print(nrow(df_indexes))
  mld <- lapply(1:nrow(df_indexes), function(x) {
    # compute mean and pcs at pressure values
    index <- df_indexes[x, 'index']
    year <- df_indexes[x, 'year']
    mean_temp_one <- predict(temp_funs[have_results_both][[
      which(indexes_have_both == index)]],p_vals = pressure_vec,index = year - 2007)
    mean_psal_one <- predict(psal_funs[have_results_both][[
      which(indexes_have_both == index)]],p_vals = pressure_vec,index = year - 2007)
    long <- df_indexes[x,'long']
    lat <- df_indexes[x,'lat']
    pc <- cbind(coefs_pc_temp[[which(grid_pc_temp$long == long &
                                       grid_pc_temp$lat == lat)]][[1]][,1:K1],
                coefs_pc_psal[[which(grid_pc_psal$long == long &
                                       grid_pc_psal$lat == lat)]][[1]][,1:K2])
    pred_pc <- sapply(as.data.frame(as.matrix(pc)), function(z) {
      predict(make_fit_object(coefficients = z, knots = knots_pc),
              pressure_vec)$y
    })


    # derivatives
    mu_temp_deriv <- predict(temp_funs[have_results_both][[
      which(indexes_have_both == index)]],p_vals = pressure_vec,index = year - 2007, deriv = 1)
    mu_psal_deriv <- predict(temp_funs[have_results_both][[
      which(indexes_have_both == index)]],p_vals = pressure_vec,index = year - 2007, deriv = 1)
    phi_deriv <- sapply(as.data.frame(as.matrix(pc)), function(z) {
      predict(make_fit_object(coefficients = z, knots = knots_pc),
              pressure_vec, deriv = 1)$y
    })

    # set up conditional distribution
    scores <- unlist(df_indexes[x, c(paste0('T', 1:K1), paste0('S', 1:K2))])
    var_scores <- matrix(nrow = K1 + K2, ncol = K1 + K2,
                         unlist(df_indexes[x, paste0('condvar',
                                                     paste0(expand.grid(1:(K1+K2), 1:(K1+K2))[,1], '_',
                                                            expand.grid(1:(K1+K2), 1:(K1+K2))[,2]))]))
    var_scores_eigen <- eigen(var_scores)
    var_scores_lambda <- ifelse(var_scores_eigen$values > 0 ,var_scores_eigen$values,
                                0)
    var_scores_half <- var_scores_eigen$vectors %*% diag(sqrt(var_scores_lambda))

    B <- 1000
    # print progress
    if (x %% 100 ==0) {
      print(x/nrow(df_indexes))
    }
    # run MLD algorithms B times
    bs_results <- t(sapply(1:B, bs_fun, scores = scores, phi_deriv = phi_deriv, pred_pc = pred_pc, var_scores_half =
                             var_scores_half, long = long, lat = lat, mean_psal_one = mean_psal_one, mean_temp_one = mean_temp_one,
                           pressure_vec = pressure_vec, mu_temp_deriv = mu_temp_deriv, mu_psal_deriv = mu_psal_deriv, K1 = K1, K2 = K2,
                           delta = delta))
  })
  return(list(df_indexes[, c('long', 'lat', 'index', 'year')],mld))
}
library(parallel)

final_list <- pb_mld(array_id,
                       mean_files_temp = mean_files_temp, mean_files_psal = mean_files_psal,
                       pressure_vec = pressure_vec, K1 = 10, K2 = 10)
save(final_list, file = paste0('analysis/results/mld/mld_', array_id,
                                        '.RData'))
q()

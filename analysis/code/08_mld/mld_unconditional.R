
library(gsw)
library(fda)
library(dplyr)
library(argofda)
setwd('/home/dyarger/argofda/')
load('analysis/results/joint_TS_20/params.RData') # conditional predictions
#source('analysis/code/08_mld/mld_source.R') # load mld funs

bs_fun <- function(z, scores, phi_deriv, pred_pc, var_scores_half, long, lat, mean_psal_one, mean_temp_one,pressure_vec,
                   mu_temp_deriv, mu_psal_deriv, K1, K2, delta, pred_interval = FALSE,
                   nugget_var_temp = rep(0, length(pressure_vec)),
                   nugget_var_psal = rep(0, length(pressure_vec))) {
  coef_boot <- as.double(scores + var_scores_half %*% rnorm(K1 + K2))
  pred_temp <- .rowSums(pred_pc[,1:K1] %*% diag(coef_boot[1:K1]), m = length(pressure_vec), n = K1)
  pred_psal <- .rowSums(pred_pc[,(K1+1):(K1+K2)] %*% diag(coef_boot[(K1+1):(K1+K2)]), m = length(pressure_vec), n = K2)
  full_psal <- mean_psal_one + pred_psal
  full_temp <- mean_temp_one + pred_temp
  bs_fun_preds(full_temp, full_psal, pressure_vec, long, lat, delta)
}
bs_fun_preds <- function(full_temp, full_psal, pressure_vec, long, lat, delta) {
  SA <- gsw_SA_from_SP(SP = full_psal,p = pressure_vec, longitude = long, latitude = lat)
  pdens <- gsw::gsw_pot_rho_t_exact(SA = SA, t = full_temp,
                                    p = pressure_vec, p_ref = 0)
  # potential density - variable
  pdens_thresh_try <- gsw::gsw_pot_rho_t_exact(SA = SA[1], t = full_temp[1] - .2,
                                               p = pressure_vec[1], p_ref = 0) - pdens[1]
  thresh_var_boot_d <- pressure_vec[which(pdens - pdens[1]  >= pdens_thresh_try)[1]]

  # potential density - constant
  pdens_sub <- pdens - pdens[1]
  thresh_boot_d <-  pressure_vec[which(pdens_sub  >= .03)[1]]
  drho_dp_all <- (pdens[2:length(pdens)] - pdens[1:(length(pdens)-1)])/delta
  drho_dp_all <- c(drho_dp_all[1] , (drho_dp_all[2:(length(drho_dp_all))] + drho_dp_all[1:(length(drho_dp_all)-1)])/2,
                   drho_dp_all[length(drho_dp_all)])
  # 0.0005 to 0.05 kg m^4 are generally used, we use .001
  deriv_boot_d <- pressure_vec[which(drho_dp_all >= 0.001)[1]]

  # Temperature
  thresh_boot <-  pressure_vec[which(full_temp - full_temp[1] <= -.2)[1]]
  #deriv_boot <- pressure_vec[which(full_deriv_temp <= -.025)[1]]

  c(thresh_var_boot_d, thresh_boot_d, deriv_boot_d, thresh_boot, NA)
}

# set seed for each job array
rseed <- 348
set.seed(rseed)
all_seeds <- sample(1:(10^6), 116)
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#array_id <- 2
set.seed(all_seeds[array_id])


df_params_only <- df_params_only[df_params_only$index %in% (500*(array_id-1)+ 1):(500*array_id),]
gc()
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

grid_to_compute <- expand.grid('long' = longvals,'lat' = latvals)
grid_to_compute$index <- 1:nrow(grid_to_compute)
grid_to_compute <- grid_to_compute[(500*(array_id-1)+ 1):(500*array_id),]
### MLD ###

# load pcs and only get the ones we need
load('analysis/results/psal_cov_pca.RData')
grid_merge <- left_join(grid_pc_psal, grid_to_compute)
coefs_pc_psal <- coefs_pc_psal[!is.na(grid_merge$index)]
grid_pc_psal <- grid_pc_psal[!is.na(grid_merge$index),]
load('analysis/results/temp_cov_pca.RData')
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
  indexes_have_both <- unique(df_params_only$index[df_params_only$index %in% indexes_have_mean])
  have_results_both <- which((1:500 + (j-1)*500) %in% indexes_have_both)
  df_indexes <- df_params_only[df_params_only$index %in% indexes_have_both,]
  print(nrow(df_indexes))
  mld <- lapply(1:nrow(df_indexes), function(x) {
    # compute mean and pcs at pressure values
    index <- df_indexes[x, 'index']
    t_funs <- temp_funs[have_results_both][[
      which(indexes_have_both == index)]]
    mean_temp_one <- sapply(1:10, function(x)
      {predict.local_function(t_funs, p_vals = pressure_vec,index = x-1)})
    mean_temp_one <- rowMeans(mean_temp_one)
    p_funs <- psal_funs[have_results_both][[
      which(indexes_have_both == index)]]
    mean_psal_one <- sapply(1:10, function(x)
    {predict.local_function(p_funs, p_vals = pressure_vec,index = x-1)})
    mean_psal_one <- rowMeans(mean_psal_one)

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

    # set up conditional distribution
    scores <- rep(0, K1+ K2)
    var_scores <- matrix(nrow = K1 + K2, ncol = K1 + K2,
                         unlist(df_indexes[x, paste0('var',
                                                     paste0(expand.grid(1:(K1+K2), 1:(K1+K2))[,1], '_',
                                                            expand.grid(1:(K1+K2), 1:(K1+K2))[,2]))]))
    var_scores_eigen <- eigen(var_scores)
    var_scores_lambda <- ifelse(var_scores_eigen$values > 0 ,var_scores_eigen$values,
                                0)
    var_scores_half <- var_scores_eigen$vectors %*% diag(sqrt(var_scores_lambda))

    B <- 1000
    # print progress
    #print(x)
    if (x %% 100 ==0) {
      print(x/nrow(df_indexes))
    }
    # run MLD algorithms B times
    bs_results <- t(sapply(1:B, bs_fun, scores = scores, phi_deriv = phi_deriv, pred_pc = pred_pc, var_scores_half =
                             var_scores_half, long = long, lat = lat, mean_psal_one = mean_psal_one, mean_temp_one = mean_temp_one,
                           pressure_vec = pressure_vec, mu_temp_deriv = mu_temp_deriv, mu_psal_deriv = mu_psal_deriv, K1 = K1, K2 = K2,
                           delta = delta))
  })
  return(list(df_indexes[, c('long', 'lat', 'index')],mld))
}
library(parallel)
final_list <- pb_mld(array_id,
                     mean_files_temp = mean_files_temp, mean_files_psal = mean_files_psal,
                     pressure_vec = pressure_vec, K1 = 10, K2 = 10)
save(final_list, file = paste0('analysis/results/mld/mld_uncond_', array_id,
                               '.RData'))
q()

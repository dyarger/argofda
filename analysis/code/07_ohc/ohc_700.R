
# this function takes a location and year, returns ohc information
ohc_fun <- function(x) {
  index_val <- which(grid_to_compute_all_year[,'index'] == df_indexes[x,'index'])
  # get mean
  pred_vals_temp <- sapply(0:9, function(x) predict.local_function(temp_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec))
  avg_pred_temp <- rowMeans(pred_vals_temp)
  pred_vals_psal <- sapply(0:9, function(x) predict.local_function(psal_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec))
  avg_pred_psal <- rowMeans(pred_vals_psal)

  mean_temp_one <- pred_vals_temp[,df_indexes[x, 'year'] - 2006]
  mean_psal_one <- pred_vals_psal[,df_indexes[x, 'year'] - 2006]
  long <- df_indexes[x,'long']
  lat <- df_indexes[x,'lat']

  # Predictions
  scores_temp <- unlist(df_indexes[x, paste0('T', 1:K1)])
  scores_psal <- unlist(df_indexes[x, paste0('S', 1:K2)])
  pc_temp <- coefs_pc_temp[[which(grid_pc_temp$long == long &
                                    grid_pc_temp$lat == lat)]][[1]][,1:K1]
  pc_psal <-  coefs_pc_psal[[which(grid_pc_psal$long == long &
                                     grid_pc_psal$lat == lat)]][[1]][,1:K2]
  pred_pc_temp <- sapply(as.data.frame(as.matrix(pc_temp)), function(z) {
    predict(make_fit_object(coefficients = z, knots = knots_pc),
            pressure_vec)$y
  })
  pred_pc_psal <- sapply(as.data.frame(as.matrix(pc_psal)), function(z) {
    predict(make_fit_object(coefficients = z, knots = knots_pc),
            pressure_vec)$y
  })
  pred_temp <- rowSums(sapply(1:K1, function(x) { scores_temp[x] * pred_pc_temp[,x]}))
  pred_psal <- rowSums(sapply(1:K2, function(x) { scores_psal[x] * pred_pc_psal[,x]}))

  # compute conservative temperature (for mean, year-avg mean, and prediction)
  absolute_salinity_pred <- gsw::gsw_SA_from_SP(SP = mean_psal_one+pred_psal, p = pressure_vec,
                                                longitude = long, latitude = lat)
  conservative_temp_pred <- gsw::gsw_CT_from_t(SA = absolute_salinity_pred, t = mean_temp_one + pred_temp,
                                               p = pressure_vec)
  absolute_salinity_mean <- gsw::gsw_SA_from_SP(SP = mean_psal_one, p = pressure_vec,
                                                longitude = long, latitude = lat)
  conservative_temp_mean <- gsw::gsw_CT_from_t(SA = absolute_salinity_mean, t = mean_temp_one, p = pressure_vec)
  absolute_salinity_avg_mean <- gsw::gsw_SA_from_SP(SP = avg_pred_psal, p = pressure_vec,
                                                    longitude = long, latitude = lat)
  conservative_temp_avg_mean <- gsw::gsw_CT_from_t(SA = absolute_salinity_avg_mean, t = avg_pred_temp, p = pressure_vec)

  conservative_temp_anomaly <- conservative_temp_pred +  273.15 - conservative_temp_mean -273.15
  conservative_temp_avg_anomaly <- conservative_temp_pred +  273.15 - conservative_temp_avg_mean -273.15

  # Kelvin for OHC predictions
  ohc_mean <- c_p/g * sum((conservative_temp_mean + 273.15) * delta)
  ohc_pred <- c_p/g * sum((conservative_temp_pred + 273.15) * delta) #mean + anomaly
  ohc_anomaly <- c_p/g * sum(conservative_temp_anomaly * delta)
  ohc_avg_anomaly <- c_p/g * sum(conservative_temp_avg_anomaly * delta)

  # Variance
  var_all <- matrix(nrow = K1+K2, ncol = K1+K2, unlist(df_indexes[x, paste0('condvar', rep(1:(K1+K2), each = K1+K2),
                                                                            '_', rep(1:(K1+K2), times = K1+K1))]))
  deriv_pred <- gsw_CT_first_derivatives_wrt_t_exact(SA = absolute_salinity_pred,
                                                     t = mean_temp_one + pred_temp,
                                                     p = pressure_vec)
  SA_deriv <- deriv_pred$CT_SA_wrt_t #absolute salinity, not practical salinity
  T_deriv <- deriv_pred$CT_t_wrt_t
  derivs <- cbind(T_deriv, SA_deriv)

  sigma11 <- pred_pc_temp %*% var_all[1:K1,1:K1] %*% t(pred_pc_temp)
  sigma12 <- pred_pc_temp %*% var_all[1:K1, (K1 + 1):(K1+K2)] %*% t(pred_pc_psal)
  sigma22 <- pred_pc_psal %*% var_all[(K1 + 1):(K1+K2), (K1 + 1):(K1+K2)] %*% t(pred_pc_psal)

  # prepare to sum over the 2-d grid in pressure for T and S
  ut <- upper.tri(sigma11, diag = T)
  sigma11ut <- sigma11[ut]
  sigma12ut <- sigma12[ut]
  sigma22ut <- sigma22[ut]
  sigma_diag_vals <- c(sigma11ut, sigma22ut)[selection_vector]
  sigma_band_vals <- c(sigma12ut, rep(0, n_grid * (n_grid+1)/2))[selection_vector]
  sigma_sparse@x <- c(sigma_diag_vals, sigma_band_vals[-length(sigma_band_vals)])

  # derivatives for delta method
  derivs_1 <- derivs[grid_temporary[,1],][selection_vector]
  derivs_2 <- derivs[grid_temporary[,2],][selection_vector]
  var_est <- t(derivs_1) %*% sigma_sparse %*% derivs_2 # sum over 2-d grid

  # We do the diagonal again, since we only did the upper triangle
  derivs_diag <- derivs[diag_indices,][selection_vector_diag]
  sigma_diag_diag_vals <- c(diag(sigma11), diag(sigma22))[selection_vector_diag]
  sigma_diag_band_vals <- c(diag(sigma12), rep(0, nrow(sigma11)))[selection_vector_diag]
  sigma_sparse_diag@x <- c(sigma_diag_diag_vals, sigma_diag_band_vals[-length(sigma_diag_band_vals)])
  var_diag_est <- t(derivs_diag) %*% sigma_sparse_diag %*% derivs_diag

  # combine the two and compute
  var_est_new <- 2*var_est[1] - var_diag_est[1] # 2 * upper tri - diag
  var_ohc <- c_p^2/g^2 * var_est_new[1] * delta * delta
  c(ohc_mean,ohc_pred, ohc_anomaly, ohc_avg_anomaly, sqrt(var_ohc))
}

library(ggplot2)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))
library(fda)
library(argofda)
library(gsw)
library(parallel)

load('analysis/results/joint_TS_20/preds.RData')

load('analysis/results/temp_one_stage_cov_pca.RData')
load('analysis/results/psal_one_stage_cov_pca.RData')
basis <- create.bspline.basis(breaks =knots_pc)


### get mean functions
calc_files_temp <- list.files('analysis/results/mean_one_stage/', pattern = 'temp')
calc_files_psal <- list.files('analysis/results/mean_one_stage/', pattern = 'psal')
final_list <- list()

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)

delta <- 3 # approximation grid in pressure
pressure_vec <- seq(0, 700, by = delta)
n_grid <- length(pressure_vec)
pressure_grid <- expand.grid('p1' = pressure_vec, 'p2' = pressure_vec)
grid_temporary <- expand.grid(1:length(pressure_vec), 1:length(pressure_vec))
grid_temporary <- as.matrix(grid_temporary[pressure_grid$p1 <= pressure_grid$p2,])
pressure_grid <- as.matrix(pressure_grid[pressure_grid$p1 <= pressure_grid$p2,])
diag_indices <- grid_temporary[grid_temporary[,1] == grid_temporary[,2],1]

c_p <- 3850 # Joules per kg per degrees celcius  - specific heat of seawater
g <- (1025)^(-1) # 1/density of seawater in (kg/m^3)^-1

for (j in 1:length(calc_files_temp)) {
  # get mean functions
  load(paste0('analysis/results/mean_one_stage/', calc_files_temp[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  temp_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][1] == 'no profiles in range') {return(NULL)} else {return(x[[1]])}
  })
  eval(parse(text=paste0('rm(', file, ')')))
  load(paste0('analysis/results/mean_one_stage/', calc_files_psal[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  psal_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][1] == 'no profiles in range') {return(NULL)} else { return(x[[1]]) }
  })
  eval(parse(text=paste0('rm(', file, ')')))
  have_results <- which(sapply(temp_funs, length) != 0  & sapply(psal_funs, length) != 0)
  indexes_have_mean <- (1:500 + (j-1)*500)[have_results]
  indexes_have_both <- unique(df_preds$index[df_preds$index %in% indexes_have_mean])
  have_results_both <- which((1:500 + (j-1)*500) %in% indexes_have_both)
  df_indexes <- df_preds[df_preds$index %in% indexes_have_both,]
  # set up sparse matrices so we don't have to recreate them
  sigma_sparse <- bandSparse(n = n_grid* (n_grid+1), m = n_grid* (n_grid+1), k = c(0,1),
                             diagonals = list(rep(0,n_grid* (n_grid+1)), rep(0, n_grid* (n_grid+1)- 1)),
                             symmetric = T, giveCsparse = F)
  sigma_sparse_diag <-  bandSparse(n = n_grid*2, m = n_grid * 2, k = c(0,1),
                                   diagonals = list(rep(0,n_grid*2), rep(0, n_grid*2- 1)),
                                   symmetric = T, giveCsparse = F)
  selection_vector <- rep(1:(n_grid * (n_grid+1)/2), each = 2) +
    rep(c(0,n_grid * (n_grid+1)/2), times = n_grid * (n_grid+1)/2)
  selection_vector_diag <- rep(1:n_grid, each = 2) +
    rep(c(0,n_grid), times = n_grid)
  # compute for each location and year
  a <- proc.time()
  ohc <- lapply(1:nrow(df_indexes), ohc_fun)
  b <- proc.time()
  ohc_df <- cbind(df_indexes[,c('index', 'long', 'lat', 'year')], do.call(rbind, ohc))
  colnames(ohc_df) <- c('index', 'long', 'lat', 'year', 'ohc_mean', 'ohc_pred',
                        'ohc_anomaly', 'ohc_avg_anomaly', 'ohc_sd')
  final_list[[j]] <- ohc_df
  print(j)
  print(b-a)
}
ohc <- do.call(rbind, final_list)
save(ohc, file = 'analysis/results/ohc/ohc_700.RData')
q()

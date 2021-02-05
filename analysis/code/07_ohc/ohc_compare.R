setwd('/home/dyarger/argofda/')

# levels <- c(2.5,   10.0,  20.0 ,  30.0,   40.0  , 50.0   ,
#                   60.0,   70.0,   80.0,   90.0,  100.0,  110.0,  120.0,  130.0,
#                   140.0,  150.0,  160.0,  170.0,  182.5,
#                   200.0,  220.0,  240.0,  260.0,  280.0,  300.0,  320.0,  340.0,  360.0,
#                   380.0,  400.0,  420.0,  440.0,  462.5,
#                   500.0,  550.0,  600.0,  650.0,  700.0)
# levels <- seq(0, 700, by = 10)

# this function takes a location and year, returns ohc information
ohc_fun <- function(x) {
  #index_val <- which(indexes_have_mean == df_indexes[x,'index'])
  #index_val <- df_indexes[x,'index']
  index_val <- which(grid_to_compute_year == df_indexes[x,'index'])
  #index_nugget <- which(rownames(df_nugget) == df_indexes[x,'index'])
  long <- df_indexes[x,'long']
  lat <- df_indexes[x,'lat']
  scores_temp <- unlist(df_indexes[x, paste0('T', 1:K1)])
  scores_psal <- unlist(df_indexes[x, paste0('S', 1:K2)])
  pc_temp <- coefs_pc_temp[[which(grid_pc_temp$long == long &
                                    grid_pc_temp$lat == lat)]][[1]][,1:K1]
  pc_psal <-  coefs_pc_psal[[which(grid_pc_psal$long == long &
                                     grid_pc_psal$lat == lat)]][[1]][,1:K2]
  # get mean
  pred_vals_temp <- sapply(0:9, function(x) predict.local_function(temp_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec))
  avg_pred_temp <- rowMeans(pred_vals_temp)
  pred_vals_psal <- sapply(0:9, function(x) predict.local_function(psal_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec))
  avg_pred_psal <- rowMeans(pred_vals_psal)

  mean_temp_one <- pred_vals_temp[,df_indexes[x, 'year'] - 2006]
  mean_psal_one <- pred_vals_psal[,df_indexes[x, 'year'] - 2006]

  pred_vals_temp_red <- sapply(0:9, function(x) predict.local_function(temp_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec_red))
  avg_pred_temp_red <- rowMeans(pred_vals_temp_red)
  pred_vals_psal_red <- sapply(0:9, function(x) predict.local_function(psal_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec_red))
  avg_pred_psal_red <- rowMeans(pred_vals_psal_red)

  mean_temp_one <- pred_vals_temp[,df_indexes[x, 'year'] - 2006]
  mean_psal_one <- pred_vals_psal[,df_indexes[x, 'year'] - 2006]

  mean_temp_one_red <- pred_vals_temp_red[,df_indexes[x, 'year'] - 2006]
  mean_psal_one_red <- pred_vals_psal_red[,df_indexes[x, 'year'] - 2006]

  # Predictions
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

  pred_pc_temp_red <- sapply(as.data.frame(as.matrix(pc_temp)), function(z) {
    predict(make_fit_object(coefficients = z, knots = knots_pc),
            pressure_vec_red)$y
  })
  pred_pc_psal_red <- sapply(as.data.frame(as.matrix(pc_psal)), function(z) {
    predict(make_fit_object(coefficients = z, knots = knots_pc),
            pressure_vec_red)$y
  })
  pred_temp_red <- rowSums(sapply(1:K1, function(x) { scores_temp[x] * pred_pc_temp_red[,x]}))
  pred_psal_red <- rowSums(sapply(1:K2, function(x) { scores_psal[x] * pred_pc_psal_red[,x]}))

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

  # reduced versions
  absolute_salinity_pred_red <- gsw::gsw_SA_from_SP(SP = mean_psal_one_red+pred_psal_red,
                                                    p = pressure_vec_red,
                                                     longitude = long, latitude = lat)
  conservative_temp_pred_red <- gsw::gsw_CT_from_t(SA = absolute_salinity_pred_red,
                                                   t = mean_temp_one_red+ pred_temp_red,
                                                    p = pressure_vec_red)
  absolute_salinity_mean_red  <- gsw::gsw_SA_from_SP(SP = mean_psal_one_red, p = pressure_vec_red,
                                                     longitude = long, latitude = lat)
  conservative_temp_mean_red  <- gsw::gsw_CT_from_t(SA = absolute_salinity_mean_red, t = mean_temp_one_red,
                                                    p = pressure_vec_red)
  absolute_salinity_avg_mean_red <- gsw::gsw_SA_from_SP(SP = avg_pred_psal_red, p = pressure_vec_red,
                                                         longitude = long, latitude = lat)
  conservative_temp_avg_mean_red <- gsw::gsw_CT_from_t(SA = absolute_salinity_avg_mean_red,
                                                        t = avg_pred_temp_red, p = pressure_vec_red)

  conservative_temp_anomaly_red <- conservative_temp_pred_red +  273.15 - conservative_temp_mean_red -273.15
  conservative_temp_avg_anomaly_red <- conservative_temp_pred_red +  273.15 - conservative_temp_avg_mean_red -273.15



  # Kelvin for OHC predictions
  ohc_mean <- c_p/g * sum((conservative_temp_mean + 273.15) * delta)
  ohc_pred <- c_p/g * sum((conservative_temp_pred + 273.15) * delta) #mean + anomaly
  ohc_anomaly <- c_p/g * sum(conservative_temp_anomaly * delta)
  ohc_avg_anomaly <- c_p/g * sum(conservative_temp_avg_anomaly * delta)

  ohc_mean_red <- c_p/g * sum((conservative_temp_mean_red + 273.15) * pressure_gaps_red)
  ohc_pred_red <- c_p/g * sum((conservative_temp_pred_red + 273.15) * pressure_gaps_red) #mean + anomaly
  ohc_anomaly_red <- c_p/g * sum(conservative_temp_anomaly_red * pressure_gaps_red)
  ohc_avg_anomaly_red <- c_p/g * sum(conservative_temp_avg_anomaly_red * pressure_gaps_red)

  # Variance
  var_all <- matrix(nrow = K1+K2, ncol = K1+K2, unlist(df_indexes[x, paste0('condvar', rep(1:(K1+K2), each = K1+K2),
                                                                            '_', rep(1:(K1+K2), times = K1+K1))]))
  #AA_est <- matrix(nrow = K1+K2, ncol = K1+K2, unlist(df_params_only[x, paste0('var', rep(1:(K1+K2), each = K1+K2),
  #                                                                         '_', rep(1:(K1+K2), times = K1+K1))]))
  #AA_eigen <- eigen(AA_est)
  #nugget_vars <- unlist(df_params_only[x, paste0('sigma2_', 1:(K1+K2))])
  #var_all_nugget <- var_all + AA_eigen$vectors %*% diag(nugget_vars) %*% t( AA_eigen$vectors)
  #diag(var_all) <- diag(var_all) + nugget_vars
  var_all_nugget <- var_all
  # if (length(index_nugget) >0 ) {
  #   psal_nugget_use <- psal_nugget_funs_use[[index_nugget]]
  #   temp_nugget_use <- temp_nugget_funs_use[[index_nugget]]
  #   temp_nugget_eval_all <- 2*exp(predict.local_function(temp_nugget_use, pressure_vec) + digamma(1))
  #   temp_nugget_eval_red <- 2*exp(predict.local_function(temp_nugget_use, pressure_vec_red) + digamma(1))
  #   if (!is.null(psal_nugget_funs_use[[index_nugget]])) {
  #     psal_nugget_eval_all <- 2*exp(predict.local_function(psal_nugget_use, pressure_vec) + digamma(1))
  #     psal_nugget_eval_red <- 2*exp(predict.local_function(psal_nugget_use, pressure_vec_red) + digamma(1))
  #   } else {
  #     psal_nugget_eval_all <- rep(0, length(pressure_vec))
  #     psal_nugget_eval_red <- rep(0, length(pressure_vec_red))
  #   }
  # } else {
    temp_nugget_eval_all <- rep(0, length(pressure_vec))
    temp_nugget_eval_red <- rep(0, length(pressure_vec_red))
    psal_nugget_eval_all <- rep(0, length(pressure_vec))
    psal_nugget_eval_red <- rep(0, length(pressure_vec_red))
  #}

  deriv_pred <- gsw_CT_first_derivatives_wrt_t_exact(SA = absolute_salinity_pred,
                                                     t = mean_temp_one + pred_temp,
                                                     p = pressure_vec)
  SA_deriv <- deriv_pred$CT_SA_wrt_t #absolute salinity, not practical salinity
  T_deriv <- deriv_pred$CT_t_wrt_t
  derivs <- cbind(T_deriv, SA_deriv)

  deriv_pred_red <- gsw_CT_first_derivatives_wrt_t_exact(SA = absolute_salinity_pred_red,
                                                     t = mean_temp_one_red + pred_temp_red,
                                                     p = pressure_vec_red)
  SA_deriv <- deriv_pred_red$CT_SA_wrt_t #absolute salinity, not practical salinity
  T_deriv <- deriv_pred_red$CT_t_wrt_t
  derivs_red <- cbind(T_deriv, SA_deriv)

  sigma11 <- pred_pc_temp %*% var_all_nugget[1:K1,1:K1] %*% t(pred_pc_temp)
  #test <- eigen(sigma11)
  sigma12 <- pred_pc_temp %*% var_all_nugget[1:K1, (K1 + 1):(K1+K2)] %*% t(pred_pc_psal)
  sigma22 <- pred_pc_psal %*% var_all_nugget[(K1 + 1):(K1+K2), (K1 + 1):(K1+K2)] %*% t(pred_pc_psal)

  sigma11_cross <- pred_pc_temp_red %*%
    var_all_nugget[1:K1,1:K1] %*% t(pred_pc_temp)
  sigma12_cross <- pred_pc_temp_red %*% var_all_nugget[1:K1, (K1 + 1):(K1+K2)] %*% t(pred_pc_psal)
  sigma22_cross <- pred_pc_psal_red %*% var_all_nugget[(K1 + 1):(K1+K2), (K1 + 1):(K1+K2)] %*% t(pred_pc_psal)

  sigma11_red <- pred_pc_temp_red %*%
    var_all_nugget[1:K1,1:K1] %*% t(pred_pc_temp_red)
  sigma12_red <- pred_pc_temp_red %*%
    var_all_nugget[1:K1, (K1 + 1):(K1+K2)] %*% t(pred_pc_psal_red)
  sigma22_red <- pred_pc_psal_red %*%
    var_all_nugget[(K1 + 1):(K1+K2), (K1 + 1):(K1+K2)] %*% t(pred_pc_psal_red)

  # prepare to sum over the 2-d grid in pressure for T and S
  ut <- upper.tri(sigma11, diag = T)
  sigma11ut <- sigma11[ut]
  sigma12ut <- sigma12[ut]
  sigma22ut <- sigma22[ut]
  sigma_diag_vals <- c(sigma11ut, sigma22ut)[selection_vector]
  sigma_band_vals <- c(sigma12ut, rep(0, n_grid * (n_grid+1)/2))[selection_vector]
  sigma_sparse@x <- c(sigma_diag_vals, sigma_band_vals[-length(sigma_band_vals)])

  ut_red <- upper.tri(sigma11_red, diag = T)
  sigma11ut_red <- sigma11_red[ut_red]
  sigma12ut_red <- sigma12_red[ut_red]
  sigma22ut_red <- sigma22_red[ut_red]
  sigma_diag_vals_red <- c(sigma11ut_red, sigma22ut_red)[selection_vector_red]
  sigma_band_vals_red <- c(sigma12ut_red, rep(0, n_grid_red * (n_grid_red+1)/2))[selection_vector_red]
  sigma_sparse_red@x <- c(sigma_diag_vals_red, sigma_band_vals_red[-length(sigma_band_vals_red)])

  #ut_cross <- upper.tri(sigma11_cross, diag = T)
  #sigma11ut_cross <- sigma11_cross[ut_cross]
  #sigma12ut_cross <- sigma12_cross[ut_cross]
  #sigma22ut_cross <- sigma22_cross[ut_cross]
  sigma_diag_vals_cross <- c(sigma11_cross, sigma22_cross)[selection_vector_cross] # groups of 38
  sigma_band_vals_cross <- as.vector(matrix(c(sigma12_cross, rep(0, length(sigma12_cross))), nrow = 2, byrow = T))#[selection_vector_cross]
  sigma_sparse_cross@x <- c(sigma_diag_vals_cross, sigma_band_vals_cross[-length(sigma_band_vals_cross)])



  # derivatives for delta method
  derivs_diag1@x <- derivs[grid_temporary[,1],][selection_vector]
  derivs_diag2@x <- derivs[grid_temporary[,2],][selection_vector]
  #var_est <- diag(t(derivs_diag1) %*% sigma_sparse %*% derivs_diag2) # sum over 2-d gd
  #var_est2 <- colSums(derivs_diag1* (sigma_sparse %*% derivs_diag2)) # sum over 2-d gd
  sigma_sparse2 <- sigma_sparse
  index_remove <- sigma_sparse2@x !=0
  sigma_sparse2@i <- sigma_sparse2@i[index_remove]
  sigma_sparse2@j <- sigma_sparse2@j[index_remove]
  sigma_sparse2@x <- sigma_sparse2@x[index_remove]
  var_est <- diag(t(derivs_diag1) %*% sigma_sparse2 %*% derivs_diag2)
  # rbenchmark::benchmark(one = {diag(t(derivs_diag1) %*% sigma_sparse %*% derivs_diag2)},
  #                       two = {diag(t(derivs_diag1) %*% sigma_sparse2 %*% derivs_diag2)},
  #                       #two = {colSums(derivs_diag1* (sigma_sparse %*% derivs_diag2)) },
  #                       #three = {colSums((t(derivs_diag1) %*% sigma_sparse) * t(derivs_diag2))},
  #                       replications = 10)
  mat_final <- matrix(nrow = length(pressure_vec), ncol = length(pressure_vec), 0)
  mat_final[upper.tri(mat_final, diag = T)] <- var_est
  mat_final <- forceSymmetric(mat_final, 'U')
  diag(mat_final) <- diag(mat_final) +
    temp_nugget_eval_all * (derivs[,1])^2 +
    psal_nugget_eval_all * (derivs[,2])^2

  derivs_diag1_red@x <- derivs_red[grid_temporary_red[,1],][selection_vector_red]
  derivs_diag2_red@x <- derivs_red[grid_temporary_red[,2],][selection_vector_red]
  var_est_red <- diag(t(derivs_diag1_red) %*% sigma_sparse_red %*% derivs_diag2_red) # sum over 2-d grid
  mat_final_red <- matrix(nrow = length(pressure_vec_red), ncol = length(pressure_vec_red), 0)
  mat_final_red[upper.tri(mat_final_red, diag = T)] <- var_est_red
  mat_final_red <- forceSymmetric(mat_final_red, 'U')
  diag(mat_final_red) <- diag(mat_final_red) +
    temp_nugget_eval_red * (derivs_red[,1])^2 +
    psal_nugget_eval_red * (derivs_red[,2])^2
  # %*% diag(temp_nugget_eval_red) %*% t(derivs_red[,1])
  #mat_final_red_eigen <- eigen(mat_final_red)

  derivs_diag1_cross@x <- derivs_red[grid_temporary_cross[,1],][selection_vector_cross] # reduced goes on left
  derivs_diag2_cross@x <- derivs[grid_temporary_cross[,2],][selection_vector_cross]
  var_est_cross <- diag(t(derivs_diag1_cross) %*% sigma_sparse_cross %*% derivs_diag2_cross) # sum over 2-d grid

  mat_final_cross <- matrix(var_est_cross,byrow =T,
                            ncol = length(pressure_vec_red), nrow = length(pressure_vec))
  # diag(mat_final) <- diag(mat_final) + .1
  # diag(mat_final_red) <- diag(mat_final_red) + .1
  var_ohc <- c_p^2/g^2 *t(rep(delta, nrow(mat_final_cross))) %*%
    mat_final %*%  rep(delta, nrow(mat_final_cross))
  #var_ohc_red <- c_p^2/g^2 * sum(mat_final_red * pressure_gaps_red * pressure_gaps_red)
  var_ohc_red <- (c_p^2/g^2 * t(pressure_gaps_red) %*%
    mat_final_red %*% pressure_gaps_red)[1]
  var_ohc_cross <- c_p^2/g^2 * t(rep(delta, nrow(mat_final_cross))) %*%
                                   mat_final_cross %*% pressure_gaps_red

  # mat_all <- cbind(rbind(mat_final, t(mat_final_cross)),
  #                  rbind(mat_final_cross, mat_final_red))
  # mat_all <- cbind(rbind(mat_final_red, mat_final_cross),
  #                  rbind(t(mat_final_cross), mat_final))
  # isSymmetric(mat_all)
  # corpcor::is.positive.definite(mat_all)
  # corpcor::is.positive.definite(mat_final_red)
  # corpcor::is.positive.definite(mat_final)
  # test <- eigen(mat_all)
  # var_ohc  + var_ohc_red - 2 * var_ohc_cross
  # sqrt(var_ohc  + var_ohc_red - 2 * var_ohc_cross)
  # sqrt(var_ohc  + var_ohc_red - 2 * var_ohc_cross)/sqrt(var_ohc  )
  #sqrt(var_ohc  + var_ohc_red - 2 * var_ohc_cross)/sqrt(var_ohc_red  )

  # combine the two and compute
  c(ohc_mean,ohc_pred, ohc_anomaly, ohc_avg_anomaly,
    ohc_mean_red,ohc_pred_red, ohc_anomaly_red, ohc_avg_anomaly_red,
    var_ohc[1], var_ohc_red, var_ohc_cross)
}

library(ggplot2)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))
library(fda)
library(argofda)
library(gsw)
library(parallel)

load('analysis/results/joint_TS_20/preds.RData')
load('analysis/results/joint_TS_20/params.RData')

load('analysis/results/temp_cov_pca.RData')
load('analysis/results/psal_cov_pca.RData')
basis <- create.bspline.basis(breaks =knots_pc)

source('analysis/code/05_nugget_variance/load_nugget.R')

### get mean functions
calc_files_temp <- list.files('analysis/results/mean_estimation/', pattern = 'temp')
calc_files_psal <- list.files('analysis/results/mean_estimation/', pattern = 'psal')
final_list <- list()

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$index <- 1:nrow(grid_to_compute)

delta <- .5 # approximation grid in pressure
pressure_vec <- seq(0.25, 700.0, by = delta)
n_grid <- length(pressure_vec)
pressure_grid <- expand.grid('p1' = pressure_vec, 'p2' = pressure_vec)
grid_temporary <- expand.grid(1:length(pressure_vec), 1:length(pressure_vec))
grid_temporary <- as.matrix(grid_temporary[pressure_grid$p1 <= pressure_grid$p2,])
pressure_grid <- as.matrix(pressure_grid[pressure_grid$p1 <= pressure_grid$p2,])
diag_indices <- grid_temporary[grid_temporary[,1] == grid_temporary[,2],1]

pressure_vec_red <- c(2.5,   10.0,  20.0 ,  30.0,   40.0  , 50.0   ,
            60.0,   70.0,   80.0,   90.0,  100.0,  110.0,  120.0,  130.0,
            140.0,  150.0,  160.0,  170.0,  182.5,
            200.0,  220.0,  240.0,  260.0,  280.0,  300.0,  320.0,  340.0,  360.0,
            380.0,  400.0,  420.0,  440.0,  462.5,
            500.0,  550.0,  600.0,  650.0,  700.0)
delta2 <- 10# approximation grid in pressure
pressure_vec_red <- seq(5, 700.0, by = delta2)
# pressure_gaps_red <- pressure_vec_red[2:length(pressure_vec_red)] -
#   pressure_vec_red[1:(length(pressure_vec_red) - 1)]
# pressure_gaps_red_under <- c(20, pressure_gaps_red)
# pressure_gaps_red_over <- c(pressure_gaps_red, 20)
pressure_gaps_red <- rep(delta2, length(pressure_vec_red))


n_grid_red <- length(pressure_gaps_red)
pressure_grid_red <- expand.grid('p1' = pressure_vec_red, 'p2' = pressure_vec_red)
grid_temporary_red <- expand.grid(1:length(pressure_vec_red), 1:length(pressure_vec_red))
grid_temporary_red <- as.matrix(grid_temporary_red[pressure_grid_red$p1 <= pressure_grid_red$p2,])
pressure_grid_red <- as.matrix(pressure_grid_red[pressure_grid_red$p1 <= pressure_grid_red$p2,])
diag_indices_red <- grid_temporary_red[grid_temporary_red[,1] == grid_temporary_red[,2],1]
list_of_vals_red <- split(rep(1, length(pressure_grid_red)), rep(1:(length(pressure_grid_red)/2), each = 2))
derivs_diag1_red <- do.call(bdiag, list_of_vals_red)
derivs_diag2_red <- derivs_diag1_red

# cross
pressure_grid_cross <- expand.grid('p1' = pressure_vec, 'p2' = pressure_vec_red)
#grid_temporary_cross <- expand.grid(1:length(pressure_vec), 1:length(pressure_vec_red))
grid_temporary_cross <- expand.grid(1:length(pressure_vec_red), 1:length(pressure_vec))
#grid_temporary_cross <- as.matrix(grid_temporary_cross[pressure_grid_cross$p1 <= pressure_grid_cross$p2,])
#pressure_grid_cross <- as.matrix(pressure_grid_cross[pressure_grid_cross$p1 <= pressure_grid_cross$p2,])
#diag_indices_cross <- grid_temporary_cross[grid_temporary_cross[,1] == grid_temporary_cross[,2],1]
list_of_vals_cross <- split(rep(1, nrow(pressure_grid_cross)*2), rep(1:(nrow(pressure_grid_cross)), each = 2))
derivs_diag1_cross <- do.call(bdiag, list_of_vals_cross)
derivs_diag2_cross <- derivs_diag1_cross


list_of_vals <- split(rep(1, length(pressure_grid)), rep(1:(length(pressure_grid)/2), each = 2))
derivs_diag1 <- do.call(bdiag, list_of_vals)
derivs_diag2 <- derivs_diag1
derivs_diag3 <- derivs_diag1

list_of_vals_red <- split(rep(1, length(pressure_grid_red)),
                          rep(1:(length(pressure_grid_red)/2), each = 2))
derivs_diag1_red <- do.call(bdiag, list_of_vals_red)
derivs_diag2_red <- derivs_diag1_red
derivs_diag3_red <- derivs_diag1_red


c_p <- 3850 # Joules per kg per degrees celcius  - specific heat of seawater
g <- (1025)^(-1) # 1/density of seawater in (kg/m^3)^-1

#for (j in 1:length(calc_files_temp)) {
j <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  # get mean functions
  load(paste0('analysis/results/mean_estimation/', calc_files_temp[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  temp_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][1] == 'no profiles in range') {return(NULL)} else {return(x[[1]])}
  })
  eval(parse(text=paste0('rm(', file, ')')))
  load(paste0('analysis/results/mean_estimation/', calc_files_psal[j]))
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
  grid_to_compute_year <- df_indexes[!duplicated(df_indexes$index), 'index']
  temp_nugget_funs_use <- temp_nugget_funs[rownames(grid_nugget) %in% indexes_have_both]
  psal_nugget_funs_use <- psal_nugget_funs[rownames(grid_nugget) %in% indexes_have_both]
  df_nugget <- grid_nugget[rownames(grid_nugget) %in% indexes_have_both,]

  # set up sparse matrices so we don't have to recreate them
  sigma_sparse <- bandSparse(n = n_grid* (n_grid+1), m = n_grid* (n_grid+1), k = c(0,1),
                             diagonals = list(rep(0,n_grid* (n_grid+1)), rep(0, n_grid* (n_grid+1)- 1)),
                             symmetric = T, giveCsparse = F)
  sigma_sparse_diag <- bandSparse(n = n_grid*2, m = n_grid * 2, k = c(0,1),
                                  diagonals = list(rep(0,n_grid*2), rep(0, n_grid*2- 1)),
                                  symmetric = T, giveCsparse = F)
  selection_vector <- rep(1:(n_grid * (n_grid+1)/2), each = 2) +
    rep(c(0,n_grid * (n_grid+1)/2), times = n_grid * (n_grid+1)/2)
  selection_vector_diag <- rep(1:n_grid, each = 2) +
    rep(c(0,n_grid), times = n_grid)

  # reduced
  sigma_sparse_red <- bandSparse(n = n_grid_red* (n_grid_red+1), m = n_grid_red* (n_grid_red+1), k = c(0,1),
                             diagonals = list(rep(0,n_grid_red* (n_grid_red+1)), rep(0, n_grid_red* (n_grid_red+1)- 1)),
                             symmetric = T, giveCsparse = F)
  sigma_sparse_diag_red <- bandSparse(n = n_grid_red*2, m = n_grid_red * 2, k = c(0,1),
                                  diagonals = list(rep(0,n_grid_red*2), rep(0, n_grid_red*2- 1)),
                                  symmetric = T, giveCsparse = F)
  selection_vector_red <- rep(1:(n_grid_red * (n_grid_red+1)/2), each = 2) +
    rep(c(0,n_grid_red * (n_grid_red+1)/2), times = n_grid_red * (n_grid_red+1)/2)
  selection_vector_diag_red <- rep(1:n_grid_red, each = 2) +
    rep(c(0,n_grid_red), times = n_grid_red)

  sigma_sparse_cross <- bandSparse(n = n_grid * n_grid_red*2, m =n_grid * n_grid_red*2, k = c(0,1),
                                 diagonals = list(rep(0,n_grid * n_grid_red*2), rep(0, n_grid * n_grid_red*2- 1)),
                                 symmetric = T, giveCsparse = F)
  selection_vector_cross <- rep(1:(n_grid * (n_grid_red+1)/2), each = 2) +
    rep(c(0,n_grid * (n_grid_red+1)/2), times = n_grid * (n_grid_red+1)/2)
  selection_vector_cross <- rep(1:(n_grid * (n_grid_red)), each = 2) +
    rep(c(0,n_grid * (n_grid_red)), times = n_grid * (n_grid_red))

  K1 <- 10;K2 <- 10
  a <- proc.time()
  ohc <- mclapply(1:nrow(df_indexes), ohc_fun, mc.cores = 3, mc.preschedule = F)
  b <- proc.time()
  ohc_df <- cbind(df_indexes[,c('index', 'long', 'lat', 'year')], do.call(rbind, ohc))
  #colnames(ohc_df) <- c('index', 'long', 'lat', 'year', 'ohc_mean', 'ohc_pred',
                        # 'ohc_anomaly', 'ohc_avg_anomaly', 'ohc_sd')
  final_list[[j]] <- ohc_df
  print(j)
  print(b-a)
#}
ohc <- do.call(rbind, final_list)
save(ohc, file = paste0( 'analysis/results/ohc/ohc_700_cross_no_nugget', j, '.RData'))
q()

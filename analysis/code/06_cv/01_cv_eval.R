
# cross validation

# Expectation step for one year
e_step_cv <- function(params_use, delayed, not_delayed, df_use,
                      dist_long_del, dist_lat_del, dist_tim_del,
                      pred_dist_long, pred_dist_lat, pred_dist_tim,
                      eig_vecs_AA,
                      nu, K1, K2) {

  pred <- list()
  for (k in 1:(K1+K2)) {
    theta <- params_use[k,]
    scale_long <- theta[1]
    scale_lat <- theta[2]
    scale_tim <- theta[3]
    sigma2 <- theta[4]
    var_score <- theta[5]

    pred_C <- Matern(sqrt( (pred_dist_long/scale_long)^2+
                             (pred_dist_lat/scale_lat)^2+
                             (pred_dist_tim/scale_tim)^2),
                     nu = nu)
    C_full <- as.matrix(Exponential(sqrt( (dist_long_del/scale_long)^2  + (dist_lat_del/scale_lat)^2 +
                                            (dist_tim_del/scale_tim)^2)))
    response_delayed <- df_use[delayed, paste0('N', k, '_new') ]
    pred[[k]] <-  t(var_score * pred_C) %*% solve(var_score*C_full +
                                              diag(sigma2, nrow(C_full)), response_delayed)
  }
  df_scores_new <- df_use
  pred_year <- do.call(cbind, pred)

  pred_old_scores <- t(eig_vecs_AA %*% t(pred_year))

  df_scores_new[not_delayed, paste0('S', 1:K2, '_old')] <-
    as.matrix(pred_old_scores[,(K1+1):(K1+K2)])

  old_scores <- df_scores_new[,c(paste0('T', 1:K1, '_old'),paste0('S', 1:K2, '_old'))]
  new_scores <- t(t(eig_vecs_AA) %*% t(old_scores))
  colnames(new_scores) <- paste0('N', 1:(K1 + K2), '_new')

  df_scores_new[paste0('N', 1:(K1 + K2), '_new')] <- new_scores
  return(df_scores_new)
}
setwd('/home/dyarger/argofda/')
load('analysis/results/joint_TS_20/params.RData')
source('analysis/code/05_nugget_variance/load_nugget.R')

library(argofda)
library(Matrix)
library(fda)
library(dplyr)
library(fields)
library(CompQuadForm)
library(parallel)

load('analysis/results/psal_cov_pca.RData')
load('analysis/results/temp_cov_pca.RData')

load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)


# for each profile in February, leave it out and predict
cv <- function(profile_num, h_space, K1, K2) {
  a <- proc.time()

  lat <- profLatAggr[profile_num]
  long <- profLongAggr[profile_num]
  day <- profJulFracofYearAggr[profile_num]
  year <- profYearAggr[profile_num]
  lat_eval <- latvals[which.min(abs(lat - latvals))]
  long_eval <- longvals[which.min(abs(long - longvals))]

  # Get nearby data
  df <- get_profile_data_subset(lat = lat, long = long, day = 45.25, h_space = h_space,
                                  RG_def = T,h_time = 45.25,mode = 'all',min_prof = 1,
                                  exclusion = T,years = year)

  # leave profile out in df_out
  df_out <- df[df$profile_num == profile_num,]
  df <- df[df$profile_num != profile_num,]
  mode_prof <- as.character(df_out$mode[1])
  df$profile <- as.numeric(as.factor(df$profile))

  # If no pca was computed at the location, skip
  if (!(lat_eval %in% grid_pc_temp$lat[grid_pc_temp$long == long_eval] &
        lat_eval %in% grid_pc_psal$lat[grid_pc_psal$long == long_eval]  ) |
      is.null(df) | is.null(nrow(df))) {
    b <- proc.time()
    final_list <- c(profile_num, 'no nearby data', (b-a)[3])
    residual_info <- df_out[, c('pressure', 'temperature', 'temperature','temperature', 'temperature',
                                'temperature', 'salinity', 'salinity', 'salinity', 'salinity', 'salinity')]
    colnames(residual_info) <- c('pressure', 'temp_mean_residual',
                                 'temp_residual', 'var_values_temp', 'nugget_var_temp', 'var_values_temp_uc',
                                 'psal_mean_residual',
                                 'psal_residual', 'var_values_psal', 'nugget_var_psal', 'var_values_psal_uc')
    final_list <- list(
      NULL,
      NULL,
      residual_info,
      NULL,
      mode_prof,
      NULL,
      c(profile_num, 'no nearby data', max(df$profile)))
    print(final_list[[6]])
    return(final_list)
  }

  # remove profiles that we can't get an accurate estimate of scores
  # these are either too short or have too few measurements
  pressure_min <- df %>% group_by(profile) %>% summarise(min(pressure), .groups = "drop")
  pressure_max <- df %>% group_by(profile) %>% summarise(max(pressure), .groups = "drop")
  short_profiles <- (1:nrow(pressure_max))[which(pressure_max[,2] - pressure_min[,2] < 800)]
  prof_lengths <- df %>% group_by(profile) %>% summarise(n(), .groups = "drop")
  small_profiles <-(1:nrow(prof_lengths))[which(prof_lengths[,2] < 15)]
  df <- df[!(df$profile %in% short_profiles) & !(df$profile %in% small_profiles),]
  df$profile <- as.numeric(as.factor(df$profile))

  if (is.null(df) | is.null(nrow(df)) ) {
    b <- proc.time()
    residual_info <- df_out[, c('pressure', 'temperature', 'temperature','temperature', 'temperature',
                                'temperature', 'salinity', 'salinity', 'salinity', 'salinity', 'salinity')]
    colnames(residual_info) <- c('pressure', 'temp_mean_residual',
                                 'temp_residual', 'var_values_temp', 'nugget_var_temp', 'var_values_temp_uc',
                                 'psal_mean_residual',
                                 'psal_residual', 'var_values_psal', 'nugget_var_psal', 'var_values_psal_uc')
    final_list <- list(
      NULL,
      NULL,
      residual_info,
      NULL,
      mode_prof,
      NULL,
      c(profile_num, 'no nearby data', max(df$profile)))
    print(final_list[[6]])
    return(final_list)
  }

  # Calculate scores
  temp_pc <- coefs_pc_temp[[which(grid_pc_temp$long == long_eval &
                                    grid_pc_temp$lat == lat_eval)]][[1]][,1:K1]
  psal_pc <- coefs_pc_psal[[which(grid_pc_psal$long == long_eval &
                                    grid_pc_psal$lat == lat_eval)]][[1]][,1:K2]

  # estimate scores
  scores_out_temp <- tryCatch(get_scores(df_out, temp_pc, K1, knots = knots_pc)[[1]],
                              error = function(q) {rep(NA, K1)})
  scores_out_psal <- tryCatch(get_scores(df_out, psal_pc, K2, 'salinity',knots = knots_pc)[[1]],
                              error = function(q) {rep(NA, K2)})
  scores_temp <- get_scores(df, temp_pc, K1, knots = knots_pc)[[1]]
  scores_psal <- get_scores(df, psal_pc, K2, 'salinity', knots = knots_pc)[[1]]
  small_profiles <- which(is.na(scores_temp[,1]))
  if (length(small_profiles) == 0) {small_profiles =0}
  scores_temp <- scores_temp[!is.na(scores_temp[,1]),]
  scores_psal <- scores_psal[!is.na(scores_psal[,1]),]

  # get longitudes and latitudes and break into year
  max_prof <- max(df$profile)
  df_scores <- df %>%
    filter(profile != small_profiles) %>%
    group_by(profile) %>%
    filter(pressure== pressure[1])

  # Get AA^\top estimate the results
  if (nrow(df_params_only[df_params_only$long == long_eval & df_params_only$lat == lat_eval,]) == 0) {
    b <- proc.time()
    residual_info <- df_out[, c('pressure', 'temperature', 'temperature','temperature', 'temperature',
                                'temperature', 'salinity', 'salinity', 'salinity', 'salinity', 'salinity')]
    colnames(residual_info) <- c('pressure', 'temp_mean_residual',
                                 'temp_residual', 'var_values_temp', 'nugget_var_temp', 'var_values_temp_uc',
                                 'psal_mean_residual',
                                 'psal_residual', 'var_values_psal', 'nugget_var_psal', 'var_values_psal_uc')
    final_list <- list(
      NULL,
      NULL,
      residual_info,
      NULL,
      mode_prof,
      NULL,
      c(profile_num, 'no nearby data', max(df$profile)))
    print(final_list[[6]])
    return(final_list)
  }

  # nugget variances
  if (length(which(grid_nugget$long == long_eval &
                   grid_nugget$lat == lat_eval)) == 0) {
    var_pred_temp_out <- rep(0, nrow(df_out))
    var_pred_psal_out <- rep(0, nrow(df_out))
  } else {
    nugget_temp_out <- temp_nugget_funs[[which(grid_nugget$long == long_eval &
                                                 grid_nugget$lat == lat_eval)]]
    nugget_psal_out <- psal_nugget_funs[[which(grid_nugget$long == long_eval &
                                                 grid_nugget$lat == lat_eval)]]
    if (is.null(nugget_temp_out)) {
      var_pred_temp_out <- rep(0, nrow(df_out))
      var_pred_psal_out <- rep(0, nrow(df_out))
    } else if (is.null(nugget_psal_out)) {
      var_pred_psal_out <- rep(0, nrow(df_out))
      var_pred_temp_out <- exp(predict(nugget_temp_out, df_out$pressure))
    } else {
      var_pred_temp_out <- exp(predict(nugget_temp_out, df_out$pressure))
      var_pred_psal_out <- exp(predict(nugget_psal_out, df_out$pressure))
    }
  }


  # AA^T estimate and form decorrelated scores
  AA_est <- matrix(nrow = K1 + K2, unlist(df_params_only[df_params_only$long == long_eval & df_params_only$lat == lat_eval,
                      paste0('var', rep(1:(K1+K2), each = K1+K2), '_', rep(1:(K1+K2), times = K1+K2))]))
  eig_AA <- eigen(AA_est)
  eig_vecs_AA <- eig_AA$vectors
  colnames(scores_temp) <- paste0('T', 1:K1, '_old')# original
  colnames(scores_psal) <- paste0('S', 1:K2, '_old')# original
  new_scores <- t(t(eig_AA$vectors) %*% t(cbind(scores_temp, scores_psal)))
  colnames(new_scores) <- paste0('N', 1:(K1+ K2), '_new') # decorrelated

  df_scores <- data.frame(df_scores[,c('longitude', 'latitude',
                                       'fracofyear', 'mode')], scores_temp, scores_psal,new_scores)
  df_scores <- df_scores[!(duplicated(df_scores$latitude) &  duplicated(df_scores$longitude)),]

  # form matrices
  delayed <- which(df_scores$mode == 'D')
  not_delayed <- which(df_scores$mode != 'D')

  dist_long <- rdist.earth(x1 = cbind(df_scores[c('longitude')], lat),
                          cbind(df_scores[c('longitude')], lat),miles = F)
  dist_lat <- rdist.earth(x1 = cbind(long, df_scores[c('latitude')]),
                          cbind(long, df_scores[c('latitude')]), miles = F)
  dist_tim <- rdist(df_scores[,'fracofyear'],df_scores[,'fracofyear'])

  pred_dist_long_all <- rdist.earth(x1 = cbind(df_scores[c('longitude')], lat), cbind(long, lat),miles = F)
  pred_dist_lat_all <- rdist.earth(x1 = cbind(long, df_scores[c('latitude')]),
                                   cbind(long, lat), miles = F)
  pred_dist_tim_all <- abs(df_scores[,'fracofyear'] - day)

  len_0_delayed <- sum(delayed) == 0 | sum(not_delayed) == 0

  # get parameters
  params <- matrix(nrow = K1 +K2, ncol = 5)
  for (k in 1:(K1+K2)) {
    params[k,] <-  as.double(df_params_only[df_params_only$long == long_eval & df_params_only$lat == lat_eval,
                                       c(paste0('scale_long_', k) ,
                                         paste0('scale_lat_', k),
                                         paste0('scale_tim_', k),
                                         paste0('sigma2_', k),
                                         paste0('var_score_', k))])
  }
  nu <- .5

  # go through e step with estimated parameters
  if (!len_0_delayed) {
    pred_dist_long <- rdist.earth(x1 = cbind(df_scores[delayed,c('longitude')], lat),
                                cbind(df_scores[not_delayed,c('longitude')], lat), miles = F)
    pred_dist_lat <- rdist.earth(x1 = cbind(long, df_scores[delayed,c('latitude')]),
                                 cbind(long, df_scores[not_delayed,c('latitude')]), miles = F)
    pred_dist_tim <- rdist(df_scores[delayed,'fracofyear'] ,
                           df_scores[not_delayed,'fracofyear'])
    df_new <- e_step_cv(params_use = params,
                     delayed = delayed,
                     not_delayed = not_delayed,df_use = df_scores,
                     dist_long_del = dist_long[delayed, delayed],
                     dist_lat_del = dist_lat[delayed, delayed],
                     dist_tim_del = dist_tim[delayed, delayed],
                     pred_dist_long = pred_dist_long,
                     pred_dist_lat = pred_dist_lat,
                     pred_dist_tim = pred_dist_tim,
                     eig_vecs_AA,nu =nu, K1 = K1, K2 = K2)
  } else {
    df_new <- df_scores
  }
  # predict
  uncond_var <- list()
  cond_var <- list()
  pred <- list()
  for (k in 1:(K1 + K2)) {
    response  <- df_new[,paste0('N', k, '_new')]
    scale_long <- params[k,1]
    scale_lat <- params[k,2]
    scale_tim <- params[k,3]
    sigma2    <- params[k,4]
    var_score <- params[k,5]

    # Predict

    pred_C <- Matern(sqrt( (pred_dist_long_all/scale_long)^2+
                             (pred_dist_lat_all/scale_lat)^2+
                             (pred_dist_tim_all/scale_tim)^2),
                     nu = nu)
    C_full <-  Matern(sqrt( (dist_long/scale_long)^2  + (dist_lat/scale_lat)^2 +
                                      (dist_tim/scale_tim)^2) ,
                              nu = nu)
    pred[[k]] <- (t(var_score * pred_C) %*% solve(var_score*C_full +
                                                              diag(sigma2, nrow(C_full)),
                                                       response))[1]
    uncond_var[[k]] <- var_score + sigma2
    cond_var[[k]] <- (var_score+ sigma2 - t(var_score*pred_C) %*%
                        solve(var_score*C_full + diag(sigma2, nrow(C_full)),
                              var_score*pred_C))[1]
  }

  # compute predictions and variances at each pressure of df_out
  pred_old_scores <- as.double(t(eig_vecs_AA %*% unlist(pred) ))
  pred_old_scores_temp <- pred_old_scores[1:K1]
  pred_old_scores_psal <- pred_old_scores[(K1+1):(K1+K2)]
  uncond_var_old_scores <- eig_vecs_AA %*% Diagonal(n = K1 + K2, unlist(uncond_var)) %*% t(eig_vecs_AA)
  uncond_var_old_scores_temp <- uncond_var_old_scores[1:K1, 1:K1]
  uncond_var_old_scores_psal <- uncond_var_old_scores[(K1+1):(K1+K2), (K1+1):(K1+K2)]
  cond_var_old_scores <- eig_vecs_AA %*% Diagonal(n = K1 + K2, unlist(cond_var)) %*% t(eig_vecs_AA)
  cond_var_old_scores_temp <- cond_var_old_scores[1:K1, 1:K1]
  cond_var_old_scores_psal <- cond_var_old_scores[(K1+1):(K1+K2), (K1+1):(K1+K2)]
  pc_values_temp <- sapply(1:K1, function(x) {
    predict(make_fit_object(knots = knots_pc, coefficients = temp_pc[,x] ), df_out$pressure)$y})
  pc_values_psal <- sapply(1:K2, function(x) {
    predict(make_fit_object(knots = knots_pc, coefficients = psal_pc[,x] ), df_out$pressure)$y})
  pc_values <- cbind(pc_values_temp, pc_values_psal)
  pred_values_temp <- rowSums(sapply(1:K1, function(x) { pred_old_scores_temp[x] * pc_values_temp[,x]}))
  pred_values_psal <- rowSums(sapply(1:K2, function(x) { pred_old_scores_psal[x] * pc_values_psal[,x]}))

  # Pointwise intervals
  var_values_uc_temp <- diag(pc_values_temp %*% uncond_var_old_scores_temp %*% t(pc_values_temp))
  var_values_uc_psal <- diag(pc_values_psal %*% uncond_var_old_scores_psal %*% t(pc_values_psal))
  var_values_temp <- diag(pc_values_temp %*% cond_var_old_scores_temp %*% t(pc_values_temp))
  var_values_psal <- diag(pc_values_psal %*% cond_var_old_scores_psal %*% t(pc_values_psal))

  lower_bound_temp <- pred_values_temp - 2*sqrt(var_values_temp+var_pred_temp_out)
  upper_bound_temp <- pred_values_temp + 2*sqrt(var_values_temp+var_pred_temp_out)

  lower_bound_psal <- pred_values_psal - 2*sqrt(var_values_psal + var_pred_psal_out)
  upper_bound_psal <- pred_values_psal + 2*sqrt(var_values_psal+ var_pred_psal_out)

  in_interval_temp <- ifelse(df_out$temperature <= upper_bound_temp & df_out$temperature >= lower_bound_temp,
                             T, F)
  in_interval_psal <- ifelse(df_out$salinity <= upper_bound_psal & df_out$salinity >= lower_bound_psal,
                             T, F)

  # Confidence bands - based on Choi and Reimherr
  e_vals <- diag(cond_var_old_scores_temp)
  c_vec <- sqrt(sqrt(e_vals))
  cond_var_norm <- cond_var_old_scores_temp * 1/sqrt(e_vals)
  cond_var_norm <- t(cond_var_norm) * 1/sqrt(e_vals)
  lambda_in <- eigen(diag(e_vals/c_vec^2) %*% cond_var_norm)$values
  quantile_all <- pnorm(q = 2) - pnorm(q = -2)
  quantile <- 1 - (1-quantile_all)/2
  max_val <- 1000
  # find quantile of W_theta
  quantile_approx <- optimise(f = function(x) {
    y <- imhof(q = x,
               lambda = lambda_in)$Qq
    abs(y - (1-quantile))
  }, lower = 0, upper = max_val)$minimum
  while ((imhof(quantile_approx, lambda_in)$Qq < .02 | imhof(quantile_approx, lambda_in)$Qq > .05)  & max_val > .2) {
    max_val <- max_val/5
    quantile_approx <- optimise(f = function(x) {
      y <- imhof(q = x,
                 lambda = lambda_in)$Qq
      abs(y - (1-quantile))
    }, lower = 0, upper = max_val)$minimum
  }
  quantile_noise <- qnorm(p = (1-quantile)/nrow(df_out)/2 )
  r_x = sqrt( quantile_approx* rowSums(sapply(1:K1, function(x)
    { c_vec[x]^2 * pc_values_temp[,x]^2})) +
      abs(quantile_noise)^2 * var_pred_temp_out)

  # confidence band for salinity
  e_vals_psal <- diag(cond_var_old_scores_psal)
  c_vec_psal <- sqrt(sqrt(e_vals_psal))
  cond_var_norm_psal <- cond_var_old_scores_psal * 1/sqrt(e_vals_psal)
  cond_var_norm_psal <- t(cond_var_norm_psal) * 1/sqrt(e_vals_psal)
  lambda_in_psal <- eigen(diag(e_vals_psal/c_vec_psal^2) %*% cond_var_norm_psal)$values
  max_val <- 1000
  # find quantile of W_theta
  quantile_approx_psal <- optimise(f = function(x) {
    max_val <- max_val/5
    y <- imhof(q = x,
               lambda = lambda_in_psal)$Qq
    abs(y - (1-quantile))
  }, lower = 0, upper = max_val)$minimum
  while ((imhof(quantile_approx_psal, lambda_in_psal)$Qq < .02 |
          imhof(quantile_approx_psal, lambda_in_psal)$Qq > .05) & max_val > .2) {
    quantile_approx_psal <- optimise(f = function(x) {
      y <- imhof(q = x,
                 lambda = lambda_in_psal)$Qq
      abs(y - (1-quantile))
    }, lower = 0, upper = max_val)$minimum
    max_val <- max_val/5
  }
  r_x_psal = sqrt(quantile_approx_psal* rowSums(sapply(1:K2, function(x) {
    c_vec_psal[x]^2 * pc_values_psal[,x]^2})) +
      abs(quantile_noise)^2 *var_pred_psal_out)

  # check bands and make residuals
  band_in_temp <- ifelse(df_out$temperature <= pred_values_temp+ r_x &
                                  df_out$temperature >= pred_values_temp- r_x,
                                T, F)
  band_in_psal <- ifelse(df_out$salinity <= pred_values_psal+ r_x_psal &
                           df_out$salinity >= pred_values_psal- r_x_psal,
                         T, F)

  temp_residual <- df_out$temperature - pred_values_temp
  psal_residual <- df_out$salinity - pred_values_psal

  residual_info <- cbind(df_out[, c('pressure', 'temperature')],
                         temp_residual, var_values_temp, var_pred_temp_out,
                         var_values_uc_temp,
                         df_out[, c('salinity')], psal_residual,
                         var_values_psal, var_pred_psal_out,var_values_uc_psal)
  colnames(residual_info) <- c('pressure', 'temp_mean_residual',
                               'temp_residual', 'var_values_temp', 'nugget_var_temp', 'var_values_temp_uc',
                               'psal_mean_residual',
                               'psal_residual', 'var_values_psal', 'nugget_var_psal', 'var_values_psal_uc')

  b <- proc.time()
  final_list <- list(
    list(pred_old_scores_temp, scores_out_temp, quantile_approx,
         c_vec, band_in_temp, in_interval_temp),
    list(pred_old_scores_psal, scores_out_psal, quantile_approx_psal,
         c_vec_psal,band_in_psal, in_interval_psal),
    residual_info,
    cond_var_old_scores,
    mode_prof,
    length(delayed),
    c(profile_num, (b-a)[3], max(df$profile)))
  print(final_list[[7]])
  return(final_list)
}

indexes <- (1:length(profMonthAggr))[profMonthAggr == 2] # February Profiles
indexes_list <- split(indexes, rep(1:40, length(profMonthAggr)))
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

a <- proc.time()
cv_results <- lapply(indexes_list[[array_id]],
                       cv, h_space = 1100,# mc.cores = 1,mc.preschedule = FALSE,
                       K1 = 10, K2 = 10)
b <- proc.time()
b-a

save(cv_results, file = paste0('analysis/results/cv/cv_results_', array_id,
                               '.RData'))

library(argofda)
library(Matrix)
library(fda)
library(sparseinv)
library(dplyr)
library(fields)
library(parallel)
setwd('/pylon5/msz3a1p/d7yar/argofda/')
#source('analysis/code/04_joint_TS_20/ml_source_vec.R')

# vectorized versions of maximum likelihood

collect_matrices_vec <- function(df_scores, long, lat) {
  delayed <- lapply(df_scores, function(x) {
    which(x$mode == 'D')
  })

  not_delayed <- lapply(df_scores, function(x) {
    which(x$mode != 'D')
  })

  dist_long <- lapply(df_scores, function(y) {
    rdist.earth(x1 = cbind(y[c('longitude')], lat), cbind(y[c('longitude')], lat),
                miles = F)})
  dist_lat <- lapply(df_scores, function(y)
    rdist.earth(x1 = cbind(long, y[c('latitude')]), cbind(long, y[c('latitude')]),
                miles = F))
  dist_tim <- lapply(df_scores, function(y) rdist(y[,'fracofyear'],y[,'fracofyear']))

  pred_dist_long <- lapply(df_scores, function(y) {
    rdist.earth(x1 = cbind(y[c('longitude')], lat), cbind(long, lat),
                miles = F)})
  pred_dist_lat <- lapply(df_scores, function(y)
    rdist.earth(x1 = cbind(long, y[c('latitude')]), cbind(long, lat),
                miles = F))
  pred_dist_tim <- lapply(df_scores, function(y) {
    abs(y[,'fracofyear'] - 45.25)
  })

  len_0_delayed <- sapply(delayed, length) == 0 | sapply(not_delayed, length) == 0

  pred_dist_long_del <- lapply((1:length(df_scores))[!len_0_delayed], function(y) {
    rdist.earth(x1 = cbind(df_scores[[y]][delayed[[y]],c('longitude')], lat),
                cbind(df_scores[[y]][not_delayed[[y]],c('longitude')], lat))})
  pred_dist_lat_del <- lapply((1:length(df_scores))[!len_0_delayed], function(y) {
    rdist.earth(x1 = cbind(long, df_scores[[y]][delayed[[y]],c('latitude')]),
                cbind(long, df_scores[[y]][not_delayed[[y]],c('latitude')]))})
  pred_dist_tim_del <- lapply((1:length(df_scores))[!len_0_delayed], function(y) {
    rdist(df_scores[[y]][delayed[[y]],'fracofyear'] , df_scores[[y]][not_delayed[[y]],'fracofyear'])})

  dist_long_del <- lapply((1:length(len_0_delayed))[!len_0_delayed], function(x) {
    dist_long[[x]][delayed[[x]], delayed[[x]]]
  })
  dist_lat_del <- lapply((1:length(len_0_delayed))[!len_0_delayed], function(x) {
    dist_lat[[x]][delayed[[x]], delayed[[x]]]
  })
  dist_tim_del <- lapply((1:length(len_0_delayed))[!len_0_delayed], function(x) {
    dist_tim[[x]][delayed[[x]], delayed[[x]]]
  })

  ut <- lapply(1:length(dist_long), function(y) upper.tri(dist_long[[y]], diag = T))
  dist_long_upper <-  lapply(1:length(dist_long), function(y) dist_long[[y]][ut[[y]]])
  dist_lat_upper <- lapply(1:length(dist_long), function(y) dist_lat[[y]][ut[[y]]])
  dist_tim_upper <- lapply(1:length(dist_long), function(y) dist_tim[[y]][ut[[y]]])
  dist_mat_final <- lapply(1:length(dist_long), function(y) matrix(0, nrow = nrow(dist_long[[y]]), ncol = ncol(dist_long[[y]])))
  list(list(dist_long, dist_lat, dist_tim),
       list(pred_dist_long, pred_dist_lat, pred_dist_tim),
       list(dist_long_del, dist_lat_del, dist_tim_del),
       list(pred_dist_long_del, pred_dist_lat_del, pred_dist_tim_del),
       len_0_delayed,
       delayed,
       not_delayed,
       list(dist_long_upper, dist_lat_upper, dist_tim_upper),
       ut,
       dist_mat_final)
}

m_step_vec <- function(df_use, dist_mats,
                       eig_vecs_AA, nu, n_year,
                       initial_params, K_total) {
  params_opt_temp <- list()
  cond_var <- list()
  pred <- list()
  for (k in 1:K_total) {
    response <- lapply(1:length(df_use), function(x) {
      df_use[[x]][, paste0('N', k, '_new') ]
    })
    indices <- (1:max(n_year)) * (1:max(n_year)+1)/2
    our_dpp_mat_all <- lapply(1:length(n_year), function(y) {
      new('dppMatrix', x  = rep(0, length(response[[y]])* (length(response[[y]])+1)/2), Dim = c(length(response[[y]]), length(response[[y]])))
    })
    our_full_matrix_all <- lapply(1:length(n_year), function(y) matrix(0, nrow = n_year[y], ncol = n_year[y]))
    params_opt_temp[[k]] <- optim(par = log(initial_params[k,]),
                                  fn = ll_joint_est_vec, method = 'L-BFGS-B',
                                  dist_long_upper = dist_mats[[8]][[1]],
                                  dist_lat_upper = dist_mats[[8]][[2]],
                                  dist_tim_upper = dist_mats[[8]][[3]],
                                  response = response,
                                  our_dpp_mat_all = our_dpp_mat_all,
                                  our_full_matrix_all = our_full_matrix_all,
                                  ut = dist_mats[[9]], indices = indices,
                                  nu = nu,
                                  n_year = n_year,
                                  lower = log(c(5,5,5, .1, .1)),
                                  upper = log(c(5000, 3500, 200, 100, 5700)),
                                  control = list(factr = 1e10),
                                  hessian = F)
    theta <- exp(params_opt_temp[[k]]$par)
    scale_long <- theta[1]
    scale_lat <- theta[2]
    scale_tim <- theta[3]
    sigma2 <- theta[4]
    var_score <- theta[5]

    # Pred
    pred_C <- lapply(1:length(dist_mats[[2]][[1]]), function(y) {
      Matern(sqrt( (dist_mats[[2]][[1]][[y]]/scale_long)^2+
                     (dist_mats[[2]][[2]][[y]]/scale_lat)^2+
                     (dist_mats[[2]][[3]][[y]]/scale_tim)^2),
             nu = nu)
    })
    C_full <- create_mats(dist_lon = dist_mats[[1]][[1]], dist_lat = dist_mats[[1]][[2]],
                          dist_tim = dist_mats[[1]][[3]],
                          scale_long = scale_long, scale_lat = scale_lat, scale_tim = scale_tim, nu = nu)
    C_full <- lapply(C_full, as.matrix)
    pred[[k]] <- t(sapply(1:length(C_full), function(y) {
      (t(var_score * pred_C[[y]]) %*% solve(var_score*C_full[[y]] +
                                              diag(sigma2, nrow(C_full[[y]])), response[[y]]))[1]
    }))

    cond_var[[k]] <- t(sapply(1:length(C_full), function(y) {
      (var_score+ sigma2 - t(var_score*pred_C[[y]]) %*%
         solve(var_score*C_full[[y]] + diag(sigma2, nrow(C_full[[y]])),
               var_score*pred_C[[y]]))[1]
    }))
  }
  return(list(params_opt_temp, pred, cond_var))
}

e_step <- function(params_use, df_use,
                   dist_mats,
                   eig_vecs_AA,
                   nu, K1, K2) {
  # predict at locations of real-time data
  pred <- list()
  for (k in 1:(K1+ K2)) {
    theta <- params_use[k,]
    scale_long <- theta[1]
    scale_lat <- theta[2]
    scale_tim <- theta[3]
    sigma2 <- theta[4]
    var_score <- theta[5]

    pred_C <- lapply(1:length(dist_mats[[4]][[1]]), function(y) {
      cov_mat <- Matern(sqrt( (dist_mats[[4]][[1]][[y]]/scale_long)^2+
                                (dist_mats[[4]][[2]][[y]]/scale_lat)^2+
                                (dist_mats[[4]][[3]][[y]]/scale_tim)^2),
                        nu = nu)
    })
    C_full <- create_mats(dist_lon = dist_mats[[3]][[1]], dist_lat = dist_mats[[3]][[2]],
                          dist_tim = dist_mats[[3]][[3]],
                          scale_long = scale_long, scale_lat = scale_lat, scale_tim = scale_tim,
                          nu = nu)
    C_full <- lapply(C_full, as.matrix)
    response_delayed <- lapply((1:length(dist_mats[[5]]))[!dist_mats[[5]]], function(x) {
      df_use[[x]][dist_mats[[6]][[x]], paste0('N', k, '_new') ]
    })
    pred[[k]] <-  lapply(1:length(C_full), function(y) {
      t(var_score * pred_C[[y]]) %*% solve(var_score*C_full[[y]] +
                                             diag(sigma2, nrow(C_full[[y]])), response_delayed[[y]])
    })
  }

  df_scores_new <- df_use
  # replace with expectations
  for (q in 1:length(pred[[1]])) {
    pred_year <- do.call(cbind, lapply(pred, function(y) {
      return(y[[q]])
    }))
    pred_old_scores <- t(eig_vecs_AA %*% t(pred_year))
    df_scores_new[[(1:length(dist_mats[[5]]))[!dist_mats[[5]]][q]]][
      dist_mats[[7]][[(1:length(dist_mats[[5]]))[!dist_mats[[5]]][q]]],
      paste0('S', 1:K2, '_old')] <- as.matrix(pred_old_scores[,(K1+ 1):(K1 + K2)])
  }

  df_scores_new <- do.call(rbind, df_scores_new)
  old_scores <- df_scores_new[, c(paste0('T', 1:K1, '_old'), paste0('S', 1:K2, '_old'))]
  new_scores <- t(t(eig_vecs_AA) %*% t(old_scores))
  colnames(new_scores) <- paste0('N', 1:(K1+K2), '_new')

  df_scores_new[paste0('N', 1:(K1+K2), '_new')] <- new_scores
  df_scores_new <- split(df_scores_new, df_scores_new$year)
  return(df_scores_new)
}




load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$order <- 1:nrow(grid_to_compute)

h_space <- 1100

ml_est <- function(index, h_space, K1, K2) {
  a <- proc.time()

  lat <- grid_to_compute[index,2]
  long <- grid_to_compute[index,1]

  df <- get_profile_data_subset(long = long, lat = lat, day = 45.25, h_time = 45.25, RG_def = FALSE,
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

  pressure_min <- df %>% group_by(profile) %>% summarise(min(pressure))
  pressure_max <- df %>% group_by(profile) %>% summarise(max(pressure))
  short_profiles <- (1:nrow(pressure_max))[which(pressure_max[,2] - pressure_min[,2] < 800)]
  prof_lengths <- df %>% group_by(profile) %>% summarise(n())
  small_profiles <-(1:nrow(prof_lengths))[which(prof_lengths[,2] < 15)]
  df <- df[!(df$profile %in% short_profiles) & !(df$profile %in% small_profiles),]
  df$profile <- as.numeric(as.factor(df$profile))

  if (is.null(df) | is.null(nrow(df)) | sum(df$mode =='D') ==0) {
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
  scores_temp <- t_t[[1]]
  s_s <- get_scores(df, psal_pc, K2, type = 'salinity', knots = knots_pc)
  scores_psal <- s_s[[1]]

  small_profiles <- which(is.na(scores_temp[,1]) | is.na(scores_psal[,1]))
  if (length(small_profiles) == 0) {small_profiles =0}
  scores_temp <- scores_temp[!is.na(scores_psal[,1]) & !is.na(scores_temp[,1]),]
  scores_psal <- scores_psal[!is.na(scores_psal[,1]) & !is.na(scores_temp[,1]),]

  # get longitudes and latitudes and break into year
  max_prof <- max(df$profile)
  df_scores <- df %>%
    filter(profile != small_profiles) %>%
    group_by(profile) %>%
    filter(pressure== pressure[1])
  df_years <- unique(df$year)[order(unique(df$year))]
  if (sum(df_scores$mode =='D') <2) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }

  # Use delayed mode scores to estimate A %*% t(A)
  AA_est <- cov(cbind(scores_temp[df_scores$mode == 'D',], scores_psal[df_scores$mode == 'D',]))
  if (nrow(AA_est) == 2) {
    AA_est <- cov(cbind(matrix(nrow = 1,scores_temp[df_scores$mode == 'D',]),
                        matrix(scores_psal[df_scores$mode == 'D',], nrow = 1)))

  }
  eig_AA <- eigen(AA_est)
  colnames(scores_temp) <- paste0('T', 1:K1, '_old')# original
  colnames(scores_psal) <- paste0('S', 1:K2, '_old')# original
  scores_psal[df_scores$mode != 'D', ] <- 0
  new_scores <- t(t(eig_AA$vectors) %*% t(cbind(scores_temp, scores_psal)))
  colnames(new_scores) <- paste0('N', 1:(K1+ K2), '_new') # decorrelated

  df_scores <- data.frame(df_scores[,c('longitude', 'latitude', 'fracofyear', 'mode', 'year')],
                          scores_temp, scores_psal, new_scores)
  df_scores <- df_scores[!(duplicated(df_scores$year) &  duplicated(df_scores$latitude) &
                             duplicated(df_scores$longitude)),]
  df_scores <- split(df_scores, df_scores$year)

  dist_mats <- collect_matrices_vec(df_scores, long, lat)

  n_year <- sapply(df_scores, nrow)

  m_step_results <- list()
  e_step_results <- list()
  e_step_results[[1]] <- df_scores
  params <- list()
  start_params <- c(150, 150, 30, 3, 10)
  params1 <- matrix(nrow = K1+K2, ncol = 5, 1)
  params2 <- matrix(nrow = K1+K2, ncol = 5, start_params, byrow = T)
  q <- 1

  if (sum(dist_mats[[5]]) ==length(dist_mats[[5]])) { # if all delayed
    m_step_results[[q]] <- m_step_vec(df_use = e_step_results[[q]],
                                  dist_mats,eig_vecs_AA = eig_AA$vectors, nu = .5,
                                  n_year, initial_params = params2, K_total = K1 + K2)
    params[[q]] <- t(sapply(m_step_results[[q]][[1]], function(x) exp(x$par)))
  } else {
    while (q < 7 & sum(abs(params1 - params2)) > .2) {
      m_step_results[[q]] <- m_step_vec(df_use = e_step_results[[q]],
                                    dist_mats, eig_vecs_AA = eig_AA$vectors, nu = .5,
                                    n_year = n_year, initial_params = params2, K_total = K1 + K2)
      params[[q]] <- t(sapply(m_step_results[[q]][[1]], function(x) exp(x$par)))
      e_step_results[[q+1]] <- e_step(params[[q]], df_scores,  dist_mats,
                                      eig_vecs_AA = eig_AA$vectors, nu = .5, K1 = K1, K2 = K2)
      params1 <- params2
      params2 <- params[[q]]
      q = q+1
    }
  }

  pred_df <- cbind(df_years, t(do.call(rbind, m_step_results[[length(m_step_results)]][[2]])))

  pred_old_scores <- t(eig_AA$vectors %*% t(pred_df[,2:ncol(pred_df)] ))
  pred_old_scores <- cbind(df_years, pred_old_scores)
  cond_var <- t(do.call(rbind,m_step_results[[length(m_step_results)]][[3]]))
  cond_var_old_scores <- sapply(1:nrow(cond_var), function(y) {
    eig_AA$vectors %*% Diagonal(n = K1+ K2, cond_var[y, ]) %*% t(eig_AA$vectors)
  })
  b <- proc.time()

  final_list <- list(m_step_results,
                     pred_df,
                     pred_old_scores, cond_var_old_scores,
                     cbind(df_years, n_year),
                     AA_est,params, q-1,
                     NULL,
                     c(index, (b-a)[3],
                       max_prof))
  print(final_list[[10]])
  return(final_list)
}
sequence <- as.list(seq(1, nrow(grid_to_compute), by = 1))
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
increments <- 50

profile_start <- 1 + (array_id - 1) * increments
indexes <- profile_start:(profile_start+increments-1)
# get PCs
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

# Run for a large number of grid points
if (array_id == ceiling(nrow(grid_to_compute)/increments)) {
  value <- sequence[[indexes[length(indexes)]]]-1
} else {
  value <- sequence[[indexes[length(indexes)]+1]]-1
}
object_name <-  paste('Grid_Pred', value, sep = '')
command <- paste(object_name, ' <- mclapply(indexes, ml_est,',
                 'mc.preschedule = F, mc.cores = 4,K1 = 10, K2 = 10,',
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
file.to.save <- paste0('analysis/results/joint_TS_20/joint_TS_20_',
                       value, '.RData')
eval(parse(text=paste('save(', object_name, ', file = file.to.save)', sep = '')))

q()

library(argofda)
library(Matrix)
library(fda)
library(sparseinv)
library(dplyr)
library(fields)
library(parallel)

source('analysis/code/04_joint_TS_20/ml_source_vec.R')
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
increments <- 500

profile_start <- 1 + (array_id - 1) * increments
indexes <- profile_start:(profile_start+increments-1)
# get PCs
grid_to_compute_only_ours <- grid_to_compute[grid_to_compute$order %in% indexes,]
load('analysis/results/psal_one_stage_cov_pca.RData')
coefs_pc_psal <- coefs_pc_psal[grid_pc_psal$long %in% grid_to_compute_only_ours$long &
                                 grid_pc_psal$lat %in% grid_to_compute_only_ours$lat]
grid_pc_psal <- grid_pc_psal[grid_pc_psal$long %in% grid_to_compute_only_ours$long &
                               grid_pc_psal$lat %in% grid_to_compute_only_ours$lat,]
load('analysis/results/temp_one_stage_cov_pca.RData')
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
                 'mc.preschedule = F, mc.cores = 12,K1 = 10, K2 = 10,',
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

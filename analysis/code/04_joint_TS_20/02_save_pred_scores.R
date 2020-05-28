calc_files <- list.files('analysis/results/joint_TS_20/',
                         pattern = 'joint_TS')

# Get results from score modelling
final_list <- list()
K1 <- 10
K2 <- 10
K <- K1 + K2
for (j in 1:length(calc_files)) {
  load(paste0('analysis/results/joint_TS_20/', calc_files[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  eval(parse(text=paste0('rm(', file, ')')))

  indexes_used <- as.numeric(sapply(now, function(x) {
    if (length(x) == 3) {
      x[1]
    } else {
      x[[10]][1]
    }
  }))
  df_return <- lapply(1:length(now), function(x) {
    if (length(now[[x]]) == 1 | length(now[[x]]) == 0 | length(now[[x]]) == 3) {
      mat <- matrix(c(indexes_used[x],
                      rep(NA, K + K^2+ 1)),nrow = 1, ncol = 1 + K + K^2 +1)
    } else {
      n_steps <- length(now[[x]][[1]])
      m_res <- now[[x]][[1]][[n_steps]][[1]]
      cond_var <- t(sapply(now[[x]][[4]], as.vector))
      mat <- cbind(indexes_used[x], now[[x]][[3]] ,cond_var)
    }
    colnames(mat) <- c('index',
                       'year', paste0('T', 1:K1), paste0('S', 1:K2),
                                paste0('condvar', rep(1:K, each = K), '_', rep(1:K, times = K)))
    mat
  })

  df_return <- do.call(rbind, df_return)
  final_list[[j]] <- list(df_return, indexes_used)
}

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$index <- 1:nrow(grid_to_compute)
df_params <- lapply(final_list, function(x) x[[1]])
df_params <- do.call(rbind, df_params)
df_params <- df_params[!is.na(df_params[,2]),]
df_preds <- merge(grid_to_compute,df_params, by = 'index')
save(df_preds, file = 'analysis/results/joint_TS_20/preds.RData')

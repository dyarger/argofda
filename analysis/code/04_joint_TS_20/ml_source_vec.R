
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



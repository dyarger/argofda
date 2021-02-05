

cov_est_current <-  function(index,h_space = 550, knots, lambda, basis, penalty_mat,
                             min_prof) {
  a <- proc.time()
  lat <- grid_to_compute_merge[index,2]
  long <- grid_to_compute_merge[index,1]
  RG_check = TRUE
  if (RG_check == TRUE) {
    RG_def <- grid_to_compute_merge[index,3]
  } else {
    RG_def <- T
  }
  df <- get_profile_data_subset(lat = lat, long = long, day = 45.25, RG_def = RG_def, h_time = 45.25,
                                h_space = h_space, mode = 'all',
                                min_prof =  min_prof, years = 2007:2016,
                                exclusion= TRUE)
  df$pressure <- round(as.numeric(df$pressure), digits = 2)

  if (is.null(df) | is.null(nrow(df)) ) {
    b <- proc.time()
    print(c(index, 'no profiles in range', (b-a)[3]))
    return(c(index, 'no profiles in range', (b-a)[3]))
  }
  h_space_evaluated <- df$h_space[1]
  max_prof <- max(df$profile)
  df <- df[, c('pressure', 'temperature', 'profile', 'fracofyear', 'wt')]

  df <- df[df$pressure <2000 & df$pressure > 0,]
  gdf = split(df, df$profile)
  press = lapply(gdf, function(x){return(x$pressure)})
  temp = lapply(gdf, function(x){return(x$temperature)})

  weights1 = as.numeric(lapply(gdf, function(x){return(x$wt[[1]])}))
  weights2 = as.numeric(lapply(gdf, function(x){return(x$wt[[1]]/nrow(x))}))
  weights <- weights2
  cv_info <-   capture.output(
    coefs_temp <- unlist(penalizedCovariance(indexes = press, values = temp, lower = 0.0,
                                             upper = 2000.0, nbasis = 102,
                                             order = 4,
                                             weights= weights, use_cv = T,
                                             mixedOnly = T, lambda = 1000,splinePenalty = as.matrix(penalty_mat), #knots = knots,
                                             lower_optim = 10^3, upper_optim = 10^10))
  )
  cv_info <- as.numeric(strsplit(cv_info, split = "\\s+")[[1]])
  cv_info <- data.frame(matrix(ncol = 8, byrow = T,cv_info)[,c(2,4,6,8)])
  colnames(cv_info) <- c("SSE", "Trace", "lambda", "gcv")
  n <- nrow(df)
  nprof <- length(gdf)
  m_total <- sum(table(df$profile) *(table(df$profile)-1))

  b <- proc.time()
  print(c(index, n, (b-a)[3],nprof,m_total,
          max_prof, h_space_evaluated))
  return(list(coefs_temp,cv_info,
              c(index, n, (b-a)[3],nprof,m_total,
                max_prof, h_space_evaluated)))
  # basis <- create.bspline.basis(breaks = knots)
  # omega_0 <- as(bsplinepen(basis,Lfdobj = 0), 'sparseMatrix')
  # omega_0_sqrt <- as(chol(omega_0), 'sparseMatrix')
  #
  # A_mat = matrix(0, length(knots)+2, length(knots)+2)
  # A_mat[lower.tri(A_mat, diag=TRUE)] = coefs_temp
  # A_mat = as.matrix(forceSymmetric(A_mat, 'L'))
  # #A_mat <- matrix(unlist(x), nrow = length(knots) + 2)
  # pca_mat <- omega_0_sqrt %*% A_mat %*% t(omega_0_sqrt)
  # standard_pca <- prcomp(pca_mat, center = FALSE, scale = FALSE)
  # vectors <- solve(omega_0_sqrt, standard_pca$rotation)

}

cov_est_salinity <-  function(index,h_space = 550, knots, lambda, basis, penalty_mat,
                              min_prof) {
  a <- proc.time()
  lat <- grid_to_compute_merge[index,2]
  long <- grid_to_compute_merge[index,1]
  RG_check = TRUE
  if (RG_check == TRUE) {
    RG_def <- grid_to_compute_merge[index,3]
  } else {
    RG_def <- T
  }
  df <-  get_profile_data_subset(lat = lat, long = long, day = 45.25, RG_def = RG_def, h_time = 45.25,
                                 h_space = h_space, mode = 'delayed',
                                 min_prof =  min_prof, years = 2007:2016,
                                 exclusion= TRUE)
  df$pressure <- round(as.numeric(df$pressure), digits = 2)

  if (is.null(df) | is.null(nrow(df)) ) {
    b <- proc.time()
    print(c(index, 'no profiles in range', (b-a)[3]))
    return(c(index, 'no profiles in range', (b-a)[3]))
  }
  h_space_evaluated <- df$h_space[1]
  max_prof <- max(df$profile)
  df <- df[, c('pressure', 'salinity', 'profile', 'fracofyear', 'wt')]

  df <- df[df$pressure <2000 & df$pressure > 0,]
  gdf = split(df, df$profile)
  press = lapply(gdf, function(x){return(x$pressure)})
  temp = lapply(gdf, function(x){return(x$salinity)})
  weights1 = as.numeric(lapply(gdf, function(x){return(x$wt[[1]])}))
  weights2 = as.numeric(lapply(gdf, function(x){return(x$wt[[1]]/nrow(x))}))
  weights <- weights2
  cv_info <-   capture.output(
    coefs_temp <- unlist(penalizedCovariance(indexes = press, values = temp, lower = 0.0,
                                             upper = 2000.0, nbasis = 102,
                                             order = 4,
                                             weights= weights, use_cv = T,
                                             mixedOnly = T, lambda = 1000,splinePenalty = as.matrix(penalty_mat), #knots = knots,
                                             lower_optim = 10^3, upper_optim = 10^10))
  )
  cv_info <- as.numeric(strsplit(cv_info, split = "\\s+")[[1]])
  cv_info <- data.frame(matrix(ncol = 8, byrow = T,cv_info)[,c(2,4,6,8)])
  colnames(cv_info) <- c("SSE", "Trace", "lambda", "gcv")
  n <- nrow(df)
  nprof <- length(gdf)
  m_total <- sum(table(df$profile) *(table(df$profile)-1))

  b <- proc.time()
  print(c(index, n, (b-a)[3],nprof,m_total,
          max_prof, h_space_evaluated))
  return(list(coefs_temp,cv_info,
              c(index, n, (b-a)[3],nprof,m_total,
                max_prof, h_space_evaluated)))
}

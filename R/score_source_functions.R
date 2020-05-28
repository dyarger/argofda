

# Covariance Source functions

#' from a data frame and principal components, computes scores for each profile
#'
#' @param df data.frame with columns profile, pressure, temperature | salinity
#' @param pc matrix of nbasis*K with coefficients for the first K principal components
#' @param K number of principal component scores to compute
#' @param type temperature or salinity
#' @param knots knots used for principal component functions
#'
#' @return list of 3: scores, residuals after subtracting scores, and score representation of profile
#'
#' @export
get_scores <- function(df, pc, K, type = 'temperature', knots) {
  scores <- matrix(nrow = length(unique(df$profile)), ncol = K)
  residuals <- list(); score_rep <- list()
  unique_prof <-  unique(df$profile)
  prof_list <- split(df, df$profile)
  fits <- lapply(1:K, function(x) make_fit_object(knots = knots, coefficients =pc[,x]))
  for (k  in 1:length(unique_prof)) {
    prof <- prof_list[[k]]
    pc1atprof1 <- sapply(1:K, function(x) {
      predict(fits[[x]], prof$pressure)$y
    })
    if (type == 'temperature') { # least squares fit
      scores[k,] <- .lm.fit(pc1atprof1,as.matrix(prof[,'temperature']) )$coefficients
      score_rep[[k]] <- rowSums(sapply(1:K, function(z) scores[k,z] * pc1atprof1[,z]))
      residuals[[k]] <- prof$temperature - score_rep[[k]]
    } else if (type == 'salinity') {
      scores[k,] <- .lm.fit(pc1atprof1,as.matrix(prof[,'salinity']) )$coefficients
      score_rep[[k]] <- rowSums(sapply(1:K, function(z) scores[k,z] * pc1atprof1[,z]))
      residuals[[k]] <- prof$salinity - score_rep[[k]]
    }
  }
  return(list(scores, residuals,score_rep))
}

# maximum likelihood estimation

#' computes likelihood of matern covariance
#'
#' @param theta logged parameters vector of length 5
#' @param dist_long list of length n_year, each is a matrix of distances for longitude
#' @param dist_lat list of length n_year, each is a matrix of distances for latitude
#' @param dist_tim list of length n_year, each is a matrix of distances for time
#' @param response list of length n_year, each is a vector of residuals
#' @param nu chosen value of nu for matern smoothness
#' @param n_year vector of length n_year, the number of profiles used for each year
#'
#' @return value of likelihood
#'
#' @export
ll_joint_est <- function(theta, dist_long, dist_lat, dist_tim,
                         response, nu, n_year) {
  scale_long <- exp(theta[1])
  scale_lat <- exp(theta[2])
  scale_tim <- exp(theta[3])
  sigma2 <- exp(theta[4])
  var_score <- exp(theta[5])
  
  #print(c(scale_long,scale_lat, scale_tim, sigma2, var_score))
  C_full <- create_mats_upper(dist_long, dist_lat , dist_tim, # covariance matrix
                        scale_long, scale_lat, scale_tim,
                        nu = nu, n_year = n_year) 
  # calculate log(det(Sigma)) and t(Y) %*% solve(Sigma,Y) for each year
  # if chol runs into an issue, return NA
  dets <-   tryCatch( {sapply(1:length(C_full), function(x) {
    mat1 <- var_score*C_full[[x]] 
    diag(mat1) <- diag(mat1) + sigma2
    L <- chol(mat1) 
    det <- sum(log(diag(L)))
    v2 <- .Internal(backsolve(r = L, x = response[[x]], k = ncol(L), upper.tri = T,transpose = T))
    v1 <- sum(v2^2)
    return(c(2*sum(det), v1))
  })}, error = function(z) {NA})
  if (length(dets) == 1) {return(NA)}
  (as.double(.5 * sum(dets[1,] + dets[2,])+ sum(n_year) * log(2*pi))) #sum over all years
}


#' create covariance matrix for each year
#'
#' @param dist_long list of length n_year, each is a matrix of distances for longitude
#' @param dist_lat list of length n_year, each is a matrix of distances for latitude
#' @param dist_tim list of length n_year, each is a matrix of distances for time
#' @param scale_long value of longitude scale parameter
#' @param scale_lat value of latitude scale parameter
#' @param scale_tim value of time scale parameter
#' @param nu chosen value of nu for matern smoothness
#'
#' @return list of length n_year, each is a covariance matrix specified
#'
#' @export
create_mats <- function(dist_long, dist_lat, dist_tim,
                        scale_long, scale_lat, scale_tim, nu) {
  if (nu != .5) {
    lapply(1:length(dist_long), function(y) {
      Matern(sqrt( (dist_long[[y]]/scale_long)^2  + (dist_lat[[y]]/scale_lat)^2 + 
                          (dist_tim[[y]]/scale_tim)^2), nu = nu)
    })
  } else {
    lapply(1:length(dist_long), function(y) {
      Exponential(sqrt( (dist_long[[y]]/scale_long)^2  + (dist_lat[[y]]/scale_lat)^2 + 
                          (dist_tim[[y]]/scale_tim)^2))
    })
  }
}

#' create upper triangle of covariance matrix for each year (faster for ML)
#'
#' @param dist_long list of length n_year, each is a matrix of distances for longitude
#' @param dist_lat list of length n_year, each is a matrix of distances for latitude
#' @param dist_tim list of length n_year, each is a matrix of distances for time
#' @param scale_long value of longitude scale parameter
#' @param scale_lat value of latitude scale parameter
#' @param scale_tim value of time scale parameter
#' @param nu chosen value of nu for matern smoothness
#' @param n_year vector of length n_year, the number of profiles used for each year
#'
#' @return list of length n_year, each is the upper triangle of covariance matrix specified
create_mats_upper <- function(dist_long, dist_lat, dist_tim,
                              scale_long, scale_lat, scale_tim, nu, n_year) {
  lapply(1:length(dist_long), function(y) {
    .Call("ExponentialUpperC", sqrt( (dist_long[[y]]/scale_long)^2  + (dist_lat[[y]]/scale_lat)^2 +
                                       (dist_tim[[y]]/scale_tim)^2), n_year[y],
          1, PACKAGE = "fields")
  })
}

####################################################################################
# the following has versions of the above functions that do not deal with the entire matrix
# but a vector that contains the elements of the upper triangle of the matrix
# more efficient memory-wise

#' computes likelihood of matern covariance
#'
#' @param theta logged parameters vector of length 5
#' @param dist_long_upper list of length length(n_year), each is a vector of distances for longitude
#' @param dist_lat_upper list of length length(n_year), each is a vector of distances for latitude
#' @param dist_tim_upper list of length length(n_year), each is a vector of distances for time
#' @param response list of length length(n_year), each is a vector of residuals
#' @param nu value of nu for matern smoothness - only .5 works for now
#' @param n_year vector of length length(n_year), the number of profiles used for each year
#' @param ut matrix that gives the upper triangular elements of matrix of dimensions n_year(y)^2  for y=1,...,10
#' @param indices vector that gives the diagonal elements of the vector of upper triangular elements
#' @param our_dpp_mat_all list with matrices, each with only upper triangular elements non-sparse, that will be filled with covariances
#' @param our_full_matrix list with matrices, each with only upper triangular elements non-sparse, that will be filled with cholesky
#'
#' @return value of likelihood
#'
#' @export
ll_joint_est_vec <- function(theta, dist_long_upper, dist_lat_upper, dist_tim_upper, response, nu, 
                             n_year, ut, indices, our_dpp_mat_all, our_full_matrix_all) {
  scale_long <- exp(theta[1])
  scale_lat <- exp(theta[2])
  scale_tim <- exp(theta[3])
  sigma2 <- exp(theta[4])
  var_score <- exp(theta[5])
  C_full <- lapply(1:length(dist_long_upper), create_mats_upper_vec,
                   var_score = var_score, dist_long_upper= dist_long_upper, dist_lat_upper=dist_lat_upper, 
                   dist_tim_upper = dist_tim_upper,
                   scale_long = scale_long, scale_lat = scale_lat, scale_tim = scale_tim,indices = indices, sigma2 = sigma2)
  dets <- sapply(1:length(C_full), get_dets, response = response,
                 our_dpp_mat = our_dpp_mat_all, our_full_matrix = our_full_matrix_all, C_full = C_full,  sigma2 = sigma2,ut = ut)
  (as.double(0.5 * sum(dets[1, ] + dets[2, ])+sum(n_year) * log(2 * pi)))
}

#' create vector for upper triangle of covariance matrix for each year
#'
#' @param y index of matrices dist_lon_upper, etc to use, generally in 1:10
#' @param var_score value of varaince parameter (gamma in paper)
#' @param dist_long_upper list of length n_year, each is a vector of distances for longitude
#' @param dist_lat_upper list of length n_year, each is a vector of distances for latitude
#' @param dist_tim_upper list of length n_year, each is a vector of distances for time
#' @param scale_long value of longitude scale parameter
#' @param scale_lat value of latitude scale parameter
#' @param scale_tim value of time scale parameter
#' @param indices indices of vector that denote the diagonal elements
#' @param sigma2 value of nugget variance parameter
#'
#' @return vector of upper triangle covariances based on the supplied parameters and distance matrices
create_mats_upper_vec <- function(y, var_score, dist_long_upper, dist_lat_upper, dist_tim_upper,
                                  scale_long, scale_lat, scale_tim, indices,sigma2) {
  cov_vec <- var_score*exp(-sqrt( (dist_long_upper[[y]]/scale_long)^2  + (dist_lat_upper[[y]]/scale_lat)^2 +
                                    (dist_tim_upper[[y]]/scale_tim)^2)) 
  # diagonal
  cov_vec[indices[indices <= length(cov_vec)]] <- cov_vec[indices[indices <= length(cov_vec)]]+ sigma2
  cov_vec
}

#' computes cholesky and evaluates parts of likelihood for each year
#'
#' @param y index to use, generally in 1:10
#' @param response list of length length(n_year), each is a vector of residuals
#' @param C_full vector with upper triangular elements of covariance matrix
#' @param ut matrix that gives the upper triangular elements of matrix of dimensions n_year(y)^2  for y=1,...,10
#' @param our_dpp_mat_all list with matrices, each with only upper triangular elements non-sparse, that will be filled with covariances
#' @param our_full_matrix list with matrices, each with only upper triangular elements non-sparse, that will be filled with cholesky
#'
#' @return c(log(det(Sigma)), y^T Sigma^(-1) y) where Sigma is a covariance matrix
#'
#' @export
get_dets <- function(y, response,our_dpp_mat, our_full_matrix, C_full, sigma2,ut) {
  our_dpp_mat[[y]]@x <- C_full[[y]]
  L_chol <- Matrix::chol(our_dpp_mat[[y]])
  det <- Matrix::determinant(L_chol, logarithm = T)$modulus # log(det(sigma))
  our_full_matrix[[y]][ut[[y]]] <- L_chol@x
  v1 <- sum(.Internal(backsolve(r = our_full_matrix[[y]], x = response[[y]],
                                k = ncol(L_chol), upper.tri = T, transpose = T))^2) # y^T sigma^(-1) y
  return(c(2 * det, v1))
}



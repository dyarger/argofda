
#' Creates weight matrix
#'
#' @param p vector of pressure from all profiles included in the fit
#' @param weights vector of kernel weights for the measurements
#' @param profile vector that groups measurements into profiles
#' @param rho correlation parameter
#' @param specify_cor if FALSE, rho is not used and gives diagonal matrix
#'
#' @return tridiagonal sparse matrix of weights to use in the fit
#' 
#' @export
construct_working_cov <- function(p, weights, profile, rho, specify_cor) {
  prof_lengths <- aggregate(weights, list(profile), length)
  W <- weights/(rep(prof_lengths$x, times = prof_lengths$x)) #normalize by profile length
  if (specify_cor == TRUE) {
    cor_mats <- lapply(split(p, profile), invert_continuous, rho = rho)
    if (length(cor_mats) > 1) {
      sigma_inv <- do.call(Matrix::bdiag, cor_mats)
    } else{
      sigma_inv <- cor_mats[[1]]
    }
    W_final <- W * sigma_inv
  } else {
    W_final <- Matrix::Diagonal(n = length(W), x = W)
  }
  return(W_final)
}

#' Creates preliminary matrices for a fixed lambda
#'
#' @param covariates matrix of covariates used in the fit; a single column of 1s will give a simple smoothing spline estimate
#' @param lambdas vector of length ncol(covariates) that specifies initial lambda for each covariate
#' @param Phi basis matrix for data used in fit (should have nrow = nrow(covariates))
#' @param Psi penalty matrix
#' @param W weight matrix created by construct_working_cov
#' @param Y response vector
#'
#' @return list of matrices
#' U - Phi_all^T W Phi_all  
#' V - Phi_all^T W Y
#' Omega - penalty matrix for all covariates
#' the solution is (U + Omega)^(-1) V
#' reorder - signifies how the basis functions were reordered to yield sparse matrices
#' Phi_all - basis matrix for all covariates
#' 
#' @export
model_set_up <- function(covariates, lambdas, Phi, Psi, W, Y) {
  covariates <- as.matrix(covariates)
  n <- nrow(covariates)
  Phi_mats <- lapply(1:ncol(covariates), function(x) {
    return(Matrix::Diagonal(n = n, covariates[,x]) %*% Phi)
  })
  Phi_all  <- do.call(cbind, Phi_mats)
  Omega_mats <- lapply(1:ncol(covariates), function(x) {
    return(lambdas[x] * Psi)
  })
  Omega <- do.call(Matrix::bdiag, Omega_mats)
  nfunc <- ncol(Phi_all)/ncol(Phi)
  selection_vector <- as.vector(sapply(1:ncol(Phi),
                                       function(x) {
                                         return(x + 0:(nfunc - 1) * ncol(Phi))
                                       }))
  Phi_all <- Phi_all[,selection_vector]
  Omega <- as(Omega[selection_vector, selection_vector], 'symmetricMatrix')
  Phi_t <- Matrix::t(Phi_all)
  U <- Phi_t %*% W %*% Phi_all
  V <- Phi_t %*% (W %*% Y)
  return(list('U' = U, 'V' = V, 'Omega' = Omega, 'reorder' = selection_vector,'Phi_all' =  Phi_all))
}

#' Uses golden section/quadratic search to choose amount of smoothing
#'
#' @param U matrix of covariates used in the fit; a single column of 1s will give a simple smoothing spline estimate
#' @param V vector of length ncol(covariates) that specifies initial lambda for each covariate
#' @param W weight matrix created by construct_working_cov
#' @param Omega penalty matrix
#' @param Y response vector
#' @param Phi_all matrix from model_set_up that contains all basis functions evaluated
#' @param lower_optim lower limit for cross validation search
#' @param upper_optim upper limit for cross validation search
#' @param tol tolerance for optimization - lower tolerance means more computation
#' @param selection_vector vector from model_set_up that reorders basis functions
#' @param knots knots used in basis
#'
#' @return local_function - list of
#' 1: beta - list of coefficients for each function
#' 2: knots - knots used for basis
#' 3: cv_info - cross validation summary information
#' 4: lambda_chosen - lambda_chosen by cross validation
#' 
#' @export
complete_optim <- function(U, V, W, Omega, Y, Phi_all,
                           lower_optim, upper_optim, tol = 1e-03, selection_vector, 
                           knots, n_basis) {
  the_trace <- sum(Matrix::diag(Omega)) / sum(Matrix::diag(U))
  cv_info <- capture.output(l_opt <- optimise(f = gcv_fun,
                                              lower = (log(base = 16, the_trace*lower_optim) + 2)/6,
                                              upper = (log(base = 16, the_trace*upper_optim) + 2)/6,
                                              tol = tol,
                                              U = U, V = V, W = W,Phi_all = Phi_all, Omega = Omega,Y = Y))
  # Compute solution at the minimizer
  the_chol_LL <- Cholesky(forceSymmetric(U + 1/the_trace * 16^(l_opt$minimum * 6 - 2) * Omega), perm = F)
  beta <- Matrix::solve(the_chol_LL, V, system = 'A')
  beta <- beta[order(selection_vector)]
  beta <- split(beta, rep(1:(ncol(Phi_all)/n_basis), each = n_basis))
  
  # organize cv information
  cv_info <- data.frame(do.call(rbind, strsplit(cv_info, split = '\\s+'))[,2:5],
                        stringsAsFactors = FALSE)
  colnames(cv_info) <- c('lambda', 'SSE', 'df', 'gcv')
  lambda_chosen <- 1/the_trace * 16^(l_opt$minimum  * 6 - 2)
  ret_done <- list('beta' = beta, 'knots' = knots, 'cv_info' = cv_info, 'lambda_chosen' = lambda_chosen)
  class(ret_done) <- 'local_function'
  return(ret_done)
}
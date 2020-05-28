
#' calculates gcv score for a given lambda
#'
#' @param lambda value of lambda to multiply Omega, but on log scale
#' @param U matrix made from model_set_up
#' @param V matrix made from model_set_up
#' @param W matrix made from construct_working_cov
#' @param Phi_all matrix made from model_set_up
#' @param Omega matrix made from model_set_up
#' @param Y response variable
#'
#' @return value of gcv for specific lambda
#'
#' @export
gcv_fun <- function(lambda, U, V, W, Phi_all, Omega, Y) {
  the_trace <- sum(Matrix::diag(Omega)) / sum(Matrix::diag(U))
  lambda_eval <- 1/the_trace * 16^(lambda * 6 - 2)
  B <- forceSymmetric(U + lambda_eval * Omega)
  the_chol_LL <- Cholesky(B, perm = F)
  chol_final <- expand(the_chol_LL)$L
  beta <- Matrix::solve(the_chol_LL, V, system = 'A')
  B_inv <- TD_adj(Q = B, cholQp = chol_final, 
                  P = expand(the_chol_LL)$P, gc = 0)
  A_trace_test <- sum(U * B_inv)
  A_trace <- ifelse(A_trace_test > 0, min(c(A_trace_test, ncol(U))),ncol(U))
  pred <- Phi_all %*% beta
  residuals <- as.matrix(Y - pred)
  SSE <- t(residuals) %*% W %*% (residuals)
  gcv <- (1/sum(W))*SSE/( (1/nrow(residuals)) * (nrow(residuals) - A_trace))^2
  results <- c(lambda_eval, as.double(SSE), A_trace, as.double(gcv))
  print(results)
  return(gcv[1])
}

#' Based on sparseinv::Takahashi_Davis, with slight changes for speed/memory
#' Computes certain elements of the inverse of a sparse matrix
#'
#' @param Q sparse matrix for which elements of the inverse are wanted
#' @param cholQp cholesky decomposition of Q with permutation matrix P
#' @param return_perm_chol should the function return the cholesky decomposition as well?
#' @param P permutation matrix for cholesky decomposition
#' @param gc should garbage collection be run during function
#'
#' @return sparse matrix that contains certain elements of Q^{-1}
#'
#' @export
TD_adj <- function (Q = NULL, cholQp = NULL, return_perm_chol = 0, P, 
                    gc = 0) {
  # sparseinv:::.check_args_Takahashi(Q = Q, cholQp = cholQp, return_perm_chol = return_perm_chol, 
  #                                   P = P, gc = gc)
  n <- nrow(Q)
  if (is.null(cholQp)) {
    X <- cholPermute(Q = Q)
    L <- X$Qpermchol
    P <- X$P
  } else {
    L <- cholQp
    P <- P
  }
  if (return_perm_chol == 0) 
    rm(cholQp)
  d <- Matrix::diag(L)
  L <- tril(L %*% Diagonal(x = 1/d), -1)
  d <- d^2
  D <- Diagonal(x = d)
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp)
  if (gc) 
    gc()
  Zpattern <- sparseMatrix(c(L@i + 1, 1:n), c(jj, 1:n), symmetric = T, giveCsparse = F)
  Zpattern <- as(Zpattern, 'nsCMatrix')
  
  
  rm(dp, jj)
  if (gc) 
    gc()
  Z <- sparseinv:::.sparseinv_wrapper(L, d, L, Zpattern)
  if (return_perm_chol == 0) {
    return(P %*% Z %*% Matrix::t(P))
  } else {
    return(list(S = P %*% Z %*% Matrix::t(P), Lp = cholQp, P = P))
  }
}

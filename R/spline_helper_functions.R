#' Makes a smooth.spline object that can be easily used to predict
#'
#' @param knots vector of knots of splines
#' @param coefficients vector of coefficients of the spline
#'
#' @return smooth.spline.fit object, one can use predict() with it
#'
#' @export
make_fit_object <- function(knots, coefficients) {
  min.knot <- min(knots)
  range.knot <- range(knots)[2] - range(knots)[1]
  fit <- list(knot = c(0,0,0, (knots - min.knot)/range.knot, 1, 1, 1),
              nk = length(knots)+2,
              min = min.knot,
              range = range.knot,
              coef = as.vector(coefficients))
  fit <- structure(fit, class = "smooth.spline.fit")
  return(fit)
}


#' creates an OU(1) matrix given pressures and rho
#'
#' @param pressure_vec observed pressures of a single profile
#' @param rho correlation between measurements 1 decibar apart
#'
#' @return sparse tridiagonal precision matrix
#'
#' @export 
invert_continuous <- function(pressure_vec, rho) {
  ni <- length(pressure_vec)
  
  lag_calc <- pressure_vec[1:(length(pressure_vec)-1)] - pressure_vec[2:length(pressure_vec)]
  cor_band1 <- exp(log(rho) * abs(lag_calc))
  cor_band2 <- cor_band1[1:(length(cor_band1)-1)] * cor_band1[2:length(cor_band1)]
  
  dets <- 1/ ((1-cor_band1[2:length(cor_band1)]^2) -
                cor_band1[1:(length(cor_band1)-1)] *
                (cor_band1[1:(length(cor_band1)-1)] - cor_band1[2:length(cor_band1)] * cor_band2))
  prelim_diag <- dets * (1-cor_band2^2)
  prelim_band <- dets * (-cor_band1[1:(length(cor_band1)-1)] + cor_band1[2:length(cor_band1)] * cor_band2)
  
  final_diag <- c(1/(1-cor_band1[1]^2), prelim_diag, 1/(1-cor_band1[length(cor_band1)]^2))
  final_band <- c(prelim_band, -cor_band1[length(cor_band1)]/(1-cor_band1[length(cor_band1)]^2))
  
  inverse_matrix <- Matrix::bandSparse(n = ni, k = c(0,1),
                                       diagonals = list(final_diag, final_band),
                                       symmetric = T)
  
  if (ni == 2) {
    inverse_matrix <- matrix(nrow = 2, ncol = 2,data =
                               1/(1-cor_band1[1]^2)* c(1, -cor_band1[1],
                                                       -cor_band1[1], 1))
  }
  return(inverse_matrix)
}


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### The rest of the functions here are not currently used, but are faster than fda::eval.basis
###### Based on smooth.spline function in fortran
#' creates a matrix of B-spline evaluations
#'
#' @param x observed values for the functions to be evaluated
#' @param knots knots of B-splines
#'
#' @return sparse matrix of B-spline evaluations
#'
#' @export 
basis_mat <- function(x,knots) {
  ileft <- vector(mode = 'integer', length = length(x))
  t <- c(min(knots), min(knots), min(knots), knots, max(knots), max(knots), max(knots))

  for (i in 1:length(x)) {
    ileft[i] <- interv(x[i], t)
    if(ileft[i] == length(t)) { ileft[i] = length(t)-4}
  }

  vnikx <- vector(mode = 'double', length = length(x)*4)
  for (i in 1:length(x)) {
    vnikx[(i*4 - 3):(i*4)] <- bsplvd(t,4,x[i],ileft[i],1)
  }
  jval <- vector(mode = 'integer', length = 4*length(x))
  jval[1:length(jval) %% 4 == 1] = ileft-3
  jval[1:length(jval) %% 4 == 2] = ileft-2
  jval[1:length(jval) %% 4 == 3] = ileft-1
  jval[1:length(jval) %% 4 == 0] = ileft

  G1 <- sparseMatrix(i = rep(1:length(x), each = 4), j =jval, x = vnikx, 
                     dims = c(length(x),length(knots)+2))
  return(G1)
}

interv <- function(xi,t) {
  left <- findInterval(xi, t, rightmost.closed = FALSE, all.inside = FALSE,
                       left.open = FALSE)
  return(left)
}


bsplvd <- function(t, k = 4, xp, left, nderiv =1) {
  mhigh = max(min(nderiv,k),1)
  kp1 = k+1
  return(bsplvb(t, kp1-mhigh,1,xp,left, matrix(0,nrow = 4, ncol = 1)))
}

bsplvb <- function(t, jhigh = 4, index, xq, left, biatx = matrix(0,nrow = 4, ncol = 1)) {
  jmax=20
  deltar= c()
  deltal= c()

  if (index ==1) {
    j=1
    biatx[1,1] =1
    if (j >= jhigh) {return(biatx)}
  }
  while (j < jhigh) {
    jp1 = j + 1
    deltar[j] = t[left+j] - xq
    deltal[j] = xq - t[left+1-j]
    saved= 0
    for (i in 1:j) {
      term = biatx[i,1]/(deltar[i] + deltal[jp1-i])
      biatx[i,1] = saved + deltar[i]*term
      saved = deltal[jp1-i]*term
    }
    biatx[jp1,1] = saved
    j = jp1
  }
  return(biatx)
}



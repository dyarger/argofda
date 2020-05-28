

#' Given a 'local_function' object, predicts for one of the functions
#'
#' @param x local_function object that contains coefficients, knots, etc
#' @param p_vals pressure values to predict at
#' @param index which function to predict - 0 = first covariate, 1 = second, etc
#' @param deriv derivative in pressure of spline to compute
#'
#' @return vector of estimates at p_vals
#'
#' @export
predict.local_function <- function(x, p_vals, index = 0, deriv = 0) {
  predict(make_fit_object(knots = x[[2]], coefficients = x[[1]][[index+1]]),
                 x = p_vals, deriv = deriv)$y
}

#' Given a 'local_function' object, plots for one of the functions, or cross validation info
#'
#' @param x local_function object that contains coefficients, knots, etc
#' @param index which function to plot - 0 = first covariate, 1 = second, etc
#' @param cv if true, plots log10(lambda) versus gcv
#' @param deriv derivative in pressure of spline to compute
#'
#' @return a plot
#'
#' @export
plot.local_function <- function(x,
                                index = 0, # which function: 0 = mean, 1 = first local linear
                                cv = FALSE,
                                deriv = 0) {# which response variable, by column of Y
  if (cv == TRUE) {
    results <- data.frame(x[[3]], stringsAsFactors = FALSE)
    results$lambda <- as.numeric(results$lambda)
    results$gcv <- as.numeric(results$gcv)
    plot(log10(results[,1]), results[,4], main = 'GCV Results', xlab = 'log10(lambda_vec)',
         ylab = 'gcv')
  } else {
    pred <- predict(make_fit_object(knots = x[[2]], coefficients = x[[1]][[index+1]]),
                    seq(from = min(x[[2]]), to = max(x[[2]]), length.out = 300), deriv = deriv)
    plot(pred$x, pred$y, type = 'l', main = paste0('Spline Estimate, index ', index), xlab = 'p', ylab = 'Y')
  }
}


#' Given a 'local_function' mean as described in the paper, compute estimated mean
#'
#' @param x local_function object that contains coefficients, knots, etc
#' @param p_vals pressure values to compute mean
#' @param deriv derivative in pressure to compute
#' @param day day of the year to use
#' @param dist_EW distance from grid point in East-West direction
#' @param dist_NS distance from grid point in North-South direction
#' @param year year to evaluate mean
#'
#' @return vector of length = length(p_vals) with the mean evaluated with the above covariates
#'
#' @export
predict_time_year <-  function(x, p_vals, deriv = 0, day = 0, dist_EW, dist_NS, year)  {
  preds <- sapply(c(0:16), predict.local_function,
                  x = x, p_vals = p_vals, deriv = deriv)
  t_s <- day - 45.25
  coef <- matrix(c(ifelse(2007:2016 == year, 1, 0), dist_EW,dist_NS,dist_EW^2,
                   dist_NS^2,dist_EW*dist_NS, t_s, t_s^2),
                 nrow = nrow(preds), ncol = 17, byrow = TRUE)
  rowSums(coef * preds)
}

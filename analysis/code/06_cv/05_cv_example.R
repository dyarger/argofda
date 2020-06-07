library(argofda)
load('analysis/data/jan_march_residuals.RData')
load('analysis/results/cv/cv_results_8.RData')
profile_nums <- sapply(cv_results, function(x) x[[7]][1])
single_results <- cv_results[[574]]
residuals <- single_results[[3]]
cond_var <- single_results[[4]]
single_results[[5]]
K1 <- 10
K2 <- 10

profLatAggr[7332]
profLongAggr[7332]
profFloatIDAggr[7332]
profCycleNumberAggr[7332]

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$order <- 1:nrow(grid_to_compute)

lat_eval <- round(profLatAggr[73372]-.5) + .5 # grid point 44788, job 90, 288
long_eval <- round(profLongAggr[73372]-.5) + .5

grid_to_compute[grid_to_compute$long == long_eval & grid_to_compute$lat == lat_eval,]

load('analysis/results/psal_one_stage_cov_pca.RData')
coefs_pc_psal <- coefs_pc_psal[[which(grid_pc_psal$long == long_eval & grid_pc_psal$lat == lat_eval)]][[1]][,1:K2]
load('analysis/results/temp_one_stage_cov_pca.RData')
coefs_pc_temp <- coefs_pc_temp[[which(grid_pc_temp$long == long_eval & grid_pc_temp$lat == lat_eval)]][[1]][,1:K1]
p <- 0:2000
pc_vals_temp <- sapply(1:K1, function(x) predict(make_fit_object(knots_pc, coefs_pc_temp[,x]), p)$y)
pc_vals_psal <- sapply(1:K2, function(x) predict(make_fit_object(knots_pc, coefs_pc_psal[,x]), p)$y)

pred_values_temp <- rowSums(pc_vals_temp %*% diag(single_results[[1]][[1]]))
pred_values_psal <- rowSums(pc_vals_psal %*% diag(single_results[[2]][[1]]))

var_values_temp <- diag(as.matrix(pc_vals_temp %*% cond_var[1:K1, 1:K1] %*% t(pc_vals_temp)))
var_values_psal <- diag(as.matrix(pc_vals_psal %*% cond_var[(K1+1):(K1 + K2), (K1+1):(K1 + K2)] %*% t(pc_vals_psal)))

load('analysis/results/nugget_variance/nugget_var_40000.RData')
nugget_single <- Grid_Pred40000[[267]]
var_pred_temp_out <- exp(predict(nugget_single[[1]], 0:2000))/(1/2*exp(digamma(1)))
var_pred_psal_out <- exp(predict(nugget_single[[2]], 0:2000))/(1/2*exp(digamma(1)))

lower_bound_temp <- pred_values_temp -2*sqrt(var_values_temp+var_pred_temp_out)
upper_bound_temp <- pred_values_temp +2*sqrt(var_values_temp+var_pred_temp_out)

lower_bound_psal <- pred_values_psal -2*sqrt(var_values_psal+var_pred_psal_out)
upper_bound_psal <- pred_values_psal +2*sqrt(var_values_psal+var_pred_psal_out)

quantile_all <- .95

quantile_all <- pnorm(q = 2) - pnorm(q = -2)
quantile <- 1 - (1-quantile_all)/2

r_x = sqrt(single_results[[1]][[3]]* rowSums(sapply(1:K1, function(x) {
  single_results[[1]][[4]][x]^2 * pc_vals_temp[,x]^2})) +
    abs(qnorm(p = (1-quantile)/nrow(residuals)/2 ) )^2 *var_pred_temp_out)

r_x_psal = sqrt(single_results[[2]][[3]]* rowSums(sapply(1:K2, function(x) {
  single_results[[2]][[4]][x]^2 * pc_vals_psal[,x]^2})) +
    abs(qnorm(p = (1-quantile)/nrow(residuals)/2 ) )^2 * var_pred_psal_out)

png('analysis/images/cv/cv_example_temp.png',  width = 400/72, height = 600/72,
   units = 'in', res = 144)
min_val <- min(-r_x)
max_val <- max(r_x)
plot(residuals$temp_mean_residual,residuals$pressure, xlim = c(min_val- .5, max_val),
     ylab = 'Pressure', xlab = 'Temperature Residual (°C)', cex.lab = 1.5,cex.axis = 1.5,
     ylim=c(2000, 0))
lines(y=0:2000,pred_values_temp + r_x,lty = 2)
lines(y=0:2000,pred_values_temp - r_x, lty = 2)
lines(y=0:2000,lower_bound_temp, col = 2)
lines(y=0:2000,upper_bound_temp, col = 2)
lines(y=0:2000, pred_values_temp, cex = .5)
legend('bottomleft',c('Observed', 'Prediction', '95.4% point\nprediction\nbound', '95.4%\nprediction\nband'),
       lty = c(NA, 1, 1, 2), col = c(1, 1,2,1), pch = c(1, NA, NA, NA),
       pt.cex = c(1, NA, NA, NA), cex = 1.2, y.intersp = 2)
dev.off()

png('analysis/images/cv/cv_example_psal.png', width = 400/72, height = 600/72,
    units = 'in', res = 144)
min_val <- min(-r_x_psal)
max_val <- max(r_x_psal)
plot(residuals$psal_mean_residual,residuals$pressure, xlim = c(min_val, max_val),
     ylab = 'Pressure', xlab = 'Salinity Residual (PSU)', cex.lab = 1.5,cex.axis = 1.5,
     ylim=c(2000, 0))
lines(y=0:2000,pred_values_psal + r_x_psal,lty = 2)
lines(y=0:2000,pred_values_psal - r_x_psal, lty = 2)
lines(y=0:2000,lower_bound_psal, col = 2)
lines(y=0:2000,upper_bound_psal, col = 2)
lines(y=0:2000, pred_values_psal, cex = .5)
# legend('bottomleft',c('Observed', 'Prediction', '95.4% point prediction bound', '95.4% point prediction band'),
#        lty = c(NA, 1, 1, 2), col = c(1, 1,2,1), pch = c(1, NA, NA, NA),
#        pt.cex = c(.5, NA, NA, NA), cex = 1.2)
dev.off()


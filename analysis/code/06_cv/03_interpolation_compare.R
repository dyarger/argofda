
library(dplyr)
library(argofda)
library(Matrix)

calc_files <- list.files('analysis/results/cv/', pattern = 'cv_resu')
load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

load('analysis/products/KS_cv_floats_resid_final.RData')
profile_nums_10 <- (1:length(profFloatIDAggr))[1:length(profFloatIDAggr) %in% data_10$profile_num[data_10$month == 2]]
profile_nums_300 <- (1:length(profFloatIDAggr))[1:length(profFloatIDAggr) %in% data_300$profile_num[data_300$month == 2]]
profile_nums_1500 <- (1:length(profFloatIDAggr))[1:length(profFloatIDAggr) %in% data_1500$profile_num[data_1500$month == 2]]

KS_profiles <- unique(c(profile_nums_10, profile_nums_300, profile_nums_1500))

rg_levels <- c(2.5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 182.5, 200,
               220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 462.5, 500, 550, 600, 650, 700, 750,
               800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1412.5, 1500, 1600, 1700, 1800, 1900, 1975)
final_list <- list()

for (j in 1:length(calc_files)) {
  load(paste0('analysis/results/cv/', calc_files[j]))
  n_delayed <- sapply(cv_results, function(x) {
    if (length(x) == 7) {
      if (is.null(x[[6]])) {
        return(0)
      }
      return(x[[6]])
    } else {
      return(0)
    }
  })
  no_pred <- sapply(cv_results, function(x) is.null(x[[1]]))
  mode <- sapply(cv_results, function(x) x[[5]])
  residuals <- lapply(1:length(cv_results), function(x) data.frame(cv_results[[x]][[3]], mode = mode[x]))
  temp_results <- lapply(cv_results, function(x) x[[1]])
  psal_results <- lapply(cv_results, function(x) x[[2]])
  profile_nums <- sapply(cv_results, function(x) x[[7]][1])
  info <- lapply(cv_results, function(x) x[[7]][2])

  cond_var_old_scores <- lapply(cv_results, function(x) x[[4]])
  final_list[[j]] <- list(mode, residuals,
                          temp_results, psal_results,
                          cond_var_old_scores,
                          profile_nums,
                          info,n_delayed,
                          no_pred)
  print(j)
}

mode <- unlist(sapply(final_list, function(x) x[[1]]))
temp_results <- unlist(lapply(final_list, function(x) x[[3]]), recursive = F)
psal_results <- unlist(lapply(final_list, function(x) x[[4]]), recursive = F)
profile_nums <- sapply(final_list, function(x) x[[6]])
profile_nums <- as.numeric(unlist(profile_nums, recursive = T))
profile_nums_use <- profile_nums[which(profile_nums %in% KS_profiles)]
temp_results_use <- temp_results[which(profile_nums %in% KS_profiles)]
psal_results_use <- psal_results[which(profile_nums %in% KS_profiles)]


latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)
load('analysis/results/psal_cov_pca.RData')
load('analysis/results/temp_cov_pca.RData')

temp_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
psal_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
temp_mean_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
psal_mean_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
temp_representation_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
temp_rep_error_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
psal_representation_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
psal_rep_error_summary <- matrix(ncol = length(rg_levels), nrow = length(profile_nums_use))
K1 <- 10
K2 <- 10
for (q in 1:length(profile_nums_use)) {
  lat <- profLatAggr[profile_nums_use[q]]
  long <- profLongAggr[profile_nums_use[q]]
  lat_eval <- latvals[which.min(abs(lat - latvals))]
  long_eval <- longvals[which.min(abs(long - longvals))]

  if (!(lat_eval %in% grid_pc_temp$lat[grid_pc_temp$long == long_eval] &
        lat_eval %in% grid_pc_psal$lat[grid_pc_psal$long == long_eval]  )) {
    temp_summary[q,] <- rep(NA, 3)
    psal_summary[q,] <- rep(NA, 3)
    temp_mean_summary[q,] <- rep(NA, 3)
    psal_mean_summary[q,] <- rep(NA, 3)
    next
  }
  temp_pc <- coefs_pc_temp[[which(grid_pc_temp$long == long_eval &
                                    grid_pc_temp$lat == lat_eval)]][[1]][,1:K1]
  psal_pc <- coefs_pc_psal[[which(grid_pc_psal$long == long_eval &
                                    grid_pc_psal$lat == lat_eval)]][[1]][,1:K2]

  pc_values_temp <- sapply(as.data.frame(as.matrix(temp_pc)), function(x) {
    predict(make_fit_object(knots = knots_pc,
                            coefficients = x), rg_levels)$y
  })
  pc_values_psal <- sapply(as.data.frame(as.matrix(psal_pc)), function(x) {
    predict(make_fit_object(knots = knots_pc,
                            coefficients = x), rg_levels)$y
  })

  temp_pred_scores <- as.double(temp_results_use[[q]][[1]])
  psal_pred_scores <- as.double(psal_results_use[[q]][[1]])

  scores_out_temp <- as.double(temp_results_use[[q]][[2]])
  scores_out_psal <- as.double(psal_results_use[[q]][[2]])

  pc_representation_temp <- rowSums(sapply(1:K1, function(x) { scores_out_temp[x] * pc_values_temp[,x]}))
  pc_representation_psal <- rowSums(sapply(1:K2, function(x) { scores_out_psal[x] * pc_values_psal[,x]}))

  pred_values_temp <- rowSums(sapply(1:K1, function(x) { temp_pred_scores[x] * pc_values_temp[,x]}))
  pred_values_psal <- rowSums(sapply(1:K2, function(x) { psal_pred_scores[x] * pc_values_psal[,x]}))

  temp_interp <- approx(x = profPresAggr[[profile_nums_use[q]]][[1]], y = profTempAggr[[profile_nums_use[q]]][[1]], xout = rg_levels)$y
  psal_interp <- approx(x = profPresAggr[[profile_nums_use[q]]][[1]], y = profPsalAggr[[profile_nums_use[q]]][[1]], xout = rg_levels)$y

  temp_summary[q,] <- temp_interp - pred_values_temp
  psal_summary[q,] <- psal_interp - pred_values_psal

  temp_representation_summary[q,] <- pc_representation_temp - pred_values_temp
  temp_rep_error_summary[q,] <- temp_interp - pc_representation_temp

  psal_representation_summary[q,] <- pc_representation_psal - pred_values_psal
  psal_rep_error_summary[q,] <- psal_interp - pc_representation_psal

  temp_mean_summary[q,] <- temp_interp
  psal_mean_summary[q,] <- psal_interp
  if (q %% 1000 == 0) {print(q)}
}
save(temp_summary, psal_summary,
     temp_mean_summary, psal_mean_summary,
     profile_nums_use,temp_representation_summary,temp_rep_error_summary,
     psal_representation_summary,psal_rep_error_summary,
     file = 'analysis/results/compare_joint_TS_20.RData')
load('analysis/results/compare_joint_TS_20.RData')
apply(temp_summary, 2, function(x) sum(!is.na(x)))
apply(temp_summary, 2, function(x) sum(!is.na(x)))[rg_levels %in% c(10, 300, 1500)]

apply(temp_mean_summary, 2, function(x) median(abs(x), na.rm = T))[rg_levels %in% c(10, 300, 1500)] # 0.31280451 0.33065208 0.05690722
apply(psal_mean_summary, 2, function(x) median(abs(x), na.rm = T))[rg_levels %in% c(10, 300, 1500)] # 0.089872655 0.041585264 0.008162833
apply(temp_mean_summary, 2, function(x) quantile(abs(x), na.rm = T, probs = .75))[rg_levels %in% c(10, 300, 1500)] # 0.6068204 0.6634698 0.1120474
apply(psal_mean_summary, 2, function(x) quantile(abs(x), na.rm = T, probs = .75))[rg_levels %in% c(10, 300, 1500)] # 0.17876025 0.08694514 0.01805677
apply(temp_mean_summary, 2, function(x) sqrt(mean(x^2, na.rm = T)))[rg_levels %in% c(10, 300, 1500)] # 0.7266018 0.8414801 0.1461088
apply(psal_mean_summary, 2, function(x) sqrt(mean(x^2, na.rm = T)))[rg_levels %in% c(10, 300, 1500)] # 0.2233626 0.1143802 0.0307123

apply(temp_summary, 2, function(x) median(abs(x), na.rm = T))[rg_levels %in% c(10, 300, 1500)] # 0.20034225 0.17738107 0.03814318
apply(psal_summary, 2, function(x) median(abs(x), na.rm = T))[rg_levels %in% c(10, 300, 1500)] # 0.059170526 0.027154927 0.006357115
apply(temp_summary, 2, function(x) quantile(abs(x), na.rm = T, probs = .75))[rg_levels %in% c(10, 300, 1500)] # 0.40162190 0.37126605 0.07636077
apply(psal_summary, 2, function(x) quantile(abs(x), na.rm = T, probs = .75))[rg_levels %in% c(10, 300, 1500)] # 0.13222699 0.06023890 0.01395079
apply(temp_summary, 2, function(x) sqrt(mean(x^2, na.rm = T)))[rg_levels %in% c(10, 300, 1500)] # 0.5234630 0.5176045 0.1095574
apply(psal_summary, 2, function(x) sqrt(mean(x^2, na.rm = T)))[rg_levels %in% c(10, 300, 1500)] # 0.19541214 0.10880223 0.02930242

# Median Absolute residual
cex <- 1
png('analysis/images/cv/temp_interp_comparison_median.png',
   width =800, height = 1200, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = apply(temp_mean_summary, 2, function(x) median(abs(x), na.rm = T)), cex = cex - .5,ylim = c(2000, 0),
     cex.lab = 1.8,cex.axis = 1.8, ylab = 'Pressure (decibars)', xlab = 'Median Absolute Residual' , xlim = c(0, max(apply(temp_mean_summary, 2, function(x) median(abs(x), na.rm = T)))))
points(y=rg_levels,x=apply(temp_summary, 2, function(x) median(abs(x), na.rm = T)), pch = 2, cex =cex - .5)
points(y=c(10, 300, 1500), x=c(.1801, .1740, .0311),pch = 3, cex = cex, col = 2)
points(y=c(10, 300, 1500), x=c(.4750, .3062, .0530),pch = 4, cex = cex, col = 3)
points(y=c(10, 300, 1500), x=c(.2556, .1991, .0356),pch = 5, cex = cex, col = 4)
legend('bottomright', c('Fun. Mean', 'Fun. Model', 'KS', 'RG mean', 'RG ref'), pt.cex = c(cex - .5, cex - .5, cex,cex,cex),
       col = c(1,1,2,3,4), cex = 1.8, pch = 1:5)
dev.off()

salinity_median_mean_summary <- apply(psal_mean_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                                 function(x) median(abs(x), na.rm = T))
salinity_median_summary <- apply(psal_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                                 function(x) median(abs(x), na.rm = T))
png('analysis/images/cv/psal_interp_comparison_median.png',
    width =800, height = 1200, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y= rg_levels, x = salinity_median_mean_summary, cex = cex - .5,ylim = c(2000, 0),
     cex.lab = 1.5,cex.axis = 1.5, ylab = 'Pressure (decibars)', xlab = 'Median Absolute Residual' , xlim = c(0, max(salinity_median_mean_summary)))
points(y=rg_levels,x=salinity_median_summary, pch = 2, cex =cex - .5)
legend('bottomright', c('Fun. Mean', 'Fun. Model'), pt.cex = rep(cex - .5, 2), col = 1,
       cex = 1.8,
       pch = 1:2)
dev.off()

# Q3
png('analysis/images/cv/temp_interp_comparison_Q3.png',
    width =800, height = 1200, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = apply(temp_mean_summary, 2, function(x) quantile(abs(x),probs = .75, na.rm = T)), cex = cex - .5,ylim = c(2000,0),
     cex.lab = 1.8,cex.axis = 1.8, ylab = 'Pressure (decibars)', xlab = 'Q3 of Absolute Residuals' , xlim = c(0, max(apply(temp_mean_summary, 2, function(x) quantile(abs(x),probs = .75, na.rm = T)))))
points(rg_levels,apply(temp_summary, 2, function(x) quantile(abs(x),probs = .75, na.rm = T)), pch = 2, cex =cex - .5)
points(y=c(10, 300, 1500), x=c(.3735, .3684, .0641),pch = 3, cex = cex, col = 2)
points(y=c(10, 300, 1500),x=c(.8670, .6320, .1043),pch = 4, cex = cex, col = 3)
points(y=c(10, 300, 1500), x=c(.5026, .4213, .0736),pch = 5, cex = cex, col = 4)
legend('bottomright', c('Fun. Mean', 'Fun. Model', 'KS', 'RG mean', 'RG ref'), pt.cex = c(cex - .5, cex -.5, cex, cex, cex),
       col = 1, cex = 1.8, pch = 1:5)
dev.off()

salinity_mean_summary <- apply(psal_mean_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                                      function(x) quantile(abs(x), probs = .75, na.rm = T))
salinity_summary <- apply(psal_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                                 function(x) quantile(abs(x), probs = .75, na.rm = T))
png('analysis/images/cv/psal_interp_comparison_Q3.png',
    width =800, height = 1200, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = salinity_mean_summary, cex = cex - .5,ylim = c(2000, 0),
     cex.lab = 1.8,cex.axis = 1.8, ylab = 'Pressure (decibars)', xlab = 'Q3 of Absolute Residuals' , xlim = c(0, max(salinity_mean_summary)))
points(y=rg_levels,x=salinity_summary, pch = 2, cex =cex - .5)
legend('bottomright', c('Fun. Mean', 'Fun. Model'), pt.cex = rep(cex - .5, 2), col = 1,
       cex = 1.8,
       pch = 1:2)
dev.off()

# RMSE
png('analysis/images/cv/temp_interp_comparison_RMSE.png',
    width =800, height = 1200, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = apply(temp_mean_summary, 2, function(x) sqrt(mean(x^2, na.rm = T))), cex = cex - .5,ylim = c(2000,0),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'RMSE Residuals' , xlim = c(0, max(apply(temp_mean_summary, 2, function(x) sqrt(mean(x^2, na.rm = T))))))
points(y=rg_levels,x=apply(temp_summary, 2, function(x) sqrt(mean(x^2, na.rm = T))), pch = 2, cex =cex - .5)
points(y=c(10, 300, 1500), c(.5072, .5124, .0883),pch = 3, cex = cex, col = 2)
points(y=c(10, 300, 1500), c(.8889, .8149, .1337),pch = 4, cex = cex, col = 3)
points(y=c(10, 300, 1500), c(.6135, .5782, .1014),pch = 5, cex = cex, col = 4)
legend('bottomright', c('Fun. Mean', 'Fun. Model', 'KS', 'RG mean', 'RG ref'), pt.cex = c(cex - .5, cex - .5, cex, cex, cex),
       col = c(1,1,2,3,4),
       cex = 1.8,
       pch = 1:5)
dev.off()


salinity_mean_summary <- apply(psal_mean_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                               function(x) sqrt(mean(x^2, na.rm = T)))
salinity_summary <- apply(psal_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                          function(x) sqrt(mean(x^2, na.rm = T)))
png('analysis/images/cv/psal_interp_comparison_RMSE.png',
    width =800, height = 1200, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = salinity_mean_summary, cex = cex - .5,ylim = c(2000, 0),
     cex.lab = 1.8,cex.axis = 1.8, ylab = 'Pressure (decibars)', xlab = 'RMSE Residuals' , xlim = c(0, max(salinity_mean_summary)))
points(y=rg_levels,salinity_summary, pch = 2, cex =cex -.5)
legend('bottomright', c('Fun. Mean', 'Fun. Model'), pt.cex = rep(cex - .5, 2), pch = 1:2,
       cex = 1.8,
       col = rep(1, 2))
dev.off()

## pc representation error

png('analysis/images/cv/temp_score_rep.png',
    height = 1200, width = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = apply(temp_mean_summary, 2, function(x) median(abs(x), na.rm = T)), cex = cex - .5,ylim = c(2000, 0),
     cex.lab = 1.8,cex.axis = 1.8, ylab = 'Pressure (decibars)', xlab = 'Median Absolute Residual' ,
     xlim = c(0, max(apply(temp_mean_summary, 2, function(x) median(abs(x), na.rm = T)))))
points(y=rg_levels,x=apply(temp_summary, 2, function(x) median(abs(x), na.rm = T)), col = 2, cex =cex - .5)
points(y=rg_levels,x=apply(temp_representation_summary, 2, function(x) median(abs(x), na.rm = T)), col = 3, cex =cex - .5)
points(y=rg_levels,x=apply(temp_rep_error_summary, 2, function(x) median(abs(x), na.rm = T)), col = 4, cex =cex - .5)
legend('bottomright', c('Mean', 'Pred', 'Score Rep', 'Score Error'), pt.cex = c(cex - .5, cex - .5, cex -.5, cex - .5),
       col = c(1:4),
       cex = 1.6,
       pch = c(1,1,1,1))
dev.off()

salinity_mean_summary <- apply(psal_mean_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                               function(x) median(abs(x),  na.rm = T))
salinity_summary <- apply(psal_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                          function(x) median(abs(x),  na.rm = T))
salinity_respresentation_summary <- apply(psal_representation_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                               function(x) median(abs(x),na.rm = T))
salinity_rep_error_summary <- apply(psal_rep_error_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2,
                          function(x) median(abs(x), na.rm = T))

png('analysis/images/cv/psal_score_rep.png',
    width = 1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(y = rg_levels, x = salinity_mean_summary, cex = cex - .5,ylim = c(2000, 0),
     cex.lab = 1.8,cex.axis = 1.8, ylab = 'Pressure (decibars)', xlab = 'Median Absolute Residual' ,
     xlim = c(0, max(salinity_mean_summary)))
points(y=rg_levels,salinity_summary, col = 2, cex =cex - .5)
points(y=rg_levels,salinity_respresentation_summary, col = 3, cex =cex - .5)
points(y=rg_levels,salinity_rep_error_summary, col = 4, cex =cex - .5)
legend('bottomright', c('Mean', 'Pred', 'Score Rep', 'Score Error'), pt.cex = rep(cex - .5, 4),
       col = c(1:4),
       cex = 1.6,
       pch = c(1,1,1,1))
dev.off()








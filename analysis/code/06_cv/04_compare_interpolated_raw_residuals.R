
# comparison of interpolation and raw residuals
library(dplyr)

calc_files <- list.files('analysis/results/cv/', pattern = 'cv_resu')
load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

load('analysis/products/KS_cv_floats_resid_final.RData')
profile_nums_10 <- (1:length(profFloatIDAggr))[1:length(profFloatIDAggr) %in% data_10$profile_num[data_10$month == 2]]
profile_nums_300 <- (1:length(profFloatIDAggr))[1:length(profFloatIDAggr) %in% data_300$profile_num[data_300$month == 2]]
profile_nums_1500 <- (1:length(profFloatIDAggr))[1:length(profFloatIDAggr) %in% data_1500$profile_num[data_1500$month == 2]]

KS_profiles <- unique(c(profile_nums_10, profile_nums_300, profile_nums_1500))
length(KS_profiles)

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
  profile_nums <- sapply(cv_results, function(x) x[[7]][1])
  residuals <- lapply(1:length(cv_results), function(x) data.frame(cv_results[[x]][[3]], mode = mode[x],
                                                                   profile_num =profile_nums[x]  ))
  temp_results <- lapply(cv_results, function(x) x[[1]])
  psal_results <- lapply(cv_results, function(x) x[[2]])
  info <- lapply(cv_results, function(x) x[[7]][2])

  cond_var_old_scores <- lapply(cv_results, function(x) x[[4]])
  final_list[[j]] <- list(mode[profile_nums %in%KS_profiles], residuals[profile_nums %in%KS_profiles],
                          temp_results[profile_nums %in%KS_profiles], psal_results[profile_nums %in%KS_profiles],
                          cond_var_old_scores[profile_nums %in%KS_profiles],
                          profile_nums[profile_nums %in%KS_profiles],
                          info[profile_nums %in%KS_profiles],n_delayed[profile_nums %in%KS_profiles],
                          no_pred[profile_nums %in%KS_profiles])
  print(j)
}
residuals <- lapply(final_list, function(x) x[[2]])
residuals <- unlist(residuals, recursive = F)
residuals <- bind_rows(residuals)



load('analysis/results/compare_joint_TS_20.RData')

midpoints_RG <- (rg_levels[1:(length(rg_levels)-1)] + rg_levels[2:length(rg_levels)])/2
residual_summary <- residuals %>%
  select(-mode, - profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) sqrt(mean(z^2, na.rm = T)))) %>%
  mutate(pgroup = as.numeric(as.character(pgroup)))
residual_summary_delayed <- residuals %>%
  select(-profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) sqrt(mean(z^2, na.rm = T)))) %>%
  filter(mode == 'D') %>%
  ungroup() %>%
  mutate(pgroup = as.numeric(as.character(pgroup)))

cex <- .8
png('analysis/images/cv/temp_interpolation_comparison.png',
    width = 1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(residual_summary$pgroup, residual_summary$temp_mean_residual, cex = cex,
     xlab = 'Pressure', ylab = 'RMSE of Residuals', ylim = c(0, 1.3),
     cex.lab = 1.8,cex.axis = 1.8)
points(rg_levels, apply(temp_mean_summary, 2, function(x) sqrt(mean(x^2, na.rm = T))),
       col = 2, cex = cex, pch = 20)
points(residual_summary$pgroup, residual_summary$temp_residual,
       col = 3, cex = cex, pch = 3)
points(rg_levels, apply(temp_summary, 2, function(x)sqrt(mean(x^2, na.rm = T))),
       col = 4, cex = cex, pch = 4)
legend('topright', c('Raw Mean Values', 'Interp Mean Values', 'Raw Pred Values ', 'Interp Pred Values'),
       pt.cex = rep(cex, 5),
       col = c(1:4),
       cex = 1.8,
       pch = c(1,20,3,4),
       pt.bg = c(NA, 2, NA, NA))
dev.off()

png('analysis/images/cv/psal_interpolation_comparison.png',
    width = 1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(residual_summary_delayed$pgroup, residual_summary_delayed$psal_mean_residual, cex = cex,
     xlab = 'Pressure', ylab = 'RMSE of Residuals', ylim = c(0, .23),
     cex.lab = 1.8,cex.axis = 1.8)
points(rg_levels, apply(psal_mean_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2, function(x) sqrt(mean(x^2, na.rm = T))),
       col = 2, cex = cex, pch = 20)
points(residual_summary_delayed$pgroup, residual_summary_delayed$psal_residual,
       col = 3, cex = cex, pch = 3)
points(rg_levels, apply(psal_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2, function(x)sqrt(mean(x^2, na.rm = T))),
       col = 4, cex = cex, pch = 4)
legend('topright', c('Raw Mean Values', 'Interp Mean Values', 'Raw Pred Values ', 'Interp Pred Values'), pt.cex = rep(cex, 5),
       col = c(1:4),
       cex = 1.8,
       pch = c(1,20,3,4),
       pt.bg = c(NA, 2, NA, NA))
dev.off()

##########

cex <- 1
# RMSE
png('analysis/images/cv/temp_pred_error_comparison_RMSE.png',
    width =1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(x = residual_summary$pgroup, y = residual_summary$temp_mean_residual, cex = cex - .5,xlim = c(0,2000),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'RMSE Residuals' , ylim = c(0, max(residual_summary$temp_mean_residual)))
points(residual_summary$pgroup,residual_summary$temp_residual, pch = 2, cex =cex - .5)
points(c(10, 300, 1500), c(.5072, .5124, .0883),pch = 3, cex = cex, col = 2)
points(c(10, 300, 1500), c(.8889, .8149, .1337),pch = 4, cex = cex, col = 3)
points(c(10, 300, 1500), c(.6135, .5782, .1014),pch = 5, cex = cex, col = 4)
legend('topright', c('Fun. Mean', 'Fun. Model', 'KS', 'RG mean', 'RG ref'), pt.cex = c(cex - .5, cex - .5, cex, cex, cex),
       col = c(1,1,2,3,4),
       cex = 1.8,
       pch = 1:5)
dev.off()
residual_summary$temp_mean_residual[residual_summary$pgroup %in% c(10, 300, 1500)]
residual_summary$temp_residual[residual_summary$pgroup %in% c(10, 300, 1500)]

png('analysis/images/cv/psal_pred_error_comparison_RMSE.png',
    width =1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(x = residual_summary_delayed$pgroup, y = residual_summary_delayed$psal_mean_residual, cex = cex - .5,xlim = c(0, 2000),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'RMSE Residuals' ,
     ylim = c(0, max(residual_summary_delayed$psal_mean_residual)))
points(residual_summary_delayed$pgroup,residual_summary_delayed$psal_residual, pch = 2, cex =cex -.5)
legend('topright', c('Fun. Mean', 'Fun. Model'), pt.cex = rep(cex - .5, 2), pch = 1:2,
       cex = 1.8,
       col = rep(1, 2))
dev.off()

# median
residual_summary <- residuals %>%
  select(-mode, - profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) median(abs(z), na.rm = T))) %>%
  mutate(pgroup = as.numeric(as.character(pgroup)))
residual_summary_delayed <- residuals %>%
  select(-profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) median(abs(z), na.rm = T))) %>%
  filter(mode == 'D') %>%
  ungroup() %>%
  mutate(pgroup = as.numeric(as.character(pgroup)))


png('analysis/images/cv/temp_pred_error_comparison_median.png',
    width =1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(x = residual_summary$pgroup, y = residual_summary$temp_mean_residual, cex = cex - .5,xlim = c(0,2000),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'Median Absolute Residual' , ylim = c(0, max(residual_summary$temp_mean_residual)))
points(residual_summary$pgroup,residual_summary$temp_residual, pch = 2, cex =cex - .5)
points(c(10, 300, 1500), c(.1801, .1740, .0311),pch = 3, cex = cex, col = 2)
points(c(10, 300, 1500), c(.4750, .3062, .0530),pch = 4, cex = cex, col = 3)
points(c(10, 300, 1500), c(.2556, .1991, .0356),pch = 5, cex = cex, col = 4)
legend('topright', c('Fun. Mean', 'Fun. Model', 'KS', 'RG mean', 'RG ref'), pt.cex = c(cex - .5, cex - .5, cex, cex, cex),
       col = c(1,1,2,3,4),
       cex = 1.8,
       pch = 1:5)
dev.off()
residual_summary$temp_mean_residual[residual_summary$pgroup %in% c(10, 300, 1500)] # 0.31928370 0.34936718 0.06200045
residual_summary$temp_residual[residual_summary$pgroup %in% c(10, 300, 1500)] # 0.19777220 0.17268008 0.03760884

png('analysis/images/cv/psal_pred_error_comparison_median.png',
    width =1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(x = residual_summary_delayed$pgroup, y = residual_summary_delayed$psal_mean_residual, cex = cex - .5,xlim = c(0, 2000),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'Median Absolute Residual' ,
     ylim = c(0, max(residual_summary_delayed$psal_mean_residual)))
points(residual_summary_delayed$pgroup,residual_summary_delayed$psal_residual, pch = 2, cex =cex -.5)
legend('topright', c('Fun. Mean', 'Fun. Model'), pt.cex = rep(cex - .5, 2), pch = 1:2,
       cex = 1.8,
       col = rep(1, 2))
dev.off()


# Q3
residual_summary <- residuals %>%
  select(-mode, - profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) quantile(abs(z), probs = .75, na.rm = T))) %>%
  mutate(pgroup = as.numeric(as.character(pgroup)))
residual_summary_delayed <- residuals %>%
  select(-profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) quantile(abs(z), probs = .75, na.rm = T))) %>%
  filter(mode == 'D') %>%
  ungroup() %>%
  mutate(pgroup = as.numeric(as.character(pgroup)))


png('analysis/images/cv/temp_pred_error_comparison_Q3.png',
    width =1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(x = residual_summary$pgroup, y = residual_summary$temp_mean_residual, cex = cex - .5,xlim = c(0,2000),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'Q3 of Absolute Residuals' , ylim = c(0, max(residual_summary$temp_mean_residual)))
points(residual_summary$pgroup,residual_summary$temp_residual, pch = 2, cex =cex - .5)
points(c(10, 300, 1500), c(.3735, .3684, .0641),pch = 3, cex = cex, col = 2)
points(c(10, 300, 1500), c(.8670, .6320, .1043),pch = 4, cex = cex, col = 3)
points(c(10, 300, 1500), c(.5026, .4213, .0736),pch = 5, cex = cex, col = 4)
legend('topright', c('Fun. Mean', 'Fun. Model', 'KS', 'RG mean', 'RG ref'), pt.cex = c(cex - .5, cex - .5, cex, cex, cex),
       col = c(1,1,2,3,4),
       cex = 1.8,
       pch = 1:5)
dev.off()
residual_summary$temp_mean_residual[residual_summary$pgroup %in% c(10, 300, 1500)] #  0.6247356 0.6845067 0.1160059
residual_summary$temp_residual[residual_summary$pgroup %in% c(10, 300, 1500)] # 0.39716038 0.36660264 0.07236247

png('analysis/images/cv/psal_pred_error_comparison_Q3.png',
    width =1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(x = residual_summary_delayed$pgroup, y = residual_summary_delayed$psal_mean_residual, cex = cex - .5,xlim = c(0, 2000),
     cex.lab = 1.8,cex.axis = 1.8, xlab = 'Pressure (decibars)', ylab = 'Q3 of Absolute Residuals' ,
     ylim = c(0, max(residual_summary_delayed$psal_mean_residual)))
points(residual_summary_delayed$pgroup,residual_summary_delayed$psal_residual, pch = 2, cex =cex -.5)
legend('topright', c('Fun. Mean', 'Fun. Model'), pt.cex = rep(cex - .5, 2), pch = 1:2,
       cex = 1.8,
       col = rep(1, 2))
dev.off()





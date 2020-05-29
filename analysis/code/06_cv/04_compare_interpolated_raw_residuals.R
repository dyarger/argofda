
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
  summarise_all(list(function(z) sqrt(mean(z^2, na.rm = T))))
residual_summary_delayed <- residuals %>%
  select(-profile_num) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) sqrt(mean(z^2, na.rm = T)))) %>%
  filter(mode == 'D')

cex <- .8
png('analysis/images/cv/temp_interpolation_comparison.png',
    width = 1200, height = 800, res = 144)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
plot(as.numeric(as.character(residual_summary$pgroup)), residual_summary$temp_mean_residual, cex = cex,
     xlab = 'Pressure', ylab = 'RMSE of Residuals', ylim = c(0, 1.3),
     cex.lab = 1.8,cex.axis = 1.8)
points(rg_levels, apply(temp_mean_summary, 2, function(x) sqrt(mean(x^2, na.rm = T))),
       col = 2, cex = cex, pch = 20)
points(as.numeric(as.character(residual_summary$pgroup)), residual_summary$temp_residual,
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
plot(as.numeric(as.character(residual_summary_delayed$pgroup)), residual_summary_delayed$psal_mean_residual, cex = cex,
     xlab = 'Pressure', ylab = 'RMSE of Residuals', ylim = c(0, .23),
     cex.lab = 1.8,cex.axis = 1.8)
points(rg_levels, apply(psal_mean_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2, function(x) sqrt(mean(x^2, na.rm = T))),
       col = 2, cex = cex, pch = 20)
points(as.numeric(as.character(residual_summary_delayed$pgroup)), residual_summary_delayed$psal_residual,
       col = 3, cex = cex, pch = 3)
points(rg_levels, apply(psal_summary[profile_nums_use %in% which(profModeAggr == 'D'),], 2, function(x)sqrt(mean(x^2, na.rm = T))),
       col = 4, cex = cex, pch = 4)
legend('topright', c('Raw Mean Values', 'Interp Mean Values', 'Raw Pred Values ', 'Interp Pred Values'), pt.cex = rep(cex, 5),
       col = c(1:4),
       cex = 1.8,
       pch = c(1,20,3,4),
       pt.bg = c(NA, 2, NA, NA))
dev.off()

residual_summary <- residuals %>%
  select(pressure, temp_residual) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) median(abs(z), na.rm = T)))

as.data.frame(residual_summary)

residual_summary <- residuals %>%
  select(pressure, temp_residual) %>%
  mutate(pgroup = base::cut(x = pressure, breaks = c(0,midpoints_RG, Inf), labels = rg_levels, include.lowest = T)) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) quantile(abs(z), probs = .75, na.rm = T)))

as.data.frame(residual_summary[, c('pressure', 'temp_residual')])





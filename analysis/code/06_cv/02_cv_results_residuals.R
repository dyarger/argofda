
library(dplyr)
library(ggplot2)
calc_files <- list.files('analysis/results/cv/', pattern = 'cv_resu')
load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

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
  profile_nums <- as.numeric(sapply(cv_results, function(x) x[[7]][1]))
  residuals <- lapply(1:length(cv_results), function(x) data.frame(cv_results[[x]][[3]], mode = mode[x],
                                                                   profile_num = profile_nums[x]))
  temp_results <- lapply(cv_results, function(x) x[[1]])
  psal_results <- lapply(cv_results, function(x) x[[2]])
  info <- lapply(cv_results, function(x) x[[7]][2])

  cond_var_old_scores <- lapply(cv_results, function(x) x[[4]])
  final_list[[j]] <- list(mode[!no_pred], residuals[!no_pred],
                          temp_results[!no_pred], psal_results[!no_pred],
                          cond_var_old_scores[!no_pred],
                          profile_nums[!no_pred],
                          info[!no_pred],n_delayed[!no_pred])
  print(j)
}
residuals <- lapply(final_list, function(x) x[[2]])
residuals <- unlist(residuals, recursive = F)
residuals <- bind_rows(residuals)

# summarize median absolute residual by pressure
residual_summary <- residuals %>%
  select(-mode) %>%
  mutate(pgroup = round(pressure/10)*10) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) median(abs(z), na.rm = T)))
residual_summary_delayed <- residuals %>%
  mutate(pgroup = round(pressure/10)*10) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) median(abs(z), na.rm = T))) %>%
  filter(mode == 'D')

png('analysis/images/cv/one_stage_compare_temp.png',
    width =1200, height = 800, res = 144)
plot(residual_summary$pgroup, residual_summary$temp_mean_residual, cex = .5,  ylim = c(0, 0.6),
     xlab = 'Pressure (decibars)', ylab = 'Median Absolute Residual (°C)')
points(residual_summary$pgroup, residual_summary$temp_residual, cex = .5, col = 2)
legend('topright', c('Mean', 'Model'), col = 1:2, pt.cex = .5, pch = c(1,1))
dev.off()

png('analysis/images/cv/one_stage_compare_psal.png',
    width =1200, height = 800, res = 144)
plot(residual_summary_delayed$pgroup, residual_summary_delayed$psal_mean_residual, cex = .5,  ylim = c(0, 0.10),
     xlab = 'Pressure (decibars)', ylab = 'Median Absolute Residual (PSU)')
points(residual_summary_delayed$pgroup, residual_summary_delayed$psal_residual, cex = .5, col = 2)
legend('topright', c('Mean', 'Model'), col = 1:2, pt.cex = .5, pch = c(1,1))
dev.off()

# summarize mean residual by pressure
residual_summary_bias <- residuals %>%
  select(-mode) %>%
  mutate(pgroup = round(pressure/10)*10) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) mean(z, na.rm = T)))
residual_summary_bias_delayed <- residuals %>%
  mutate(pgroup = round(pressure/10)*10) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) mean(z, na.rm = T))) %>%
  ungroup() %>% filter(mode == 'D')
png('analysis/images/cv/one_stage_compare_bias_temp.png', width =1200, height = 800, res = 144)
plot(residual_summary_bias$pgroup, residual_summary_bias$temp_mean_residual, cex = .5,  ylim = c(-.05, .05),
     ylab = 'Mean Residual (°C)', xlab = 'Pressure')
points(residual_summary_bias$pgroup, residual_summary_bias$temp_residual, cex = .5, col = 2)
abline(h= 0 )
legend('topright', c('Mean', 'Model'), col = 1:2, pt.cex = .5, pch = c(1,1))
dev.off()

png('analysis/images/cv/one_stage_compare_bias_psal.png', width =1200, height = 800, res = 144)
plot(residual_summary_bias_delayed$pgroup, residual_summary_bias_delayed$psal_mean_residual, cex = .5,  ylim = c(-.005, .005),
     ylab = 'Mean Residual (PSU)', xlab = 'Pressure')
points(residual_summary_bias_delayed$pgroup, residual_summary_bias_delayed$psal_residual, cex = .5, col = 2)
abline(h= 0 )
legend('topright', c('Mean', 'Model'), col = 1:2, pt.cex = .5, pch = c(1,1))
dev.off()

# summarize residual variance by pressure
residual_summary_var <- residuals %>%
  select(-mode) %>%
  mutate(pgroup = round(pressure/10)*10) %>%
  group_by(pgroup) %>%
  summarise_all(list(function(z) var(z, na.rm = T)))
residual_summary_var_delayed <- residuals %>%
  mutate(pgroup = round(pressure/10)*10) %>%
  group_by(pgroup, mode) %>%
  summarise_all(list(function(z) var(z, na.rm = T))) %>%
  ungroup() %>% filter(mode == 'D')
png('analysis/images/cv/one_stage_compare_empirical_var_temp.png', width =1200, height = 800, res = 144)
plot(residual_summary_var$pgroup, residual_summary_var$temp_mean_residual, cex = .5,  ylim = c(0, 1.6),
     ylab = 'Variance of Residuals (°C^2)', xlab = 'Pressure')
points(residual_summary_var$pgroup, residual_summary_var$temp_residual, cex = .5, col = 2)
points(residual_summary_bias$pgroup, residual_summary_bias$var_values_temp +  residual_summary_bias$nugget_var_temp,
       cex = .5, col = 4)
points(residual_summary_bias$pgroup, residual_summary_bias$var_values_temp_uc +  residual_summary_bias$nugget_var_temp,
       cex = .5, col = 3)
legend('topright', c('Mean', 'Model', 'Avg UC Model Variance',
                     'Avg C Model Variance'), col = 1:4, pt.cex = .5, pch = c(1,1,1,1))
dev.off()

png('analysis/images/cv/one_stage_compare_empirical_var_psal.png', width =1200, height = 800, res = 144)
plot(residual_summary_var_delayed$pgroup, residual_summary_var_delayed$psal_mean_residual, cex = .5,  ylim = c(0, .05),
     ylab = 'Variance of Residuals (PSU^2)', xlab = 'Pressure')
points(residual_summary_var_delayed$pgroup, residual_summary_var_delayed$psal_residual, cex = .5, col = 2)
points(residual_summary_bias$pgroup, residual_summary_bias_delayed$var_values_psal +  residual_summary_bias_delayed$nugget_var_psal,
       cex = .5, col = 4)
points(residual_summary_bias$pgroup, residual_summary_bias_delayed$var_values_psal_uc +  residual_summary_bias_delayed$nugget_var_psal,
       cex = .5, col = 3)
legend('topright', c('Mean', 'Model', 'Avg UC Model Variance',
                     'Avg C Model Variance'), col = 1:4, pt.cex = .5, pch = c(1,1,1,1))
dev.off()

# summarize coverages by pressure
# some profiles did not have neought profiles nearby, so we don't have good variance estimates for them
lower_bound_temp <- residuals$temp_mean_residual - residuals$temp_residual - 2*sqrt(residuals$var_values_temp+residuals$nugget_var_temp)
upper_bound_temp <- residuals$temp_mean_residual - residuals$temp_residual + 2*sqrt(residuals$var_values_temp+residuals$nugget_var_temp)

lower_bound_psal <- residuals$psal_mean_residual - residuals$psal_residual- 2*sqrt(residuals$var_values_psal + residuals$nugget_var_psal)
upper_bound_psal <- residuals$psal_mean_residual - residuals$psal_residual+ 2*sqrt(residuals$var_values_psal + residuals$nugget_var_psal)

residuals$in_interval_temp <- ifelse(residuals$temp_mean_residual <= upper_bound_temp &
                                       residuals$temp_mean_residual >= lower_bound_temp,
                                     1, 0)
residuals$in_interval_psal <- ifelse(residuals$psal_mean_residual <= upper_bound_psal &
                                       residuals$psal_mean_residual >= lower_bound_psal,
                                     1, 0)
mean(residuals$in_interval_temp)
mean(residuals$in_interval_psal[residuals$mode == 'D'])

residual_summary_coverage <- residuals %>%
  select(-mode) %>%
  mutate(pgroup = round(pressure/20)*20) %>%
  group_by(pgroup) %>%
  summarise(in_interval_temp = mean(in_interval_temp, na.rm = T))
residual_summary_coverage_delayed <- residuals %>%
  mutate(pgroup = round(pressure/20)*20) %>%
  group_by(pgroup, mode) %>%
  summarise(in_interval_psal = mean(in_interval_psal, na.rm = T)) %>%
  ungroup() %>% filter(mode == 'D')
png('analysis/images/cv/one_stage_coverage_both.png', width =900, height = 600, res = 144)
par(mar=c(5.1-1, 4.1, 4.1-3, 2.1-1))
plot(residual_summary_coverage$pgroup, residual_summary_coverage$in_interval_temp, cex = .7, # ylim = c(-.05, .05),
     ylab = 'Pointwise Coverage', xlab = 'Pressure (dbar)',
     ylim = c(.8, 1))
abline(h=  pnorm(q = 2) - pnorm(q = -2))
points(residual_summary_coverage_delayed$pgroup, residual_summary_coverage_delayed$in_interval_psal, cex = .7,  #ylim = c(-.005, .005),
       ylab = 'Pointwise Coverage, Salinity', xlab = 'Pressure',ylim = c(.8, 1), pch = 2)
legend('bottomright', c('Temperature', 'Salinity'), pch = c(1,2), pt.cex = c(.7, .7))
dev.off()

## QQ plots
prob_vals <- seq(.01, .99, by = .005)
residual_summary_quantiles <- residuals %>%
  select(-mode) %>%
  mutate(pgroup = round(pressure/50)*50,
                norm_temp_resid = temp_residual/sqrt(var_values_temp + nugget_var_temp)) %>%
  group_by(pgroup) %>%
  summarise(x=list(tibble::enframe(quantile(norm_temp_resid, probs=prob_vals,
                                                        na.rm = T), "quantiles", "resid"))) %>%
                tidyr::unnest(x) %>%
  mutate(quantiles = readr::parse_number(quantiles))
normal_quantiles <- data.frame(quantiles = prob_vals*100,
                               norm_values = qnorm(p = prob_vals))
residual_summary_quantiles <- merge(residual_summary_quantiles, normal_quantiles)
residual_summary_quantiles$resid_diff <- residual_summary_quantiles$resid - residual_summary_quantiles$norm_values
ggplot(data = residual_summary_quantiles,
       aes(x = norm_values, y = resid_diff, group = pgroup, color = pgroup))+
  geom_line()+
  coord_cartesian(ylim = c(-1.2, 1.2))+
  scale_color_viridis_c()+theme_bw() +
  labs(x = 'Theoretical Normal Quantile', y = 'Empirical Quantile - Theoretical Normal Quantile',
       color = 'Pressure\nGroup')
ggsave(file = 'analysis/images/cv/qq_temp.png',
       scale = .8,height = 6, width = 7.25)


# QQ salinity
residual_summary_quantiles_delayed <- residuals %>%
  mutate(pgroup = round(pressure/50)*50,
         norm_temp_resid = psal_residual/sqrt(var_values_psal + nugget_var_psal)) %>%
  group_by(pgroup, mode) %>%
  summarise(x=list(tibble::enframe(quantile(norm_temp_resid, probs=prob_vals,
                                            na.rm = T), "quantiles", "resid"))) %>%
  tidyr::unnest(x) %>%
  mutate(quantiles = readr::parse_number(quantiles)) %>%
  ungroup() %>% filter(mode == 'D')
normal_quantiles <- data.frame(quantiles = prob_vals*100,
                               norm_values = qnorm(p = prob_vals))
residual_summary_quantiles_delayed <- merge(residual_summary_quantiles_delayed, normal_quantiles)
residual_summary_quantiles_delayed$resid_diff <- residual_summary_quantiles_delayed$resid - residual_summary_quantiles_delayed$norm_values
ggplot(data = residual_summary_quantiles_delayed,
       aes(x = norm_values, y = resid_diff, group = pgroup, color = pgroup))+
  geom_line()+
  coord_cartesian(ylim = c(-1.2, 1.2))+
  scale_color_viridis_c()+theme_bw() +
  labs(x = 'Theoretical Normal Quantile', y = 'Empirical Quantile - Theoretical Normal Quantile',
       color = 'Pressure\nGroup')
ggsave(file = 'analysis/images/cv/qq_psal.png',
       scale = .8,height = 6, width = 7.25)

### space comparison 0-20 dbar
residual_loc_summary <- residuals %>%
  filter(pressure <20) %>%
  left_join(data.frame(profile_num = as.numeric(1:length(profFloatIDAggr)),
                       long = profLongAggr,
                       lat = profLatAggr)) %>%
  select(-mode) %>%
  mutate(long_group = findInterval(x = long, vec =seq(-180,180, by = 3)),
         lat_group = findInterval(x = lat, vec = seq(-80,80, by = 3)),
         long = round(long/3)*3,
         lat = round(lat/3)*3) %>%
  select(-profile_num) %>%
  group_by(long_group, lat_group) %>%
  summarise(temp_residual = sqrt(mean(temp_residual^2, na.rm = T)) ,
            long = mean(long, na.rm = T), lat = mean(lat, na.rm = T))
residual_loc_summary$long_test <- round(residual_loc_summary$long/3)*3
residual_loc_summary$lat_test <- round(residual_loc_summary$lat/3)*3
ggplot()+
  geom_raster(data = residual_loc_summary,
              aes(x = ifelse(long_test < 0, long_test+ 360, long_test),
                                               y = lat_test, fill = temp_residual))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  scale_fill_gradientn(limits = c(0, 6), values = scales::rescale(c(0,.3,1,2, 6)),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  theme_bw()+
  labs(fill = 'RMSE (°C)', x= 'Longitude', y = 'Latitude')+
  theme(panel.background = element_blank() )
ggsave(file = 'analysis/images/cv/location_summary_0_20.png',
       scale = .8,height = 4, width = 7.25)

### space comparison 300-500 dbar
residual_loc_summary <- residuals %>%
  filter(pressure <500 & pressure > 300) %>%
  left_join(data.frame(profile_num = as.numeric(1:length(profFloatIDAggr)),
                       long = profLongAggr,
                       lat = profLatAggr)) %>%
  select(-mode) %>%
  mutate(long_group = findInterval(x = long, vec =seq(-180,180, by = 3)),
         lat_group = findInterval(x = lat, vec = seq(-80,80, by = 3)),
         long = round(long/3)*3,
         lat = round(lat/3)*3) %>%
  select(-profile_num) %>%
  group_by(long_group, lat_group) %>%
  summarise(temp_residual = sqrt(mean(temp_residual^2, na.rm = T)) ,
            long = mean(long, na.rm = T), lat = mean(lat, na.rm = T))
residual_loc_summary$long_test <- round(residual_loc_summary$long/3)*3
residual_loc_summary$lat_test <- round(residual_loc_summary$lat/3)*3
ggplot()+
  geom_raster(data = residual_loc_summary,
              aes(x = ifelse(long_test < 0, long_test+ 360, long_test),
                  y = lat_test, fill = temp_residual))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  scale_fill_gradientn(limits = c(0, 8.7), values = scales::rescale(c(0,.3,1,2, 6)),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  theme_bw()+
  labs(fill = 'RMSE (°C)', x= 'Longitude', y = 'Latitude')+
  theme(panel.background = element_blank() )
ggsave(file = 'analysis/images/cv/location_summary_300_500.png',
       scale = .8,height = 4, width = 7.25)


### space comparison 1500-2000 dbar
residual_loc_summary <- residuals %>%
  filter(pressure >1500) %>%
  left_join(data.frame(profile_num = as.numeric(1:length(profFloatIDAggr)),
                       long = profLongAggr,
                       lat = profLatAggr)) %>%
  select(-mode) %>%
  mutate(long_group = findInterval(x = long, vec =seq(-180,180, by = 3)),
         lat_group = findInterval(x = lat, vec = seq(-80,80, by = 3)),
         long = round(long/3)*3,
         lat = round(lat/3)*3) %>%
  select(-profile_num) %>%
  group_by(long_group, lat_group) %>%
  summarise(temp_residual = sqrt(mean(temp_residual^2, na.rm = T)) ,
            long = mean(long, na.rm = T), lat = mean(lat, na.rm = T))
residual_loc_summary$long_test <- round(residual_loc_summary$long/3)*3
residual_loc_summary$lat_test <- round(residual_loc_summary$lat/3)*3

ggplot()+
  geom_raster(data = residual_loc_summary,
              aes(x = ifelse(long_test < 0, long_test+ 360, long_test),
                  y = lat_test, fill = temp_residual))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  scale_fill_gradientn(limits = c(0, 1.75), values = scales::rescale(c(0,.3,1,2, 6)),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  theme_bw()+
  labs(fill = 'RMSE (°C)', x= 'Longitude', y = 'Latitude')+
  theme(panel.background = element_blank())
ggsave(file = 'analysis/images/cv/location_summary_1500_2000.png',
       scale = .8,height = 4, width = 7.25)

## overall coverages
temp_results <- unlist(sapply(final_list, function(x) x[[3]]), recursive = F)
temp_in_band <- sapply(temp_results, function(x) all(x[[5]]))
temp_in_ci <- t(sapply(temp_results, function(x) c(sum(x[[6]]), length(x[[6]]))))
temp_prop <- apply(temp_in_ci,2, sum, na.rm = T)
temp_prop[1]/temp_prop[2]
mean(temp_in_band)

# Salinity coverages - delayed mode only
delayed <- unlist(sapply(final_list, function(x) x[[1]]), recursive = F)
psal_results <- unlist(sapply(final_list, function(x) x[[4]]), recursive = F)
psal_results <- psal_results[delayed == 'D']
psal_in_band <- sapply(psal_results, function(x) all(x[[5]]))
psal_in_ci <- t(sapply(psal_results, function(x) c(sum(x[[6]]), length(x[[6]]))))
psal_prop <- apply(psal_in_ci, 2, sum, na.rm = T)
psal_prop[1]/psal_prop[2]
table(psal_in_band)[2]/length(psal_in_band)



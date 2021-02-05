library(ggplot2)
library(dplyr)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))

load('analysis/results/joint_TS_20/mean_anomalies_10_300_1500.RData')
mean_anomalies <- mean_anomalies %>%
  mutate(full_pred_10_T = pred_10_T + mean_10_T,
         full_pred_300_T = pred_300_T + mean_300_T,
         full_pred_1500_T = pred_1500_T + mean_1500_T,
         full_pred_10_S = pred_10_S + mean_10_S,
         full_pred_300_S = pred_300_S + mean_300_S,
         full_pred_1500_S = pred_1500_S + mean_1500_S)
mean_anomalies <- mean_anomalies[mean_anomalies$year == 2012,]

load('analysis/results/compare_mean.RData')

# 10 decibars
anomaly_kuusela <- R.matlab::readMat('analysis/products/anomalySpaceTimeExp_10_02_2007_2016_2012.mat')
pred_KS <- data.frame(long = rep(20:380, times = 181),
                      lat = rep(-90:90, each = 361),
                      pred = as.vector(anomaly_kuusela[[3]]))
compare_single_df <- compare_df[compare_df$pressure == 10,]
pred_KS$long <- ifelse(pred_KS$long > 360, pred_KS$long - 360, pred_KS$long)
pred_KS$long <- pred_KS$long-.5
pred_KS$lat <- pred_KS$lat-.5
RG_KS <- full_join(pred_KS, compare_single_df)
RG_KS$full_pred <- RG_KS$pred + RG_KS$temperature_RG
mean_anomalies$long <- ifelse(mean_anomalies$long < 0, mean_anomalies$long + 360, mean_anomalies$long)
compare_anom <- inner_join(RG_KS, mean_anomalies)
ggplot()+
  geom_raster(data = compare_anom, aes(x = long, y = lat, fill = full_pred_10_T- full_pred))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  scale_fill_gradientn(limits = c(-2.7,2.7), colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-2.5, -1, 0,1, 2.5)))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Fun-KS (°C)')
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_10_compare_KS.png', scale = .8,height = 4, width = 7.25,dpi = 150)
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_10_compare_KS_large.png', scale = .8,height = 4, width = 7.25)
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_10_compare_KS.eps', scale = .8,height = 4, width = 7.25)

# 300 decibars
anomaly_kuusela <- R.matlab::readMat('analysis/products/anomalySpaceTimeExp_300_02_2007_2016_2012.mat')
pred_KS <- data.frame(long = rep(20:380, times = 181),
                      lat = rep(-90:90, each = 361),
                      pred = as.vector(anomaly_kuusela[[3]]))
compare_single_df <- compare_df[compare_df$pressure == 300,]
pred_KS$long <- ifelse(pred_KS$long > 360, pred_KS$long - 360, pred_KS$long)
pred_KS$long <- pred_KS$long-.5
pred_KS$lat <- pred_KS$lat-.5
RG_KS <- full_join(pred_KS, compare_single_df)
RG_KS$full_pred <- RG_KS$pred + RG_KS$temperature_RG
mean_anomalies$long <- ifelse(mean_anomalies$long < 0, mean_anomalies$long + 360, mean_anomalies$long)
compare_anom <- inner_join(RG_KS, mean_anomalies)

ggplot()+
  geom_raster(data = compare_anom, aes(x = long, y = lat, fill = full_pred_300_T- full_pred))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  scale_fill_gradientn(limits = c(-3,3), colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-2.5, -1, 0,1, 2.5)))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Fun-KS (°C)')
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_300_compare_KS.png', scale = .8,height = 4, width = 7.25,dpi = 150)
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_300_compare_KS_large.png', scale = .8,height = 4, width = 7.25)
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_300_compare_KS.eps', scale = .8,height = 4, width = 7.25)


# 1500  decibars
anomaly_kuusela <- R.matlab::readMat('analysis/products/anomalySpaceTimeExp_1500_02_2007_2016_2012.mat')
pred_KS <- data.frame(long = rep(20:380, times = 181),
                              lat = rep(-90:90, each = 361),
                              pred = as.vector(anomaly_kuusela[[3]]))
compare_single_df <- compare_df[compare_df$pressure == 1500,]
pred_KS$long <- ifelse(pred_KS$long > 360, pred_KS$long - 360, pred_KS$long)
pred_KS$long <- pred_KS$long-.5
pred_KS$lat <- pred_KS$lat-.5
RG_KS <- full_join(pred_KS, compare_single_df)
RG_KS$full_pred <- RG_KS$pred + RG_KS$temperature_RG
mean_anomalies$long <- ifelse(mean_anomalies$long < 0, mean_anomalies$long + 360, mean_anomalies$long)
compare_anom <- inner_join(RG_KS, mean_anomalies)

ggplot()+
  geom_raster(data = compare_anom, aes(x = long, y = lat, fill = full_pred_1500_T- full_pred))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  scale_fill_gradientn(limits = c(-1,1), colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-1, -.4, 0, .4, 1)))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Fun-KS (°C)')
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_1500_compare_KS.png', scale = .8,height = 4, width = 7.25,dpi = 150)
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_1500_compare_KS_large.png', scale = .8,height = 4, width = 7.25)
ggsave(filename = 'analysis/images/joint_TS_20/anomaly_1500_compare_KS.eps', scale = .8,height = 4, width = 7.25)

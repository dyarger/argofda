

# MLD results
files <- list.files('analysis/results/ohc/', pattern = 'cross_no_nugget')
ohc_list <- list()
for (i in 1:length(files)){
  load(paste0('analysis/results/ohc/', files[i]))
  ohc_list[[i]] <- ohc
}
library(ggplot2)
library(dplyr)
# 394386 49664  3.5 57.5 2007
ohc <- do.call(rbind, ohc_list)
summary(nchar(ohc$ohc_anomaly))
ohc$ohc_anomaly[nchar(ohc$ohc_anomaly) > 20]
head()
ohc[6485,]
which(nchar(ohc$ohc_anomaly) > 20)
colnames(ohc) <-   c('index', 'long', 'lat', 'year',
  'ohc_mean','ohc_pred', 'ohc_anomaly', 'ohc_avg_anomaly',
                     'ohc_mean_red','ohc_pred_red', 'ohc_anomaly_red', 'ohc_avg_anomaly_red',
                     'var_ohc', 'var_ohc_red', 'var_ohc_cross')
ohc$var_diff <- ohc$var_ohc + ohc$var_ohc_red - 2 * ohc$var_ohc_cross
ohc$var_diff_percent <- ohc$var_diff/ohc$var_ohc

ohc2 <- ohc[!is.na(ohc$var_diff_percent),]
ohc2[ohc2$var_diff_percent > 1,] %>%
  dplyr::select(starts_with('var'))
ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat, fill = ohc_pred_red))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(7 * 10^11, .85*10^12))

ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat, fill = ohc_pred))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(7 * 10^11, .85*10^12))


ggplot(data = ohc[ohc$year == 2016,], aes(x = long, y = lat, fill = ohc_avg_anomaly))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-7 * 10^9, 7*10^9))
ggplot(data = ohc[ohc$year == 2012,], aes(x = long, y = lat, fill = var_diff > 0))+
  geom_raster()
ggplot(data = ohc[ohc$year == 2011,], aes(x = long, y = lat,
                                          fill = ifelse(var_diff > 0, var_diff, 0)))+
    geom_raster()+labs(fill = 'diff')+
  scale_fill_viridis_c(limits = c(0, 8*10^16))
ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat, fill = var_ohc_cross))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, 10^20))
ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat, fill = var_ohc))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, 10^20))
ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat, fill = var_ohc_red))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, 10^20))

ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat,
                                          fill = ifelse(var_diff_percent > 0, var_diff_percent, 0)))+
  geom_raster()+labs(fill = 'diff')+
  scale_fill_viridis_c(limits = c(0, .02))
ggplot(data = ohc[ohc$year == 2015,], aes(x = long, y = lat, fill = ohc_pred - ohc_pred_red))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-.001*10^10, .001*10^10)) +
  labs(fill = 'Bias Int')
load('analysis/data/RG_Defined_mask.RData')
ggplot(data = ohc[ohc$year == 2015,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = long, y = lat,
                                          fill = ifelse(var_diff >0,
                                                        (ohc_pred - ohc_pred_red)^2 + var_diff,
                                                        (ohc_pred - ohc_pred_red)^2)))+
  geom_raster()+
  scale_fill_gradient2(limits = c(0*10^10, 5*10^16)) +
  labs(fill = 'MSE')

ggplot(data = ohc[ohc$year == 2015,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = long, y = lat,
                                   fill = ifelse(var_diff >0,
                                                 (ohc_pred - ohc_pred_red)^2 + var_diff,
                                                 (ohc_pred - ohc_pred_red)^2)/var_ohc))+
  geom_raster()+
  scale_fill_gradient2(limits = c(0*10^10, .02)) +
  labs(fill = 'MSE')

ggplot(data = ohc[ohc$year == 2015,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = long, y = lat,
                                   fill = ifelse(var_diff >0,
                                                 (ohc_pred - ohc_pred_red)^2 ,
                                                 (ohc_pred - ohc_pred_red)^2)))+
  geom_raster()+
  scale_fill_gradient2(limits = c(0*10^10, 1*10^12)) +
  labs(fill = 'bias2')

ggplot(data = ohc[ohc$year == 2015,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = long, y = lat,
                                   fill = ifelse(var_diff >0,
                                                 (ohc_pred - ohc_pred_red)^2,
                                                 (ohc_pred - ohc_pred_red)^2)/var_ohc))+
  geom_raster()+
  scale_fill_gradient2(limits = c(0*10^10, .00001)) +
  labs(fill = 'Bias2')

ggplot(data = ohc[ohc$year == 2015,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = long, y = lat,
                                   fill = var_ohc_red- var_ohc))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-10^15, 10^16)) +
  labs(fill = 'variance difference')
year <- 2016
h <- 6.5
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))
ggplot(data = ohc[ohc$year == year,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = ifelse(long < 0, long + 360, long), y = lat,
                                   fill = (var_ohc_red- var_ohc)/var_ohc_red))+
  geom_raster()+
  scale_fill_gradient2(low = 'blue',limits = c(-.004, .004), # limits = c(-.004, .004)
                       mid = 'white', high = 'red',
                       guide = guide_colorbar(barwidth = 12)) +
  geom_polygon(data = map_data('world2'), aes(x = long , y = lat, group = group),
                color = 'black', fill = 'white', size = .2)+
  labs(fill =  expression(frac(sigma[coarse]^2 - sigma[fine]^2, sigma[fine]^2)),
     #  expression(frac(sigma[coarse]^2 - sigma[fine]^2, sigma[fine]^2))
       x = 'Longitude', y = 'Latitude')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_var_diff_', year, '.png'),
       scale = .8,height = h, width = 7.25)

ggplot(data = ohc[ohc$year == year,] %>% inner_join(RG_defined_long) %>%
         filter(value > 1999), aes(x = ifelse(long < 0, long + 360, long), y = lat,
                                   fill = (sqrt(var_ohc_red)- sqrt(var_ohc))/sqrt(var_ohc)))+
  geom_raster()+
  scale_fill_gradient2(low = 'blue',limits = c(-.004, .004), # limits = c(-.004, .004)
                       mid = 'white', high = 'red',
                       guide = guide_colorbar(barwidth = 12)) +
  geom_polygon(data = map_data('world2'), aes(x = long , y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(fill = expression(frac(sigma[coarse] - sigma[fine], sigma[fine])),
       #  expression(frac(sigma[coarse]^2 - sigma[fine]^2, sigma[fine]^2))
       x = 'Longitude', y = 'Latitude')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_var_diff_', year, '_sd2.png'),
       scale = .8,height = h, width = 7.25)


summary(ohc$var_ohc_cross - ohc$var_ohc_red
        )
summary(ohc$var_ohc_cross - ohc$var_ohc)
library(ncdf4)
library(raster)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))
# Adjust for area of 1 degree by 1 degree grid
r <- raster()
ar <- area(r)
m <- as.matrix(ar) * 1000^2 # area of each grid point in m^2
latitudes <- seq(-89.5, 89.5, by = 1)
adjust_area <- data.frame(lat = latitudes, area = m[,1])

ohc <- merge(ohc, adjust_area)
ohc$ohc_anomaly_adjust <- ohc$ohc_anomaly* ohc$area
ohc$ohc_sd_adjust <- ohc$ohc_sd* ohc$area
ohc$ohc_avg_anomaly_adjust <- ohc$ohc_avg_anomaly* ohc$area
ohc$long_p <- ifelse(ohc$long < 0, ohc$long+360, ohc$long)
# plot parameters
h <- 6.5
rescaling <- 10^(21)
x_min <- -14*10^(19)
x_max <- 13*10^(19)
scale <- c(x_min, -3*10^(19), 0, 3*10^(19), x_max)/rescaling

map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
                           fill = 'white', color = 'black',size=  .2)
year <- 2016
colors_plot <-  c(scales::muted("blue"), "white", scales::muted("red"))
colors_plot <-  c("blue", "white", "red")

ggplot(data = ohc[ohc$year == year,], aes(x = long_p, y  =lat,
                                          fill = ohc_anomaly_adjust/rescaling))+
  geom_raster()+
  map_plot+
  scale_fill_gradientn(limits = c(scale[1], scale[length(scale)]), values = scales::rescale(scale),
                       colours =colors_plot,
                       breaks = c(-.13, 0, .13))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'OHC Pred (ZJ)')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_pred_', year, '_large.png'),
       scale = .8,height = h, width = 7.25)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_pred_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_pred_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

ggplot(data = ohc[ohc$year == year,], aes(x = long_p, y  =lat,
                                          fill = ohc_avg_anomaly_adjust/rescaling))+
  geom_raster()+
  map_plot+
  scale_fill_gradientn(limits = c(scale[1], scale[length(scale)]), values = scales::rescale(scale),
                       colours = colors_plot,
                       breaks = c(-.13, 0, .13))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'OHC Anomaly (ZJ)')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_anom_', year, '_large.png'),
       scale = .8,height = h, width = 7.25)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_anom_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_anom_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

# Compare with NOAA estimates available at https://www.nodc.noaa.gov/OC5/3M_HEAT_CONTENT/
ohc_noaa <- nc_open('analysis/products/heat_content_anomaly_0-700_seasonal.nc')
ohc_noaa_vals <- ncvar_get(ohc_noaa, 'h18_hc')
time <- ohc_noaa$dim$time$vals
ohc_noaa_single <- ohc_noaa_vals[,,which(1955 + time/12 == year + .125)] # Jan-Mar year

# Adjust for area of 1 degree by 1 degree grid
r <- raster()
ar <- area(r)
m <- as.matrix(ar) * 1000^2 # area of each grid point in m^2
latitude <- ohc_noaa$dim$lat$vals
longitude <- ohc_noaa$dim$lon$vals
noaa_df <- data.frame(expand.grid('long' = longitude, 'lat' = latitude), value = as.vector(ohc_noaa_single))
noaa_df$long_p <- ifelse(noaa_df$long < 0, noaa_df$long+360, noaa_df$long)

rescaling_noaa <- rescaling^(-1)

ggplot(data = noaa_df, aes(x = long_p, y  =lat, fill = value * 10^18/10^21))+
  geom_raster()+
  map_plot +
  scale_fill_gradientn(limits = c(scale[1], scale[length(scale)]), values = scales::rescale(scale),
                       colours = colors_plot,
                       breaks = c(-.13, 0, .13))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'OHC Anomaly (ZJ)')
ggsave(filename = paste0('analysis/images/ohc/noaa_ohc_700_', year, '_large.png'),
       scale = .8,height = h, width = 7.25)
ggsave(filename = paste0('analysis/images/ohc/noaa_ohc_700_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(filename = paste0('analysis/images/ohc/noaa_ohc_cons_temp_700_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

rescaling <- 10^(21)
x_min <- -14*10^(19)
x_max <- 13*10^(19)
scale <- c(x_min, -3*10^(19), 0, 3*10^(19), x_max)/rescaling

ggplot(data = ohc[ohc$year == year,], aes(x = long_p,
                                          y  =lat, fill = ohc_sd_adjust/rescaling))+
  geom_raster()+
  scale_fill_gradientn(colours = colorRamps::matlab.like(10), limits = c(0, .06),
                       breaks = c(0, .025, .05)) +
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'OHC SD (ZJ)')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_sd_', year, '_large.png'),
       scale = .8,height = h, width = 7.25)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_sd_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_sd_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

ggplot(data = ohc[ohc$year == year & !is.na(ohc$ohc_sd),],
       aes(x = long_p, y  =lat,fill =
             factor(ifelse(ohc_anomaly > 2 * ohc_sd, 1,
                           ifelse(ohc_anomaly < - 2 * ohc_sd, -1, 0)))))+
  geom_raster()+
  scale_fill_manual(values = c('blue', 'white', 'darkred'))+
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'OHC within 2 SD of 0')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_pred_sig_', year, '_large.png'),
       scale = .8,height = h, width = 7.25)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_pred_sig_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_pred_sig_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

ggplot(data = ohc[ohc$year == year & !is.na(ohc$ohc_sd),],
       aes(x = long_p, y  =lat,fill =
             factor(ifelse(ohc_avg_anomaly > 2 * ohc_sd, 1,
                           ifelse(ohc_avg_anomaly < - 2 * ohc_sd, -1, 0)))))+
  geom_raster()+
  scale_fill_manual(values = c('blue', 'white', 'darkred'))+
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'OHC within 2 SD of 0')
ggsave(filename = paste0('analysis/images/ohc/ohc_700_anom_sig_', year, '_large.png'),
       scale = .8,height = h, width = 7.25)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_anom_sig_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(filename = paste0('analysis/images/ohc/ohc_700_anom_sig_', year, '.eps'),
       scale = .8,height = h, width = 7.25)



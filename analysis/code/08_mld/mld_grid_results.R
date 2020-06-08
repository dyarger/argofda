# MLD results
get_mld <- list()
for (j in 1:116) {
  load(paste0('analysis/results/mld/mld_',j ,'.RData'))
  summaries <- lapply(final_list[[2]], function(x) {
    apply(x, 2, quantile, na.rm = T, probs = c(.025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .975))})
  mean_summaries <- lapply(final_list[[2]], function(x) {
    apply(x, 2, mean, na.rm = T)})
  sd_summaries <- lapply(final_list[[2]], function(x) {
    apply(x, 2, sd, na.rm = T)})
  get_mld[[j]] <- list(final_list[[1]], summaries,mean_summaries,sd_summaries)
  print(j)
}

save(get_mld, file = 'analysis/results/mld/mld_grid_summaries.RData')
q()

library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(ncdf4)
library(dplyr)
load('analysis/results/mld/mld_grid_summaries.RData')
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))

locations <- do.call(rbind, lapply(get_mld, function(x) x[[1]]))
mld_quantiles <- unlist(lapply(get_mld, function(x) x[[2]]), recursive = F)
mld_var_dens_quan <- do.call(rbind, lapply(mld_quantiles, function(x) t(x)[1,]))
colnames(mld_var_dens_quan) <- paste0('quant_', c(.025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .975))
mld_var_dens_mean <- unlist(lapply(get_mld, function(x) x[[3]]), recursive = F)
mld_var_dens_mean <- do.call(c, lapply(mld_var_dens_mean, function(x) x[1]))
mld_var_dens_sd <- unlist(lapply(get_mld, function(x) x[[4]]), recursive = F)
mld_var_dens_sd <- do.call(c, lapply(mld_var_dens_sd, function(x) x[1]))
mld_df <- data.frame(locations, mld_var_dens_quan, mean = mld_var_dens_mean, sd = mld_var_dens_sd)
mld_df$long_p <- ifelse(mld_df$long < 0, mld_df$long + 360, mld_df$long)

mld_avg <- mld_df %>%
  group_by(long, lat) %>%
  summarise_all(.funs = mean)

scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)

a <- ggplot(data = mld_avg, aes(x = long_p, y = lat, fill = mean))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_var_dens_mean_mean.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_var_dens_mean_mean.eps',
       scale = .8,height = h, width = 7.25)

### Median and quantiles
a <- ggplot(data = mld_avg,
            aes(x = long_p, y = lat, fill = quant_0.5))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_median_mean.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_median_mean.eps'),
       scale = .8,height = h, width = 7.25)
a <- ggplot(data = mld_avg,
            aes(x = long_p, y = lat, fill = quant_0.8))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_80_mean.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_80_mean.eps'),
       scale = .8,height = h, width = 7.25)
a <- ggplot(data = mld_avg,
            aes(x = long_p, y = lat, fill = quant_0.7))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_70_mean.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_70_mean.eps'),
       scale = .8,height = h, width = 7.25)

mld_avg$quant_len_0.9 <- mld_avg$quant_0.95 - mld_avg$quant_0.05
scale <- c(seq(0, 100, by  =20),120, 140, 160, 180, 200, 300, 400, 500, 1000, 1500)
limits <- c(0, 2000)
a <- ggplot(data = mld_avg,
            aes(x = long_p, y = lat, fill = quant_len_0.9))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')+
  theme(legend.text = element_text(size = 8))
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_90_range_mean.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_90_range_mean.eps'),
       scale = .8,height = h, width = 7.25)

##### Specifically for each year
year <- 2014
scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5

a <- ggplot(data =mld_df[mld_df$year == year,],
            aes(x = long_p, y = lat, fill = mean))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_mean_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_mean_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

a <- ggplot(data =mld_df[mld_df$year == year,],
            aes(x = long_p, y = lat, fill = quant_0.5))+
  geom_raster() +
  map_plot+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_median_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_median_', year, '.eps'),
       scale = .8,height = h, width = 7.25)

mld_df$quant_len_0.9 <- mld_df$quant_0.95 - mld_df$quant_0.05
scale <- c(seq(0, 100, by  =20),120, 140, 160, 180, 200, 300, 400, 500, 1000, 1500)
limits <- c(0, 2000)
a <- ggplot(data = mld_df[mld_df$year == mld_df$year,],
            aes(x = long_p, y = lat, fill = quant_len_0.9))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')+
  theme(legend.text = element_text(size = 8))
a
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_90_range_', year, '.png'),
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = paste0('analysis/images/mld/mld_var_dens_90_range_', year, '.eps'),
       scale = .8,height = h, width = 7.25)


######## Compare with ucsd estimates from http://mixedlayer.ucsd.edu/
mld <- nc_open('analysis/products/Argo_mixedlayers_monthlyclim_05092018.nc')
mld_da <- ncvar_get(mld,'mld_da_mean')[2,,] # February density algorithm
mld_dt <- ncvar_get(mld,'mld_dt_mean')[2,,] # February density threshold
mld_da_std <- ncvar_get(mld,'mld_da_std')[2,,]
mld_dt_std <- ncvar_get(mld,'mld_dt_std')[2,,]
long <- ncvar_get(mld,'lon')
lat <- ncvar_get(mld,'lat')

mld_da_df <- cbind(expand.grid('long' = long, 'lat' = lat), mld_var = as.vector(mld_da))
mld_dt_df <- cbind(expand.grid('long' = long, 'lat' = lat), mld_var = as.vector(mld_dt))
mld_dt_std <- cbind(expand.grid('long' = long, 'lat' = lat), mld_var = as.vector(mld_dt_std))
mld_da_std <- cbind(expand.grid('long' = long, 'lat' = lat), mld_var = as.vector(mld_da_std))

scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5

a <- ggplot(data = mld_da_df, aes(x = ifelse(long < 0, long + 360, long), y = lat,  fill = mld_var))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits)+
  geom_polygon(data = map_data('world2'), aes(x  = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_algo.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_algo.eps',
       scale = .8,height = h, width = 7.25)

a <- ggplot(data = mld_dt_df, aes(x = ifelse(long < 0, long + 360, long), y = lat,  fill = mld_var))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits)+
  geom_polygon(data = map_data('world2'), aes(x  = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a
ggsave(plot = a,
       filename =  'analysis/images/mld/scr_mld_density_thresh.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_thresh.eps',
       scale = .8,height = h, width = 7.25)

a <- ggplot(data = mld_dt_std, aes(x = ifelse(long < 0, long + 360, long), y = lat,  fill = mld_var))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits)+
  geom_polygon(data = map_data('world2'), aes(x  = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'SD')
a
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_algo_sd.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_algo_sd.eps',
       scale = .8,height = h, width = 7.25)
a <- ggplot(data = mld_da_std, aes(x = ifelse(long < 0, long + 360, long), y = lat,  fill = mld_var))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits)+
  geom_polygon(data = map_data('world2'), aes(x  = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'SD')
a
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_thresh_sd.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
a
ggsave(plot = a,
       filename = 'analysis/images/mld/scr_mld_density_thresh_sd.eps',
       scale = .8,height = h, width = 7.25)


# Compare
mld_all <- full_join(mld_avg, mld_dt_df, by = c('long', 'lat'))
a <- ggplot(data = mld_all, aes(x = long_p, y = lat,  fill = (mld_var - mean)/mld_var))+
  geom_raster()+
  scale_fill_gradientn(limits = c(-2, 2), values = scales::rescale(c(-2,-1,-.5, 0,.5,1, 2)),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  geom_polygon(data = map_data('world2'), aes(x  = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nDifference')
a
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_comparison_relative_thresh.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_comparison_relative_thresh.eps',
       scale = .8,height = h, width = 7.25)

mld_all <- full_join(mld_avg, mld_da_df, by = c('long', 'lat'))
a <- ggplot(data = mld_all, aes(x = long_p, y = lat,  fill = (mld_var - mean)/mld_var))+
  geom_raster()+
  scale_fill_gradientn(limits = c(-2, 2), values = scales::rescale(c(-2,-1,-.5, 0,.5,1, 2)),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  geom_polygon(data = map_data('world2'), aes(x  = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Relative\nDifference')
a
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_comparison_relative_algo.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_comparison_relative_algo.eps',
       scale = .8,height = h, width = 7.25)

a <- ggplot(data = mld_all, aes(x = mld_var, y = mean))+
  geom_point(size = .2, alpha = .2)+
  labs(x = 'Holte et al. (2017) Density\nAlgorithm Climatology (decibars)', y = 'Functional Bootstrap Estimate (decibars)')+
  scale_x_continuous(limits = c(0, 500))+
  scale_y_continuous(limits = c(0, 500))+
  scale_color_viridis_c()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw() + theme(text = element_text(size = 15),
                     legend.position = 'bottom',legend.text = element_text(size = 10))
a
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_comparison_line_algo.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_comparison_line_also.eps',
       scale = .8,height = h, width = 7.25)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))
# variance estimates
mld_all <- full_join(mld_avg, mld_da_std, by = c('long', 'lat'))

a <- ggplot(data = mld_all, aes(x = long_p, y = lat, fill = sd))+
  geom_raster() +
  map_plot+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'SD')
a
ggsave(plot = a,
       filename = 'analysis/images/mld/mld_var_dens_sd_mean.png',
       scale = .8,height = h, width = 7.25, dpi = 200)
ggsave(plot = a,
       filename ='analysis/images/mld/mld_var_dens_sd_mean.eps',
       scale = .8,height = h, width = 7.25)



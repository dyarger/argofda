
# RG_comparison
load('analysis/results/mean_one_stage/mean_review.RData')

avg_mean <- cbind(avg_mean, grid_used)
avg_mean$long <- rep(ifelse(grid_used$long< 0, grid_used$long + 360, grid_used$long), each = 206)
avg_mean$lat <- rep(grid_used$lat, each = 206)

load('products/RG_mean_mid_feb.RData')

mean_field$longitude <- round(mean_field$longitude*2)/2
mean_field$latitude <- round(mean_field$latitude*2)/2

library(dplyr)
library(ggplot2)

mean_field <- mean_field %>%
  rename(temperature_RG = temperature, salinity_RG = salinity,
         dens_RG = density, pdens_RG = pdensity, long = longitude, lat = latitude) %>%
  mutate(long = ifelse(long > 360, long - 360, long))

compare_df <- full_join(x = avg_mean, mean_field, by = c('long', 'lat', 'pressure'))
save(compare_df, file = 'analysis/results/compare_mean.RData',
     version = 2)

load('analysis/results/compare_mean.RData')
load('analysis/data/RG_Defined.RData')

compare_df$diff <- compare_df$temperature_RG - compare_df$temperature

compare_df$diff_salinity <- compare_df$salinity_RG - compare_df$salinity

compare_df <- full_join(compare_df, RG_defined_long, by = c('long', 'lat'))
compare_df <- compare_df[!is.na(compare_df$value),] # limit to RG mask
p_val <- 1500

quants <- max(quantile(compare_df$diff[compare_df$pressure == p_val], probs = c(.01, .99),na.rm = T))

a <- ggplot()+
  stat_ecdf(data = compare_df[compare_df$pressure == p_val,], aes(x = diff))+
  coord_cartesian(xlim = c(-quants, quants))+
  theme_bw()+
  theme( text = element_text(size = 15))+
  labs(x = 'Difference between RG and Functional Means (°C)', y = 'Proportion of Differences')
a
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_temp_ecdf_', p_val, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_temp_ecdf_', p_val, '.eps'),
       scale = .8,height = 4, width = 7.25)


quants <- max(quantile(compare_df$diff_salinity[compare_df$pressure == p_val], probs = c(.01, .99),na.rm = T))


a <- ggplot()+
  stat_ecdf(data = compare_df[compare_df$pressure == p_val,], aes(x = diff_salinity))+
  coord_cartesian(xlim = c(-quants, quants))+
  theme_bw()+
  theme( text = element_text(size = 15))+
  labs(x = 'Difference between RG and Functional Means (PSU)', y = 'Proportion of Differences')
a
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_psal_ecdf_', p_val, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_psal_ecdf_', p_val, '_eps.eps'),
       scale = .8,height = 4, width = 7.25)

# compare by locations

colors_plot <- RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1]
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white', color = 'black',
                           size = .2)

a <- ggplot()+
  geom_raster(data = compare_df[compare_df$pressure == p_val,],
              aes(x = long, y = lat, fill = as.numeric(diff))) +
  colors_plot+
  scale_fill_gradientn(limits = c(-2, 2), values = scales::rescale(c(-1,-.5,-.25, 0, .25,.5,1)),
                       colours = colors_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 15))+
  labs(x = 'Longitude', y = 'Latitude')
a
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_temp_', p_val, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_temp_', p_val, '.eps'),
       scale = .8,height = 4, width = 7.25)


a <- ggplot()+
  geom_raster(data = compare_df[compare_df$pressure == p_val,],
              aes(x = long, y = lat, fill = as.numeric(diff_salinity))) +
  map_plot+
  scale_fill_gradientn(limits = c(-.4, .4), values = scales::rescale(c(-.2,-.1, 0,.1,.2)),
                       colours = colors_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), text = element_text(size = 15))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'PSU')
a
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_psal_', p_val, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  paste0('analysis/images/mean_one_stage/rg_compare_psal_', p_val, '_eps.eps'),
       scale = .8,height = 4, width = 7.25)

# MLD results
load('analysis/results/mld/mld_grid_summaries.RData')
locations <- do.call(rbind, lapply(get_mld, function(x) x[[1]]))
mld_quantiles <- unlist(lapply(get_mld, function(x) x[[2]]), recursive = F)
mld_var_dens_quan <- do.call(rbind, lapply(mld_quantiles, function(x) t(x)[1,]))
colnames(mld_var_dens_quan) <- paste0('quant_', c(.025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .975))
mld_var_dens_mean <- unlist(lapply(get_mld, function(x) x[[3]]), recursive = F)
mld_var_dens_mean <- do.call(c, lapply(mld_var_dens_mean, function(x) x[1]))
mld_var_dens_sd <- unlist(lapply(get_mld, function(x) x[[4]]), recursive = F)
mld_var_dens_sd <- do.call(c, lapply(mld_var_dens_sd, function(x) x[1]))
mld_df <- data.frame(locations, mld_var_dens_quan, mean = mld_var_dens_mean,
                     sd = mld_var_dens_sd)
get_mld <- list()
for (j in 1:116) {
  load(paste0('analysis/results/mld/mld_uncond_', j,'.RData'))
  locations_uncond <- final_list[[1]]
  for (l in 1:nrow(locations_uncond)) {
    indexes <- mld_df$long == locations_uncond$long[[l]] &
      mld_df$lat == locations_uncond$lat[[l]]
    yearly_medians <- mld_df[indexes, 'quant_0.5']
    mld_df$summaries_above[indexes] <- sapply(yearly_medians, function(x) {mean(
      final_list[[2]][[l]][,1] > x, na.rm = T)})
    mld_df$summaries_below[indexes] <- sapply(yearly_medians, function(x) {mean(
      final_list[[2]][[l]][,1] < x, na.rm = T)})
  }
  # summaries <- lapply(final_list[[2]], function(x) {
  #   apply(x, 2, quantile, na.rm = T, probs = c(.025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .975))})
  # mean_summaries <- lapply(final_list[[2]], function(x) {
  #   apply(x, 2, mean, na.rm = T)})
  # sd_summaries <- lapply(final_list[[2]], function(x) {
  #   apply(x, 2, sd, na.rm = T)})


  get_mld[[j]] <- list(final_list[[1]], summaries,mean_summaries,sd_summaries)
  print(j)
}

save(get_mld, file = 'analysis/results/mld/mld_uncond_grid_summaries.RData')
q()

library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(ncdf4)
library(dplyr)
load('analysis/results/mld/mld_uncond_grid_summaries.RData')
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

scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
load('analysis/data/RG_Defined_mask.RData')
a <- ggplot(data = mld_df %>% inner_join(RG_defined_long) %>% filter(value > 1999), aes(x = long_p, y = lat, fill = mean))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a

a <- ggplot(data = mld_df %>% inner_join(RG_defined_long) %>% filter(value > 1999),
            aes(x = long_p, y = lat, fill = sd))+
  geom_raster() +
  map_plot +
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a

mld_df_year <- mld_df
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
mld_df <- data.frame(locations,mld_var_dens_quan_year =  mld_var_dens_quan,
                     mean_year = mld_var_dens_mean, sd_year = mld_var_dens_sd)
mld_df$long_p <- ifelse(mld_df$long < 0, mld_df$long + 360, mld_df$long)

mld_df_all <- full_join(mld_df_year, mld_df, by = c('long', 'lat', 'long_p'))
scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
load('analysis/data/RG_Defined_mask.RData')
a <- ggplot(data = mld_df_all %>%
              filter(year == 2014) %>% inner_join(RG_defined_long) %>%
              filter(value > 1999), aes(x = long_p, y = lat,
                                        fill = (mean_year - mean)/sd))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2(limits = c(-2.5, 2.5))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a

b <- ggplot(data = mld_df_all %>%
              filter(year == 2014) %>% inner_join(RG_defined_long) %>%
              filter(value > 1999), aes(x = long_p, y = lat, fill = (mean_year - mean)/(sd+sd_year)))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2(limits = c(-2.5, 2.5))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
b
library(patchwork)
a+b
ggsave(a/b, file = 'analysis/images/mld/mld_uncon_compare.png',
       height = 9, width = 5)
mld_df_all$z_scores <- (mld_df_all$mean_year - mld_df_all$mean)/
  #(mld_df_all$sd + mld_df_all$sd_year)
  (mld_df_all$sd)
mld_p_val <- mld_df_all %>%
  group_by(long, lat, long_p) %>%
  summarise(t_stat = sum(z_scores^2),
            nyear = n())
mld_p_val <- mld_p_val %>%
  mutate(p_val = pchisq(lower.tail = FALSE, q = t_stat,   df = nyear - 1))
head(mld_p_val)
summary(mld_p_val$p_val)

a <- ggplot(data = mld_p_val  %>% inner_join(RG_defined_long) %>%
              filter(value > 1999), aes(x = long_p, y = lat, fill = p_val))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2(limits = c(0,1), midpoint = .5)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a


a <- ggplot(data = mld_df_all %>%
              filter(year == 2012) %>% inner_join(RG_defined_long) %>%
              filter(value > 1999), aes(x = long_p, y = lat, fill = abs((mean_year - mean)/sd) > 2))+
  geom_raster() +
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a

a <- ggplot(data = mld_df_all %>%
              filter(year == 2014) %>% inner_join(RG_defined_long) %>%
              filter(value > 1999), aes(x = long_p, y = lat, fill = sd - sd_year))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2(limits = c(-100, 100))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a


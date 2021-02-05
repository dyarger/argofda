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
  #get_mld[[j]] <- list(final_list[[1]], summaries,mean_summaries,sd_summaries)
  print(j)
}

save(mld_df, file = 'analysis/results/mld/mld_uncond_summaries_above_below.RData')
q()
load('analysis/results/mld/mld_uncond_summaries_above_below.RData')
scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
mld_df$long_p <- ifelse(mld_df$long < 0, mld_df$long + 360, mld_df$long)

load('analysis/data/RG_Defined_mask.RData')
a <- ggplot(data = mld_df %>% filter(year == 2014) %>%
              inner_join(RG_defined_long) %>% filter(value > 1999),
            aes(x = long_p, y = lat, fill = summaries_above))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2(midpoint = .5)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Prob')
# P (uncondition > median(conditional_2014))
a
ggsave(plot = a, filename = 'analysis/images/mld/exceedance_2014.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)

a <- ggplot(data = mld_df %>% filter(year == 2014) %>%
              inner_join(RG_defined_long) %>% filter(value > 1999),
            aes(x = long_p, y = lat, fill = summaries_above < .1))+
  geom_raster() +
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
# P (uncondition > median(conditional_2014))
a

a <- ggplot(data = mld_df %>% filter(year == 2014) %>%
              inner_join(RG_defined_long) %>% filter(value > 1999),
            aes(x = long_p, y = lat, fill = summaries_below))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2()+
  labs(x = 'Longitude', y = 'Latitude', fill = 'MLD')
a


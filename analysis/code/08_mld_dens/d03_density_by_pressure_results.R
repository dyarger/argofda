
locations <- list()
proportions <- list()

for(i in 1:116) {
  load(paste0('analysis/results/dens/dens_p_', i,
       '.RData'))
  locations[[i]] <- final_list[[1]]
  proportions[[i]] <- final_list[[2]]
}
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))

locations <- do.call(rbind, locations)
proportions <- unlist(proportions, recursive = F)
proportions <- do.call(rbind, proportions)

delta <- 2
pressure_vec <- seq(0, 2000, by = delta)
proportions_100 <- proportions[,400]
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
                           color = 'black', fill = 'white', size = .2)
library(dplyr)
load('analysis/data/RG_Defined_mask.RData')
ggplot()+
  geom_raster(data = data.frame(locations, proportions_100)[locations$year == 2015,] %>%
                inner_join(RG_defined_long) %>% filter(value > 1999),
              aes(x = ifelse(long < 0, long + 360, long), y =lat, fill = proportions_100))+
  map_plot +
  scale_fill_gradient2(midpoint = .5, low = 'blue', mid = 'white', high= 'red')+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Prop Invert')

proportions_100 <- proportions[,875]

ggplot()+
  geom_raster(data = data.frame(locations, proportions_100)[locations$year == 2015,] %>%
                inner_join(RG_defined_long) %>% filter(value > 1999),
              aes(x = ifelse(long < 0, long + 360, long), y =lat, fill = proportions_100))+
  map_plot +
  scale_fill_gradient2(midpoint = .5, low = 'blue', mid = 'white', high= 'red')+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Prop Invert')

proportions_100 <- proportions[,276]
ggplot()+
  geom_raster(data = data.frame(locations, proportions_100)[locations$year == 2015,] %>%
                inner_join(RG_defined_long) %>% filter(value > 1999),
              aes(x = ifelse(long < 0, long + 360, long), y =lat,
                  fill = 1- proportions_100))+
  map_plot +
  scale_fill_gradient2(midpoint = .5, low = 'blue', mid = 'white', high= 'red')+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion')
ggsave('analysis/images/misc/dens_monotone_pressure_550.png',
       height = 4, width = 7.25)

ggplot()+
  geom_raster(data = data.frame(locations, proportions_100)[locations$year == 2015,] %>%
                inner_join(RG_defined_long) %>% filter(value > 1999),
              aes(x = ifelse(long < 0, long + 360, long), y =lat,
                  fill = 1- proportions_100))+
  map_plot +
  scale_fill_gradient2(midpoint = .5, low = 'darkred', mid = 'gray', high= 'yellow')+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion')
ggsave('analysis/images/misc/dens_monotone_pressure_550_bw.png',
       height = 4, width = 7.25)

proportions_100 <- proportions[,276]
ggplot()+
  geom_raster(data = data.frame(locations, proportions_100)[locations$year == 2015,] %>%
                inner_join(RG_defined_long) %>% filter(value > 1999),
              aes(x = ifelse(long < 0, long + 360, long), y =lat, fill = proportions_100 > .95))+
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Prop Invert')
ggsave('analysis/images/misc/dens_monotone_pressure_550_ind.png',
       height = 4, width = 7.25)




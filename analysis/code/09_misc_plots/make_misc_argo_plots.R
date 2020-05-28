
library(argofda)
library(gridExtra)
library(ggplot2)
library(viridis)

load('analysis/data/jan_march_data.RData')
relevant_indices <- which(profYearAggr == 2016 & profJulFracofYearAggr > 31
                          & profJulFracofYearAggr < 60)
lat <- profLatAggr[relevant_indices]
long <- profLongAggr[relevant_indices]
long_p <- ifelse(long > 360, long - 360, long)
temp <- sapply(profTempAggr[relevant_indices], function(x) x[[1]][1])
g <- ggplot() +
  geom_polygon(data = map_data('world2'), aes(x = long, y =lat, group = group),fill = 'white', color = 'black', size = .2)+
  geom_point(data = data.frame(lat, long_p, temp),
             aes(x = long_p, y = lat, color = temp),size = .05) +
  labs(color = 'Surface temperature measurement for profile (°C)',
       x = 'Longitude', y = 'Latitude') +
  #scale_color_viridis(option = 'D')  +
  scale_color_gradientn(colours = colorRamps::matlab.like(10))  +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(size = 13))
g
ggsave(filename = 'analysis/images/misc/profile_locations.png', g, height = 4, width = 6)
ggsave(filename = 'analysis/images/misc/profile_locations.eps',
       g, height = 4, width = 6)

# Observations per profile
load('analysis/data/Argo_data_aggr.RData')
prof_lengths <- sapply(profPresAggr, function(x) length(x[[1]]))
a <- ggplot(data = data.frame(len = prof_lengths))+
  geom_histogram(aes(x = len), binwidth = 25)+
  labs(x = 'Observations in Profile', y = 'Number of Profiles')
ggsave(filename = 'analysis/images/misc/obs_per_prof.png', a,
       scale = 1.1,height = 4, width = 7.25)
ggsave(filename = 'analysis/images/misc/obs_per_prof_eps.eps', a,
       scale = 1.1,height = 4, width = 7.25)



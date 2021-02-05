

# check argo profiles

load('analysis/data/jan_march_data.RData')
profDensAggr <- list()
for (i in 1:length(profLatAggr)) {
  temp <- profTempAggr[[i]][[1]]
  psal <- profPsalAggr[[i]][[1]]
  pressure <- profPsalAggr[[i]][[1]]
  profDensAggr[[i]] <- gsw::gsw_pot_rho_t_exact(psal, temp, pressure, p_ref = 0)
  if (i %% 2000 == 0) {
    print(i)
  }
}
library(dplyr)

mean_decreasing <- sapply(1:length(profDensAggr), function(x)  {
  dens <- profDensAggr[[x]][,1]
  pres <- profPresAggr[[x]][[1]]
  lagged <- lag(dens)
  lagged_pressure <- lag(pres)
  weights <- pres - lagged_pressure
  sum(weights * ((dens - lagged) >= 0), na.rm = T)/sum(weights, na.rm = T)
})

mean_decreasing_550 <- sapply(1:length(profDensAggr), function(x)  {
  dens <- profDensAggr[[x]][,1]
  pres <- profPresAggr[[x]][[1]]
  lagged <- lag(dens)
  #lagged_pressure <- lag(pres)
  #weights <- pres - lagged_pressure
  #findInterval(550, pres, all.inside = T)
  (dens - lagged)[findInterval(550, pres, all.inside = T)+1] >= 0
})

mean_decreasing_550_val <- sapply(1:length(profDensAggr), function(x)  {
  dens <- profDensAggr[[x]][,1]
  pres <- profPresAggr[[x]][[1]]
  lagged <- lag(dens)
  lagged_pressure <- lag(pres)
  #weights <- pres - lagged_pressure
  #findInterval(550, pres, all.inside = T)
  interv <- findInterval(550, pres, all.inside = T)+1
  (dens[interv] - lagged[interv])/(pres[interv] - lagged_pressure[interv])
})

load('analysis/data/RG_Defined_mask.RData')

df <- data.frame(profLongAggr, profLatAggr, mean_decreasing,mean_decreasing_550,mean_decreasing_550_val,
                 profYearAggr, profMonthAggr)
df$long_grid <- round(ifelse(df$profLongAggr > 180,df$profLongAggr - 360,
                             df$profLongAggr)+ .5)- .5
df$lat_grid <- round(df$profLatAggr + .5)- .5
df_comb <- inner_join(df, RG_defined_long, by = c('long_grid' = 'long',
                                                  'lat_grid'  = 'lat'))
df_comb <- df_comb[df_comb$value > 1999,]
df_comb <- df_comb[!is.na(df_comb$value),]
library(ggplot2)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))
ggplot()+
  geom_point(data = df_comb,
             aes(x = ifelse(profLongAggr >360, profLongAggr - 360, profLongAggr) ,
                 y = profLatAggr, color = mean_decreasing),
             size = .1)+
  scale_color_gradient2(low = 'blue', mid = 'white', name = 'Proportion',high = 'red',
                        midpoint = .6,
                        limits = c(.18, 1))+
  geom_polygon(data=  map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2) +
  labs(x='Longitude', y = 'Latitude') +
  theme_gray()+ theme(panel.grid = element_blank(), text = element_text(size = 15))
ggsave('analysis/images/misc/dens_monotone.png',
       height = 4, width = 7.25)

ggplot()+
  geom_point(data = df_comb[df_comb$profYearAggr == 2015 & df_comb$profMonthAggr==2,],
             aes(x = ifelse(profLongAggr >360, profLongAggr - 360, profLongAggr) ,
                 y = profLatAggr, color = mean_decreasing_550),
             size = .01)+
  scale_color_discrete()+
  geom_polygon(data=  map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2) +
  labs(x='Longitude', y = 'Latitude') +
  theme_gray()+ theme(panel.grid = element_blank(), text = element_text(size = 15))
library(viridis)
ggplot()+
  geom_point(data = df_comb[df_comb$profYearAggr == 2015 & df_comb$profMonthAggr==2,],
             aes(x = ifelse(profLongAggr >360, profLongAggr - 360, profLongAggr) ,
                 y = profLatAggr, color = mean_decreasing_550_val),
             size = .2)+
  # scale_color_gradientn(colors = rev(viridis_pal()(50)),values = scales::rescale(c(-.0005,  0, .001, .003, .006, .008)),
  #                       limits = quantile(df_comb$mean_decreasing_550_val[df_comb$profYearAggr == 2015 & df_comb$profMonthAggr==2],
  #                                         probs = c(.01, .99)))+
  scale_color_gradient2(limits = quantile(df_comb$mean_decreasing_550_val[df_comb$profYearAggr == 2015 & df_comb$profMonthAggr==2],
                                          probs = c(.0001, .995)),
                        low = 'blue', mid = 'white', high = 'red')+
  geom_polygon(data=  map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2) +
  #g + xlab( expression(Value~is~sigma~R^{2}==0.6))
  # labs(x='Longitude', y = 'Latitude',
  #      color = expression(Density'\n'~is)) +
  labs(x='Longitude', y = 'Latitude',
       #color = expression(paste('Density\nGradient\n', 'a (', kg, '/',m^{3},'/p)'))) +
  color =expression(atop("Density Gradient",
                         "(kg/"*m^{3}*"/p)"))) +
  theme_gray()+ theme(legend.title = element_text(size = 12),
                      panel.grid = element_blank(), text = element_text(size = 15))
ggsave('analysis/images/misc/dens_monotone_550_profiles.png',
       height = 4, width = 7.25)

b <- ggplot()+
  geom_polygon(data=  map_data('world'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2) +
  geom_point(data = df_comb,
             aes(x = ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)  ,
                 y = profLatAggr, color = mean_decreasing),
             size = .1)+
  scale_color_gradient2(low = 'blue', mid = 'white',
                        name = 'Prop Monotone',high = 'red', midpoint = .6,
                        limits = c(.18, 1))+
  labs(x='Longitude', y = 'Latitude') +
  coord_cartesian(xlim = c(-100, 20), ylim = c(0, 60))+
  #facet_wrap(~profYearAggr, ncol = 3) +
  theme_gray()
ggsave('analysis/images/misc/dens_monotone_atl.png',
       height = 4, width = 7.25)


load('analysis/results/density_check_pred.RData')
avg_prop_nonnegative <- sapply(density_check, function(x) x[[1]])
library(ggplot2)

df_preds_use_summary <- data.frame(df_preds_use, avg_prop_nonnegative) %>%
  group_by(long, lat) %>%
  summarise(dens = mean(avg_prop_nonnegative, na.rm = T))

a <- ggplot(data = df_preds_use_summary %>%
              inner_join( RG_defined_long) %>% filter(value > 1999),
       aes(x =ifelse(long < 0, long + 360, long), y = lat, fill = dens))+
  geom_raster()+
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',midpoint = .6,
                       limits = c(.18, 1))+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  #coord_cartesian(xlim = c(-100, 20), ylim = c(0, 60))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion')
a
ggsave('analysis/images/misc/dens_prop_ours.png',
       height = 4, width = 7.25)
library(patchwork)
a / b

a <- ggplot(data = df_preds_use_summary #%>%
       #       inner_join( RG_defined_long) %>% filter(value > 1999)
            ,
            aes(x =ifelse(long < 0, long + 360, long), y = lat, fill = dens))+
  geom_raster()+
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',midpoint = .6,
                       limits = c(.18, 1))+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  coord_cartesian(xlim = c(160, 250), ylim = c(-75, -55))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'February Predictions, reference pressure 0 dbar')
a
b <- ggplot()+
  geom_polygon(data=  map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2) +
  geom_point(data = df,
             aes(x = profLongAggr  ,
                 y = profLatAggr, color = mean_decreasing),
             size = .1)+
  scale_color_gradient2(low = 'blue', mid = 'white',
                        name = 'Prop Monotone',high = 'red', midpoint = .6,
                        limits = c(.18, 1))+
  labs(x='Longitude', y = 'Latitude') +
  coord_cartesian(xlim = c(160, 250), ylim = c(-75, -55))+
  theme_gray()
b
a/b

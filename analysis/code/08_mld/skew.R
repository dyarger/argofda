# MLD results
get_mld <- list()
library(dplyr)
for (j in 1:116) {
  load(paste0('analysis/results/mld/mld_', j ,'.RData'))
  locations <- final_list[[1]]
  groups <- ((paste0(locations$long, '_', locations$lat)))
  # summaries <- sapply(final_list[[2]], function(x) {
  #   median(x[[1]][,1], na.rm = T)})
  #variances <- sapply(1:length(final_list[[2]]), function(x) {
  #  median(abs(final_list[[2]][[x]][[1]][,1] - summaries[x]), na.rm = T)})
  mean_summaries <- sapply(final_list[[2]], function(x) {
    mean(x[[1]][,1], na.rm = T)})
  variance_summarise <- sapply(1:length(final_list[[2]]), function(x) {
    var(final_list[[2]][[x]][[1]][,1], na.rm = T)})

  skew_measure <-  sapply(1:length(final_list[[2]]), function(x) {
    mean(((final_list[[2]][[x]][[1]][,1] - mean_summaries[x])/sqrt(variance_summarise[x]))^3, na.rm = T)})

  get_mld[[j]] <- data.frame(locations, skew_measure)
  print(j)
}
get_mld_try <- do.call(rbind, get_mld)
library(ggplot2)
ggplot(data = get_mld_try[get_mld_try$year == 2014,], aes(x = long, y = lat, fill = skew_measure))+
  geom_raster()+
  scale_fill_gradient2(limits = c(-5, 10)) +
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))
limits <- c(0, 1000)
h <- 6.5
scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400, 500, 1000)
limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
get_mld_try$long_p <- ifelse(get_mld_try$long < 0, get_mld_try$long + 360, get_mld_try$long)


load('analysis/data/RG_Defined_mask.RData')
a <- ggplot(data = get_mld_try %>% filter(year == 2012) %>%
              inner_join(RG_defined_long) %>% filter(value > 1999),
            aes(x = long_p, y = lat,
                fill = ifelse(skew_measure > 23 | is.na(skew_measure), 23, skew_measure)))+
  geom_raster() +
  map_plot +
  scale_fill_gradient2(limits = c(-2, 23.1), low = 'blue',
                       mid = 'white', high = 'red')+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Skew')
# P (uncondition > median(conditional_2014))
a
ggsave(plot = a, filename = 'analysis/images/mld/skew_2012.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)


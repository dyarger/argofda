# MLD results
get_mld <- list()
library(dplyr)
for (j in 1:116) {
  load(paste0('analysis/results/mld/mld_uncond_', j ,'.RData'))
  locations_uncond <- final_list[[1]]
  summaries_uncond <- sapply(final_list[[2]], function(x) {
    median(x[,1], na.rm = T)})
  load(paste0('analysis/results/mld/mld_', j ,'.RData'))
  locations <- final_list[[1]]
  groups <- ((paste0(locations$long, '_', locations$lat)))
  # summaries <- sapply(final_list[[2]], function(x) {
  #   median(x[[1]][,1], na.rm = T)})
  #variances <- sapply(1:length(final_list[[2]]), function(x) {
  #  median(abs(final_list[[2]][[x]][[1]][,1] - summaries[x]), na.rm = T)})
  summaries <- sapply(final_list[[2]], function(x) {
    median(x[[1]][,1], na.rm = T)})
  variances <- sapply(1:length(final_list[[2]]), function(x) {
    median(abs(final_list[[2]][[x]][[1]][,1] - summaries[x]), na.rm = T)})

  variances_vals <- t(sapply(1:length(final_list[[2]]), function(x) {
    abs(final_list[[2]][[x]][[1]][,1] - summaries[x])}))
  mat_groups <- matrix(groups, nrow = length(groups), ncol = 1000, byrow = F)
  splits <- split(variances_vals, mat_groups)
  split_names <- names(splits)
  split_summaries <- sapply(splits, median, na.rm =T)
  split_summaries <- data.frame('location' = names(split_summaries),
                                med = split_summaries)
  split_summaries <- tidyr::separate(split_summaries,
                                              'location', into = c('long','lat'),
                                              sep = '_')

  variances_uncond <- t(sapply(1:length(final_list[[2]]), function(x) {
    summary_use <- summaries_uncond[locations_uncond$long == locations$long[x] &
                                      locations_uncond$lat == locations$lat[x]]
    abs(final_list[[2]][[x]][[1]][,1] - summary_use)}))
  mat_groups_uncond <- matrix(groups, nrow = length(groups), ncol = 1000, byrow = F)
  split_summaries_uncond <- split(variances_uncond, mat_groups_uncond)
  split_summaries_uncond <- sapply(splits_uncond, median, na.rm =T)
  split_summaries_uncond <- data.frame('location' = names(split_summaries_uncond),
                                med_uncond = split_summaries_uncond)
  split_summaries_uncond <- tidyr::separate(split_summaries_uncond,
                                     'location', into = c('long','lat'),
                                     sep = '_')
  variance_loc_all <- full_join(split_summaries,
                                split_summaries_uncond)
  get_mld[[j]] <- variance_loc_all
  print(j)
}
get_mld_try <- do.call(rbind, get_mld)
get_mld_try$long <- as.numeric(as.character(get_mld_try$long))
get_mld_try$lat <- as.numeric(as.character(get_mld_try$lat))
library(ggplot2)

scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400)
limits <- c(0, 400)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
ggplot(data = get_mld_try  %>% inner_join(RG_defined_long) %>% filter(value > 1999),
       aes(x = ifelse(long < 0, long + 360, long), y = lat, fill = med))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion')
ggsave(filename = 'analysis/images/mld/MAbsoluteE_y.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)

theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))

limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)

load('analysis/data/RG_Defined_mask.RData')
ggplot(data = get_mld_try  %>% inner_join(RG_defined_long) %>% filter(value > 1999),
       aes(x = ifelse(long < 0, long + 360, long), y = lat,
           fill = med/(med_uncond+med)))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, 1))+
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion')
ggsave(filename = 'analysis/images/mld/mld_prop_var.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)

scale <- c(seq(0, 100, by  =10),120, 140, 160, 180, 200, 300, 400)
limits <- c(0, 400)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
load('analysis/results/mld/mld_grid_summaries_median_anova_mean_real.RData')
get_mld_try <- do.call(rbind, get_mld)
get_mld_try$long <- as.numeric(as.character(get_mld_try$long))
get_mld_try$lat <- as.numeric(as.character(get_mld_try$lat))
ggplot(data = get_mld_try  %>% inner_join(RG_defined_long) %>% filter(value > 1999),
       aes(x = ifelse(long < 0, long + 360, long), y = lat, fill = variance_sum))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = colorRamps::matlab.like(10), limits = limits) +
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Value')
ggsave(filename = 'analysis/images/mld/MAE_y.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)
scale <- c(seq(0, 50, by  =3),90, 130, 170, 210, 300, 400)
ggplot(data = get_mld_try  %>% inner_join(RG_defined_long) %>% filter(value > 1999),
       aes(x = ifelse(long < 0, long + 360, long), y = lat, fill = variance_sum))+
  geom_raster()+
  scale_fill_gradientn(values = scales::rescale(scale),
                       colours = viridis::magma(10), limits = limits) +
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Value')
ggsave(filename = 'analysis/images/mld/MAE_y_bw.png',
       scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)


#save(get_mld, file = 'analysis/results/mld/mld_grid_summaries_median_anova_mean_real.RData')
q()

# MLD results
get_mld <- list()
library(dplyr)
for (j in 1:116) {
  load(paste0('analysis/results/mld/mld_', j ,'.RData'))
  locations <- final_list[[1]]
  locations_unique <- locations[!duplicated(paste0(locations$long, '_', locations$lat)),]
  pvals <- list()
  for (y in 1:nrow(locations_unique)) {
    locations_one <- which(locations$long == locations_unique$long[y] &
                             locations$lat == locations_unique$lat[y])
    gridded_vals <- expand.grid(year1 = locations$year[locations_one],
                                year2 = locations$year[locations_one])
    gridded_vals <- gridded_vals[gridded_vals$year1 < gridded_vals$year2,]
    gridded_vals$p_val <- 0
    gridded_vals$year1 <- as.character(gridded_vals$year1)
    gridded_vals$year2 <- as.character(gridded_vals$year2)
    samples <- lapply(final_list[[2]][locations_one],
                      function(x) { x[[1]][,1]})
    names(samples) <- locations$year[locations_one]#;q
    for(m in 1:nrow(gridded_vals)) {
      s1  <- samples[[gridded_vals[m,1]]]
      s2 <- samples[[gridded_vals[m,2]]]
      if (sum(!is.na(s1)) < 2  | sum(!is.na(s2)) < 2) {
        gridded_vals$p_val[m] <- NA
      } else {
        gridded_vals$p_val[m] <- wilcox.test(samples[[gridded_vals[m,1]]],
                                             samples[[gridded_vals[m,2]]])$p.value
      }
    }
    pvals[[y]] <- gridded_vals
  }
  get_mld[[j]] <- list(locations_unique, pvals)
  print(j)
}

locations <- lapply(get_mld, function(x) x[[1]])
locations <- do.call(rbind, locations)
pvals <- lapply(get_mld, function(x) x[[2]])
pvals <- unlist(pvals, recursive = F)
pvals <- sapply(pvals, function(x) {mean(x$p_val, na.rm = T)})

#get_mld_try <- do.call(rbind, get_mld)
library(ggplot2)
ggplot(data = cbind(locations, pvals),
       aes(x = long, y = lat, fill = pvals))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, .35))

ggplot(data = cbind(locations, pvals),
       aes(x = long, y = lat, fill = pvals > .05))+
  geom_raster()

#save(get_mld, file = 'analysis/results/mld/wilcox_test.RData')
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15),
                             legend.position = 'bottom',legend.text = element_text(size = 10)))
limits <- c(0, 1000)
h <- 6.5
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white',
                           color = 'black', size = .2)
load('analysis/data/RG_Defined_mask.RData')
ggplot(data = get_mld_try  %>% inner_join(RG_defined_long) %>% filter(value > 1999),
       aes(x = ifelse(long < 0, long + 360, long), y = lat, fill = variance_sum/sum_all))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, 1.01))+
  map_plot +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion')
# ggsave(filename = 'analysis/images/mld/mld_prop_var.png',
#        scale = .8,height = h, width = 7.25, units = 'in', dpi = 200)

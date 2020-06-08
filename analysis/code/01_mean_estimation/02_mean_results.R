# Make plots from mean
library(ggplot2)
library(viridis)
library(argofda)
library(dplyr)
type <- 'psal'
label <- c('Temperature', 'Salinity')[c('temp', 'psal') == type]
label_units <- c('Temp (°C)', 'PSU')[c('temp', 'psal') == type]
color_scale <- list(c(-2, 32),c(30, 38))[[which(c('temp', 'psal') == type)]]
calc_files <- list.files('analysis/results/mean_estimation/', pattern = type)

# Load in mean funcitons
final_list <- list()
for (j in 1:length(calc_files)) {
  load(paste0('analysis/results/mean_estimation/', calc_files[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  eval(parse(text=paste0('rm(', file, ')')))
  funs <- lapply(now, function(x) {
    if (x[[3]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[1]])
    }
  })
  samples <- lapply(now, function(x) {
    if (x[[3]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[2]])
    }
  })
  prof_info <- lapply(now, function(x) {
    if (x[[3]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[3]])
    }
  })
  final_list[[j]] <- list(funs,prof_info, samples)
  print(j)
}

funs <- lapply(final_list, function(x) return(x[[1]]))
funs <- unlist(funs, recursive = FALSE)

have_results <- which(sapply(funs, length) != 1)
funs <- funs[have_results]
cv_info <- lapply(funs, function(x) x[[3]])
cv_niter <- sapply(cv_info, nrow)
summary(cv_niter)


prof_info <- unlist(lapply(final_list, function(x) return(x[[2]])), recursive = F)
samples <- unlist(lapply(final_list, function(x) return(x[[3]])), recursive = F)
samples <- samples[have_results]
samples <- do.call(rbind, samples)
prof_info <- prof_info[have_results]

prof_info <- do.call(rbind, prof_info)
colnames(prof_info) <- c('index', 'seconds', 'measurements', 'profiles', 'h_space')
summary(prof_info)

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute_only <- grid_to_compute[have_results,]
grid_to_compute_only$long_p <- ifelse(grid_to_compute_only$long  < 0,
                                    grid_to_compute_only$long + 360,
                                    grid_to_compute_only$long)

## single example
long_example <- -170.5
lat_example <- 0.5

mean_single_location <- function(long_example, lat_example, pressure_vals = seq(0, 2000, by = 1)) {
  one_func <- funs[[which(grid_to_compute_only$long == long_example & grid_to_compute_only$lat == lat_example)]]
  pred_vals <- sapply(0:9, function(x) predict.local_function(one_func, index = x, deriv = 0, p_vals = pressure_vals))
  colnames(pred_vals) <- 2007:2016
  pred_vals_long <- tidyr::gather(as.data.frame(pred_vals), year, value)
  pred_vals_long$year <- readr::parse_number(pred_vals_long$year)
  pred_vals_long$pressure <- pressure_vals
  pred_vals_long
}

ll_single_location <- function(long_example, lat_example) {
  one_func <- funs[[which(grid_to_compute_only$long == long_example & grid_to_compute_only$lat == lat_example)]]
  pressure_vals <- seq(0, 2000, by = 1)
  pred_vals <- sapply(10:16, function(x) predict.local_function(one_func, index = x, deriv = 0, p_vals = pressure_vals))
  colnames(pred_vals) <- c('ll_long', 'll_lat', 'lq_long', 'lq_lat', 'lq_long_lat', 'll_tim', 'lq_tim')
  data.frame(pressure = pressure_vals, pred_vals)
}

# For one location, get different functions
pred_vals_long <- mean_single_location(long_example, lat_example)
ll_vals <- ll_single_location(long_example, lat_example)
time_15 <- pred_vals_long$value[pred_vals_long$year == 2010] - 200 * ll_vals[, 'll_lat'] +(-200)^2/2 * ll_vals[, 'lq_lat']
pred_vals_wide <- tidyr::spread(pred_vals_long, year, value)
avg_fun <- apply(pred_vals_wide[,-1], 1, mean)
pred_vals_long <- rbind(pred_vals_long,
                        data.frame(year = '2010, 200 km Lat', value = time_15, pressure =pred_vals_long$pressure[pred_vals_long$year == 2010]),
                        data.frame(year = 'Avg', value = avg_fun, pressure =pred_vals_long$pressure[pred_vals_long$year == 2010]))
theme_set(theme_bw())

# reproduce plot from paper
load('analysis/data/jan_march_data.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

data_mode <- c('all', 'delayed')[c('temp', 'psal') == type]
data <- get_profile_data_subset(lat = lat_example, long = long_example, day = 45.25,RG_def = T,h_space = 900, h_time = 45.25,
                                mode = data_mode,
                                  exclusion = T,min_prof = 30, years = c(2010))
ylim_example <- list(c(10, 31),c(34, 36.5))[[which(c('temp', 'psal') == type)]]
data$response <- list(data$temperature,data$salinity)[[which(c('temp', 'psal') == type)]]


ggplot(data = pred_vals_long[pred_vals_long$year %in% c('2010', '2013', 'Avg'),], aes(x = pressure, y = value))+
  geom_line(data= data[data$year == 2010,], mapping = aes(x = pressure, y = response,group = profile),
            alpha = .1)+
  geom_point(data= data[data$year ==  2010,],
             aes(x = pressure, y = response),size = .2, alpha = .1)+
  geom_line(size = 1.05, aes(color = as.factor(year)))+
  coord_flip(xlim = c(350, 0), ylim = ylim_example)+
  scale_x_continuous(trans = 'reverse')+
  scale_color_viridis_d()+
  labs(x = 'Pressure (decibars)', y = label_units,color = 'Function')+
  theme(legend.position = 'bottom')
ggsave(filename = paste0('analysis/images/mean_estimation/mean_', type, '_example_', long_example, '_', lat_example, '.png'),
       scale = 1.2, height = 5, width = 3.75, dpi = 150)
ggsave(filename = paste0('analysis/images/mean_estimation/mean_', type, '_example_', long_example, '_', lat_example, '_large.png'),
       scale = 1.2, height = 5, width = 3.75, dpi = 600)

# extrapolate local linear day function
ll_vals <- ll_single_location(-19.5+360, -19.5)
time_seq <- seq(0, 90, by = 1)
deriv_time_mat <- sapply(time_seq, function(x) {
  (x- 45.25) * ll_vals[,7] + (x- 45.25)^2 * ll_vals[,8]
})
fields::image.plot(deriv_time_mat)

# Consider average mean function values along long_section of the ocean
long_section <- -179.5
pred_vals_section <- lapply(seq(-50.5, 50.5, by = 1),function(x) {
  test <- tidyr::spread(mean_single_location(long_section, x), year, value)
  apply(test[,-1], 1, mean)})
pred_vals_section <- do.call(rbind, pred_vals_section)
colnames(pred_vals_section) <- 0:2000
rownames(pred_vals_section) <- seq(-50.5, 50.5, by = 1)
png(filename =  paste0('images/mean_estimation/mean_', type, '_section_', long_section, '.png'),
    width = 600, height = 450)
fields::image.plot(pred_vals_section, ylim = c(1, 0), xlab = 'Latitude', ylab = 'Pressure (decibar)', axes=F,
                   legend.cex = 8, cex.axis = 2.2)
axis(2, at = seq(0, 2000, by = 250)/2000, seq(0, 2000, by = 250), cex = 5)
axis(1, at=(seq(-50.5, 50.5, by = 10) + 50.5)/103.5, labels=seq(-50.5, 50.5, by = 10))
dev.off()

# Consider the LL function in time along long_section of the ocean
long_section <- -179.5
pred_vals_section <- lapply(seq(-50.5, 50.5, by = 1),function(x) {
 ll_single_location(long_section, x)[,7]})
pred_vals_section <- do.call(rbind, pred_vals_section)
colnames(pred_vals_section) <- 0:2000
rownames(pred_vals_section) <- seq(-50.5, 50.5, by = 1)
png(filename =  paste0('images/mean_estimation/mean_', type, '_section_ll_time_', long_section, '.png'),
    width = 600, height = 450)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
fields::image.plot(pred_vals_section, ylim = c(1, 0), xlab = 'Latitude', ylab = 'Pressure (decibar)', axes=F,
                   legend.cex = 1.8, cex.lab = 1.8, legend.lab = '°C/day', legend.line = 4,
                   legend.mar = 7, axis.args=list(cex.axis=1.6))
axis(2, at = seq(0, 2000, by = 250)/2000, seq(0, 2000, by = 250), cex.axis = 1.6)
axis(1, at=(seq(-50.5, 50.5, by = 10) + 50.5)/103.5, labels=seq(-50.5, 50.5, by = 10),
     cex.axis = 1.6)
dev.off()

# Consider the LL function in latitude along long_section of the ocean
long_section <- -179.5
pred_vals_section <- lapply(seq(-50.5, 50.5, by = 1),function(x) {
  ll_single_location(long_section, x)[,3]})
pred_vals_section <- do.call(rbind, pred_vals_section)
colnames(pred_vals_section) <- 0:2000
rownames(pred_vals_section) <- seq(-50.5, 50.5, by = 1)
png(filename =  paste0('images/mean_estimation/mean_', type, '_section_ll_lat_', long_section - .5, '.png'),
    width = 600, height = 450)
fields::image.plot(pred_vals_section, ylim = c(1, 0), xlab = 'Latitude', ylab = 'Pressure (decibar)', axes=F,
                   legend.cex = 8, cex.axis = 2.2)
axis(2, at = seq(0, 2000, by = 250)/2000, seq(0, 2000, by = 250), cex = 5)
axis(1, at=(seq(-50.5, 50.5, by = 10) + 50.5)/103.5, labels=seq(-50.5, 50.5, by = 10))
dev.off()

png(filename =  paste0('images/mean_estimation/mean_', type, '_section_ll_lat_', long_section - .5, '.png'),
    width = 600, height = 450)
par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
fields::image.plot(pred_vals_section*100, ylim = c(1, 0), xlab = 'Latitude', ylab = 'Pressure (decibar)', axes=F,
                   legend.cex = 1.8, cex.lab = 1.8, legend.lab = '°C/100km', legend.line = 4,
                   legend.mar = 7, axis.args=list(cex.axis=1.6))
axis(2, at = seq(0, 2000, by = 250)/2000, seq(0, 2000, by = 250), cex.axis = 1.6)
axis(1, at=(seq(-50.5, 50.5, by = 10) + 50.5)/103.5, labels=seq(-50.5, 50.5, by = 10),
     cex.axis = 1.6)
dev.off()


# Mean at fixed pressure for all locations

theme_set(theme_bw()+theme(panel.grid = element_blank()))
map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group), fill = 'white', color = 'black',
                           size = .2)
colors_plot <- RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1]
# Yearly mean at 10 decibars
make_mean_at_pressure <- function(p_val) {
  all_years <- t(sapply(funs, function(y) {
    sapply(0:9, function(z) argofda::predict.local_function(y, p_vals = p_val, index = z))
  }))
  colnames(all_years) <- 2007:2016
  pred_df <- tidyr::gather(data.frame(grid_to_compute_only, all_years), year, value, -c(1,2,3))
  pred_df$year <- readr::parse_number(pred_df$year)
  pred_df
}

p_val <- 10
pred_df <- make_mean_at_pressure(p_val = p_val)
ggplot()+
  geom_raster(data = pred_df, aes(x = long_p, y = lat, fill = value))+
  scale_fill_gradientn(limits = color_scale,
                       colours = colors_plot)+
  map_plot+
  facet_wrap(~year) +
  labs(x = 'Longitude', y = 'Latitude',
       fill = label_units,title = paste0('Yearly Estimated ', label, ', ',  p_val,' decibars'))
ggsave(filename = paste0('analysis/images/mean_estimation/mean_', type, '_', p_val,'.png'), scale = 1.2,
       height = 5, width = 7.25)

# average over years
pred_avg_df <- pred_df %>% group_by(long_p, lat) %>% summarise(value = mean(value))
ggplot()+
  geom_raster(data = pred_avg_df, aes(x = long_p, y = lat, fill = value))+
  map_plot +
  scale_fill_gradientn(limits = color_scale,colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude',
       fill = label_units, title = paste0('Average Estimated ', label, ', ', p_val,' decibars'))
ggsave(filename = paste0('analysis/images/mean_estimation/avg_mean_', type, '_', p_val, '.png'), scale = 1.2,
       height = 5, width = 7.25)


# Lambdas chosen
lambdas <- sapply(funs, function(x) {return(x[[4]])})
ggplot()+
  geom_raster(data = cbind(grid_to_compute_only, lambdas = as.numeric(lambdas)),
              aes(x = long_p, y = lat, fill = log10(lambdas)))+
  map_plot +
  scale_fill_viridis()+
  labs(x = 'Longitude', y = 'Latitude', fill = 'log10(lambda)',
       title = paste0('Lambda Chosen for ', label))
ggsave(filename = paste0('analysis/images/mean_estimation/lambdas_', type, '.png'), scale = 1.2,
       height = 5, width = 7.25)


samples <- tidyr::gather(cbind(grid_to_compute_only, samples), year, n_prof , -long, -lat, -long_p)
ggplot(data = samples, aes(x =long_p, y = lat, fill = n_prof))+
  geom_raster()+
  facet_wrap(~year)+
  map_plot +
  scale_fill_viridis_c()+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Profiles\nUsed',
       title = paste0('Profiles Used for Each Year, ', label))
ggsave(filename = paste0('analysis/images/mean_estimation/profiles_year_', type, '.png'), scale = 1.2,
       height = 5, width = 7.25)

# amount of time
ggplot(data = cbind(grid_to_compute_only, seconds = prof_info[,2]), aes(x =long_p, y = lat, fill = seconds))+
  geom_raster()+
  map_plot +
  scale_fill_viridis_c()+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Seconds',
       title = paste0('Time to compute, ', label))
ggsave(filename = paste0('analysis/images/mean_estimation/time_', type, '.png'), scale = 1.2,
       height = 5, width = 7.25)

# number of measurements
ggplot(data = cbind(grid_to_compute_only, measurements = prof_info[,3]), aes(x =long_p, y = lat, fill = measurements))+
  geom_raster()+
  map_plot +
  scale_fill_viridis_c()+
  labs(x = 'Longitude', y = 'Latitude', fill = '#',
       title = paste0('Number of Measurements Used, ', label))
ggsave(filename = paste0('analysis/images/mean_estimation/measurements_', type, '.png'), scale = 1.2,
       height = 5, width = 7.25)
# number of profiles
ggplot(data = cbind(grid_to_compute_only, profiles = prof_info[,4]), aes(x =long_p, y = lat, fill = profiles))+
  geom_raster()+
  map_plot
  scale_fill_viridis_c()+
  labs(x = 'Longitude', y = 'Latitude', fill = '#',
       title = paste0('Number of Profiles Used, ', label))
ggsave(filename = paste0('analysis/images/mean_estimation/profiles_', type, '.png'), scale = 1.2,
       height = 5, width = 7.25)

###### Derivative functions #####
p_val <- 300

# in time
pred_ll_time <- sapply(funs, argofda::predict.local_function, p_vals = p_val, index = 15)

limit_time <- quantile(abs(pred_ll_time), probs = .98)

ggplot()+
  geom_raster(data =data.frame(grid_to_compute_only, pred_ll_time),
              aes(x = long_p, y = lat, fill = pred_ll_time))+
  map_plot
  scale_fill_gradientn(limits = c(-limit_time, limit_time),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill =  paste0(label_units,  '\nper day'))
ggsave(filename = paste0('analysis/images/mean_estimation/gradient_time_', type, '_', p_val, 'new.png'), scale = .8,
       height = 5, width = 7.25)

#in space
pred_ll_long <- sapply(funs, argofda::predict.local_function, p_vals = p_val, index = 10)
pred_ll_lat <- sapply(funs, argofda::predict.local_function, p_vals = p_val, index = 11)

ll_df <- data.frame(grid_to_compute_only, long_grad = pred_ll_long, lat_grad = pred_ll_lat)
ll_df$angle <- atan2(x = ll_df$long_grad, y = ll_df$lat_grad)
ll_df$length <- sqrt(ll_df$long_grad^2+ ll_df$lat_grad^2)

color_scale_deriv <- list(c(0, .25, .5, 1, 2, 3, 5),c(0, .05,.15, .25, .5, 1, 2, 3))[[which(c('temp', 'psal') == type)]]
color_limits_deriv <- list(c(0, 5),c(0, 1))[[which(c('temp', 'psal') == type)]]

h <- 4.5
# gradient length
ggplot()+
  geom_raster(data =ll_df,
            aes(x = long_p, y = lat, fill = as.numeric(length)*100))+
  map_data +
  scale_fill_gradientn(limits = color_limits_deriv, values = color_scale_deriv,
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude',
    fill = paste0('Gradient\nLength\n(', label_units,  '/100km)'))
ggsave(filename = paste0('analysis/images/mean_estimation/gradient_length_', type, '_', p_val, '.png'), scale = 1.2,
       height = 5, width = 7.25)

# gradient length and arrows
ggplot()+
  geom_tile(data =ll_df[ll_df$lat > -15,],
            aes(x = long, y = lat, fill = as.numeric(length)*100,
                color = as.numeric(length)*100))+
  geom_spoke(data = ll_df[ll_df$lat > -15,],
             aes(x = long, y = lat, angle = as.numeric(angle)),
             color = 'black',
             radius = as.numeric(1.4), size = .3,
             arrow = arrow(length = unit(0.07, "cm"))
  )+
  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size =.1)+
  # coord_map("ortho",orientation = c(-85/5, 70, 0), #indian ocean
  #           xlim = c(40, 120), ylim = c(-40,20))+
  # coord_map("ortho",orientation = c(85/7, -150, 0), # north pacific
  #           xlim = c(-220, -80), ylim  = c(-10,45))+
  coord_map("ortho",orientation = c(85/7, -25, 0), #atl
            xlim = c(-40, -10), ylim  = c(10,55))+
  # coord_map("ortho",orientation = c(85/7, -25, 0), #atl
  #           xlim = c(-80, 20), ylim  = c(-10,55))+
  scale_fill_gradientn(values = scales::rescale(color_scale_deriv),
                       colours = c('white', 'red'), limits = color_limits_deriv) +
  scale_color_gradientn(values = scales::rescale(color_scale_deriv),
                        colours = c('white', 'red'), limits = color_limits_deriv) +
  labs(x = 'Longitude', y = 'Latitude',
       color = paste0('Gradient\nLength\n(', label_units,  '/100km)'),
       fill = paste0('Gradient\nLength\n(', label_units,  '/100km)'))
ggsave(filename = paste0('analysis/images/mean_estimation/deriv_', type, '_', p_val, '.png'),
       scale = .8,height = h*1.3, width = h*1.4, dpi = 100)
ggsave(filename = paste0('analysis/images/mean_estimation/deriv_', type, '_', p_val, '_large.png'),
       scale = .8,height = h*1.3, width = h*1.4, dpi = 400)
ggsave(filename = paste0('analysis/images/mean_estimation/deriv_', type, '_', p_val, '.eps'),
       scale = .8,height = h*1.3, width = h*1.4)






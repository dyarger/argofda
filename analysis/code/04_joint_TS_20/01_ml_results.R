library(fda)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library(argofda)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))

calc_files <- list.files('analysis/results/joint_TS_20/',
                         pattern = 'joint')
# Get results from score modelling
final_list <- list()
prof_info <- list()
K1 <- 10
K2 <- 10
K <- K1 + K2
for (j in 1:length(calc_files)) {
  load(paste0('analysis/results/joint_TS_20/', calc_files[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  if (sum(sapply(now, length)==0) > 0 | sum(sapply(now, length)==1) > 0) {
     print(j)
  }
  eval(parse(text=paste0('rm(', file, ')')))

  indexes_used <- as.numeric(sapply(now, function(x) {
    if (length(x) == 3) {
        x[1]
      } else {
        x[[10]][1]
      }
  }))
  df_return <- lapply(1:length(now), function(x) {
    if (length(now[[x]]) == 1 | length(now[[x]]) == 0 | length(now[[x]]) == 3) {
      mat <- matrix(c(indexes_used[x],
                      rep(NA, K^2 + 5*K + K + 1)),nrow = 1, ncol = 1 +K^2 + 5*K + K + 1)
    } else {
      n_steps <- length(now[[x]][[1]])
      m_res <- now[[x]][[1]][[n_steps]][[1]]
      mat <- cbind(matrix(c(indexes_used[x],as.vector(now[[x]][[6]]),
                            as.vector(sapply(1:K, function(x) exp(m_res[[x]]$par)))
                            ),
                          nrow = nrow(now[[x]][[3]]), ncol = 1 + K^2 + K*5 ,byrow = T),
                   now[[x]][[3]])
    }
    colnames(mat) <- c('index',
                       paste0('var', rep(1:K, each = K), '_', rep(1:K, times = K)),
                       paste0(rep(c('scale_long_', 'scale_lat_', 'scale_tim_', 'sigma2_', 'var_score_'), times = K),
                              rep(1:K, each = 5)),
                       'year', paste0('T', 1:K1), paste0('S', 1:K2))
    mat
  })
  final_list[[j]] <- list(do.call(rbind, df_return), indexes_used)
  prof_info[[j]] <- lapply(1:length(now), function(x) {
    if (length(now[[x]]) == 1 | length(now[[x]]) == 0 | length(now[[x]]) == 3) {
      NULL
    } else {
      now[[x]][[10]]
    }
  })
}
# add longitude/latitude information
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$index <- 1:nrow(grid_to_compute)

df_params <- lapply(final_list, function(x) {return(x[[1]])})
df_params <- do.call(rbind, df_params)
df_params <- df_params[!is.na(df_params[,2]),]
df_params <- merge(df_params, grid_to_compute, by = 'index')
df <- df_params
df$long_p <- ifelse(df$long <0 , df$long+360, df$long)

# save parameters
df_params_only <- df_params %>%
  dplyr::filter(!duplicated(df_params$index)) %>%
  dplyr::select(-year, -starts_with('T', ignore.case = 'F'), -starts_with('S', ignore.case = 'F'))
save(df_params_only, file = 'analysis/results/joint_TS_20/params.RData')

# information on calculations
prof_info <- unlist(prof_info, recursive = F)
prof_info <- prof_info[!sapply(prof_info, is.null)]
prof_info <- do.call(rbind, prof_info)
summary(prof_info)

### Predictions at 10, 300, and 1500 decibars###
load('analysis/results/psal_one_stage_cov_pca.RData')
load('analysis/results/temp_one_stage_cov_pca.RData')
year <- 2012

# Salinity
df <- merge(df, grid_pc_psal)
p_vals <- c(10, 300, 1500)
anomaly_pred_psal <- t(sapply(which(df$year == year), function(x) {
  coef <- as.double(df[x, paste0('S', rep(1:K2))])
  pc_coef <- coefs_pc_psal[[which(grid_pc_psal$long == df[x,'long'] &
                               grid_pc_psal$lat == df[x,'lat'])]][[1]][,1:K2]
  vals <- rowSums(sapply(1:K2, function(y) coef[y] *pc_coef[,y] ))
  predict(make_fit_object(coefficients = vals, knots = knots_pc),
          p_vals)$y
}))
colnames(anomaly_pred_psal) <- paste0('pred', p_vals)

map_plot <-   geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
                           color = 'black', fill = 'white', size = .2)
colors_plot <-  c("blue", "white", "red")

### 10 decibars ###
ggplot()+
  geom_raster(data = cbind(df[df$year ==year,], pred = anomaly_pred_psal[,1]),
              aes(x = long_p, y = lat, fill = pred))+
  map_plot +
  scale_fill_gradientn(limits = c(-.9, .9), colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Pred (PSU)')
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[1], '_psal_', year, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[1], '_psal_', year, '.eps'),
       scale = .8,height = 4, width = 7.25)

### 300 decibars ###
ggplot()+
  geom_raster(data = cbind(df[df$year ==year,], pred = anomaly_pred_psal[,2]),
              aes(x = long_p, y = lat, fill = pred))+
  map_plot +
  scale_fill_gradientn(limits = c(-.55, .55), values = scales::rescale(c(-.4,-.25, -.1, 0, .1,.25,.4)),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Pred (PSU)')
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[2], '_psal_', year, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[2], '_psal_', year, '.eps'),
       scale = .8,height = 4, width = 7.25)

### 1500 decibars ###
a <- ggplot()+
  geom_raster(data = cbind(df[df$year ==year,], pred = anomaly_pred_psal[,3]),
              aes(x = long_p, y = lat, fill = pred))+
  map_plot +
  scale_fill_gradientn(limits = c(-.2, .2), values = scales::rescale(c(-.12, -.05,-.02, -.01, 0,.01, .02, .05,.12)),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Pred (PSU)')
a
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[3], '_psal_', year, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[3], '_psal_', year, '.eps'),
       scale = .8,height = 4, width = 7.25)

## temperature
# calculate temperature anomalies
anomaly_pred <- t(sapply(which(df$year == year), function(x) {
  coef <- as.double(df[x, paste0('T', rep(1:K1))])
  pc_coef <- coefs_pc_temp[[which(grid_pc_temp$long == df[x,'long'] &
                                    grid_pc_temp$lat == df[x,'lat'])]][[1]][,1:K1]
  vals <- rowSums(sapply(1:K1, function(y) coef[y] *pc_coef[,y] ))
  predict(make_fit_object(coefficients = vals, knots = knots_pc),
          c(10, 300, 1500))$y
}))
colnames(anomaly_pred) <- paste0('pred', p_vals)

### 10 decibars ###
x_min <- -4.657153
x_max <- 3.716673
ggplot()+
  geom_raster(data = cbind(df[df$year ==year,], pred = anomaly_pred[,1]),
              aes(x = long_p, y = lat, fill = pred))+
  map_plot +
  scale_fill_gradientn(limits = c(x_min, x_max), values = scales::rescale(c(x_min,-2, 0,2, x_max)),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Pred (°C)')
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[1], '_temp_', year, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[1], '_temp_', year, '.eps'),
       scale = .8,height = 4, width = 7.25)



### 300 decibars ###
x_min <- -5.079887
x_max <-6.767269
ggplot()+
  geom_raster(data = cbind(df[df$year ==year,], pred = anomaly_pred[,2]),
              aes(x = long_p, y = lat, fill = pred))+
  map_plot +
  scale_fill_gradientn(limits = c(x_min, x_max), values = scales::rescale(c(x_min,-2, 0,2, x_max)),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Pred (°C)')
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[2], '_temp_', year, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[2], '_temp_', year, '.eps'),
       scale = .8,height = 4, width = 7.25)

## 1500 decibars ###
x_min <- -0.6243476
x_max <- 1.276535
ggplot()+
  geom_raster(data = cbind(df[df$year ==2012,], pred = anomaly_pred[,3]),
              aes(x = long_p, y = lat, fill = pred))+
  map_plot +
  scale_fill_gradientn(limits = c(x_min, x_max), values = scales::rescale(c(x_min, -.5, 0,.5, x_max)),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Pred (°C)')
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[3], '_temp_', year, '.png'),
       scale = .8,height = 4, width = 7.25)
ggsave(filename = paste0('analysis/images/joint_TS_20/anomaly_', p_vals[3], '_temp_', year, '.eps'),
       scale = .8,height = 4, width = 7.25)

########################
# Parameters - don't need more than 1 per year
# plots of scale and variance parameters, perhaps not necessarily needed for all of them
# there are more plots than needed here

df_params_only <- df[!duplicated(df$index),]

a <- ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                      y = lat, fill = sqrt(var1_1)))+
  map_plot +
  scale_fill_gradientn(colours = colorRamps::matlab.like(10), limits = c(0, 70)) +
  labs(x = 'Longitude', y = 'Latitude', fill = expression('SD\n°C/' * phi[1] * '(p)'))
a
ggsave(plot = a,
       filename = 'analysis/images/joint_TS_20/sd_score_1.png',
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_1.eps',
       scale = .8,height = 4, width = 7.25)

a <- ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                      y = lat, fill = sqrt(var2_2)))+
  geom_polygon(data = map_data('world2'), aes(x = long, y =lat, group = group),
               color = 'black', fill = 'white', size  = .2)+
  scale_fill_gradientn(colours = colorRamps::matlab.like(10), limits = c(0, sqrt(1200))) +
  labs(x = 'Longitude', y = 'Latitude', fill = expression('SD\n°C/' * phi[2] * '(p)'))
a
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_2.png',
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_2.eps',
       scale = .8,height = 4, width = 7.25)
a <- ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = sqrt(var11_11)))+
  map_plot +
  scale_fill_gradientn(colours = colorRamps::matlab.like(10), limits = c(0, 13)) +
  labs(x = 'Longitude', y = 'Latitude', fill = expression('SD\nPSU/' * psi[1] * '(p)'))
a
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_1_psal.png',
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_1_psal.eps',
       scale = .8,height = 4, width = 7.25)
a <- ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = sqrt(var12_12)))+
  map_plot +
  scale_fill_gradientn(colours = colorRamps::matlab.like(10), limits = c(0, 7)) +
  labs(x = 'Longitude', y = 'Latitude', fill = expression('SD\nPSU/' * psi[2] * '(p)'))
a
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_2_psal.png',
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/sd_score_2_psal.eps',
       scale = .8,height = 4, width = 7.25)

a <- ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = abs(var1_11)/sqrt(var1_1*var11_11)))+
  map_plot +
  scale_fill_gradient2(limits = c(-1, 1)) +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')
a
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/cor_score_1_11.png',
       scale = .8,height = 4, width = 7.25)
ggsave(plot = a,
       filename =  'analysis/images/joint_TS_20/cor_score_1_11.eps',
       scale = .8,height = 4, width = 7.25)

a <- ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = scale_lon_1))+
  map_plot +
  scale_fill_gradientn(limits = c(0,4000),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')
a

ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = scale_lat_1))+
  map_plot +
  scale_fill_gradientn(limits = c(0,600),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')

ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = scale_tim_2))+
  map_plot +
  scale_fill_gradientn(limits = c(0, 200),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')

ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = var_score_1))+
  map_plot +
  scale_fill_gradientn(limits = c(0, 2000),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')
ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = var_score_2))+
  map_plot +
  scale_fill_gradientn(limits = c(0,600),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')

ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = sigma2_1))+
  map_plot +
  scale_fill_gradientn(limits = c(0,50),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')

ggplot()+
  geom_raster(data = df_params_only, aes(x = long_p,
                                         y = lat, fill = sigma2_2))+
  map_plot +
  scale_fill_gradientn(limits = c(0,10),
                       colours = colors_plot)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')

# Make plots from mean
setwd('/home/dyarger/argofda/')
library(ggplot2)
library(viridis)
library(argofda)
library(dplyr)
type <- 'temp'
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

temp_funs <- funs
grid_temp <- grid_to_compute_only
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


#
psal_funs <- funs
grid_psal <- grid_to_compute_only


# load('/Users/dyarger/Documents/argofda/analysis/results/joint_TS_20/preds.RData')
# load('/Users/dyarger/Documents/argofda/analysis/results/temp_complete_cov_pca.RData')
# load('/Users/dyarger/Documents/argofda/analysis/results/psal_complete_cov_pca.RData')
load('analysis/results/joint_TS_20/preds.RData')
load('analysis/results/temp_cov_pca.RData')
load('analysis/results/psal_cov_pca.RData')

temp_mean <- which(rownames(grid_temp)  %in% rownames(grid_pc_psal))#$long &
psal_mean <- which(rownames(grid_psal)  %in% rownames(grid_pc_psal))
temp_pcs_red <- which(rownames(grid_pc_temp)  %in% rownames(grid_pc_psal))
psal_funs <- psal_funs[psal_mean]
temp_funs <- temp_funs[temp_mean]
coefs_pc_temp <- coefs_pc_temp[temp_pcs_red]


delta = 1
p_vals = seq(0, 2000, by = delta)

library(gsw)
library(Matrix)
density_check <- list()
for (i in 1:nrow(df_preds)) {
  long <- df_preds$long[i]
  lat <- df_preds$lat[i]
  year <- df_preds$year[i]
  index <- which(grid_pc_psal$long == long &
                   grid_pc_psal$lat == lat)
  psal_fun_single <- argofda::predict.local_function(psal_funs[[index]], p_vals = p_vals, index = year - 2007)
  temp_fun_single <- argofda::predict.local_function(temp_funs[[index]], p_vals = p_vals, index = year - 2007)


  temp_pc_funs <- coefs_pc_temp[[index]][[1]]
  psal_pc_funs <- coefs_pc_psal[[index]][[1]]
  temp_pc <- sapply(1:10, function(x) predict(make_fit_object(knots = knots_pc,
                                                      coefficients = temp_pc_funs[,x]), p_vals)$y)
  temp_pred <- temp_fun_single + rowSums(temp_pc %*% Diagonal(n = 10, as.vector(unlist(df_preds[i, paste0('T', 1:10)]))))

  psal_pc <- sapply(1:10, function(x) predict(make_fit_object(knots = knots_pc,
                                                              coefficients = psal_pc_funs[,x]), p_vals)$y)
  psal_pred <- psal_fun_single + rowSums(psal_pc %*% Diagonal(n = 10, as.vector(unlist(df_preds[i, paste0('S', 1:10)]))))


  psal_pc_deriv <- sapply(1:10, function(x) predict(make_fit_object(knots = knots_pc,
                                                              coefficients = psal_pc_funs[,x]), p_vals,deriv = 1)$y)

  temp_pc_deriv <- sapply(1:10, function(x) predict(make_fit_object(knots = knots_pc,
                                                                    coefficients = temp_pc_funs[,x]), p_vals,deriv = 1)$y)

  psal_pred_deriv <-  argofda::predict.local_function(psal_funs[[index]], p_vals = p_vals, index = year - 2007, deriv = 1)+
    rowSums(psal_pc_deriv %*%Diagonal(n = 10, as.vector(unlist(df_preds[i, paste0('S', 1:10)]))))
  temp_pred_deriv <-  argofda::predict.local_function(temp_funs[[index]], p_vals = p_vals, index = year - 2007, deriv = 1)+
    rowSums(temp_pc_deriv %*%Diagonal(n = 10, as.vector(unlist(df_preds[i, paste0('T', 1:10)]))))
  SA <- gsw::gsw_SA_from_SP(psal_pred, p_vals, latitude = lat, longitude = long)
  CT <- gsw::gsw_CT_from_t(SA, t = temp_pred,  p_vals)
  N2_list <- gsw::gsw_Nsquared(SA = SA, CT = CT, p = p_vals, latitude = lat)
  N2 <- N2_list$N2
  pdens_P <- gsw::gsw_pot_rho_t_exact(psal_pred, temp_pred, p_vals, p_ref = 0)

  # plot(N2, type = 'l')
  # plot(pdens_P)

  nc_file <- ncdf4::nc_open('~/Downloads/CZ16_1_2000m_PD_year_2018_month_05.nc')
  ncdf4::ncvar_get(nc_file, 'depth_std')
  levels <- c(1,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,
              200,250,300,350,400,450,500,550,600,650,700,750,800,
              850,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,2000)
  plot(p_vals[-1], N2, type = 'l')
  points(levels, rep(0, length(levels)), col = 2)


  pdens_P <- gsw::gsw_pot_rho_t_exact(psal_pred, temp_pred, p_vals, p_ref = 0)
  pdens_T <- gsw::gsw_pot_rho_t_exact(psal_pred, temp_pred + rep(.01, times = length(temp_pred)),
                                      p_vals, p_ref = 0)
  pdens_S <- gsw::gsw_pot_rho_t_exact(psal_pred +  rep(.01, times = length(temp_pred)), temp_pred, p_vals, p_ref = 0)
  f_pdens <- ((lead(pdens_P) - lag(pdens_P))/delta)[c(-1, -length(pdens_P))]
  f_T <- (pdens_T - pdens_P)/(.01)
  f_S <- (pdens_S - pdens_P)/(.01)

  deriv_pdens <- f_T * temp_pred_deriv + f_S * psal_pred_deriv + c(NA, f_pdens, NA)
  deriv_pdens <- deriv_pdens[-c(1, length(deriv_pdens))]
  pred_prop_nonnegative <- mean(deriv_pdens >=0)

  pred_integral_nonnegative <- sum((deriv_pdens * (deriv_pdens < 0))) * delta

  density_check[[i]] <- list(pred_prop_nonnegative, pred_integral_nonnegative)
  if (i %%500 == 0) {
    print(i)
  }
}

df_preds_use <- df_preds[, c('long', 'lat', 'year')]

save(density_check, grid_psal, df_preds_use,
     file = 'analysis/results/density_check_pred.RData')
q()
load('analysis/results/density_check_pred.RData')
avg_prop_nonnegative <- sapply(density_check, function(x) x[[1]])
library(ggplot2)
theme_set(theme_bw())

df_preds_use_summary <- data.frame(df_preds_use, avg_prop_nonnegative) %>%
  group_by(long, lat) %>%
  summarise(dens = mean(avg_prop_nonnegative, na.rm = T))

ggplot(data = df_preds_use_summary,
       aes(x = ifelse(long <0, long + 360, long), y = lat, fill = dens))+
  geom_raster()+
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = .5)+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'February Predictions, reference pressure 0 dbar')
ggsave(file = 'analysis/images/misc/dens_avg.png', height = 4, width = 6)


library(ncdf4)
RG <- nc_open('/Volumes/DrewHardDrive/Roemmich-Gilson/RG_ArgoClim_Temperature_2017.nc')
BM <- ncvar_get(RG, 'BATHYMETRY_MASK')
MM <- ncvar_get(RG, 'MAPPING_MASK')
image.plot(MM[,,1])

longitude <- ncvar_get(RG, "LONGITUDE")
longitude <- ifelse(longitude > 360, longitude - 360, longitude)
longitude <- ifelse(longitude > 180, longitude - 360, longitude)
latitude <- ncvar_get(RG,"LATITUDE")
pressure <- ncvar_get(RG,"PRESSURE")
time <- ncvar_get(RG,"TIME")

mean_field_temp <- array(ncvar_get(RG,"MAPPING_MASK"),dim = c(360, 145, 58),
                         dimnames = list(longitude, latitude, pressure))
library(reshape2)
mean_field_long_temp <- melt(mean_field_temp, variable.name = "key",
                             value.names = "value")
RG_defined <- mean_field_temp[,,1]
RG_defined_long <- melt(RG_defined, variable.name = "key",
                        value.names = "value")
colnames(RG_defined_long) <- c('long', 'lat', 'value')
#save(RG_defined_long, file = 'analysis/data/RG_Defined_mask.RData')
ggplot(data = df_preds_use_summary %>% inner_join(RG_defined_long) %>%
         filter(value > 1999),
       aes(x = ifelse(long <0, long + 360, long), y = lat, fill = dens))+
  geom_raster()+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'February Predictions, reference pressure 0 dbar')


image.plot(BM[,,1])

ggplot(data = cbind(df_preds_use[1:length(avg_prop_nonnegative),], avg_prop_nonnegative),
       aes(x = ifelse(long <0, long + 360, long), y = lat, fill = avg_prop_nonnegative))+
  geom_raster()+
  facet_wrap(~year, ncol =5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = .5)+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'February Predictions, reference pressure 0 dbar')
ggsave(file = '/Users/dyarger/Documents/pred_prop_avg.png', height = 4, width = 6)

load('analysis/data/RG_Defined.RData')
ggplot(data = cbind(df_preds_use, avg_prop_nonnegative)[df_preds_use$year == 2015,] %>%
         inner_join(RG_defined_long),
       aes(x = ifelse(long <0, long + 360, long), y = lat, fill = avg_prop_nonnegative))+
  geom_raster()+
  #scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = .5)+
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = .5)+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'February Predictions, reference pressure 0 dbar')
ggsave(file = '/Users/dyarger/Documents/pred_prop_avg.png', height = 4, width = 6)




ggplot(data = cbind(df_preds_use, avg_prop_nonnegative),
       aes(x = ifelse(long <0, long + 360, long), y = lat, fill = avg_prop_nonnegative>.95))+
  geom_raster()+
  facet_wrap(~year, ncol =5) +
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion > .95',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'February Predictions, reference pressure 0 dbar')
ggsave(file = '/Users/dyarger/Documents/pred_prop_avg_95.png', height = 4, width = 6)

prop_nonnegative <- sapply(density_check, function(x) x[[2]])
#df_prop <- tidyr::gather(data.frame(grid_psal, prop_nonnegative), year, value, -c(1:3))

ggplot(data = cbind(df_preds_use[1:length(avg_prop_nonnegative),], prop_nonnegative),
       aes(x = ifelse(long < 0, long+ 180, long), y = lat, fill = -prop_nonnegative))+
  geom_raster()+
  facet_wrap(~year, ncol = 5) +
  #scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = .5)+
  scale_fill_viridis_c(limits = c(0, 14))+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Integral',
       title = 'Integral of Pressure Dimension with nonnegative derivative',
       subtitle = 'Yearly Mean Functions, reference pressure 0 dbar')
ggsave(file = '/Users/drewyarger/Documents/prop.png', height = 4, width = 8)
ggplot(data = cbind(df_preds_use[1:length(avg_prop_nonnegative),], avg_prop_nonnegative),
       aes(x = ifelse(long <0, long + 360, long), y = lat, fill = avg_prop_nonnegative > .95))+
  geom_raster()+
  facet_wrap(~year, ncol = 5) +
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion > .95',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'Yearly Mean Functions, reference pressure 0 dbar')
ggsave(file = '/Users/drewyarger/Documents/prop_95.png', height = 4, width = 8)

ggplot(data = cbind(df_preds_use[1:length(avg_prop_nonnegative),], prop_nonnegative),
       aes(x = long_p, y = lat, fill = -prop_nonnegative))+
  geom_raster()+
  scale_fill_viridis_c(limits = c(0, 5))+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Integral',
       title = 'Integral of Pressure Dimension with nonnegative derivative',
       subtitle = 'Year-Averaged Mean Function, reference pressure 0 dbar')
ggsave(file = '/Users/drewyarger/Documents/integral_avg.png', height = 4, width = 6)

ggplot(data = cbind(grid_psal, avg_integral_nonnegative),
       aes(x = long_p, y = lat, fill = -avg_integral_nonnegative > .1))+
  geom_raster()+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Integral > .1',
       title = 'Integral of Pressure Dimension with nonnegative derivative',
       subtitle = 'Year-Averaged Mean Function, reference pressure 0 dbar')
ggsave(file = '/Users/drewyarger/Documents/integral_avg_0_1.png', height = 4, width = 6)


int_nonnegative <- t(sapply(density_check, function(x) x[[3]]))
colnames(int_nonnegative) <- 2007:2016
df_int <- tidyr::gather(data.frame(grid_psal, int_nonnegative), year, value, -c(1:3))

ggplot(data = df_int,
       aes(x = long_p, y = lat, fill = -value))+
  geom_raster()+
  facet_wrap(~year, ncol = 5) +
  scale_fill_viridis_c(limits = c(0, 5))+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Integral',
       title = 'Integral of Pressure Dimension with nonnegative derivative',
       subtitle = 'Yearly Mean Functions, reference pressure 0 dbar')
ggsave(file = '/Users/drewyarger/Documents/integral.png', height = 4, width = 8)
ggplot(data = df_prop,
       aes(x = long_p, y = lat, fill = value > .95))+
  geom_raster()+
  facet_wrap(~year, ncol = 5) +
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Integral > .1',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'Yearly Mean Functions, reference pressure 0 dbar')
ggsave(file = '/Users/dyarger/Documents/integral_0_1.png', height = 4, width = 8)






load('/Users/dyarger/Documents/argofda/analysis/results/density_RG_mean.RData')
density_prop <- as.data.frame(density_prop)
colnames(density_prop) <- month.name
avg_prop_nonnegative <- tidyr::gather(cbind(grid_to_compute_only,density_prop),
                                      month, int, c(-1,-2))

avg_prop_nonnegative$month <- factor(avg_prop_nonnegative$month,
                                     levels = month.name)

library(ggplot2)
theme_set(theme_bw())

ggplot(data = avg_prop_nonnegative,
       aes(x = ifelse(Var1 <0, Var1 + 360, Var1), y = Var2, fill = int))+
  geom_raster()+
  facet_wrap(~month, ncol =6) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = .5)+
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'RG-type mean, reference pressure 0 dbar')
ggsave(file = '/Users/dyarger/Documents/RG_mean_prop_avg.png', height = 4, width = 6)
ggplot(data =avg_prop_nonnegative,
       aes(x =  ifelse(Var1 <0, Var1 + 360, Var1), y = Var2, fill = int>.95))+
  geom_raster()+
  facet_wrap(~month, ncol =6) +
  geom_polygon(data = map_data('world2'), aes(x =long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Proportion > .95',
       title = 'Proportion of Pressure Dimension with nonnegative derivative',
       subtitle = 'RG-type mean, reference pressure 0 dbar')
 ggsave(file = '/Users/dyarger/Documents/RG_mean_avg_95.png', height = 4, width = 6)







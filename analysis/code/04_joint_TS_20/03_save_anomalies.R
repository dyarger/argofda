
get_mean_anomalies <- function(x) {
  index_val <- which(grid_to_compute_all_year[,'index'] == df_indexes[x,'index'])
  pred_vals_temp <- sapply(0:9, function(x) predict.local_function(temp_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec))
  avg_pred_temp <- rowMeans(pred_vals_temp)
  pred_vals_psal <- sapply(0:9, function(x) predict.local_function(psal_funs[have_results_both][[index_val]],
                                                                   index = x, deriv = 0, p_vals = pressure_vec))
  avg_pred_psal <- rowMeans(pred_vals_psal)

  mean_temp_one <- pred_vals_temp[,df_indexes[x, 'year'] - 2006]
  mean_psal_one <- pred_vals_psal[,df_indexes[x, 'year'] - 2006]
  long <- df_indexes[x,'long']
  lat <- df_indexes[x,'lat']

  # Predictions
  scores_temp <- unlist(df_indexes[x, paste0('T', 1:K1)])
  scores_psal <- unlist(df_indexes[x, paste0('S', 1:K2)])
  pc_temp <- coefs_pc_temp[[which(grid_pc_temp$long == long &
                                    grid_pc_temp$lat == lat)]][[1]][,1:K1]
  pc_psal <-  coefs_pc_psal[[which(grid_pc_psal$long == long &
                                     grid_pc_psal$lat == lat)]][[1]][,1:K2]
  pred_pc_temp <- sapply(as.data.frame(as.matrix(pc_temp)), function(z) {
    predict(make_fit_object(coefficients = z, knots = knots_pc),
            pressure_vec)$y
  })
  pred_pc_psal <- sapply(as.data.frame(as.matrix(pc_psal)), function(z) {
    predict(make_fit_object(coefficients = z, knots = knots_pc),
            pressure_vec)$y
  })
  pred_temp <- rowSums(sapply(1:K1, function(x) { scores_temp[x] * pred_pc_temp[,x]}))
  pred_psal <- rowSums(sapply(1:K2, function(x) { scores_psal[x] * pred_pc_psal[,x]}))

  c(mean_temp_one,mean_psal_one,avg_pred_temp, avg_pred_psal, pred_temp, pred_psal)
}

library(ggplot2)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))
library(RColorBrewer)
library(colorRamps)
library(fda)
library(argofda)
library(gsw)
library(parallel)

load('analysis/results/joint_TS_20/preds.RData')

# get PCs
load('analysis/results/temp_one_stage_cov_pca.RData')
coefs_pc_temp <- coefs_pc
grid_pc_temp <- grid_pc
basis <- create.bspline.basis(breaks =knots_pc)
load('analysis/results/psal_one_stage_cov_pca.RData')
grid_pc_psal <- grid_pc
coefs_pc_psal <- coefs_pc

### get mean functions
calc_files_temp <- list.files('analysis/results/mean_estimation/', pattern = 'temp')
calc_files_psal <- list.files('analysis/results/mean_estimation/', pattern = 'psal')
final_list <- list()

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)

pressure_vec <- c(10, 300, 1500)
n_grid <- length(pressure_vec)
# for each mean file, predict and save
K1 <- 10; K2 <- 10
for (j in 1:length(calc_files_temp)) {
  load(paste0('analysis/results/mean_estimation/', calc_files_temp[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  temp_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][1] == 'no profiles in range') {return(NULL)} else {return(x[[1]])}
  })
  eval(parse(text=paste0('rm(', file, ')')))
  load(paste0('analysis/results/mean_estimation/', calc_files_psal[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  psal_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][1] == 'no profiles in range') {return(NULL)} else { return(x[[1]]) }
  })
  eval(parse(text=paste0('rm(', file, ')')))
  have_results <- which(sapply(temp_funs, length) != 0  & sapply(psal_funs, length) != 0)
  indexes_have_mean <- (1:500 + (j-1)*500)[have_results]
  indexes_have_both <- unique(df_preds$index[df_preds$index %in% indexes_have_mean])
  have_results_both <- which((1:500 + (j-1)*500) %in% indexes_have_both)
  df_indexes <- df_preds[df_preds$index %in% indexes_have_both,]
  a <- proc.time()
  ohc <- lapply(1:nrow(df_indexes), get_mean_anomalies)
  b <- proc.time()
  anomaly_df <- cbind(df_indexes[,c('index', 'long', 'lat', 'year')], do.call(rbind, ohc))
  colnames(anomaly_df) <- c('index', 'long', 'lat', 'year', 'mean_10_T', 'mean_300_T', 'mean_1500_T',
                        'mean_10_S', 'mean_300_S', 'mean_1500_S',
                        'avg_mean_10_T', 'avg_mean_300_T', 'avg_mean_1500_T',
                        'avg_mean_10_S', 'avg_mean_300_S', 'avg_mean_1500_S',
                        'pred_10_T', 'pred_300_T', 'pred_1500_T',
                        'pred_10_S', 'pred_300_S', 'pred_1500_S')
  final_list[[j]] <- anomaly_df
  print(j)
  print(b-a)
}
mean_anomalies <- do.call(rbind, final_list)
save(mean_anomalies, file = 'analysis/results/joint_TS_20/mean_anomalies_10_300_1500.RData')


# plot anomalies around both means
mean_anomalies$long_p <- ifelse(mean_anomalies$long < 0, mean_anomalies$long+360,
                                mean_anomalies$long)
library(ggplot2)
ggplot(data = mean_anomalies[mean_anomalies$year == 2012,],
       aes(long_p, lat, fill = pred_300_T))+
  geom_raster()+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', size = .2, color = 'black')+
  scale_fill_gradientn(limits = c(-3.5, 3.5),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  labs(x = NULL, y = NULL, fill = 'Pred (°C)')+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
ggsave(filename = 'analysis/images/joint_TS_20/around_mean_300.png')
ggplot(data = mean_anomalies[mean_anomalies$year == 2012,],
       aes(long_p, lat, fill = pred_300_T- mean_300_T + avg_mean_300_T))+
  geom_raster()+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', size = .2, color = 'black')+
  scale_fill_gradientn(limits = c(-3.5, 3.5),
                       colours = RColorBrewer::brewer.pal(name = 'RdYlBu',n = 10)[10:1])+
  labs(x = NULL, y = NULL, fill = 'Pred (°C)')+
  theme(axis.ticks = element_blank(), axis.text = element_blank())
ggsave(filename = 'analysis/images/joint_TS_20/around_avg_mean_300.png')


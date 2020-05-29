# we compute density from the mean estimates using the gsw packages
library(gsw)
library(argofda)
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)

# Roemmich and gilson pressure levels
mfpres <- c(2.5,   10.0 ,  20.0 ,  30.0,   40.0  , 50.0   ,
            60.0,   70.0,   80.0,   90.0,  100.0,  110.0,  120.0,  130.0,
            140.0,  150.0,  160.0,  170.0,  182.5,
            200.0,  220.0,  240.0,  260.0,  280.0,  300.0,  320.0,  340.0,  360.0,
            380.0,  400.0,  420.0,  440.0,  462.5,
            500.0,  550.0,  600.0,  650.0,  700.0,  750.0,  800.0,  850.0,  900.0,
            950.0 ,1000.0, 1050.0, 1100.0, 1150.0,
            1200.0, 1250.0, 1300.0, 1350.0, 1412.5, 1500.0, 1600.0, 1700.0, 1800.0,
            1900.0 ,1975.0)
# additional pressures to evaluate
pressure<- seq(0, 2000, length.out = 201)
pressure <- c(mfpres, pressure)
pressure <- unique(pressure)
pressure <- pressure[order(pressure)]
density_list <- list()

calc_files_temp <- list.files('analysis/results/mean_estimation/', pattern = 'mean_one_stage_temp')
calc_files_psal <- list.files('analysis/results/mean_estimation/', pattern = 'mean_one_stage_psal')

for (i in 1:length(calc_files_temp)) {
  load(paste0('analysis/results/mean_estimation/', calc_files_temp[i]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  eval(parse(text=paste0('rm(', file, ')')))
  removed <- sapply(now, function(x)
    if (x[[2]][[1]][1] == 'no profiles in range') {
      T} else {F})
  now <- now[!removed]
  indexes <- sapply(now, function(x) return(x[[3]][1]))

  long_temp <- grid_to_compute[indexes, 1]
  lat_temp <- grid_to_compute[indexes, 2]
  avg_temp_funs <- lapply(now, function(y) {
    rowMeans(sapply(0:9, function(z) argofda::predict.local_function(y[[1]], p_vals = pressure, index = z)))
  })

  load(paste0('analysis/results/mean_estimation/', calc_files_psal[i]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  eval(parse(text=paste0('rm(', file, ')')))
  now <- now[!removed]
  indexes <- sapply(now, function(x) return(x[[3]][1]))

  long_psal <- grid_to_compute[indexes, 1]
  lat_psal <- grid_to_compute[indexes, 2]
  avg_psal_funs <- lapply(now, function(y) {
    if (y[[2]][[1]][1] == 'no profiles in range') {
      rep(NA, length(pressure))
    } else {
      rowMeans(sapply(0:9, function(z) argofda::predict.local_function(y[[1]], p_vals = pressure, index = z)))
    }
  })

  grid_combined <- data.frame(long = long_temp, lat = lat_temp)


  SA <- lapply(1:nrow(grid_combined), function(x) {
    gsw_SA_from_SP(SP=avg_psal_funs[[x]], p=pressure, longitude=grid_combined$long[x],
                   latitude=grid_combined$lat[x]) # absolute salinity
  })
  SA <- lapply(1:length(SA), function(x) {if (length(SA[[x]]) == 0){rep(NA, length(pressure))} else {SA[[x]]}})
  CT <- lapply(1:nrow(grid_combined), function(x) {
    gsw_CT_from_t(SA=SA[[x]], t=avg_temp_funs[[x]], p=pressure) # conservative temperature
  })

  dens <- lapply(1:nrow(grid_combined), function(x) {
    gsw_rho(SA=SA[[x]], CT=CT[[x]], p=pressure) # density
  })
  pdens <- lapply(1:nrow(grid_combined), function(x) { # potential density, reference 0 decibars
    gsw_pot_rho_t_exact(SA=SA[[x]], t=avg_temp_funs[[x]], p=pressure, p_ref=0)
  })
  density_list[[i]] <- list(grid_combined, dens, pdens, CT, SA, removed, avg_psal_funs, avg_temp_funs)
  print(i)
}

grid_used <- do.call(rbind, lapply(density_list, function(x) {return(x[[1]])}))
dens <- unlist(lapply(density_list, function(x) {return(x[[2]])}), recursive = F)
pdens <- unlist(lapply(density_list, function(x) {return(x[[3]])}), recursive = F)
temperature <- unlist(lapply(density_list, function(x) {return(x[[8]])}), recursive = F)
salinity <- unlist(lapply(density_list, function(x) {return(x[[7]])}), recursive = F)

avg_mean <- data.frame(temperature = unlist(temperature), salinity = unlist(salinity),
                       dens = unlist(dens), pdens = unlist(pdens),
                       pressure)

save(grid_used, avg_mean,
     file = 'analysis/results/mean_estimation/mean_review.RData')


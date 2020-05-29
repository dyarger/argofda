
#compute residuals
library(argofda)
library(Matrix)
load('analysis/data/jan_march_data.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)

grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$index <- 1:nrow(grid_to_compute)
# Load mean results
calc_files_psal <- list.files('analysis/results/mean_estimation/', pattern = 'psal')
calc_files_temp <- list.files('analysis/results/mean_estimation/', pattern = 'temp')
final_list <- list()
for (j in 1:length(calc_files_temp)) {
  load(paste0('analysis/results/mean_estimation/', calc_files_temp[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  temp_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[1]])
    }
  })
  eval(parse(text=paste0('rm(', file, ')')))
  load(paste0('analysis/results/mean_estimation/', calc_files_psal[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  psal_funs <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[1]])
    }
  })
  eval(parse(text=paste0('rm(', file, ')')))
  final_list[[j]] <- list(temp_funs,psal_funs)
  print(j)
}
temp_funs <- lapply(final_list, function(x) {return(x[[1]])})
temp_funs <- unlist(temp_funs, recursive = FALSE)
psal_funs <- lapply(final_list, function(x) {return(x[[2]])})
psal_funs <- unlist(psal_funs, recursive = FALSE)

longvals_const <- seq(20.5, 379.5, by = 1)

# for each profile, subtract the mean
for (i in 1:length(profLatAggr)) {
  lat <- profLatAggr[i]
  long <- profLongAggr[i]
  near_lat <- latvals[findInterval(x = lat, vec = -80:80)]
  if (is.na(near_lat) & lat > 80) {
    near_lat <- 79.5
  }
  near_long <- longvals_const[findInterval(x = ifelse(long > 20, long, long +360),
                                         vec = 20:380)]
  near_long <- ifelse(near_long > 180, near_long - 360, near_long)

  index <- which(grid_to_compute$lat == near_lat & grid_to_compute$long == near_long)
  spec_fun_temp <- temp_funs[[index]]
  spec_fun_psal <- psal_funs[[index]]
  dist_EW <- fields::rdist.earth.vec(matrix(nrow = 1, c(long, near_lat)), matrix(nrow = 1, c(near_long, near_lat)), miles = F)
  dist_NS <- fields::rdist.earth.vec(matrix(nrow = 1, c(near_long, lat)), matrix(nrow = 1, c(near_long, near_lat)), miles = F)
  if (long < near_long) {dist_EW = -dist_EW}
  if (lat < near_lat) {dist_NS = -dist_NS}
  temp_pred <- predict_time_year(x = spec_fun_temp,
                                 p_vals = as.vector(profPresAggr[[i]][[1]]), day = profJulFracofYearAggr[i],
                                 dist_EW = dist_EW,
                                 year = profYearAggr[i],
                                 dist_NS = dist_NS)
  if (is.null(spec_fun_psal)) {
  } else {
    psal_pred <- predict_time_year(x = spec_fun_psal,
                                   p_vals = as.vector(profPresAggr[[i]][[1]]), day = profJulFracofYearAggr[i],
                                   dist_EW = dist_EW,
                                   year = profYearAggr[i],
                                   dist_NS = dist_NS)
    profPsalAggr[[i]][[1]] <- profPsalAggr[[i]][[1]] - psal_pred
  }

  profTempAggr[[i]][[1]] <- profTempAggr[[i]][[1]] - temp_pred
  if (i %% 1000 == 0) {
    print(i)
  }
}
save(profCycleNumberAggr,
     profFloatIDAggr,
     profJulDayAggr,
     profLatAggr,
     profLongAggr,
     profModeAggr,
     profMonthAggr,
     profYearAggr,
     profJulDayofYearAggr,
     profJulFracofYearAggr,
     profIndexAggr,profPresAggr, profTempAggr, profPsalAggr,
     file = 'analysis/data/jan_march_residuals.RData')

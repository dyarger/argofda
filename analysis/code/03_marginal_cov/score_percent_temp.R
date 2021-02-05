library(argofda)
library(Matrix)
library(fda)
library(sparseinv)
library(Rcpp)
library(dplyr)
library(corpcor)
library(fields)
library(parallel)
setwd('/pylon5/msz3a1p/d7yar/argofda/')
source('analysis/code/04_joint_TS_20/ml_source_vec.R')
load('analysis/results/temp_cov_pca.RData')

load('analysis/data/jan_march_residuals.RData')
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)

longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$order <- 1:nrow(grid_to_compute)

cv <- function(index, h_space) {
  a <- proc.time()

  lat <- grid_to_compute[index,2]
  long <- grid_to_compute[index,1]

  df <- get_profile_data_subset(long = long, lat = lat, day = 45.25, h_time = 45.25, RG_def = FALSE,
                                h_space = h_space, mode = 'all',  min_prof=  5,
                                exclusion= TRUE)
  df$pressure <- round(as.numeric(df$pressure), digits = 2)

  if (!(lat %in% grid_pc_temp$lat[grid_pc_temp$long == long])) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }

  if (is.null(df) | is.null(nrow(df)) ) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }
  pressure_min <- df %>% group_by(profile) %>% summarise(min(pressure))
  pressure_max <- df %>% group_by(profile) %>% summarise(max(pressure))
  short_profiles <- (1:nrow(pressure_max))[which(pressure_max[,2] - pressure_min[,2] < 800)]
  prof_lengths <- df %>% group_by(profile) %>% summarise(n())
  small_profiles <-(1:nrow(prof_lengths))[which(prof_lengths[,2] < 15)]
  df <- df[!(df$profile %in% short_profiles) & !(df$profile %in% small_profiles),]
  df$profile <- as.numeric(as.factor(df$profile))

  if (is.null(df) | is.null(nrow(df))  |nrow(df) == 0) {
    b <- proc.time()
    final_list <- c(index, 'no profiles in range', (b-a)[3])
    print(final_list)
    return(final_list)
  }

  temp_pc <- coefs_pc_temp[[which(grid_pc_temp$long == long &
                                    grid_pc_temp$lat == lat)]][[1]][,1:10]

  # get first few pc at this location
  # Calculate scores for each profile
  scores_5 <- unlist(get_scores(df, temp_pc[,1:5], K = 5, type = 'temperature',
                                knots = knots_pc)[[2]])
  scores_10 <- unlist(get_scores(df, temp_pc[,1:10], K = 10, type = 'temperature',
                                 knots = knots_pc)[[2]])
  per05 <- 1- mean(scores_5^2)/mean(df$temperature^2)
  per10 <- 1- mean(scores_10^2)/mean(df$temperature^2)
  cut_list <- list(cut(df$pressure, breaks = c(seq(0, 200, by = 10), seq(300, 2000, by = 100))))
  temp_aggr <- aggregate(df$temperature^2,cut_list, mean)[,2]
  per05_break <- aggregate(scores_5^2, cut_list, mean)
  per05_break$x <- 1 - per05_break$x/(temp_aggr)

  per10_break <- aggregate(scores_10^2, cut_list, mean)
  per10_break$x <- 1 - per10_break$x/(temp_aggr)

  b <- proc.time()
  final_list <- list(c(per05, per10),
                     per05_break,
                     per10_break,
                     c(index, (b-a)[3]
                     ))
  print(final_list[[4]])
  return(final_list)
}
# sequence <- as.list(seq(1, nrow(grid_to_compute), by = 1))
# array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# increments <- 500
# h_space <- 750
#
# profile_start <- 1 + (array_id - 1) * increments
# indexes <- profile_start:(profile_start+increments-1)
#
# if (array_id == ceiling(nrow(grid_to_compute)/increments)) {
#   indexes <- indexes[sapply(sequence[indexes], function(x) {!is.null(x)})]
# }
#
# # Run for a large number of grid points
# if (array_id == ceiling(nrow(grid_to_compute)/increments)) {
#   value <- sequence[[indexes[length(indexes)]]]-1
# } else {
#   value <- sequence[[indexes[length(indexes)]+1]]-1
# }
#object_name <-  paste('Grid_Pred', value, sep = '')
indexes <- 1:nrow(grid_to_compute)
h_space <- 500
result <- mclapply(indexes, cv,
                  mc.preschedule = FALSE, mc.cores = 10,
                 h_space = h_space)
# result <-lapply(indexes, cv,
#                  # mc.preschedule = FALSE, mc.cores = 10,
#                   h_space = h_space)
save(result, file = 'analysis/results/joint_TS_20/score_percent_temp.Rdata')
#file.to.save <- paste0('analysis/results/joint_TS_20/score_percent_temp.Rdata')
#eval(parse(text=paste('save(', result, ', file = file.to.save)', sep = '')))

q()

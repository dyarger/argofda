
# implement grid
library(argofda)
library(Matrix)
library(fda)
library(sparseinv)
library(parallel)
library(dplyr)
library(functionalCovariance)
setwd('/home/dyarger/argofda/')
source('analysis/code/03_marginal_cov/covariance_source.R')

load('analysis/data/jan_march_residuals.RData')
rm(profPsalAggr)
profLongAggr <- ifelse(profLongAggr > 180, profLongAggr - 360, profLongAggr)

load('analysis/data/RG_Defined.RData')
latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute$order <- 1:nrow(grid_to_compute)
grid_to_compute_merge <- left_join(grid_to_compute, RG_defined_long)
grid_to_compute_merge$RG_defined <- ifelse(is.na(grid_to_compute_merge$value),
                                           FALSE,
                                           TRUE)
grid_to_compute_merge <- grid_to_compute_merge[order(grid_to_compute_merge$order), c(1,2,5)]

knots <- seq(0, 2000, length.out = 100)
basis <- create.bspline.basis(knots)
penalty_mat <- bsplinepen(basis)

sequence <- as.list(seq(1, nrow(grid_to_compute), by = 1))
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
increments <- 500
profile_start <- 1 + (array_id - 1) * increments
indexes <- profile_start:(profile_start+increments-1)

if (array_id == ceiling(nrow(grid_to_compute)/500)) {
  indexes <- indexes[sapply(sequence[indexes], function(x) {!is.null(x)})]
}

# Run for a large number of grid points
if (array_id == ceiling(nrow(grid_to_compute)/500)) {
  value <- sequence[[indexes[length(indexes)]]]-1
} else {
  value <- sequence[[indexes[length(indexes)]+1]]-1
}
object_name <-  paste('Grid_Pred', value, sep = '')

command <- paste(object_name, ' <- mclapply(indexes, cov_est_current,',
            #     'mi_max = 50,',
                 'h_space = 600, knots = knots,  min_prof = 1,',
                 'lambda = 3, basis = basis, penalty_mat = penalty_mat,',
                 'mc.preschedule = FALSE, mc.cores = 4',
                 ')', sep="")
a <- proc.time()
eval(parse(text=command))
b <- proc.time();b-a

if (nchar(value) == 3) {
  value <- paste0('00',value)
} else if(nchar(value) == 4) {
  value <- paste0('0',value)
}
file.to.save <- paste('analysis/results/marginal_cov/temp_',
                      value, '.Rdata', sep = '')
eval(parse(text=paste('save(', object_name, ', file = file.to.save)', sep = '')))
q()

library(fda)
library(ggplot2)
library(RColorBrewer)
library(argofda)

calc_files <- list.files('analysis/results/nugget_variance/')
# Get results from score modelling
final_list <- list()
for (j in 1:length(calc_files)) {
  load(paste0('analysis/results/nugget_variance/', calc_files[j]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file, '')))
  eval(parse(text=paste0('rm(', file, ')')))

  indexes_used <- as.numeric(sapply(now, function(x) {
    if (length(x) == 3) {
      x[1]
    } else {
      x[[4]][1]
    }
  }))
  temp_nugget_funs <- lapply(1:length(now), function(x) {
    if (length(now[[x]]) == 1 | length(now[[x]]) == 0 | length(now[[x]]) == 3) {
      NULL
    } else {
      now[[x]][[1]]
    }
  })
  psal_nugget_funs <- lapply(1:length(now), function(x) {
    if (length(now[[x]]) == 1 | length(now[[x]]) == 0 | length(now[[x]]) == 3) {
      NULL
    } else {
      now[[x]][[2]]
    }
  })
  prof_info <- lapply(1:length(now), function(x) {
    if (length(now[[x]]) == 1 | length(now[[x]]) == 0 | length(now[[x]]) == 3) {
      NULL
    } else {
      now[[x]][[4]]
    }
  })
  final_list[[j]] <- list(temp_nugget_funs, psal_nugget_funs, prof_info)
}


temp_nugget_funs <- unlist(lapply(final_list, function(x) return(x[[1]])), recursive = F)
psal_nugget_funs <- unlist(lapply(final_list, function(x) return(x[[2]])), recursive = F)
prof_info <- do.call(rbind, unlist(lapply(final_list, function(x) return(x[[3]])), recursive = F))
have_results <- which(sapply(temp_nugget_funs, length) != 0 & sapply(temp_nugget_funs, length) != 0)
summary(prof_info)

temp_nugget_funs <- temp_nugget_funs[have_results]
psal_nugget_funs <- psal_nugget_funs[have_results]

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 180, longvals - 360, longvals)
grid_nugget <- expand.grid(long = longvals, lat = latvals)
grid_nugget <- grid_nugget[have_results,]

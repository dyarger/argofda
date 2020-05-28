

calc_files <- list.files('analysis/results/joint_TS_20/', pattern = 'cv')
final_list <- list()
for (j in 1:length(calc_files)) {
  load(paste0('analysis/results/joint_TS_20/', calc_files[j]))
  file <- ls(pattern = 'cv_results')
  eval(parse(text=paste0('now <- ', file, '')))
  table(sapply(now, length))

  n_delayed <- sapply(now, function(x) {
    if (length(x) == 7) {
      if (is.null(x[[6]])) {
        return(0)
      }
      return(x[[6]])
    } else {
      return(0)
    }
  })
  no_pred <- sapply(now, function(x) is.null(x[[1]]))
  mode <- sapply(now, function(x) x[[5]])
  residuals <- lapply(now, function(x) x[[3]])
  temp_results <- lapply(now, function(x) x[[1]])
  psal_results <- lapply(now, function(x) x[[2]])
  profile_nums <- sapply(now, function(x) x[[7]][1])
  info <- lapply(now, function(x) x[[7]][2])


  cond_var_old_scores <- lapply(now, function(x) x[[4]])
  final_list[[j]] <- list(mode, NULL, temp_results, psal_results,
                          cond_var_old_scores, NULL,profile_nums,NULL,n_delayed,
                          no_pred, info)
  print(j)
}

mode <- unlist(sapply(final_list, function(x) x[[1]]))
temp_results <- unlist(lapply(final_list, function(x) x[[3]]), recursive = F)
psal_results <- unlist(lapply(final_list, function(x) x[[4]]), recursive = F)
profile_nums <- sapply(final_list, function(x) x[[7]])
profile_nums <- as.numeric(unlist(profile_nums, recursive = T))
cond_var <- unlist(sapply(final_list, function(x) x[[5]]), recursive = F)

prof_info <- unlist(lapply(final_list, function(x) return(x[[11]])), recursive = F)
prof_info <- do.call(c, prof_info)
summary(as.numeric(prof_info))


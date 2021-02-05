# Cov results

library(ggplot2)
library(viridis)
library(fda)
theme_set(theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 15)))

type <- 'psal'
label <- c('Temperature', 'Salinity')[c('temp', 'psal') == type]
label_units <- c('Temp (°C)', 'PSU')[c('temp', 'psal') == type]
cov_files <- list.files('analysis/results/marginal_cov/', pattern = paste0(type))
cov_files <- cov_files[nchar(cov_files) < 20]
if (substr(cov_files[1], 1,4) == 'lags') {
  cov_files <- cov_files[-1]
}
#cov_files <- list.files('analysis/results/marginal_cov/', pattern = paste0(type, '_complete_'))

funs <- list()
cv_info <- list()
prof_info <- list()
for (i in 1:length(cov_files)) {
  load(paste0('analysis/results/marginal_cov/', cov_files[i]))
  file <- ls(pattern = 'Grid_Pred')
  eval(parse(text=paste0('now <- ', file)))
  eval(parse(text=paste0('rm(', file, ')')))

  funs[[i]] <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[1]])
    }
  })
  cv_info[[i]] <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[2]])
    }
  })
  prof_info[[i]] <- lapply(now, function(x) {
    if (x[[2]][[1]][[1]][1] == 'no profiles in range') {
      return(NULL)
    } else {
      return(x[[3]])
    }
  })
}
funs <- unlist(funs, recursive = F)
beta <- lapply(funs, function(x) {
  if (is.null(x)) {return(NULL)} else {return(x)}
})
have_results <- !sapply(beta, is.null)
beta <- beta[have_results]
knots <- funs[have_results][[1]][[2]]
knots <- seq(0, 2000, length.out = 100)
cv_info <- unlist(cv_info, recursive = F)
cv_niter <- sapply(cv_info[have_results], nrow)
summary(cv_niter)

prof_info <- unlist(prof_info, recursive = F)
prof_info <- prof_info[have_results]
prof_info <- do.call(rbind, prof_info)
colnames(prof_info) <- c('index', 'orig_samples', 'seconds', 'profiles', 'measurements', 'profiles2',
                         'h_space')
summary(prof_info)

latvals <- seq(-79.5, 79.5, by = 1)
longvals <- seq(20.5, 379.5, by = 1)
longvals <- ifelse(longvals > 360, longvals - 360, longvals)
grid_to_compute <- expand.grid('long' = longvals, 'lat' = latvals)
grid_to_compute_only <- grid_to_compute[have_results,]

#### Plots ####
vals <- substr(cov_files, 15, 19)
plot(as.numeric(vals))
# correlation at two different pressures

save_cor_image <- function(p1, p2) {
  basis <- create.bspline.basis(breaks = knots)
  B <- eval.basis(c(p1, p2), basis)
  cor_val <- sapply(beta, function(x) {
   # beta_mat <- do.call(cbind, x)
    A_mat = matrix(0, length(knots)+2, length(knots)+2)
    A_mat[lower.tri(A_mat, diag=TRUE)] = x
    beta_mat = as.matrix(forceSymmetric(A_mat, 'L'))
    (t(B[1,]) %*% beta_mat %*% B[2,])/sqrt((t(B[1,]) %*% beta_mat %*% B[1,]) *
                                             (t(B[2,]) %*% beta_mat %*% B[2,]))
  })
  cor_val <- ifelse(cor_val < -1, -1, ifelse(cor_val > 1, 1, cor_val))
  im <- ggplot(data = cbind(grid_to_compute_only, cor_val))+
    geom_raster(aes(x =  long, y = lat,
                    fill= cor_val))+
    scale_fill_gradientn(limits = c(-1, 1), colours = c('blue', 'white', 'red'))+
    geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
                 color = 'black', fill = 'white', size = .2)+
    labs(x = 'Longitude', y = 'Latitude', fill = 'Cor')
  print(im)
  ggsave(plot = im,
         filename = paste0('analysis/images/marginal_cov/cor_',p1, '_', p2, '_', type, '.png'),
         scale = .8,height = 4, width = 7.25)
  ggsave(plot = im,
         filename = paste0('analysis/images/marginal_cov/cor_',p1, '_', p2, '_', type, '_eps.eps'),
         scale = .8,height = 4, width = 7.25)
}

save_cor_image(10, 800)
save_cor_image(800, 1200)
save_cor_image(1200, 1600)
save_cor_image(10, 200)
save_cor_image(10, 100)
save_cor_image(10, 400)

####### Make and save first 15 FPCs and all eigenvalues #####
basis <- create.bspline.basis(breaks = knots)
omega_0 <- as(bsplinepen(basis,Lfdobj = 0), 'sparseMatrix')
omega_0_sqrt <- as(chol(omega_0), 'sparseMatrix')

a <- proc.time()
coefs_pc <- lapply(beta, function(x) {
  A_mat = matrix(0, length(knots)+2, length(knots)+2)
  A_mat[lower.tri(A_mat, diag=TRUE)] = x
  A_mat = as.matrix(forceSymmetric(A_mat, 'L'))
  #A_mat <- matrix(unlist(x), nrow = length(knots) + 2)
  pca_mat <- omega_0_sqrt %*% A_mat %*% t(omega_0_sqrt)
  standard_pca <- prcomp(pca_mat, center = FALSE, scale = FALSE)
  vectors <- solve(omega_0_sqrt, standard_pca$rotation)
  return(list(vectors[,1:15], standard_pca$sdev))
})
b <- proc.time()
b-a
index <- sample(1:length(coefs_pc), 1)
plot(coefs_pc[[index]][[1]][,1])
plot(coefs_pc[[index]][[1]][,2])
plot(coefs_pc[[index]][[1]][,3])

grid_to_compute_only$long <- ifelse(grid_to_compute_only$long > 180,
                                    grid_to_compute_only$long - 360,
                                    grid_to_compute_only$long)
knots_pc <- knots
grid_pc <- grid_to_compute_only
if (type == 'temp') {
  coefs_pc_temp <- coefs_pc
  grid_pc_temp <- grid_pc
  save(coefs_pc_temp, grid_pc_temp, knots_pc, file = paste0( 'analysis/results/',
                                                   type,
                                                   '_cov_pca.RData'))
} else if (type == 'psal') {
  coefs_pc_psal <- coefs_pc
  grid_pc_psal <- grid_pc
  save(coefs_pc_psal, grid_pc_psal, knots_pc, file = paste0( 'analysis/results/',
                                                   type,
                                                   '_cov_pca.RData'))
}
save(coefs_pc_psal, grid_pc_psal, knots_pc, file = paste0( 'analysis/results/',
                                                           type,
                                                           '_cov_pca_v2.RData'),
     version = 2)
grid_to_compute_only$long <- ifelse(grid_to_compute_only$long < 0,
                                    grid_to_compute_only$long + 360,
                                    grid_to_compute_only$long)

# Number of components needed to explain variance
# or amount of variance explained by first K components
# these are just estimates

eigen_pc <- t(sapply(coefs_pc, function(x) {return(x[[2]])}))

save_eigen_image <- function(pcs_needed_prop = .95, num_pcs = NULL, max_eigen = 7) {
  if (is.null(pcs_needed_prop)) {
    first_comp <- apply(eigen_pc, 1, function(x) {sum(x[1:num_pcs]^2)/sum(x^2)})
    a <- ggplot(data = data.frame(grid_to_compute_only, first_comp))+
      geom_raster(aes(x = long, y = lat, fill = first_comp))+
      scale_fill_viridis()+
      geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),
                   fill = 'white', color = 'black', size = .2)+
      labs(x = 'Longitude', y = 'Latitude', fill = '%',
           title = paste('Percent of variance from first', num_pcs, 'components, temperature'))+
      theme(text = element_text(size= 15))
    print(a)
    ggsave(plot = a,
           filename = paste0('analysis/images/marginal_cov/eigen_', num_pcs, '_pc_temp.png'),
           scale = .8,height = 4, width = 7.25)
  } else {
    pc_needed <- apply(eigen_pc, 1, function(x) {which(cumsum(x^2)/
                                                           sum(x^2) > pcs_needed_prop)[1]})
    a <- ggplot(data = data.frame(grid_to_compute_only, pc_needed))+
      geom_raster(aes(x = ifelse(long < 0, long + 360, long),
                      y = lat, fill = as.factor(pc_needed)))+
      scale_fill_viridis_d(limits = as.character(seq(1, max_eigen)))+
      geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
                   fill = 'white', color = 'black', size = .2)+
      labs(x = 'Longitude', y = 'Latitude', fill = '# PC')+
      theme(text = element_text(size= 15))
    print(a)
    ggsave(plot = a,
           filename = paste0('analysis/images/marginal_cov/eigen_', 100*pcs_needed_prop, '_percent_temp.png'),
           scale = .8,height = 4, width = 7.25)
    ggsave(plot = a,
           filename = paste0('analysis/images/marginal_cov/eigen_', 100*pcs_needed_prop, '_percent_temp_eps.eps'),
           scale = .8,height = 4, width = 7.25)
  }
}
save_eigen_image( folder = images_folder, max_eigen = 6)

## Plot of lambda chosen
lambda_pcs <- as.numeric(unlist(sapply(funs, function(x) {
  return(x[[4]])
})))
a <- ggplot(data = data.frame(grid_to_compute_only[1:length(lambda_pcs),], lambda_pcs))+
  geom_raster(aes(x=ifelse(long < 0, long + 360, long), y = lat, fill = log10(lambda_pcs)))+
  geom_polygon(data= map_data('world2'), aes(x = long, y = lat, group = group),
               fill = 'white', color = 'black', size = .2)+
  scale_fill_viridis(limits = c(-1, 4.5))+
  theme_bw()+
  labs(x = 'Longitude', y = 'Latitude', fill = 'log(lambda)',
       title = paste0('Lambda Chosen for Covariance Estimation, ', label))
a
ggsave(plot = a, filename =paste0('analysis/images/marginal_cov/covariance_lambda_', type, '.png'), scale = 1.2,
       height = 5, width = 7.25)


## example
theme_set(theme_bw() + theme(text = element_text(size = 15)))
long_example <- 90.5
lat_example <- -10.5
basis <- create.bspline.basis(breaks = knots)
omega_0 <- as(bsplinepen(basis,Lfdobj = 0), 'sparseMatrix')
omega_0_sqrt <- as(chol(omega_0), 'sparseMatrix')


A_mat = matrix(0, length(knots)+2, length(knots)+2)
A_mat[lower.tri(A_mat, diag=TRUE)] = unlist(beta[[which(grid_to_compute_only$long == long_example &
                                                          grid_to_compute_only$lat == lat_example)]])
A_mat = as.matrix(forceSymmetric(A_mat, 'L'))

pca_mat <- omega_0_sqrt %*% A_mat %*% t(omega_0_sqrt)
standard_pca <- prcomp(pca_mat, center = FALSE, scale = FALSE)
vectors <- solve(omega_0_sqrt, standard_pca$rotation)

pressure_vec <- 0:2000

pred_pc1 <- predict(argofda::make_fit_object(knots = knots, coefficients = vectors[,1]), x = pressure_vec)$y
ggplot(data = data.frame(pressure_vec, pred_pc1))+
  geom_line(aes(x = pressure_vec, y = pred_pc1))+
  geom_hline(yintercept = 0, alpha = .5)+
  coord_flip()+
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = label_units, title = '1st PC')
ggsave(filename = paste0('analysis/images/marginal_cov/pc_example_1_', type, '.png'),
       scale = .8,height = 8, width = 4)

pred_pc2 <- predict(argofda::make_fit_object(knots = knots, coefficients = vectors[,2]), x = pressure_vec)$y
ggplot(data = data.frame(pressure_vec, pred_pc2))+
  geom_line(aes(x = pressure_vec, y = pred_pc2))+
  geom_hline(yintercept = 0, alpha = .5)+
  coord_flip()+
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = label_units, title = '2nd PC')
ggsave(filename = paste0('analysis/images/marginal_cov/pc_example_2_', type, '.png'),
       scale = .8,height = 8, width = 4)

## save correlations at 50 decibars intervals and different lags

make_lags <- function(lag) {
  pressure_check_1 <- seq(0, 2000-lag, by = 50)
  pressure_check_2 <- pressure_check_1+ lag
  basis <- create.bspline.basis(breaks = knots_pc)
  B_1 <- eval.basis(c(pressure_check_1, 2000), basis)
  B_2 <- eval.basis(c(pressure_check_2, 2000), basis)
  cor_lag <- sapply(beta, function(x) {

   # beta_mat <- do.call(cbind, x)

    A_mat = matrix(0, length(knots)+2, length(knots)+2)
    A_mat[lower.tri(A_mat, diag=TRUE)] = unlist(x)
    beta_mat = as.matrix(forceSymmetric(A_mat, 'L'))

    sapply(1:length(pressure_check_1), function(z) {
      ((t(B_1[z,]) %*% beta_mat %*% B_2[z,])/sqrt((t(B_1[z,]) %*% beta_mat %*% B_1[z,]) *
                                                    (t(B_2[z,]) %*% beta_mat %*% B_2[z,])))[1]
    })
  })
  lags <- data.frame(grid_to_compute_only, t(cor_lag))
  colnames(lags) <- c('long', 'lat', pressure_check_1)
  lags_long <- tidyr::gather(lags, pressure, value, -c(1,2))

  lags_summary_mean <- aggregate(lags_long$value[lags_long$value <=1 & lags_long$value >=-1],
                                 list(lags_long$pressure[lags_long$value <=1 & lags_long$value >=-1]),
                                 mean)
  lags_summary_med <- aggregate(lags_long$value[lags_long$value <=1 & lags_long$value >=-1],
                                list(lags_long$pressure[lags_long$value <=1 & lags_long$value >=-1]),
                                median)
  lags_summary <- cbind(lags_summary_mean, lags_summary_med[,2])
  colnames(lags_summary) <- c('pressure', 'Mean', 'Median')
  lags_summary$pressure <- as.numeric(as.character(lags_summary$pressure))
  tidyr::gather(lags_summary, Type, value, -1)
}
knots_pc <- knots
lag_010 <- make_lags(10)
lag_025 <- make_lags(25)
lag_050 <- make_lags(50)
lag_100 <- make_lags(100)
lag_200 <- make_lags(200)
lag_500 <- make_lags(500)
save(lag_010, lag_025, lag_050, lag_100, lag_200, lag_500,
     file = paste0('analysis/results/marginal_cov/lags_', type, '.RData'))
theme_set(theme_bw() + theme(text = element_text(size = 15)))

lag_010 <- data.frame(lag_010, 'lag' = 10 )
lag_025 <- data.frame(lag_025, 'lag' = 25)
lag_050 <- data.frame(lag_050, 'lag' = 50 )
lag_100 <- data.frame(lag_100, 'lag' = 100 )
lag_200 <- data.frame(lag_200, 'lag' = 200 )
lag_500 <- data.frame(lag_500, 'lag' = 500 )

lags <- rbind(lag_010, lag_025, lag_050, lag_100, lag_200, lag_500)


lag_values <- c(10, 25, 50, 100, 200, 500)
a <- ggplot(data = lags) +
  geom_jitter(aes(x = as.factor(lag), y = value), height = 0, size = .6)+
  facet_wrap(~Type) +
  geom_point(data = data.frame(lag = as.factor(lag_values),
                               cor =  exp(lag_values  * log(.999))),
             aes(x = lag, y = cor), col = 'red')+
  scale_y_continuous(limits = c(0, 1))+
  labs(x= 'Lag (decibars)', y = 'Correlation between measurements\n at p and p+lag')
a
ggsave(plot = a,
       filename = paste0('analysis/images/marginal_cov/cor_incr_', type, '_complete.png'),
       scale = 1.2,height = 4, width = 7.25)
ggsave(plot = a,
       filename = paste0('analysis/images/marginal_cov/cor_incr_', type, '_complete.eps'),
       scale = 1.2,height = 4, width = 7.25)


# check signs

load(file = paste0( 'analysis/results/temp_cov_pca.RData'))

# implement sign switching
inprod_mat <- as(fda::eval.penalty(basisobj = create.bspline.basis(knots_pc), Lfdobj = 0) ,
                 'sparseMatrix')
which(grid_pc_temp$long == 90.5 &
        grid_pc_temp$lat == -50.5)
beta_test <- list(coefs_pc_temp[[27072+45]][[1]][,1],
                  coefs_pc_temp[[27003+45]][[1]][,1],
                  coefs_pc_temp[[8735-50]][[1]][,1])
grid_use <- grid_pc_temp[c(8735 -50, 27003+45, 27072+45),]
sign_val <- sapply(1:nrow(grid_pc_temp), function(x) {
  if (length(beta_test) ==0) {
    return(1)
  }
  tests <- sapply(beta_test, function(y) {
    sign((t(coefs_pc_temp[[x]][[1]][,1]) %*% inprod_mat %*% y)[1])
  })
  sign(sum(tests))
})

pc1_val <- sapply(1:length(coefs_pc_temp), function(x) {
  coefs_pc_temp[[x]][[1]][5,1] * sign_val[x]
})

pc1_val_old <- sapply(1:length(coefs_pc_temp), function(x) {
  coefs_pc_temp[[x]][[1]][5,1]
})

im <- ggplot()+
  geom_raster(data = cbind(grid_pc_temp, pc1_val),aes(x =  ifelse(long < 0, long+  360, long), y = lat,
                  fill= pc1_val))+
  scale_fill_gradient2(low = c('blue'), high = c( 'red'))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  geom_point(data = grid_use, aes(x =  ifelse(long < 0, long+  360, long), y = lat))+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Coef')
print(im)
ggsave(plot = im,
       filename = paste0('analysis/images/marginal_cov/temp_1_5_aligned.png'),
       scale = .8,height = 4, width = 7.25)


im <- ggplot()+
  geom_raster(data = cbind(grid_pc_temp, pc1_val_old),aes(x =  ifelse(long < 0, long+  360, long), y = lat,
                                                      fill= pc1_val_old))+
  scale_fill_gradient2(low = c('blue'), high = c( 'red'))+
  geom_polygon(data = map_data('world2'), aes(x = long, y = lat, group = group),
               color = 'black', fill = 'white', size = .2)+
  labs(x = 'Longitude', y = 'Latitude', fill = 'Coef')
ggsave(plot = im,
       filename = paste0('analysis/images/marginal_cov/temp_1_5_unaligned.png'),
       scale = .8,height = 4, width = 7.25)

print(im)
grid_pc_temp[14550,]
grid_pc_temp[8750,]


#sign_val_lr <-  sign_val

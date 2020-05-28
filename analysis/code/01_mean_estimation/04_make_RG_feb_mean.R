
# Make RG residuals

files <- list.files('analysis/products/gilson_files/', pattern = 'fil_temp_pa99')

long <- seq(20 + 1/12, 380 - 1/12, by = 1/6)
lat <- seq(-90 + 1/12, 90 - 1/12, by = 1/6)
grid <- expand.grid('long' = long, 'lat' = lat)

coef <- list()
for (j in 1:length(files)) {
  file <- files[j]
  mean_file <- read.delim(paste0('analysis/products/gilson_files/',file), sep = '', header = F)
  colnames(mean_file) <- c('mean', 1:12)
  print(file)
  pres_vals <- c(substr(file, 21,24),substr(file, 26,29))
  mean_file[mean_file == '-99.999'] <- NA
  feb_mean <- mean_file[,3]
  RG_mean <- cbind(grid, feb_mean)
  RG_mean <- RG_mean[abs(RG_mean$long %% 1 - 0.41666667) < .01 & abs(RG_mean$lat %% 1 - 0.41666667) < .01,]
  coef[[j]] <- RG_mean
}
save(coef, file = 'analysis/products/RG_temp_coef.RData')

# same for salinity
files <- list.files('analysis/products/gilson_files/', pattern = 'fil_psal_pa99')

long <- seq(20 + 1/12, 380 - 1/12, by = 1/6)
lat <- seq(-90 + 1/12, 90 - 1/12, by = 1/6)
grid <- expand.grid('long' = long, 'lat' = lat)

coef <- list()
for (j in 1:length(files)) {
  file <- files[j]
  mean_file <- read.delim(paste0('panalysis/roducts/gilson_files/',file), sep = '', header = F)
  colnames(mean_file) <- c('mean', 1:12)
  print(file)
  pres_vals <- c(substr(file, 21,24),substr(file, 26,29))
  mean_file[mean_file == '-99.999'] <- NA
  feb_mean <- mean_file[,3]
  RG_mean <- cbind(grid, feb_mean)
  RG_mean <- RG_mean[abs(RG_mean$long %% 1 - 0.41666667) < .01 & abs(RG_mean$lat %% 1 - 0.41666667) < .01,]
  coef[[j]] <- RG_mean
}
save(coef, file = 'analysis/products/RG_psal_coef.RData')


# make RG mean
mfpres <- c(2.5,   10.0,  20.0 ,  30.0,   40.0  , 50.0   ,
            60.0,   70.0,   80.0,   90.0,  100.0,  110.0,  120.0,  130.0,
            140.0,  150.0,  160.0,  170.0,  182.5,
            200.0,  220.0,  240.0,  260.0,  280.0,  300.0,  320.0,  340.0,  360.0,
            380.0,  400.0,  420.0,  440.0,  462.5,
            500.0,  550.0,  600.0,  650.0,  700.0,  750.0,  800.0,  850.0,  900.0,
            950.0 ,1000.0, 1050.0, 1100.0, 1150.0,
            1200.0, 1250.0, 1300.0, 1350.0, 1412.5, 1500.0, 1600.0, 1700.0, 1800.0,
            1900.0 ,1975.0)
load('analysis/products/RG_temp_coef.RData')
coef_temp <- lapply(1:length(coef), function(x) {cbind(coef[[x]], pressure=  mfpres[x])})
coef_temp <- do.call(rbind, coef_temp)
load('analysis/products/RG_psal_coef.RData')
coef_psal <- lapply(1:length(coef), function(x) {cbind(coef[[x]], pressure = mfpres[x])})
coef_psal <- do.call(rbind, coef_psal)

RG_mean <- data.frame(coef_temp, feb_mean_psal = coef_psal$feb_mean)
colnames(RG_mean) <- c('long', 'lat', 'temp_RG', 'pressure','psal_RG')
RG_mean$long <- ifelse(RG_mean$long > 360, RG_mean$long -360, RG_mean$long)

library(gsw)
SA <- gsw_SA_from_SP(SP=RG_mean$psal_RG, p=RG_mean$pressure,
                     longitude=RG_mean$long,
                     latitude=RG_mean$lat)
CT <- gsw_CT_from_t(SA=SA, t=RG_mean$temp_RG, p=RG_mean$pressure)
dens <- gsw_rho(SA=SA, CT=CT, p=RG_mean$pressure)
pdens <-  gsw_pot_rho_t_exact(SA=SA, t = RG_mean$temp_RG,
                              p=RG_mean$pressure, p_ref=0)
mean_field <- data.frame(longitude = RG_mean$long,
                         latitude = RG_mean$lat,
                         pressure = RG_mean$pressure,
                         salinity = RG_mean$psal_RG,
                         density = dens,  pdensity = pdens,
                         temperature = RG_mean$temp_RG)
save(mean_field, file = 'analysis/products/RG_mean_mid_feb.RData')


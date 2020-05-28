
# Process Data

# In matlab, change version of Kuusela's Github data so that R.matlab package can handle

# load('Argo_data_aggr.mat')
# save('Argo_data_aggr_v6.mat','endYear','nProf','profCycleNumberAggr','profFloatIDAggr','profJulDayAggr','profLatAggr','profLongAggr','profModeAggr','profMonthAggr','profPresAggr','profPsalAggr','profTempAggr','profYearAggr','startYear','-v6');

# Load into R and change into .RData file
library(R.matlab)
argo_data <- R.matlab::readMat('analysis/data/Argo_data_aggr_v6.mat')

endYear <- argo_data$endYear
nProf <- argo_data$nProf
profCycleNumberAggr <- argo_data$profCycleNumberAggr
profFloatIDAggr <- argo_data$profFloatIDAggr
profJulDayAggr <- argo_data$profJulDayAggr
profLatAggr <- argo_data$profLatAggr
profLongAggr <- argo_data$profLongAggr
profModeAggr <- argo_data$profModeAggr
profMonthAggr <- argo_data$profMonthAggr
profPresAggr <- argo_data$profPresAggr
profPsalAggr <- argo_data$profPsalAggr
profTempAggr <- argo_data$profTempAggr
profYearAggr <- argo_data$profYearAggr
startYear <- argo_data$startYear
save(file = 'analysis/data/Argo_data_aggr.RData', endYear, nProf, profCycleNumberAggr, profFloatIDAggr, profJulDayAggr,
     profLatAggr, profLongAggr, profModeAggr, profMonthAggr, profPresAggr,
     profPsalAggr, profTempAggr, profYearAggr, startYear)

# With .RData version of data, change a few things, then resave February
load('analysis/data/Argo_data_aggr.RData')

profCycleNumberAggr <- as.vector(profCycleNumberAggr)
profFloatIDAggr <- as.vector(profFloatIDAggr)
profJulDayAggr <- as.vector(profJulDayAggr)
profLatAggr <- as.vector(profLatAggr)
profLongAggr <- as.vector(profLongAggr)
profModeAggr <- as.vector(profModeAggr)
profMonthAggr <- as.vector(profMonthAggr)
profYearAggr <- as.vector(profYearAggr)

profModeAggr <- strsplit(profModeAggr, split  =  '')[[1]]

# Create day of year variables (2 of them)

# profJulDayofYearAggr denotes the number of days since the beginning of the year
# takes values 0-366 for non leap years, 0-365 for leap years
# profJulFracofYearAggr denotes the number of days since the beginning of the year

make_day_of_year <- function(x) {
  begin_JD <- 733042
  offset_JD <- x - begin_JD
  # 2008, 2012, and 2016 are leap years
  year_boundaries <- c(0,365,366,365,365,365,366,365,365,365)
  day_boundaries <- cumsum(year_boundaries)+1
  next_min <- findInterval(x = offset_JD, vec = day_boundaries)
  profJulDayofYearAggr <- offset_JD - day_boundaries[next_min]

  day_boundaries_frac <- 365.25*(0:10)+1 -365.25*2 + 365*2
  next_min_frac <- findInterval(x = offset_JD, vec = day_boundaries_frac)
  profJulFracofYearAggr <- offset_JD - day_boundaries_frac[next_min_frac]

  return(list(profJulDayofYearAggr,profJulFracofYearAggr))
}

profJulDayofYearAggr <- make_day_of_year(profJulDayAggr)[[1]]
profJulFracofYearAggr <- make_day_of_year(profJulDayAggr)[[2]]

profCycleNumberAggr<-profCycleNumberAggr[profMonthAggr%in%1:3]
profFloatIDAggr<-profFloatIDAggr[profMonthAggr%in%1:3]
profJulDayAggr<-profJulDayAggr[profMonthAggr%in%1:3]
profLatAggr<-profLatAggr[profMonthAggr%in%1:3]
profLongAggr<-profLongAggr[profMonthAggr%in%1:3]
profModeAggr<-profModeAggr[profMonthAggr%in%1:3]
profYearAggr<-profYearAggr[profMonthAggr%in%1:3]
profJulDayofYearAggr<-profJulDayofYearAggr[profMonthAggr%in%1:3]
profJulFracofYearAggr<-profJulFracofYearAggr[profMonthAggr%in%1:3]
profPresAggr <- profPresAggr[profMonthAggr%in%1:3]
profTempAggr <- profTempAggr[profMonthAggr%in%1:3]
profPsalAggr <- profPsalAggr[profMonthAggr%in%1:3]
profIndexAggr <- (1:length(profMonthAggr))[profMonthAggr%in%1:3]
profMonthAggr<-profMonthAggr[profMonthAggr%in%1:3]


save(file = 'analysis/data/jan_march_data.RData',
     profCycleNumberAggr,
     profFloatIDAggr,
     profJulDayAggr,
     profLatAggr,
     profLongAggr,
     profModeAggr,
     profMonthAggr,
     profYearAggr,
     profJulDayofYearAggr,
     profJulFracofYearAggr,
     profIndexAggr,profPresAggr, profTempAggr, profPsalAggr)



# Save locations for where RG is defined
###
RG_temp <- nc_open('analysis/products/RG_ArgoClim_Temperature_2017.nc', write = FALSE)
longitude <- ncvar_get(RG_temp, "LONGITUDE")
longitude <- ifelse(longitude > 360, longitude - 360, longitude)
longitude <- ifelse(longitude > 180, longitude - 360, longitude)
latitude <- ncvar_get(RG_temp,"LATITUDE")
pressure <- ncvar_get(RG_temp,"PRESSURE")
time <- ncvar_get(RG_temp,"TIME")

mean_field_temp <- array(ncvar_get(RG_temp,"ARGO_TEMPERATURE_MEAN"),dim = c(360, 145, 58),
                         dimnames = list(longitude, latitude, pressure))
mean_field_long_temp <- melt(mean_field_temp, variable.name = "key",
                             value.names = "value")
RG_defined <- mean_field_temp[,,1]
RG_defined_long <- melt(RG_defined, variable.name = "key",
                        value.names = "value")
colnames(RG_defined_long) <- c('long', 'lat', 'value')
save(RG_defined_long, file = 'analysis/data/RG_Defined.RData')

# Wrangle KS residuals
###

library(R.matlab)
KS_pred <- list(R.matlab::readMat('analysis/products/CVPredsSpaceTimeExp_10_02_2007_2016.mat'),
                R.matlab::readMat('analysis/products/CVPredsSpaceTimeExp_300_02_2007_2016.mat'),
                R.matlab::readMat('analysis/products/CVPredsSpaceTimeExp_1500_02_2007_2016.mat'))
sum(sapply(KS_pred[[1]]$preds, function(x) nrow(x[[1]]))) ## 66110
sum(sapply(KS_pred[[2]]$preds, function(x) nrow(x[[1]]))) ## 69666
sum(sapply(KS_pred[[3]]$preds, function(x) nrow(x[[1]]))) ## 51064

p_levels <- c(10, 300, 1500)
year <- 2007:2016
final_list <- list()
for (j in 1:length(p_levels)) {
  p <- p_levels[j]
  year_data <- list()
  for (k in 1:length(year)) {
    KS_test <- R.matlab::readMat(paste0('analysis/products/residualsJohnSpaceTime_',
                                        p, '_02_', year[k], '_extended.mat'))
    year_data[[k]] <- lapply(KS_test, function(x) x[KS_pred[[j]]$leaveOutIdx[[k]][[1]]])
  }
  final_list[[j]] <- year_data
  print(j)
}

data <- lapply(1:length(final_list), function(x) {
  do.call(rbind, lapply(1:length(final_list[[x]]), function(y) {
    data.frame(float = as.vector(final_list[[x]][[y]]$interpFloatIDYear),
               day = as.vector(final_list[[x]][[y]]$interpJulDayYear),
               lat = as.vector(final_list[[x]][[y]]$interpLatYear),
               lon = as.vector(final_list[[x]][[y]]$interpLongYear),
               residual = as.vector(final_list[[x]][[y]]$interpResYear),
               KS_residual = as.vector(KS_pred[[x]]$preds[[y]][[1]]) - as.vector(final_list[[x]][[y]]$interpResYear),
               year = year[y])
  }))
})

data_10 <- data[[1]]
data_300 <- data[[2]]
data_1500 <- data[[3]]
data_10$id <- 1:nrow(data_10)
data_10$lon <- ifelse(data_10$lon > 180, data_10$lon - 360, data_10$lon)
data_10 <- merge(data_10, data.frame(day = profJulDayAggr,
                                     month = profMonthAggr,
                                     lon = profLongAggr,
                                     lat = profLatAggr,
                                     float = profFloatIDAggr,
                                     cycle = profCycleNumberAggr,
                                     profile_num  = 1:length(profFloatIDAggr)), all.x = T)
data_10 <- data_10[!duplicated(data_10$id),]
sum(data_10$month == 2)
# 66110

data_300$id <- 1:nrow(data_300)
data_300$lon <- ifelse(data_300$lon > 180, data_300$lon - 360, data_300$lon)
data_300 <- merge(data_300, data.frame(day = profJulDayAggr,
                                       month = profMonthAggr,
                                       lon = profLongAggr,
                                       lat = profLatAggr,
                                       float = profFloatIDAggr,
                                       cycle = profCycleNumberAggr,
                                       profile_num  = 1:length(profFloatIDAggr)), all.x = T)
data_300 <- data_300[!(duplicated(data_300$float) & duplicated(data_300$day) & duplicated(data_300$lon)
                       & duplicated(data_300$lat)& duplicated(data_300$year) & duplicated(data_300$id)),]
sum(data_300$month == 2)## 69666

data_1500$id <- 1:nrow(data_1500)
data_1500$lon <- ifelse(data_1500$lon > 180, data_1500$lon - 360, data_1500$lon)
data_1500 <- merge(data_1500, data.frame(day = profJulDayAggr,
                                         month = profMonthAggr,
                                         lon = profLongAggr,
                                         lat = profLatAggr,
                                         float = profFloatIDAggr,
                                         cycle = profCycleNumberAggr,
                                         profile_num  = 1:length(profFloatIDAggr)), all.x = T)
data_1500 <- data_1500[!duplicated(data_1500$id),]
sum(data_1500$month == 2)#51064

save(data_10, data_300, data_1500, file = 'analysis/products/KS_cv_floats_resid_final.RData')

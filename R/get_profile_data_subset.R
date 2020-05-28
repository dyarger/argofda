
#' Calculated great circle distance between two points
#'
#' based on fields::rdist.earth.vec, but checks for an error in the sqrt
#'
#' @param x1 matrix with two columns that together specify the location of point(s)
#' @param x2 matrix with two columns that together specify the location of point(s) that will be compared with x1
#' @param miles return distance in miles or kilometers
#' @param R radius of body
#'
#' @return Distance between points in x1 an dx2
#'
#' @export
rdist.earth.vec2=function(x1, x2, miles = TRUE, R = NULL) {
  if (is.null(R)) {
    if (miles)
      R = 3963.34
    else R = 6378.388
  }
  x1 = x1 * (pi/180)
  x2 = x2 * (pi/180)
  lonDist2 = (x2[, 1] - x1[, 1]) * (1/2)
  latDist2 = (x2[, 2] - x1[, 2]) * (1/2)
  a = sin(latDist2) * sin(latDist2) + cos(x1[, 2]) * cos(x2[, 2]) * sin(lonDist2) * sin(lonDist2)
  error.check=tryCatch(sqrt(1 - a), error = function(x) {sqrt(1 - round(a, 15))})
  dist=2 * atan2(sqrt(a), error.check) * R
  return(2 * atan2(sqrt(a), error.check) * R)
}



#' gets data near a location
#'
#' @param long longitude of interest
#' @param lat latitude of interest
#' @param RG_def is the RG mask defined? if so, expands bandwidth if few profiles
#' @param h_space distance in kilometers used to find nearby profiles
#' @param h_time number of days in time to find nearby profiles
#' @param mode if 'delayed', returns only delayed mode profiles
#' @param min_prof if RG_def = T, then makes sure there are min_prof from each year
#' @param exculsion if T removes one profile that has some issues
#' @param years vector with the years for that data is wanted, from 2007-2016
#'
#' @return data.frame with relevant data and information with each row a measurement, profiles identified by profile
#'
#' @export
get_profile_data_subset <- function(long, lat, day, RG_def, h_space, h_time = NULL, mode, 
                                           min_prof, exclusion = TRUE, years = 2007:2016 )  {
  if (is.null(lat)) {
    return(NULL)
  }
  # find distances between each profile and long/lat
  space_distances = rdist.earth.vec2(matrix(nrow = 1, c(long, lat)), 
                                     matrix(ncol = 2, c(profLongAggr, profLatAggr)), 
                                     miles = FALSE)
  space_distances_delayed = space_distances[profModeAggr == 'D']
  # ensure that there are min_prof profiles used for each year if RG is unmasked
  if (mode == 'delayed'){
    sd_year <- split(space_distances_delayed, profYearAggr[profModeAggr == 'D'])
  } else {
    sd_year <- split(space_distances, profYearAggr)
  }
  sd_year <- lapply(sd_year, function(x) x[order(x)])
  h_space_years <- sapply(sd_year, function(x) x[min_prof+1])
  h_space_2 <- h_space
  if (RG_def & max(h_space_years) > h_space) {
    h_space_2 <- max(c(max(h_space_years), h_space))
  }

  # calculate kernel weights
  wt_space = ifelse(space_distances/h_space_2 <= 1, (1/h_space_2) * 
                      (3/4) * (1 - (space_distances/h_space_2)^2), 0)
  wt_time = ifelse(abs(day - profJulFracofYearAggr)/h_time <= 1, (1/h_time) * 
                     (3/4) * (1 - (abs(day - profJulFracofYearAggr)/h_time)^2), 0)
  wt = wt_space * wt_time

  # pick out relevant pressure, temperature, and salinity data
  if (mode == "delayed") {
    relevant_indices = which(wt != 0 & profModeAggr == "D"  & profYearAggr %in% years)
  } else {
    relevant_indices = which(wt != 0  & profYearAggr %in% years)
  }
  if (length(relevant_indices) == 0) {
    print("no profiles in range")
    return(NULL)
  }
  Pres = lapply(profPresAggr[relevant_indices], function(x) x[[1]])
  if (!exists('profTempAggr')) {
    Temp = lapply(relevant_indices, function(x) rep(NA, length(profPresAggr[[relevant_indices[x]]][[1]])))
    Psal = lapply(profPsalAggr[relevant_indices], function(x) x[[1]])
  } else if(!exists('profPsalAggr')) {
    Temp = lapply(profTempAggr[relevant_indices], function(x) x[[1]])
    Psal = lapply(relevant_indices, function(x) rep(NA, length(profPresAggr[[relevant_indices[x]]][[1]])))
  } else {
    Temp = lapply(profTempAggr[relevant_indices], function(x) x[[1]])
    Psal = lapply(profPsalAggr[relevant_indices], function(x) x[[1]])
  }

  # calculate distance in EW and NS directions
  EW_dist_km = rdist.earth.vec2(matrix(nrow = 1, c(long,  lat)), 
                                matrix(ncol = 2, c(profLongAggr[relevant_indices], rep(lat, length(profLatAggr[relevant_indices])))), 
                                miles = FALSE)
  EW_dist_km = ifelse(abs(profLongAggr[relevant_indices] - long) > 180, 
                      ifelse(profLongAggr[relevant_indices] > long,  -EW_dist_km, EW_dist_km), 
                      ifelse(profLongAggr[relevant_indices] > long, EW_dist_km, -EW_dist_km))
  NS_dist_km = rdist.earth.vec2(matrix(nrow = 1, c(long,  lat)),
                                matrix(ncol = 2, c(rep(long, length(profLatAggr[relevant_indices])),  profLatAggr[relevant_indices])),
                                miles = FALSE)
  NS_dist_km = ifelse(profLatAggr[relevant_indices] > lat, NS_dist_km, -NS_dist_km)
  pressure = unlist(Pres)
  pres_remove = pressure <= 2000 & pressure >= 0
  pressure = pressure[pres_remove]
  temperature = unlist(Temp)[pres_remove]
  salinity = unlist(Psal)[pres_remove]
  prof_lengths = sapply(lapply(Pres, function(y) {
    return(y[y <= 2000 & y >= 0])
  }), length)
  
  # put it all in one big data frame
  data = data.frame(pressure = round(as.numeric(pressure), digits = 2), 
                    temperature, salinity, profile = rep(1:length(relevant_indices), 
                                                                   times = prof_lengths),
                    profile_num = rep(relevant_indices, times = prof_lengths),
                    latitude = rep(profLatAggr[relevant_indices], times = prof_lengths),
                    longitude = rep(profLongAggr[relevant_indices], times = prof_lengths), 
                    day = rep(profJulDayAggr[relevant_indices], times = prof_lengths), 
                    fracofyear = rep(profJulFracofYearAggr[relevant_indices],  times = prof_lengths), 
                    julian_date = rep(profJulDayAggr[relevant_indices],  times = prof_lengths), 
                    wt = rep(wt[relevant_indices], times = prof_lengths), 
                    float = rep(profFloatIDAggr[relevant_indices], times = prof_lengths), 
                    year = rep(profYearAggr[relevant_indices], times = prof_lengths),
                    total_dist_km = rep(space_distances[relevant_indices], times = prof_lengths),
                    EW_dist_km = rep(EW_dist_km, times = prof_lengths), 
                    NS_dist_km = rep(NS_dist_km, times = prof_lengths), 
                    cycle = rep(profCycleNumberAggr[relevant_indices],times = prof_lengths), 
                    mode = rep(profModeAggr[relevant_indices],  times = prof_lengths), 
                    month = rep(profMonthAggr[relevant_indices], times = prof_lengths), h_space = h_space_2)
  if (exclusion == TRUE) {
    data <- data[data$float != 2901534 | data$cycle != 188, ] # one profile has some issues
    prof_lengths <- table(data$profile)
    data$profile <- rep(1:length(prof_lengths), times = prof_lengths)
  }
  return(data)
}

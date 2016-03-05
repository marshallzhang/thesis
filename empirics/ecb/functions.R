# Break down by indicator and horizon.
years = 1999:2015
quarters = 1:4
forecasts = vector("list", length(years))
names(forecasts) = as.character(years)
for (year in years) {
  forecasts[[as.character(year)]] = yearly.forecast = vector("list", length(quarters))
  names(forecasts[[as.character(year)]]) = as.character(quarters)
  for (quarter in quarters) {
#     raw = read.csv(paste("/Users/marshall/Documents/senior/thesis/empirics/data/inflation/", year, "Q", quarter, ".csv", sep = ""), header = T, stringsAsFactors = F)[,-1]
    raw = read.csv(paste("../data/inflation/", year, "Q", quarter, ".csv", sep = ""), header = T, stringsAsFactors = F)[,-1]
    horizons = unique(raw$horizon)
    forecasts[[as.character(year)]][[as.character(quarter)]] = vector("list", length(horizons))
    names(forecasts[[as.character(year)]][[as.character(quarter)]]) = as.character(horizons)
    for (horizon.i in 1:length(horizons)) {
      forecasts[[as.character(year)]][[as.character(quarter)]][[as.character(horizons[horizon.i])]] = 
        raw[which(raw$horizon == horizons[horizon.i]), -1]
    }
  }
}

# Build list of actual economists who responded to everything.
responders = vector("list", length(years))
names(responders) = as.character(years)

for (year in years) {
  responders[[as.character(year)]] = vector("list", length(quarters))
  names(responders[[as.character(year)]]) = as.character(quarters)
  for (quarter in quarters) {
    responders[[as.character(year)]][[as.character(quarter)]] = 
      1:max(forecasts[[as.character(year)]][[as.character(quarter)]][[1]]$economist)
    for (frame in forecasts[[as.character(year)]][[as.character(quarter)]]) {
      responders[[as.character(year)]][[as.character(quarter)]]  = 
        intersect(responders[[as.character(year)]][[as.character(quarter)]], 
                  frame[which(!is.na(frame$forecast)), "economist"])  
    }
  }
}

# Build forecast matrices.

# current = read.csv(paste("/Users/marshall/Documents/senior/thesis/empirics/data/inflation/actual_inf_all_months.csv", sep = ""), header = T, stringsAsFactors = F)[,-1]
current = read.csv("../data/inflation/actual_inf_all_months.csv", header = T, stringsAsFactors = F)

forecast.mats = vector("list", length(years))
names(forecast.mats) = as.character(years)
for (year in years) {
  forecast.mats[[as.character(year)]] = vector("list", length(quarters))
  names(forecast.mats[[as.character(year)]]) = as.character(quarters)
  for (quarter in quarters) {
    forecast.mats[[as.character(year)]][[as.character(quarter)]] = 
      matrix(0, nrow = length(responders[[as.character(year)]][[as.character(quarter)]]), 
             ncol = length(forecasts[[as.character(year)]][[as.character(quarter)]]))
    counter = 1
    for (frame in forecasts[[as.character(year)]][[as.character(quarter)]]) {
      forecast.mats[[as.character(year)]][[as.character(quarter)]][, counter] = 
        frame[frame$economist %in% responders[[as.character(year)]][[as.character(quarter)]], "forecast"]  
      counter = counter + 1
    }  
    if (quarter == 1) {
      cur = paste("X", year - 1, "M", 12, sep = "")
    } else if (quarter == 2) {
      cur = paste("X", year, "M", "03", sep = "")
    } else if (quarter == 3) { 
      cur = paste("X", year, "M", "06", sep = "")
    } else {
      cur = paste("X", year, "M", "09", sep = "")
    }
    
    forecast.mats[[as.character(year)]][[as.character(quarter)]] =
      cbind(rep(as.numeric(current[cur]), length(responders[[as.character(year)]][[as.character(quarter)]])), 
            forecast.mats[[as.character(year)]][[as.character(quarter)]])
    colnames(forecast.mats[[as.character(year)]][[as.character(quarter)]]) = c("0", names(forecasts[[as.character(year)]][[as.character(quarter)]]))
  }
}

fc = function(year, quarter) {
  forecast.mats[[as.character(year)]][[as.character(quarter)]]
}

true.inf = function(year, quarter) {
    if (quarter == 1) {
      cur = paste("X", year, "M", "03", sep = "")
    } else if (quarter == 2) {
      cur = paste("X", year, "M", "06", sep = "")
    } else if (quarter == 3) { 
      cur = paste("X", year, "M", "09", sep = "")
    } else {
      cur = paste("X", year, "M", "12", sep = "")
    }
  current[cur]
}

ou.ll = function(theta, forecast) {
  beta = theta[2]
  alpha = theta[1] / theta[2]
  sigma = theta[3]
  horizons = as.numeric(colnames(forecast))
  ll = array(0, nrow(forecast))
  for (horizon.i in 2:length(horizons)) {
    dt = horizons[horizon.i] - horizons[horizon.i - 1]
    ll = ll + dnorm(forecast[,horizon.i], 
          alpha + (forecast[, horizon.i - 1] - alpha) * exp(- beta * dt),
          sqrt((sigma^2) * (1 - exp(-2 * beta * dt)) / (2 * beta)), log = T)
  }
  -mean(ll)
}
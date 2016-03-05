source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/ecb/functions.R")

quantiles = array(NA, 1000)

# Test how good economist predictions are.
counter = 1
for (year in years) {
  for (quarter in quarters) {
    if (year < 2000 | year > 2010) {
      next
    }
    forecasts = fc(year, quarter)
    horizons = as.numeric(colnames(forecasts))
    for (horizon.i in 2:length(horizons)) {
#     for (horizon.i in 2:2) {
      horizon = horizons[horizon.i]
      quarters.forward = quarter + horizon - 1
      add.year = ceiling((quarters.forward - 1) / 4) - 1
      new.quarter = 1 + (quarters.forward - 1) %% 4
      print(c(year, quarter, horizon, add.year, new.quarter))
      actual = true.inf(year + add.year, new.quarter)
      quantiles[counter] = ecdf(forecasts[,horizon.i])(actual)
      counter = counter + 1
    }
  }
}

quantiles = quantiles[!is.na(quantiles)]
plot(sort(runif(length(quantiles))), sort(quantiles))
abline(a=0,b=1)
ks.test(quantiles, runif(length(quantiles)))

# Test how good our unconstrained MLEs are.
mle.quantiles = array(NA, 1000)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    if (year < 2000 | year > 2010) {
      next
    }
    forecasts = fc(year, quarter)
    horizons = as.numeric(colnames(forecasts))
    for (horizon.i in 2:length(horizons)) {
#     for (horizon.i in 2:2) {
      horizon = horizons[horizon.i]
      
      # Simulate my process.
      if (quarter == 1) {
        mle.year = year -1
        mle.quarter = 4
      } else {
        mle.year = year
        mle.quarter = quarter
      }
      theta = gmle(mle.year, mle.quarter)
      preds = array(NA, nrow(forecasts))
      for (i in 1:nrow(forecasts)) {
#         preds[i] = ou(as.numeric(forecasts[1,1]), 100, theta, time = horizon)[100]
#         preds[i] = hyp(as.numeric(forecasts[1,1]), 100, theta, time = horizon)[100]
        preds[i] = hyp(as.numeric(forecasts[1,1]), 100, c(theta[1], theta[2], 0.15, theta[3]), time = horizon)[100]
      }
      
      # Find true.
      quarters.forward = quarter + horizon - 1
      add.year = ceiling((quarters.forward - 1) / 4) - 1
      new.quarter = 1 + (quarters.forward - 1) %% 4
      print(c(year, quarter, horizon, add.year, new.quarter))
      actual = true.inf(year + add.year, new.quarter)
      mle.quantiles[counter] = ecdf(preds)(actual)
      counter = counter + 1
    }
    
  }
}

mle.quantiles = mle.quantiles[!is.na(mle.quantiles)]
plot(sort(runif(length(mle.quantiles))), sort(mle.quantiles))
abline(a = 0, b= 1)
ks.test(mle.quantiles, runif(length(mle.quantiles)))


preds = matrix(NA, nrow = 200, ncol = length(2009:2015))
counter = 1
for (year in 2009:2015) {
    print(year)
    forecast = fc(year - 5,4)[,5]
    preds[1:length(forecast), counter] = forecast
    counter = counter + 1
}
boxplot(preds, names = 2009:2015, ylim = c(-1,4))
for (i in 1:length(2009:2015)) {
  points(i, true.inf((2009:2015)[i], 4), pch = 15)
}


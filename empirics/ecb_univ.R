source("/Users/marshall/Documents/senior/thesis/empirics/functions.R")
library(mvtnorm) 

# Compare economist distributions to mine.
r.economists = function(data, responders, N, start, end) {
  samples = matrix(0, nrow = N, ncol = 6)
  r.i = sample(1:length(responders), N, replace = T)
  for (i in 1:N) {
    samples[i, ] = c(sapply(data,function(x) x[r.i[i], start]), sapply(data, function(x) x[r.i[i], end]))
  } 
  samples  
}  

# Get parameters of inf/gdp/une process.
lik.ou = function(theta, dat, horizons, responders) {
  if (theta[2] <= 0 | theta[3] <= 0) {
    return(Inf)
  }
  mu = theta[1] / theta[2]
  lambda = theta[2]
  sigma = theta[3]
  ll = 0
  for (responder in 1:length(responders)) {
    response = dat[responder, ]
    for (horizon.i in 1:length(horizons)) {
      current = response[horizon.i + 1]
      old = response[horizon.i]
      n.mean = mu + (old - mu) * exp(-lambda * horizons[horizon.i])
      n.sig = sqrt(sigma^2 * (1 - exp(-2 * lambda * horizons[horizon.i])) / (2 * lambda))
      ll = ll + dnorm(current, mean = n.mean, sd = n.sig, log = T)
    }
  }
  print(ll)
  -ll
} 

analyze.predictions = function(raw, current.indics, samples) {
  indicators = c("inf", "gdp", "une")
  
  # Break down by indicator and horizon.
  raw.list = vector("list", length(indicators))
  names(raw.list) = indicators
  raw.list$inf = raw[which(raw$indicator == "inflation"), -1]
  raw.list$gdp = raw[which(raw$indicator == "growth"), -1]
  raw.list$une = raw[which(raw$indicator == "unemployment"), -1] 
  horizons = unique(raw$horizon)
  economists = max(raw.list$inf$source)
  
  data = vector("list", 3)
  names(data) = indicators
  for (i in 1:length(indicators)) {
    data[[i]] = vector("list", length(horizons))
    names(data[[i]]) = horizons
    for (horizon.i in 1:length(horizons)) {
      data[[i]][[horizon.i]] = raw.list[[i]][which(raw.list[[i]]$horizon == horizons[horizon.i]), - 1]  
    }
  }
  
  # Build list of actual economists who responded to everything.
  responders = 1:economists
  for (i in 1:length(indicators)) {
    for (j in 1:length(horizons)) {
      frame = data[[i]][[j]]
      responders = intersect(responders, frame[which(!is.na(frame$point)), "source"])  
    }
  }
  
  # Build forecast matrices.
  current = current.indics
  forecasts = vector("list", 3)
  names(forecasts) = indicators
  for (i in 1:length(indicators)) {
    forecasts[[i]] = matrix(0, nrow = length(responders), ncol = length(horizons))
    for (j in 1:length(horizons)) {
      frame = data[[i]][[j]]
      forecasts[[i]][,j] = frame[frame$source %in% responders, "point"]
    }
    forecasts[[i]] = cbind(rep(current[i], length(responders)), forecasts[[i]])
  }
  
#   theta = optim(c(0, 1, 1),
#         method = "CG",
#         fn = function(theta) lik.ou(theta, forecasts[[1]], horizons, responders),
# #         lower = c(-Inf, 1e-6, 1e-6),
# #         upper = c(Inf,  Inf,  Inf),
#         control = list(trace = 6))$par

  theta = DEoptim(function(theta) lik.ou(theta, forecasts[[1]], horizons, responders),
                  lower = c(-5, 0, 0),
                  upper = c(5, 10, 10),
                  control = DEoptim.control(trace = T))

  print(theta$bestval)
  print(theta$bestmem)
  
#   implied = gen.mv.nu.bridge(function(n) r.economists(forecasts, responders, n, 1, 5), theta, samples, horizons[4] + 1)
#   
#   interpolated = vector("list", 3)
#   names(interpolated) = indicators
#   for (i in 1:length(indicators)) {
#     interpolated[[i]] = matrix(0, nrow = samples, ncol = length(horizons) + 1)
#     for (j in 1:samples) {
#       interpolated[[i]][j, ] = implied[[j]][i, c(1, horizons + 1)]
#     }
#   }

#   print(theta)
#   list(nu.bridges = implied, forecasts = forecasts, theta = theta, interpolated = interpolated, forecast.years = year + horizons / 4)
}

raw = read.csv("/Users/marshall/Documents/senior/thesis/empirics/data/1999-1.csv")
current.indics = c(1.1, 2, 10)

res = analyze.predictions(raw, current.indics, 100)
source("/Users/marshall/Documents/senior/thesis/empirics/functions.R")
library(mvtnorm) 

# Make orthogonal matrix
ortho = function(a, b, c) {
  matrix(c(cos(a) * cos(c) - sin(a) * sin(b) * sin(c), -sin(a) * cos(b), -cos(a)*sin(c) - sin(a) * sin(b) * cos(c),
           cos(a)*sin(b)*sin(c)+sin(a)*cos(c), cos(a)*cos(b), cos(a)*sin(b)*cos(c)-sin(a)*sin(c),
           cos(b) * sin(c), -sin(b), cos(b)*cos(c)), nrow = 3, ncol = 3, byrow =T)
}

gen.cov = function(t) {
  a = t[1]
  b = t[2]
  c = t[3]
  d = t[4]
  e = t[5]
  f = t[6]
  G = ortho(a,b,c)
  l1 = d + e + f
  l2 = e + f
  l3 = f
  t(G) %*% diag(c(l1,l2,l3)) %*% G
}

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
  mu = theta[1:3]
  lambda = rep(theta[4], 3)
  cov = gen.cov(theta[5:10])
  ll = 0
  for (responder in 1:length(responders)) {
    response = sapply(dat, function(x) x[responder, ])
    for (horizon.i in 1:length(horizons)) {
      current = response[horizon.i + 1, ]
      old = response[horizon.i, ]
      Lambda = exp(diag(lambda) * horizons[horizon.i])
      Lambda[which(Lambda == 1)] = 0
      n.mean = mu + Lambda %*% (old - mu)
      n.cov = cov - Lambda %*% cov %*% Lambda
      ll = ll + dmvnorm(current, mean = n.mean, sigma = n.cov, log = T)
    }
  }
  print(ll)
  -ll
} 

analyze.predictions = function(raw, current.indics) {
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
  
  theta = optim(c(2, 2, 10, -0.2, -1, 0.6, 1, 0.1, 0.1, 0.1),
        method = "L-BFGS-B",
        fn = function(theta) lik.ou(theta, forecasts, horizons, responders),
        lower = c(-Inf, -Inf, 1e-6, -Inf, -Inf, -Inf, -Inf, 1e-6, 1e-6, 1e-6),
        upper = c(Inf,  Inf,  Inf,  -1e-6, Inf,  Inf,  Inf,   Inf,  Inf, Inf),
        control = list(trace = 6, factr = 5e11))$par
  
  implied = gen.mv.nu.bridge(function(n) r.economists(forecasts, responders, n, 1, 5), theta, 100, horizons[4] + 1)
  
  

  list(nu.bridges = implied, forecasts = forecasts, theta = theta)
}

raw = read.csv("/Users/marshall/Documents/senior/thesis/empirics/data/2016-1.csv")
current.indics = c(0.4, 1.2, 10.5)

res = analyze.predictions(raw, current.indics)

# for (i in 1:4){
#   inf16 = sapply(implied, function(x) x[1,c(4,8,12,20)[i]])
#   hist(forecasts$inf[,i+1], freq = F, xlim = c(0, 2.5), ylim = c(0,3))
#   hist(inf16,freq=F, add=T, xlim = c(0,2),col = rgb(1,0,0,0.4))
# }
# 
# infs = sapply(implied, function(x) x[1,])
# plot(infs[,1], type = "l", col = rgb(0,0,0,0.1))
# for (i in 2:100)  {
#   lines(infs[,i], col = rgb(0,0,0,0.1))
# }

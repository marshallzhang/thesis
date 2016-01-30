library(mvtnorm)
raw = read.csv("/Users/marshall/Documents/senior/thesis/empirics/data/2016-1.csv")
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
current = c(0.4, 1.2, 10.5)
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

# Likelihood of multivariate OU.
lik.ou = function(theta, dat) {
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

optim(c(1, 1.5, 9.5, -0.2, 0.05, 0, 0.05, 0,0,0.05),
      method = "L-BFGS-B",
      fn = function(theta) lik.ou(theta, forecasts),
      lower = c(-Inf, -Inf, 1e-6, -Inf, -Inf, -Inf, -Inf, 1e-6, 1e-6, 1e-6),
      upper = c(Inf,  Inf,  Inf,  -1e-6, Inf,  Inf,  Inf,   Inf,  Inf, Inf),
      control = list(trace = 6))

lik.ou(c(0, 0, 9.5, -0.2, 0.05, 0, 0.05, 0,0,0.05), forecasts)

theta = c(1.91364, 1.85883, 9.10795, -0.0844359, -0.0672126, -0.658418, 1.3999, 0.101771, 1e-06, 0.128031)
gen.cov(theta[5:10])

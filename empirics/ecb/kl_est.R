
est = function(x, xs) {
  xs = sort(xs)
  ec = ecdf(xs)
  extended.xs = unique(c(min(xs) - 1, xs, max(xs) + 1))
  i = findInterval(x, extended.xs)
  if (i == length(extended.xs)) {
    1
  } else if (i == 0) {
    0 
  } else {
    points = extended.xs[i:(i+1)]
    xi1 = ec(points[1])
    xi = ec(points[2])
    prop = (x - points[1]) / (points[2] - points[1])
    prop * xi + (1-prop) * xi1
  }
}

est.div = function(Ps, Qs) {
  n = length(Ps)
  Ps = sort(Ps)
  Qs = sort(Qs)
  ep = min(diff(Ps))/2
  div = 0
  for (i in 2:(length(Ps) - 1)) {
    diff =  log((est(Ps[i], Ps) - est(Ps[i] - ep, Ps)) / (est(Ps[i], Qs) - est(Ps[i] - ep, Qs))) / (n - 2)
    if (!is.finite(diff)) {
      next
    }
    div = div + diff
  }
  div
}

Ps = rnorm(50,0,1)
Qs = c(rnorm(25,2,1), rnorm(25,-2,1))

print(est.div(Ps, Qs) - 1)
est.div(trues[["1999"]][["4"]]), linears[["1999"]][["4"]])

# Turn forecast mats into linear interpolation.
linears = sapply(forecast.mats, function(x) sapply(x, function(y) {
  names = as.numeric(colnames(y))
  prop = names[2] / names[3]
  prop * y[,3] + (1-prop) * y[,1]
}), simplify = F)

trues = sapply(forecast.mats, function(x) sapply(x, function(y) {
  y[,2]
}), simplify = F)

linear.true.kl = array(NA, 68)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    print(c(year, quarter))
    l = length(trues[[as.character(year)]][[as.character(quarter)]])
    kl = est.div(trues[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01), linears[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01))
    print(kl)
    linear.true.kl[counter] = kl  - 1
    counter = counter + 1
  }
}
linear.true.kl[linear.true.kl < 0 ] = NA
mean(linear.true.kl, na.rm = T)
sd(linear.true.kl, na.rm = T)
hist(linear.true.kl)
 
ous = linears
counter = 1
for (year in years) {
  print(year)
  for (quarter in quarters) {
    print(quarter)
    if (year == 2010 & quarter == 4) {
      next
    }
    if (year == 2014 & quarter == 2) {
      next
    }
    for (i in 1:length(linears[[as.character(year)]][[as.character(quarter)]])) {
      mat = forecast.mats[[as.character(year)]][[as.character(quarter)]]
      horizons = as.numeric(colnames(mat))
      points = mat[i,c(1,3)]
      ous[[as.character(year)]][[as.character(quarter)]][i] = gen.nu.bridge(function(n) rep(0,n), ou.fast, thetas[[as.character(year)]][[as.character(quarter)]], 1, 100, set.start = points, time = horizons[3] - horizons[1])[, 100 * horizons[2] / horizons[3]]
    }
    counter = counter + 1
  }
}

ous.true.kl = array(NA, 68)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    l = length(trues[[as.character(year)]][[as.character(quarter)]])
    kl = est.div(trues[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01), ous[[as.character(year)]][[as.character(quarter)]]+ rnorm(l, 0, 0.01))
    print(kl)
    ous.true.kl[counter] = kl - 1
    counter = counter + 1
  }
}
ous.true.kl[ous.true.kl < 0 ] = NA
mean(ous.true.kl, na.rm = T)
sd(ous.true.kl, na.rm = T)
hist(ous.true.kl)

hyp.thetas = matrix(NA, ncol = 4, nrow = 68)
hyp = linears
counter = 1
for (year in years) {
  print(year)
  for (quarter in quarters) {
    print(quarter)
    if (year == 2010 & quarter == 4) {
      next
    }
    if (year == 2014 & quarter == 2) {
      next
    }
    if (year == 2015 & quarter == 4) {
      next
    }
    for (i in 1:length(linears[[as.character(year)]][[as.character(quarter)]])) {
      mat = forecast.mats[[as.character(year)]][[as.character(quarter)]]
      horizons = as.numeric(colnames(mat))
      points = mat[i,c(1,3)]
      sigma = thetas[[as.character(year)]][[as.character(quarter)]][3]
      real.mu = thetas[[as.character(year)]][[as.character(quarter)]][1]/thetas[[as.character(year)]][[as.character(quarter)]][2]
      coefs = coef(fit.hypuv(mat[,ncol(mat)], symmetric = T, mu = real.mu, opt.pars = c(lambda = T, alpha.bar = T, mu = F, sigma = T, gamma = F), silent =T), type = "alpha.delta")
      theta = c(coefs$mu, coefs$alpha, coefs$delta, sigma)
      hyp.thetas[counter, ] = theta
      hyp[[as.character(year)]][[as.character(quarter)]][i] = gen.nu.bridge(function(n) rep(0,n), hyp.fast, theta, 1, 100, set.start = points, time = horizons[3] - horizons[1])[, 100 * horizons[2] / horizons[3]]
    }
    print(theta)
    counter = counter + 1
  }
}

hyp.true.kl = array(NA, 68)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    l = length(trues[[as.character(year)]][[as.character(quarter)]])
    kl = est.div(trues[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01), hyp[[as.character(year)]][[as.character(quarter)]]+ rnorm(l, 0, 0.01))
    print(kl)
    hyp.true.kl[counter] = kl - 1
    counter = counter + 1
  }
}
hyp.true.kl[hyp.true.kl < 0 ] = NA
hyp.true.kl[which(hyp.thetas[,2] > 30)] = NA
hyp.true.kl[38] = NA
mean(hyp.true.kl, na.rm = T)
sd(hyp.true.kl, na.rm = T)
hist(hyp.true.kl)



kl.divs = data.frame(linear = linear.true.kl, ou = ous.true.kl, hyp = hyp.true.kl)
kl.divs = kl.divs[complete.cases(kl.divs), ]
kl.divs = data.frame(cond = factor(rep(c("Linear Interpolation", "Gen. Bridge (OU)", "Gen. Bridge (Hyp.)"), each = nrow(kl.divs))),
                     kl = c(kl.divs[,1], kl.divs[,2], kl.divs[,3]))
pretty.plot(kl.divs) + geom_boxplot(aes(x=cond, y= kl))

















# Turn forecast mats into f.linear interpolation.
f.linears = sapply(forecast.mats, function(x) sapply(x, function(y) {
  names = as.numeric(colnames(y))
  n = length(names)
  prop = (names[n-1] - names[n-2]) / (names[n] - names[n-2])
  prop * y[,n] + (1-prop) * y[,n-2]
}), simplify = F)

trues = sapply(forecast.mats, function(x) sapply(x, function(y) {
  n = ncol(y)
  y[,n-1]
}), simplify = F)

f.linear.true.kl = array(NA, 68)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    print(c(year, quarter))
    l = length(trues[[as.character(year)]][[as.character(quarter)]])
    kl = est.div(trues[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01), f.linears[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01))
    print(kl)
    f.linear.true.kl[counter] = kl  - 1
    counter = counter + 1
  }
}
f.linear.true.kl[f.linear.true.kl < 0 ] = NA
mean(f.linear.true.kl, na.rm = T)
sd(f.linear.true.kl, na.rm = T)
hist(f.linear.true.kl)
 
f.ous = f.linears
counter = 1
for (year in years) {
  print(year)
  for (quarter in quarters) {
    print(quarter)
    if (year == 2010 & quarter == 4) {
      next
    }
    if (year == 2014 & quarter == 2) {
      next
    }
    for (i in 1:length(f.linears[[as.character(year)]][[as.character(quarter)]])) {
      mat = forecast.mats[[as.character(year)]][[as.character(quarter)]]
      horizons = as.numeric(colnames(mat))
      n = length(horizons)
      points = mat[i,c(n-2,n)]
      f.ous[[as.character(year)]][[as.character(quarter)]][i] = gen.nu.bridge(function(n) rep(0,n), ou.fast, thetas[[as.character(year)]][[as.character(quarter)]], 1, 100, set.start = points, time = horizons[n] - horizons[n-2])[, 100 * (horizons[n-1] - horizons[n-2]) / (horizons[n] - horizons[n-2])]
    }
    counter = counter + 1
  }
}

f.ous.true.kl = array(NA, 68)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    l = length(trues[[as.character(year)]][[as.character(quarter)]])
    kl = est.div(trues[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01), f.ous[[as.character(year)]][[as.character(quarter)]]+ rnorm(l, 0, 0.01))
    print(kl)
    f.ous.true.kl[counter] = kl - 1
    counter = counter + 1
  }
}
f.ous.true.kl[f.ous.true.kl < 0 ] = NA
mean(f.ous.true.kl, na.rm = T)
sd(f.ous.true.kl, na.rm = T)
hist(f.ous.true.kl)




f.hyp = linears
counter = 1
for (year in years) {
  print(year)
  for (quarter in quarters) {
    print(quarter)
    if (year == 2010 & quarter == 4) {
      next
    }
    if (year == 2014 & quarter == 2) {
      next
    }
    if (year == 2015 & quarter == 4) {
      next
    }
    for (i in 1:length(linears[[as.character(year)]][[as.character(quarter)]])) {
      mat = forecast.mats[[as.character(year)]][[as.character(quarter)]]
      horizons = as.numeric(colnames(mat))
      n = length(horizons)
      points = mat[i,c(n-2,n)]
      theta = hyp.thetas[counter, ]
      f.hyp[[as.character(year)]][[as.character(quarter)]][i] = gen.nu.bridge(function(n) rep(0,n), hyp.fast, theta, 1, 100, set.start = points, time = horizons[n] - horizons[n-2])[, 100 * (horizons[n-1] - horizons[n-2]) / (horizons[n] - horizons[n-2])]
    }
    print(theta)
    counter = counter + 1
  }
}

f.hyp.true.kl = array(NA, 68)
counter = 1
for (year in years) {
  for (quarter in quarters) {
    l = length(trues[[as.character(year)]][[as.character(quarter)]])
    kl = est.div(trues[[as.character(year)]][[as.character(quarter)]] + rnorm(l, 0, 0.01), f.hyp[[as.character(year)]][[as.character(quarter)]]+ rnorm(l, 0, 0.01))
    print(kl)
    f.hyp.true.kl[counter] = kl - 1
    counter = counter + 1
  }
}
f.hyp.true.kl[f.hyp.true.kl < 0 ] = NA
f.hyp.true.kl[which(hyp.thetas[,2] > 30)] = NA
f.hyp.true.kl[38] = NA
mean(f.hyp.true.kl, na.rm = T)
sd(f.hyp.true.kl, na.rm = T)
hist(f.hyp.true.kl)













f.kl.divs = data.frame(linear = f.linear.true.kl, ou = f.ous.true.kl, hyp = f.hyp.true.kl)
f.kl.divs = f.kl.divs[complete.cases(f.kl.divs), ]
f.kl.divs = data.frame(cond = factor(rep(c("Linear Interpolation", "Gen. Bridge (OU)", "Gen. Bridge (Hyp.)"), each = nrow(f.kl.divs))),
                     kl = c(f.kl.divs[,1], f.kl.divs[,2], f.kl.divs[,3]))

c = pretty.plot(kl.divs) + geom_boxplot(aes(x=cond, y= kl)) + 
              xlab("") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Estimated K-L Divergence") + 
              scale_y_continuous(breaks = seq(0,6,0.5)) + 
              coord_cartesian(ylim=c(0,6))+
              theme(legend.position = "none") + 
            ggtitle("1-4 Month Horizon")

f = pretty.plot(f.kl.divs) + geom_boxplot(aes(x=cond, y= kl)) + 
              xlab("") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Estimated K-L Divergence") + 
              scale_y_continuous(breaks = seq(0,6,0.5)) + 
              theme(legend.position = "none") + 
              coord_cartesian(ylim=c(0,6))+
            ggtitle("7-10 Month Horizon")

multiplot(c, f, cols = 2)

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/est_kl.pdf", width= 13, height = 6.5, #' see how it looks at this size
    useDingbats=F)
multiplot(c, f, cols = 2)
dev.off()

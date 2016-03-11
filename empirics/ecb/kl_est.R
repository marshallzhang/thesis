Ps = rnorm(10000,0,1)
Qs = rexp(10000,1)

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
      return(div)
    }
    div = div + diff
  }
  div
}

print(est.div(Qs, Ps) - 1)
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


kl.divs = data.frame(linear = linear.true.kl, ou = ous.true.kl)
kl.divs = kl.divs[complete.cases(kl.divs), ]
kl.divs = data.frame(cond = factor(rep(c("Linear Interpolation", "Gen. OU Bridge"), each = nrow(kl.divs))),
                     kl = c(kl.divs[,1], kl.divs[,2]))
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

f.kl.divs = data.frame(linear = f.linear.true.kl, ou = f.ous.true.kl)
f.kl.divs = f.kl.divs[complete.cases(f.kl.divs), ]
f.kl.divs = data.frame(cond = factor(rep(c("Linear Interpolation", "Gen. OU Bridge"), each = nrow(f.kl.divs))),
                     kl = c(f.kl.divs[,1], f.kl.divs[,2]))

c = pretty.plot(kl.divs) + geom_boxplot(aes(x=cond, y= kl)) + 
              xlab("") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Estimated K-L Divergence") + 
              scale_y_continuous(breaks = seq(0,6,0.5)) + 
              coord_cartesian(ylim=c(-0.5,6))+
              theme(legend.position = "none") + 
            ggtitle("1-4mo Horizon")

f = pretty.plot(f.kl.divs) + geom_boxplot(aes(x=cond, y= kl)) + 
              xlab("") + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 12)) + 
              ylab("Estimated K-L Divergence") + 
              scale_y_continuous(breaks = seq(0,6,0.5)) + 
              theme(legend.position = "none") + 
              coord_cartesian(ylim=c(-0.5,6))+
            ggtitle("7-10mo Horizon")

multiplot(c, f, cols = 2)

pdf(file = "/Users/marshall/Documents/senior/thesis/figures/est_kl.pdf", width= 10, height = 5, #' see how it looks at this size
    useDingbats=F)
multiplot(c, f, cols = 2)
dev.off()

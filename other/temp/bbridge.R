n = 1000

path = array(0, dim = n)
for (i in 2:n) {
  path[i] = 
    if (runif(1) > 0.5) {
      path[i-1] + 
      ((-5 - path[i-1]) / (1 - (i/n))) * (1/n) +
      rnorm(1, 0, 1) * sqrt((1/n))
    } else {
      path[i-1] + 
      ((5 - path[i-1]) / (1 - (i/n))) * (1/n) +
      rnorm(1, 0, 1) * sqrt((1/n))
    }
}
plot(path, type = "l")
path
abline(h=0)
abline(v=0)
abline(v=1000)
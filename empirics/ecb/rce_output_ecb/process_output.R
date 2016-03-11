
ids = 0:67
ests = matrix(NA, nrow = 68, ncol = 3)
counter = 1
for (id in ids) {
  f = readLines(file(paste("/Users/marshall/Documents/senior/thesis/empirics/ecb/rce_output_ecb/out_fcast.", id, sep = ""), open="r"))
  last.par = 0
  for (i in 1:length(f)) {
    if (length(grep("par1", f[i])) > 0) {
      last.par = i
    }
  }
  print(floor(last.par/200))
  pars = as.numeric(na.omit(as.numeric(strsplit(f[last.par + 1], " ")[[1]])))
  if (floor(last.par/200) <10) {
    pars = c(NA, NA, NA)
  } 
  ests[counter, ] = pars
  counter = counter + 1
}
matrix.thetas = matrix(unlist(thetas), ncol = 3, byrow = T)
ests = ests[4:68, ]
matrix.thetas = matrix.thetas[4:68, ]
ests[,1] = ests[,1] / ests[,2]
matrix.thetas[,1] = matrix.thetas[,1] / matrix.thetas[,2]
bias = ests - matrix.thetas
bias = bias[complete.cases(bias), ]

print(colMeans(ests - matrix.thetas, na.rm = T))
cov(bias)
print((colSds(ests - matrix.thetas, na.rm = T)))

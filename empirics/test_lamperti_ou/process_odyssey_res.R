
ids = c(4) + 1:99 * 6
ests = matrix(NA, nrow = 100, ncol = 3)
counter = 1
for (id in ids) {
  f = readLines(file(paste("/Users/marshall/Documents/senior/thesis/empirics/thetas/out_exact_many.", id, sep = ""), open="r"))
  last.par = 0
  for (i in 1:length(f)) {
    if (length(grep("par1", f[i])) > 0) {
      last.par = i
    }
  }
  if (floor(last.par/200) <4 ) {
    next
  }
  pars = as.numeric(na.omit(as.numeric(strsplit(f[last.par + 1], " ")[[1]])))
  ests[counter, ] = pars
  counter = counter + 1
}

print(counter)
print(colMeans(ests, na.rm = T))
ests[,1] = ests[,1]/ests[,2]
print((colSds(ests, na.rm = T)))

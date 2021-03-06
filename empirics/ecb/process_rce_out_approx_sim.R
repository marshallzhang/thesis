
ids = 0:9*6+5
ests = matrix(NA, nrow = 100, ncol = 3)
counter = 1
for (id in ids) {
  f = readLines(file(paste("/Users/marshall/Documents/senior/thesis/empirics/ecb/rce_output_ecb/out.", id,sep = ""), open="r"))
  last.par = 0
  for (i in 1:length(f)) {
    if (length(grep("par1", f[i])) > 0) {
      last.par = i
    }
  }
  print(floor(last.par/200))
  pars = as.numeric(na.omit(as.numeric(strsplit(f[last.par + 1], " ")[[1]])))
  ests[counter, ] = pars
  counter = counter + 1
}

print(colMeans(ests, na.rm = T))
print((colSds(ests, na.rm = T)))

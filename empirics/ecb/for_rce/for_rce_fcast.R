source("../../general_functions.R")
source("functions.R")
source("../mcem_functions.R")

my_id = as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1
years = 1999:2015
quarters = 1:4

times = vector("list", length(years) * length(quarters))

for (i in 1:length(years)) {
  for (j in 1:length(quarters)) {
    times[[(i - 1) * 4 + j]] = c(years[i], quarters[j])
  }
}

year = times[[my_id]][1]
quarter = times[[my_id]][2]

forecast = fc(year, quarter)

nu = function(n) {(as.matrix(forecast[sample(1:nrow(forecast), n, replace = T),]))}

horizons = as.numeric(colnames(forecast))

steps = max(horizons) * 50
oldpar = c(1,1,1)
for (i in 1:10) {
  obj = function(theta) {
    source("/Users/marshall/Documents/senior/thesis/empirics/ecb/for_rce/to_source_parallel.R")
    q.beskos.unknown(theta, trans.ou.fast, oldpar, nu, round(i^1.2)*100, steps, horizons/max(horizons), max(horizons), exact = F)
  }
  oldpar = DEoptim2(fn = obj,
                   lower = c(-1.5, 0, 0),
                   upper = c(1.5, 2, 2),
                   control = DEoptim.control(trace = TRUE, strategy = 6, parallelType = 1, limitCores= 8,
                                             packages = c(),
                                             parVar = list("q.beskos.unknown", "eta", "nu", "i", "A",
                                                           "steps", "oldpar", "horizons",
                                                           "vec.seq", "a.2.a.prime", "forecast")
                                                        ))[["optim"]]$bestmem
  print(oldpar)
}
source("../../general_functions.R")
source("../mcem_functions.R")

process = as.integer(as.character(commandArgs(trailingOnly = TRUE))) + 1

theta = c(0, 1, 1)

steps = 100 

r.diffusions = c(ou.fast,
                 wiener,
                 function(start, steps, theta, time = 1) {as.numeric(fbm(hurst = 0.2, steps - 1))},
                 ou.fast,
                 wiener,
                 function(start, steps, theta, time = 1) {as.numeric(fbm(hurst = 0.2, steps - 1))}
)

r.diffusion = r.diffusions[[1 + ((process - 1) %% 6)]]

horizons = c(.5)
tot.time = 1
print("Making nu.")
if (1 + ((process - 1)) %% 6 >= 4) {
  nu = make.joint.many(1000, theta, r.diffusion, function(n) rep(0,n ), horizons = horizons, steps = 100, time = tot.time, all = F)
} else {
  nu = function(n) {
    sims2 = array(0, dim = c(n, steps))
    sims.start = rep(0,n)
    for (i in 1:n) {
      sims2[i, ] = r.diffusion(sims.start[i], steps, theta, tot.time)
    }
    sims2[,1] = 0
    as.matrix(sims2[,c(1/steps,horizons,1) * steps])
  }
}

print("Starting DEoptim.")
oldpar = c(4,2,2)
steps = 100
for (i in 1:10) {
  obj = function(theta) {
    source("/Users/marshall/Documents/senior/thesis/empirics/ecb/for_rce/to_source_parallel.R")
    q.beskos.unknown(theta, trans.ou.fast, oldpar, nu, round(i^1.2)*100, steps, c(0, horizons, 1), tot.time, exact = T)
  }
  oldpar = DEoptim2(fn = obj,
                   lower = c(-1.5, 0, 0),
                   upper = c(1.5, 2, 2),
                   control = DEoptim.control(trace = TRUE, strategy = 6, parallelType = 1, limitCores= 8,
                                             packages = c("somebm"),
                                             parVar = list("q.beskos.unknown", "eta", "nu", "i", "A",
                                                           "steps", "oldpar", "horizons", "tot.time",
                                                           "vec.seq", "a.2.a.prime", "r.diffusion", "theta")
                                                        ))[["optim"]]$bestmem
  print(oldpar)
}
if (Sys.getenv("SLURM_JOB_ID") != "") {
  
  source("../general_functions.R")
  source("functions.R")
  source("mcem_functions.R")
  
  my_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
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
  
  print("Starting optim.")
  set.seed(2 + my_id)
  theta.star = c(.42, .33, .76)
  for (i in 1:10) {
    theta.star = DEoptim(fn = function(theta) q.beskos.unknown(theta, trans.ou.fast, theta.star, function(n) forecast[sample(1:nrow(forecast), n, replace = T), ], 100 * i^1.2, 100, as.numeric(colnames(forecast)), exact = F),
                     lower = c(0, 0, 0),
                     upper = c(16, 8, 6),
                     control = DEoptim.control(trace = TRUE, strategy = 6, itermax = 200))[["optim"]]$bestmem
    print(theta.star)
    if (i == 5) {
      cat(paste("i5", year, quarter, theta.star[1], theta.star[2], theta.star[3], sep = ","), file="/dump/ou_optim.txt", append=TRUE, sep = "\n")
    } else if (i == 10) {
      cat(paste("i10", year, quarter, theta.star[1], theta.star[2], theta.star[3], sep = ","), file="/dump/ou_optim.txt", append=TRUE, sep = "\n")
    }
  }
}
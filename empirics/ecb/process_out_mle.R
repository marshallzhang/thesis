source("/Users/marshall/Documents/senior/thesis/empirics/general_functions.R")
source("/Users/marshall/Documents/senior/thesis/empirics/ecb/functions.R")

mle = read.csv("/Users/marshall/Documents/senior/thesis/empirics/ecb/out_mle.txt", sep = ",", header = F)

get.theta = function(year, quarter) {
  as.numeric(mle[intersect(which(year == mle[,1]), which(quarter == mle[,2])), ][,3:5])
}

thetas = vector("list", length(years))
names(thetas) = as.character(years)
for (year in years) {
  thetas[[as.character(year)]] = vector("list", length(quarters))
  names(thetas[[as.character(year)]]) = as.character(quarters)
  for (quarter in quarters) {
    thetas[[as.character(year)]][[as.character(quarter)]] = get.theta(year, quarter)
  }
}

gmle = function(year, quarter) {
  thetas[[as.character(year)]][[as.character(quarter)]]
}

plot(as.vector(unlist(sapply(thetas, function(x) sapply(x, function(z) z[1])))))



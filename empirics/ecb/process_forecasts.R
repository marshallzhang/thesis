years = 1999:2015
quarters = 1:4
for (year in years) {
  for (quarter in quarters) {
    hi = read.csv(paste("/Users/marshall/Documents/senior/thesis/empirics/data/SPF_individual_forecasts/", year, "Q", quarter, ".csv", sep = ""), header = T, stringsAsFactors = F)
    preds = hi[2:which(hi[,1] == "")[1], 1:3]
    preds = sapply(preds, as.numeric)
    preds = preds[which(!is.na(preds[,1])),]
    preds[,1] = (preds[,1] - year + 1) * 4 - (quarter - 1)
    colnames(preds) = c("horizon","economist","forecast")  
    write.csv(preds, paste("/Users/marshall/Documents/senior/thesis/empirics/data/inflation/", year, "Q", quarter, ".csv", sep = ""))
  }
}

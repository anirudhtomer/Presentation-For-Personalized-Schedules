files = list.files("/home/a_tomer/Data/demo_patients_csv", full.names = T)

result = lapply(files, function(fileName){
  data = read.csv(fileName, header = T)
  
  data$log2psaplus1 = log(data$psa + 1, base = 2)
  data$high_dre = ifelse(data$dre=="T1c", 0, 1)
  
  totalRows = nrow(data)
  
  lapply(1:totalRows, function(rowNum){
    temp = data[1:rowNum,]
    
    curVisitTime = max(temp$visitTimeYears[!is.na(temp$psa) | !is.na(temp$dre)])
    lastBiopsyTime = max(temp$visitTimeYears[!is.na(temp$gleason)])
    
    if(curVisitTime<=lastBiopsyTime){
      curVisitTime = lastBiopsyTime + 0.1
    }

    summaryPlot = summaryGraph(temp, curVisitTime, lastBiopsyTime, FONT_SIZE=20, POINT_SIZE = 4, DRE_PSA_Y_GAP = 0.2)
    psaPredictionPlot = psaPredictionGraph(temp, FONT_SIZE=20, POINT_SIZE = 4)
    psaVelocityPlot = psaVelocityGraph(temp, FONT_SIZE=20, POINT_SIZE = 4)
    return(list(summaryPlot = summaryPlot, psaPredictionPlot = psaPredictionPlot, 
                psaVelocityPlot = psaVelocityPlot))
  })
})

names(result) = as.character(sapply(files, function(fileName){
  read.csv(fileName, header = T)[1,1]
}))


source("/home/a_tomer/Google Drive/PhD/src/prias/src/decision_analytic/predictPSADRE.R")
load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_light.Rdata")
load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/DRE_PSA/thresholdsList.Rdata")
load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/Simulation/scheduleResCombined.Rdata")
levels(scheduleResCombined$methodName)[5] = "Risk: Dynamic"
scheduleResCombined = scheduleResCombined[scheduleResCombined$methodName!=levels(scheduleResCombined$methodName)[6],]
scheduleResCombined$methodName = droplevels(scheduleResCombined$methodName)

print("Loaded Rdata")

THEME_COLOR = "#04A4DC"
LEGEND_POSITION = list(orientation = "h", y=1.125, x=0.4)

getBoxplotStatsDf=function(progression_time_low, progression_time_high, attribute){
  temp = scheduleResCombined[scheduleResCombined$progression_time>=progression_time_low & scheduleResCombined$progression_time<=progression_time_high,]
  res = unlist(by(temp$methodName, data = temp[, attribute], FUN = function(x){
    boxplot.stats(x)$stats
  }))
  
  resDf = data.frame(matrix(res, ncol=5, byrow = T))
  resDf$methodName = levels(scheduleResCombined$methodName)
  return(resDf)
}

psaObsDataGraph = function(data){
  data$Date = c(sapply(1:nrow(data), function(index){
    paste0(format(as.POSIXct(data$dom[index], origin = "1582-10-14"), format = "%b %e, %Y"),
           "\nYears since diagnosis: ", round(data$visitTimeYears[index],2), 
           " years\nPSA: ", data$psa[index], " ng/mL")
  }))
  
  plot = ggplot() + 
    geom_point(data=data[!is.na(data$psa),], aes(x=visitTimeYears, y=psa, label=Date), color=THEME_COLOR, size=3) +
    geom_line(data=data[!is.na(data$psa),], aes(x=visitTimeYears, y=psa), color=THEME_COLOR) + 
    geom_vline(aes(xintercept = data$visitTimeYears[!is.na(data$gleason)], 
                   linetype="Older biopsies")) +
    scale_linetype_manual(name="", values="dotted") +
    theme_bw()+ 
    theme(legend.position = "top") + 
    xlab("Follow-up time (years)") + 
    ylab("Observed PSA (ng/mL)")
  return(plot)
}

dreObsDataGraph = function(data){
  data$Date = c(sapply(1:nrow(data), function(index){
    paste0(format(as.POSIXct(data$dom[index], origin = "1582-10-14"), format = "%b %e, %Y"),
           "\nYears since diagnosis: ", round(data$visitTimeYears[index],2), 
           " years\nDRE: ", data$dre[index])
  }))
  
  plot = ggplot() + 
    geom_point(data=data[!is.na(data$dre),], aes(x=visitTimeYears, y=dre, label=Date), color=THEME_COLOR, shape=17, size=3) +
    geom_vline(aes(xintercept = data$visitTimeYears[!is.na(data$gleason)], 
                   linetype="Older biopsies")) +
    scale_linetype_manual(name="", values="dotted") +
    theme_bw()+ 
    xlab("Follow-up time (years)") + ylab("Observed DRE")
  return(plot)
}

psaPredictionGraph = function(data, FONT_SIZE=12, POINT_SIZE = 2){
  maxBiomarkerTime = max(data$visitTimeYears[!is.na(data$psa)])
  
  futureTimes = seq(0, 10, by = 0.1)
  predictedPSA_DRE = predictPSADRE(mvJoint_dre_psa_dre_value_light, data, idVar = "P_ID", survTimes = futureTimes)
  
  meanPredictedLog2psaplus1 = apply(predictedPSA_DRE$trueLog2psaplus1, MARGIN = 1, mean)
  lower95 = apply(predictedPSA_DRE$trueLog2psaplus1, MARGIN = 1, quantile, probs=0.025)
  upper95 = apply(predictedPSA_DRE$trueLog2psaplus1, MARGIN = 1, quantile, probs=0.975)
  
  predictionData = data.frame(visitTimeYears = futureTimes, 
                              meanPredictedLog2psaplus1, lower95, upper95)
  
  plot = ggplot() + 
    geom_point(data=data[!is.na(data$psa),], aes(x=visitTimeYears, y=log2psaplus1), color=THEME_COLOR, size=POINT_SIZE) +
    geom_line(data=predictionData[predictionData$visitTimeYears<=maxBiomarkerTime,], aes(x=visitTimeYears, y=meanPredictedLog2psaplus1, linetype="Fitted")) +
    geom_vline(aes(xintercept = data$visitTimeYears[!is.na(data$gleason)], 
                   linetype="Older biopsies"), show.legend=FALSE) +
    scale_linetype_manual(name="",
                          values=c("dashed", "dotted")) +
    theme_bw()+ 
    theme(text = element_text(size=FONT_SIZE),
          legend.background = element_blank(), legend.position = "top",
          legend.text = element_text(size=FONT_SIZE-3))  +
    xlab("Follow-up time (years)") + 
    ylab(expression('log'[2]*'(PSA + 1)')) 
  return(plot)
}

psaVelocityGraph = function(data, FONT_SIZE=12, POINT_SIZE = 2){
  maxBiomarkerTime = max(data$visitTimeYears[!is.na(data$psa)])
  
  futureTimes = seq(0, 10, by = 0.1)
  predictedPSA_DRE = predictPSADRE(mvJoint_dre_psa_dre_value_light, data, idVar = "P_ID", survTimes = futureTimes)
  
  meanPredictedLog2psaplus1 = apply(predictedPSA_DRE$trueLog2psaplus1_velocity, MARGIN = 1, mean)
  lower95 = apply(predictedPSA_DRE$trueLog2psaplus1_velocity, MARGIN = 1, quantile, probs=0.025)
  upper95 = apply(predictedPSA_DRE$trueLog2psaplus1_velocity, MARGIN = 1, quantile, probs=0.975)
  
  predictionData = data.frame(visitTimeYears = futureTimes, 
                              meanPredictedLog2psaplus1, lower95, upper95)
  
  plot = ggplot() + 
    geom_line(data=predictionData[predictionData$visitTimeYears<=maxBiomarkerTime,], aes(x=visitTimeYears, y=meanPredictedLog2psaplus1, linetype="Fitted")) +
    geom_vline(aes(xintercept = data$visitTimeYears[!is.na(data$gleason)], 
                   linetype="Older biopsies"), show.legend=FALSE) +
    scale_linetype_manual(name="",
                          values=c("solid", "dotted")) +
    theme_bw()+ 
    theme(text = element_text(size=FONT_SIZE),
          legend.background = element_blank(), legend.position = "top",
          legend.text = element_text(size=FONT_SIZE-3))  +
    xlab("Follow-up time (years)") + 
    ylab(expression('Rate of change (velocity) of log'[2]*'(PSA + 1)')) 
  return(plot)
}

drePredictionGraph = function(data){
  data$Info = c(sapply(1:nrow(data), function(index){
    paste0("\nTime = ", round(data$visitTimeYears[index],2), " years\nDRE = ", data$dre[index])
  }))
  
  maxBiomarkerTime = max(data$visitTimeYears[!is.na(data$dre)])
  
  futureTimes = seq(0, 10, by=0.1)
  predictedPSA_DRE = predictPSADRE(mvJoint_dre_psa_dre_value_light, data, idVar = "P_ID", survTimes = futureTimes)
  
  meanPredictedProb = plogis(apply(predictedPSA_DRE$trueDRELogOdds, MARGIN = 1, mean))
  lower95 = plogis(apply(predictedPSA_DRE$trueDRELogOdds, MARGIN = 1, quantile, probs=0.025))
  upper95 = plogis(apply(predictedPSA_DRE$trueDRELogOdds, MARGIN = 1, quantile, probs=0.975))
  
  predictionData = data.frame(visitTimeYears = futureTimes, 
                              meanPredictedProb, lower95, upper95)
  
  plot = ggplot() + 
    geom_point(data=data[!is.na(data$dre),], aes(x=visitTimeYears, y=high_dre, label=Info), color=THEME_COLOR, shape=17, size=3) +
    geom_line(data=predictionData[predictionData$visitTimeYears<=maxBiomarkerTime,], aes(x=visitTimeYears, y=meanPredictedProb, linetype="Fitted"), color=THEME_COLOR) +
    geom_line(data=predictionData[predictionData$visitTimeYears>maxBiomarkerTime,], aes(x=visitTimeYears, y=meanPredictedProb, linetype="Future prediction"), color=THEME_COLOR) +
    geom_ribbon(data=predictionData[predictionData$visitTimeYears>maxBiomarkerTime,],
                aes(x = visitTimeYears, ymin=lower95, ymax=upper95), alpha=0.5, fill="gray")+
    geom_vline(aes(xintercept = data$visitTimeYears[!is.na(data$gleason)], 
                   linetype="Older biopsies")) +
    scale_linetype_manual(name="", 
                          values=c("solid", "dashed", "dotted")) +
    scale_y_continuous(breaks=seq(0,1, length.out = 5), labels=paste0(seq(0,1, length.out = 5)*100, "%"),
                       limits=c(0,1)) +
    theme_bw()+ 
    xlab("Follow-up time (years)") + ylab("Probability of obtaining DRE more than T1c (%)")
  return(plot)
}

riskPredictionGraph = function(data){
  lastBiopsyTime = max(data$visitTimeYears[!is.na(data$gleason)])
  futureTimes = seq(lastBiopsyTime, 10, length.out = 20)
  sfit = survfitJM(mvJoint_dre_psa_dre_value_light, data, idVar="P_ID", 
                   survTimes = futureTimes, last.time = lastBiopsyTime)
  
  riskDf = data.frame(sfit$summaries[[1]])
  riskDf = rbind(c("times"=futureTimes[1], "Mean"=1, "Median"=1, "Lower"=1, "Upper"=1), riskDf)
  # -1 because 1st column is time
  riskDf[,-1] = 1 - riskDf[,-1]
  
  plot = ggplot() + geom_line(data = riskDf, aes(x=times, y=Mean, linetype="Risk of cancer progression (%)"), color="Red") +
    geom_ribbon(data = riskDf, aes(x=times, ymin=Lower, ymax=Upper), 
                fill="lightcoral", alpha=0.5) +  
    geom_ribbon(aes(x=seq(0, lastBiopsyTime, length.out = 10), ymin=0, ymax=1, fill="Cancer progression not detected"),
                alpha=0.25) + 
    geom_vline(aes(xintercept = data$visitTimeYears[!is.na(data$gleason)], 
                   linetype="Older biopsies")) +
    scale_linetype_manual(name="", 
                          values=c("dotted", "solid")) +
    scale_fill_manual(name="", 
                      values=c("limegreen")) +
    theme_bw()+
    theme(legend.title = element_blank()) + 
    scale_y_continuous(breaks=seq(0,1, length.out = 11), labels=paste0(seq(0,1, length.out = 11)*100, "%"),
                       limits=c(0,1), position = "left") +
    xlab("Follow-up time (years)") +
    ylab("Risk of cancer progression (%)")
  return(plot)
}

getSurvThreshold = function(lastBiopsyTime, curVisitTime){
  riskMethodName = "F1score"
  #Dt = curVisitTime - lastBiopsyTime
  Dt = 0.5
  
  availableDt = as.numeric(names(thresholdsList))
  Dt = availableDt[which.min(abs(Dt-availableDt))]
  
  cutpoints = thresholdsList[[as.character(Dt)]]$cut_points
  
  lastBiopsyTime = min(max(lastBiopsyTime, cutpoints[1,1]), cutpoints[nrow(cutpoints),1])
  
  #1st column in cutpoints is the column of times
  upperTimeIndex = which(lastBiopsyTime <= cutpoints[,1])[1]
  lowerTimeIndex = tail(which(lastBiopsyTime >= cutpoints[,1]),1)
  
  probDiffPerUnitTimeDiff = (cutpoints[upperTimeIndex, riskMethodName] - cutpoints[lowerTimeIndex, riskMethodName]) / (cutpoints[upperTimeIndex, 1] - cutpoints[lowerTimeIndex,1])
  
  diffProb = probDiffPerUnitTimeDiff * (lastBiopsyTime - cutpoints[lowerTimeIndex,1])
  
  if(is.nan(diffProb)){
    diffProb = 0
  }
  
  #2nd column in cutpoints is the column for F1 score thresholds
  return(cutpoints[lowerTimeIndex, riskMethodName] + diffProb)
}

riskColumnGraph = function(data, curVisitTime, riskThreshold, meanRiskProb, firstVisitDom){
  if(meanRiskProb >= riskThreshold){
    fillColor = "red"
  }else{
    fillColor = "green"
  }
  
  xlabel = format(as.POSIXct(curVisitTime * 365 * 24 * 60 * 60 + data$dom[1], origin = "1582-10-14"),
                  format = "%b %e, %Y")
  
  plot = ggplot() + geom_col(aes(x=c(xlabel,xlabel), y=c(meanRiskProb, 1-meanRiskProb)),
                                 fill=c(fillColor, "white"), color='black', width = 0.2) + 
    scale_y_continuous(breaks=seq(0,1, length.out = 11), labels=paste0(seq(0,1, length.out = 11)*100, "%"),
                       limits=c(0,1), position = "left") + theme_bw() + 
    theme(legend.position = "top", text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, color=fillColor)) + 
    xlab("") +
    ylab("Risk of cancer progression (%)") +
    ggtitle(paste0("Current risk: ", round(meanRiskProb * 100,2), "%"))
  
  plot = plot + geom_hline(aes(yintercept = riskThreshold, 
                               linetype=paste0('Maximum acceptable risk of \ncancer progression: ', round(riskThreshold,2)*100, "%")), 
                           color='black') + 
    scale_linetype_manual(name="", values = "dashed")
  
  return(plot)
}

summaryGraph = function(data, curVisitTime=10, 
                        FONT_SIZE=12, POINT_SIZE = 2, DRE_PSA_Y_GAP=0.1){
  patientDs =  data[!(is.na(data$dre) & is.na(data$psa)),]
  
  lastBiopsyTime = max(data$visitTimeYears[!is.na(data$gleason)])
  
  meanRiskProb = NA
  
  if(curVisitTime == lastBiopsyTime){
    meanRiskProb = 0
    curVisitTime = curVisitTime + 0.1
  }
  
  sfit = survfitJM(mvJoint_dre_psa_dre_value_light, patientDs, idVar="P_ID", 
                   survTimes = curVisitTime, last.time = lastBiopsyTime)
  
  patientDs$fitted_high_dre_prob = plogis(sfit$fitted.y[[1]]$high_dre)
  patientDs$fitted_log2psaplus1 = sfit$fitted.y[[1]]$log2psaplus1
  
  #The base of axes in this plot is of DRE
  minYLeft = 0
  maxYleft = 2 + DRE_PSA_Y_GAP * 2
  
  psaDs = patientDs[, c("visitTimeYears", "log2psaplus1", "fitted_log2psaplus1")]
  maxPSA = max(psaDs[,-1], na.rm=T)
  psaDs[,-1] = psaDs[,-1] / maxPSA + maxYleft/2 + DRE_PSA_Y_GAP
  
  if(is.na(meanRiskProb)){
    meanRiskProb = 1 - sfit$summaries[[1]][, "Mean"]
  }
  maxMeanRiskScaled = meanRiskProb * (maxYleft - minYLeft) + minYLeft
  
  curVisitDate = format(as.POSIXct(curVisitTime * 365 * 24 * 60 * 60 + data$dom[1], origin = "1582-10-14"),
                        format = "%b %e, %Y")
  
  if(curVisitTime - max(data$visitTimeYears) > 1.5){
    curVisitTimeTick = curVisitTime
  }else{
    curVisitTimeTick = max(data$visitTimeYears) + 1
  }
  
  if(curVisitTimeTick<5){
    col_width = 0.5
  }else{
    col_width = 1
  }
  
  xTicks = seq(0, max(data$visitTimeYears), length.out = 6)
  xTicks = xTicks[abs(xTicks - lastBiopsyTime) >= 0.5]
  xTicks = c(xTicks, lastBiopsyTime, curVisitTimeTick)
  xLabels = round(xTicks,1)
  xLabels[length(xLabels)-1] = paste0(round(lastBiopsyTime,2), "\n(Latest biopsy)")
  xLabels[length(xLabels)] = paste0(round(curVisitTime,2), "\n(Today: ", curVisitDate, ")")
  
  xLabelColors = c(rep("black", length(xTicks)-1), "red")
  
  p=ggplot() +
    geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
    geom_line(data = psaDs, aes(x = visitTimeYears, y=fitted_log2psaplus1, linetype="Fitted PSA")) +
    geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=log2psaplus1, shape="Observed PSA", color="Observed PSA")) +
    geom_line(data = patientDs, aes(x = visitTimeYears, y=fitted_high_dre_prob, linetype="Fitted DRE")) +
    geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre, shape="Observed DRE", color="Observed DRE")) +
    geom_col(aes(x=c(curVisitTimeTick,curVisitTimeTick), y=c(maxMeanRiskScaled, maxYleft-maxMeanRiskScaled)), 
             fill=c("red", "white"),  color=c('black'), width = col_width) + 
    geom_segment(aes(x=-Inf, xend=lastBiopsyTime, y=maxYleft/2, yend=maxYleft/2), 
                 linetype="solid", color="gray", size=1) +
    geom_text(aes(x=curVisitTimeTick, y=maxYleft - 0.5, 
                  label=paste0("Current risk\n", round(meanRiskProb * 100,2), "%")),
              color = "red", size=6)+
    xlab("Follow-up time (years)") + 
    ylab(expression('Pr (DRE > T1c)            '*'log'[2]*'(PSA + 1)')) +
    scale_linetype_manual(name="",
                          labels= c("Fitted Pr (DRE > T1c)", expression('Fitted log'[2]*'(PSA + 1)')),
                          values = c("dotted", "dashed")) +       
    scale_shape_manual(name="",
                       labels=c("Observed DRE", expression('Observed log'[2]*'(PSA + 1)')),
                       values = c(17,16)) + 
    scale_color_manual(name="",
                       labels=c("Observed DRE", expression('Observed log'[2]*'(PSA + 1)')),
                       values = c("darkorchid", THEME_COLOR)) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size=FONT_SIZE, color = rep(c("darkorchid", THEME_COLOR), each=4)),
          axis.title.y = element_text(size=FONT_SIZE, color = "black"),
          axis.title.y.right = element_text(size=FONT_SIZE, color = "red"),
          axis.text.y.right = element_text(size=FONT_SIZE, color = "red"),
          axis.text.x = element_text(size=FONT_SIZE, color=xLabelColors),
          legend.background = element_blank(), legend.position = "top",
          legend.text = element_text(size=FONT_SIZE-3))  +
    scale_x_continuous(breaks=xTicks, labels = xLabels, limits = c(0, curVisitTimeTick + col_width)) +
    scale_y_continuous(limits = c(minYLeft, maxYleft), 
                       breaks = c(seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 4), seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 4)),
                       labels = c(paste0(round(seq(0, 1, length.out = 4),2) * 100, "%"), 
                                  round(seq(0, maxPSA, length.out = 4),2)), 
                       sec.axis = sec_axis(~(.-minYLeft)/(maxYleft-minYLeft),
                                           breaks = seq(0,1, length.out = 5),
                                           labels = paste0(seq(0,1, length.out = 5)*100,"%"),
                                           name = "Risk of cancer progression"))
  
  return(p)
}

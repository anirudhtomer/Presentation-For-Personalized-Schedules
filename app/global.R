# source("/home/a_tomer/Google Drive/PhD/src/prias/src/decision_analytic/predictPSADRE.R")
# load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_light.Rdata")
# load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/DRE_PSA/thresholdsList.Rdata")
# load("/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/Simulation/scheduleResCombined.Rdata")
# levels(scheduleResCombined$methodName)[6] = "Risk: Dynamic"
# scheduleResCombined = scheduleResCombined[scheduleResCombined$methodName!=levels(scheduleResCombined$methodName)[5],]
# scheduleResCombined$methodName = droplevels(scheduleResCombined$methodName)

load("appdata.Rdata")
print("Loaded data necessary to run the app")

THEME_COLOR = "dodgerblue1"
THEME_COLOR_DARK = "dodgerblue4"
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

psaObsDataGraph = function(data, FONT_SIZE=15, POINT_SIZE=4){
  data$Date = c(sapply(1:nrow(data), function(index){
    paste0(format(as.POSIXct(data$dom[index], origin = "1582-10-14"), format = "%b %e, %Y"),
           "\nYears since diagnosis: ", round(data$visitTimeYears[index],2), 
           " years\nPSA: ", data$psa[index], " ng/mL")
  }))
  
  max_psa_val = max(data$psa, na.rm = T)
  last_biomarker_time = max(data$visitTimeYears[!is.na(data$psa) | !is.na(data$dre)])
  last_psa_val = tail(data$psa[!is.na(data$psa)],1)
  
  ybreaks = pretty(seq(0, max_psa_val, length.out = 4), n = 4)
  
  curVisitLabelYPosition = max(ybreaks) * 0.1
  curVisitLabelXPosition = max(0.1 * 10, last_biomarker_time)
  
  repeat{
    #10% space at least
    if(abs(last_psa_val - curVisitLabelYPosition)/max(ybreaks) > 0.1){
      break
    }
    curVisitLabelYPosition = curVisitLabelYPosition + 0.5
  }
  
  minY = - 0.1 * max(ybreaks)
  maxy = max(max(ybreaks), curVisitLabelYPosition + 1)
  
  plot = ggplot() + 
    geom_segment(aes(x = data$visitTimeYears[!is.na(data$gleason)],
                     xend = data$visitTimeYears[!is.na(data$gleason)],
                     y = -Inf, yend = 0,
                     linetype="Older biopsies"), size=0.7) +
    geom_segment(aes(x=last_biomarker_time, xend=last_biomarker_time, y=0, yend=Inf), 
                 linetype="twodash", color=THEME_COLOR_DARK, size=0.7)+
    geom_label(aes(x=curVisitLabelXPosition, y=curVisitLabelYPosition, 
                   label=paste0("Current visit\n", round(last_biomarker_time,1), " years")), 
               size=6, color=THEME_COLOR_DARK) +
    geom_hline(yintercept = 0, color="grey60") + 
    geom_point(data=data[!is.na(data$psa),], aes(x=visitTimeYears, y=psa, label=Date), color=THEME_COLOR, size=POINT_SIZE) +
    geom_line(data=data[!is.na(data$psa),], aes(x=visitTimeYears, y=psa), color=THEME_COLOR, alpha=0.3) + 
    scale_linetype_manual(name="", values="dotted") +
    theme_bw() + 
    theme(axis.text = element_text(size = FONT_SIZE),
          axis.title = element_text(size = FONT_SIZE),
          legend.text = element_text(size = FONT_SIZE - 2),
      legend.position = "bottom", legend.direction = "horizontal") + 
    xlim(0, max(10, last_biomarker_time + 1)) + 
    scale_y_continuous(breaks = ybreaks, limits = c(minY, maxy)) + 
    xlab("Follow-up time (years)") + 
    ylab("Observed PSA (ng/mL)")
  return(plot)
}

dreObsDataGraph = function(data, FONT_SIZE=15, POINT_SIZE=4){
  data$Date = c(sapply(1:nrow(data), function(index){
    paste0(format(as.POSIXct(data$dom[index], origin = "1582-10-14"), format = "%b %e, %Y"),
           "\nYears since diagnosis: ", round(data$visitTimeYears[index],2), 
           " years\nDRE: ", data$dre[index])
  }))
  
  dreTimes = data$visitTimeYears[!is.na(data$dre)]
  dreValues = as.numeric(data$dre[!is.na(data$dre)]) + 1
  
  ybreaks = min(dreValues):max(dreValues)
  ylabels = levels(data$dre[!is.na(data$dre)])
  minY = - 0.1 * max(ybreaks)
  
  last_biomarker_time = max(data$visitTimeYears[!is.na(data$psa) | !is.na(data$dre)])
  curVisitLabelYPosition = (max(ybreaks) + 1) * 0.1
  
  plot = ggplot() + 
    scale_linetype_manual(name="", values="dotted") +
    geom_segment(aes(x=last_biomarker_time, xend=last_biomarker_time, y=0, yend=Inf), 
                 linetype="twodash", color=THEME_COLOR_DARK, size=0.7)+
    geom_label(aes(x=last_biomarker_time, y=curVisitLabelYPosition, 
                   label=paste0("Current visit\n", round(last_biomarker_time,1), " years")), 
               size=6, color=THEME_COLOR_DARK) +
    geom_hline(yintercept = 0, color="grey60") + 
    geom_point(aes(x=dreTimes, y=dreValues), 
               color="darkorchid", shape=17, size=POINT_SIZE) +
    geom_segment(aes(x = data$visitTimeYears[!is.na(data$gleason)],
                     xend = data$visitTimeYears[!is.na(data$gleason)],
                     y = -Inf, yend = 0,
                     linetype="Older biopsies"), size=0.7) +
    theme_bw()+ 
    theme(axis.text = element_text(size = FONT_SIZE),
          axis.title = element_text(size = FONT_SIZE),
          legend.text = element_text(size = FONT_SIZE - 2),
          legend.position = "bottom", legend.direction = "horizontal") + 
    scale_y_continuous(breaks = ybreaks, labels = ylabels, limits = c(minY, max(ybreaks)+1)) + 
    xlim(0,10) +
    xlab("Follow-up time (years)") + ylab("Observed DRE")
  return(plot)
}

psaPredictionGraph = function(data, FONT_SIZE=12, POINT_SIZE = 2){
  last_biomarker_time = max(data$visitTimeYears[!is.na(data$psa) | !is.na(data$dre)])
  max_x_time = max(data$visitTimeYears, na.rm = T)
  
  futureTimes = seq(0, max_x_time + 1, length.out = 40)
  predictedPSA_DRE = predictPSADRE(mvJoint_dre_psa_dre_value_light, data, idVar = "P_ID", survTimes = futureTimes)
  
  meanPredictedLog2psaplus1 = apply(predictedPSA_DRE$trueLog2psaplus1, MARGIN = 1, mean)
  lower95 = apply(predictedPSA_DRE$trueLog2psaplus1, MARGIN = 1, quantile, probs=0.025)
  upper95 = apply(predictedPSA_DRE$trueLog2psaplus1, MARGIN = 1, quantile, probs=0.975)
  
  max_psa_val = max(data$log2psaplus1, meanPredictedLog2psaplus1, na.rm = T)
  ybreaks = pretty(seq(0, max_psa_val, length.out = 4), n = 4)
  minY = - 0.1 * max(ybreaks)
  
  predictionData = data.frame(visitTimeYears = futureTimes, 
                              meanPredictedLog2psaplus1, lower95, upper95)
  
  plot = ggplot() + 
    geom_point(data=data[!is.na(data$psa),], aes(x=visitTimeYears, y=log2psaplus1), color=THEME_COLOR, size=POINT_SIZE) +
    geom_line(data=predictionData[predictionData$visitTimeYears<=last_biomarker_time,],
              aes(x=visitTimeYears, y=meanPredictedLog2psaplus1, linetype="Underlying profile"), size=0.7) +
    geom_segment(aes(x = data$visitTimeYears[!is.na(data$gleason)],
                     xend = data$visitTimeYears[!is.na(data$gleason)], y = -Inf, yend=0,
                   linetype="Older biopsies"), show.legend=FALSE, size=0.7) +
    geom_hline(yintercept = 0, color="grey60") + 
    scale_linetype_manual(name="",
                          labels=c("Older biopsies", expression('Underlying log'[2]*'(PSA + 1) profile')),
                          values=c("dotted", "dashed")) +
    theme_bw()+ 
    theme(text = element_text(size=FONT_SIZE),
          legend.background = element_blank(), 
          legend.text = element_text(size = FONT_SIZE - 2),
          legend.position = "bottom", legend.direction = "horizontal") + 
    scale_y_continuous(breaks = ybreaks, limits = c(minY,  max(ybreaks))) + 
    xlab("Follow-up time (years)") + 
    xlim(0, max_x_time) + 
    ylab(expression('log'[2]*'(PSA + 1)')) 
  return(plot)
}

psaVelocityGraph = function(data, FONT_SIZE=12, POINT_SIZE = 2){
  last_biomarker_time = max(data$visitTimeYears[!is.na(data$psa) | !is.na(data$dre)])
  max_x_time = max(data$visitTimeYears, na.rm = T)
  
  futureTimes = seq(0, max_x_time + 1, length.out = 40)
  predictedPSA_DRE = predictPSADRE(mvJoint_dre_psa_dre_value_light, data, idVar = "P_ID", survTimes = futureTimes)
  
  meanPredictedLog2psaplus1 = apply(predictedPSA_DRE$trueLog2psaplus1_velocity, MARGIN = 1, mean)
  lower95 = apply(predictedPSA_DRE$trueLog2psaplus1_velocity, MARGIN = 1, quantile, probs=0.025)
  upper95 = apply(predictedPSA_DRE$trueLog2psaplus1_velocity, MARGIN = 1, quantile, probs=0.975)
  
  predictionData = data.frame(visitTimeYears = futureTimes, 
                              meanPredictedLog2psaplus1, lower95, upper95)
  
  ybreaks = pretty(seq(min(meanPredictedLog2psaplus1), max(meanPredictedLog2psaplus1), length.out = 4), n = 4)
  minY = min(ybreaks) - 0.1 * (max(ybreaks) - min(ybreaks))
  
  plot = ggplot() + 
    geom_line(data=predictionData[predictionData$visitTimeYears<=last_biomarker_time,], 
              aes(x=visitTimeYears, y=meanPredictedLog2psaplus1, linetype="Underlying velocity")) +
    geom_segment(aes(x = data$visitTimeYears[!is.na(data$gleason)],
                     xend = data$visitTimeYears[!is.na(data$gleason)], 
                     y = -Inf, yend=min(ybreaks),
                     linetype="Older biopsies"), show.legend=FALSE, size=0.7) +
    geom_hline(yintercept = min(ybreaks), color="grey60") + 
    scale_linetype_manual(name="", 
                          labels=c("Older biopsies", expression('Underlying log'[2]*'(PSA + 1) velocity')),
                          values=c("dotted", "solid")) +
    theme_bw()+ 
    theme(text = element_text(size=FONT_SIZE),
          legend.background = element_blank(), legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size=FONT_SIZE-2))  +
    scale_y_continuous(breaks=ybreaks, limits=c(minY, max(ybreaks)))+
    xlab("Follow-up time (years)") + 
    xlim(0, max_x_time) + 
    ylab(expression('Velocity of log'[2]*'(PSA + 1)')) 
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
  Dt = curVisitTime - lastBiopsyTime
  #Dt = 0.5
  
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
    fillColor = "red3"
  }else{
    fillColor = "forestgreen"
  }
  
  xlabel = format(as.POSIXct(curVisitTime * 365 * 24 * 60 * 60 + data$dom[1], origin = "1582-10-14"),
                  format = "%b %e, %Y")
  xlabel = paste0(xlabel, "\n(", round(curVisitTime,1), " years since diagnosis)")
  
  plot = ggplot() + geom_col(aes(x=c(xlabel,xlabel), y=c(meanRiskProb, 1-meanRiskProb)),
                                 fill=c(fillColor, "white"), color=fillColor, width = 0.2) + 
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"),
                       limits=c(0,1), sec.axis = dup_axis(trans = ~.,
                                                          breaks = riskThreshold,
                                                          labels = paste0(round(riskThreshold*100,2), "%\n(Biopsy threshold)"))) + 
                         theme_bw() + 
    geom_label(aes(x=xlabel, y=0.75, 
                   label=paste0("Current risk\n", round(meanRiskProb * 100,2), "%\n")),
               color = fillColor, size=5)+
    theme(legend.position = "top", text = element_text(size = 20),
          axis.text.y.right = element_text(color = "firebrick1"),
          axis.title.y.right = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, color=fillColor)) + 
    geom_segment(aes(x=-Inf, xend=Inf, y=riskThreshold, yend=riskThreshold), 
                 linetype="dashed", color="firebrick1", size=0.5) +
    xlab("") +
    ylab("Risk of cancer progression (%)")
  return(plot)
}

summaryGraph = function(data, curVisitTime=10, lastBiopsyTime,
                        FONT_SIZE=12, POINT_SIZE = 2, DRE_PSA_Y_GAP=0.1){
  
  curVisitDate = format(as.POSIXct(curVisitTime * 365 * 24 * 60 * 60 + data$dom[1], origin = "1582-10-14"),
                        format = "%b %e, %Y")
  
  patientDs =  data[!(is.na(data$dre) & is.na(data$psa)),]
  
  sfit = survfitJM(mvJoint_dre_psa_dre_value_light, patientDs, idVar="P_ID", 
                   survTimes = curVisitTime, last.time = lastBiopsyTime)
  
  meanRiskProb = 1 - sfit$summaries[[1]][, "Mean"]
  
  patientDs$fitted_high_dre_prob = plogis(sfit$fitted.y[[1]]$high_dre)[1:nrow(patientDs)]
  patientDs$fitted_log2psaplus1 = sfit$fitted.y[[1]]$log2psaplus1[1:nrow(patientDs)]
  
  psaDs = patientDs[, c("visitTimeYears", "log2psaplus1", "fitted_log2psaplus1")]
  
  maxPSA = max(psaDs[,-1], na.rm=T)
  minPSA = min(psaDs[,-1], na.rm = T)
  
  #real PSA, and fake DRE values
  psaPlot = ggplot() + 
    geom_vline(aes(xintercept = lastBiopsyTime, linetype="Latest biopsy"), show.legend = F, size=0.7) +
    geom_vline(aes(xintercept = curVisitTime), show.legend = F, color=THEME_COLOR_DARK, linetype="twodash", size=0.7) +
    geom_point(aes(x = 0, y=900, shape="Observed DRE", color="Observed DRE")) +
    geom_line(aes(x = 0, y=900, linetype="Fitted")) +
    geom_label(aes(x=curVisitTime, y=0.5, 
                   label=paste0("Current visit\n", round(curVisitTime,1), " years")), 
               size=5, color=THEME_COLOR_DARK) +
    geom_line(data = psaDs, aes(x = visitTimeYears, y=fitted_log2psaplus1, linetype="Fitted")) +
    geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=log2psaplus1, shape="Observed PSA", color="Observed PSA")) +
    scale_linetype_manual(name="",
                          labels= c(expression(atop('Underlying log'[2]*'(PSA + 1)',  'Underlying Pr (DRE > T1c)')), "Latest biopsy"),
                          values = c("dashed", "dotted")) +       
    scale_shape_manual(name="",
                       labels=c("Observed DRE (T1c / above T1c)", expression('Observed log'[2]*'(PSA + 1)')),
                       values = c(17,16)) + 
    scale_color_manual(name="",
                       labels=c("Observed DRE (T1c / above T1c)", expression('Observed log'[2]*'(PSA + 1)')),
                       values = c("darkorchid", THEME_COLOR)) + 
    ylab(expression('log'[2]*'(PSA + 1)')) + 
    ylim(0, maxPSA) +
    xlim(0, curVisitTime + curVisitTime * 0.1) + 
    xlab("Follow-up time (years)") + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), 
          axis.line = element_line(),
          legend.text = element_text(size = FONT_SIZE-4),
          legend.position = "bottom", legend.direction = "horizontal",
          plot.margin = margin(0, 5.5, 2, 5.5, "pt")) 
  
  drePlot = ggplot() + geom_vline(aes(xintercept = lastBiopsyTime), linetype="dotted", size=0.7) +
    geom_vline(aes(xintercept = curVisitTime), show.legend = F, color=THEME_COLOR_DARK, linetype="twodash", size=0.7) +
    geom_label(aes(x=curVisitTime, y=0.5, 
                   label=paste0("Current visit\n", round(curVisitTime,1), " years")), 
               size=5, color=THEME_COLOR_DARK) +
    geom_line(data = patientDs, aes(x = visitTimeYears, y=fitted_high_dre_prob), linetype="dashed") +
    geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre), shape=17, color="darkorchid") +
    ylab("Pr (DRE > T1c)") + 
    scale_y_continuous(breaks = seq(0,1, length.out = 3), 
                       labels=paste0(seq(0,1,length.out = 3)*100, "%"), limits = c(0,1)) +
    xlim(0, curVisitTime + curVisitTime * 0.1) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), 
          axis.line = element_line(),
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          legend.position = "none",
          plot.margin = margin(0, 5.5, 0, 5.5, "pt")) 
  
  riskGaugeLabel = paste0("Risk of cancer progression on current visit \n", curVisitDate, " (", round(curVisitTime,1)," years)")
  riskGauge = ggplot(data = NULL, 
                     aes(ymax = meanRiskProb, ymin = 0, xmax = 2, xmin = 1, fill="Risk")) +
    geom_rect(aes(ymax=1, ymin=0, xmax=2, xmin=1), fill ="white", color="red3") +
    geom_rect() +
    scale_fill_manual("", values="red3") +
    coord_polar(theta = "y",start=-pi/2) + xlim(c(0, 2)) + ylim(c(0,2)) +
    geom_text(aes(x = 0, y = 0, label = paste0(round(meanRiskProb*100,2),"%")), color='red3', size=6) +
    geom_text(aes(x=0.5, y=1.5, label = riskGaugeLabel), size=6, color="red3") +
    theme_void() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    guides(fill=FALSE) +
    guides(colour=FALSE)
  
  p=ggpubr::ggarrange(drePlot, psaPlot, nrow=2, ncol=1, heights = c(0.7,1),
                      align = "v")
  p=ggpubr::ggarrange(p, riskGauge, nrow=1, ncol=2, widths = c(1.75,1))
  
  return(p)
}
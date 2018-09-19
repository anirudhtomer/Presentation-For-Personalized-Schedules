#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)
library(JMbayes)
library(splines)
library(ggplot2)
library(plotly)
library(DT)
library(ggpubr)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  output$table_obs_data <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$patientFile
    
    if (is.null(inFile))
      return(NULL)
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    
    dataToShow = data[1, ]
    dataToShow$lastBiopsyDate = format(as.POSIXct(max(data$dom[!is.na(data$gleason)]), origin = "1582-10-14"), format = "%b %e, %Y")
    dataToShow$firstVisitDom_human = format(as.POSIXct(dataToShow$firstVisitDom, origin = "1582-10-14"), format = "%b %e, %Y")
    dataToShow$progressed = ifelse(dataToShow$progressed==0, "No", "Yes")
    dataToShow$lastVisitDom_human = format(as.POSIXct(max(data$dom), origin = "1582-10-14"), format = "%b %e, %Y")
    
    dataToShow = dataToShow[, c("P_ID", "Age", "firstVisitDom_human", "lastVisitDom_human", "lastBiopsyDate", "progressed", "nr_visits")]
    
    colnames(dataToShow) = c("Patient ID", "Age (years)", "First visit", "Last visit", "Last biopsy", "Cancer progression", "Nr. of visits")
    
    return(dataToShow)
  })
  
  output$graph_obs_psa <- renderPlot({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      output$graph_obs_dre = NULL
      return(NULL)
    }
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    
    # output$graph_obs_dre = renderPlotly({ggplotly(dreObsDataGraph(data), tooltip=c("label")) %>% 
    #         config(displayModeBar = F) %>% 
    #         layout(legend = LEGEND_POSITION)})
    
    # print(ggplotly(psaObsDataGraph(data), tooltip=c("label")) %>% 
    #         config(displayModeBar = F, mathjax = 'cdn') %>% 
    #         layout(legend = LEGEND_POSITION))
    
    output$graph_obs_dre = renderPlot(dreObsDataGraph(data))
    
    return(psaObsDataGraph(data))
  })
  
  output$graph_prediction <- renderPlot({
    inFile <- input$patientFile
    
    if (is.null(inFile))
      return(NULL)
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data$log2psaplus1 = log(data$psa + 1, base = 2)
    data$high_dre = ifelse(data$dre=="T1c", 0, 1)
    
    if(input$pred_type=="SUMMARY"){
      curVisitTime = max(data$visitTimeYears[!is.na(data$psa) | !is.na(data$dre)])
      lastBiopsyTime = max(data$visitTimeYears[!is.na(data$gleason)])
      
      if(curVisitTime<=lastBiopsyTime){
        curVisitTime = lastBiopsyTime + 0.1
      }
      
      plot = summaryGraph(data, curVisitTime, lastBiopsyTime, FONT_SIZE=20, POINT_SIZE = 4, DRE_PSA_Y_GAP = 0.2)
    }else if(input$pred_type=="PSA"){
      plot = psaPredictionGraph(data, FONT_SIZE=20, POINT_SIZE = 4)
    }else{
      plot = NULL
    }
    
    return(plot)
  })
  
  output$graph_prediction_psa_velocity <- renderPlot({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$patientFile
    
    if (is.null(inFile))
      return(NULL)
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data$log2psaplus1 = log(data$psa + 1, base = 2)
    data$high_dre = ifelse(data$dre=="T1c", 0, 1)
    
    if(input$pred_type=="PSA"){
      return(psaVelocityGraph(data, FONT_SIZE=20, POINT_SIZE = 4))
    }else {
      return(NULL)
    }
  })
  
  output$biopsy_threshold_reason = renderUI({
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      return(NULL)
    }
    
    if(input$risk_choice_biopsy=="5_PERC"){
      reason = "Fixed biopsy threshold."
    }else if(input$risk_choice_biopsy=="15_PERC"){
      reason = "Fixed biopsy threshold."
    }else{
      reason = "The currently chosen biopsy threshold is a follow-up time dependent 
risk threshold. Using information from the observed cancer progression times in 
the PRIAS dataset, at a follow-up visit this threshold was chosen because it best discriminated
between patients who obtain cancer progression versus others. We chose the well known F1 score as the discrimination index for this purpose.
This is because the F1 score specifically focuses on identifying patients who may observe cancer progression by combining true positive rate, 
and positive predictive value."
    }
    
    return(HTML(reason))
  })
  
  output$graph_risk_now = renderPlot({
    inFile <- input$patientFile
    
    if (is.null(inFile)){
      output$biopsy_decision_yes = renderText("")
      output$biopsy_decision_no = renderText("")
      output$biopsy_decision_reason = renderUI("")
      return(NULL)
    }
    
    data = read.csv(inFile$datapath, header=TRUE, dec = input$dec,
                    sep = input$sep, quote = input$quote)
    data$log2psaplus1 = log(data$psa + 1, base = 2)
    data$high_dre = ifelse(data$dre=="T1c", 0, 1)
    
    curVisitTime = max(data$visitTimeYears[!is.na(data$psa) | !is.na(data$dre)])
    lastBiopsyTime = max(data$visitTimeYears[!is.na(data$gleason)])
    
    if(curVisitTime<=lastBiopsyTime){
      curVisitTime = lastBiopsyTime + 0.1
    }
    
    if(input$risk_choice_biopsy=="5_PERC"){
      riskThreshold = 0.05
    }else if(input$risk_choice_biopsy=="15_PERC"){
      riskThreshold = 0.15
    }else{
      riskThreshold = 1 - getSurvThreshold(lastBiopsyTime, curVisitTime)
    }
    
    if(curVisitTime==lastBiopsyTime){
      meanRiskProb = 0
    }else{
      sfit = survfitJM(mvJoint_dre_psa_dre_value_light, data, idVar="P_ID", 
                       survTimes = curVisitTime, last.time = lastBiopsyTime)
      meanRiskProb = 1 - sfit$summaries[[1]][, "Mean"]
    }
    
    if(meanRiskProb >= riskThreshold){
      if(curVisitTime - lastBiopsyTime >= 1 | input$year_gap_biopsy=="No"){
        output$biopsy_decision_yes = renderText("Biopsy recommended")
        output$biopsy_decision_no = renderText("")
        output$biopsy_decision_reason = renderUI(HTML(paste0("<b> Reason: </b> Patient's current risk of cancer progression is ",
                                                             round(meanRiskProb*100,2), "%, which is more than the maximum acceptable risk (biopsy threshold) of ",
                                                             round(riskThreshold*100,2), "%.")))
        
      }else{
        output$biopsy_decision_yes = renderText("")
        output$biopsy_decision_no = renderText("Biopsy not recommended")
        output$biopsy_decision_reason = renderUI(HTML(paste0("<b> Reason: </b> Although patient's current risk of cancer progression of ",
                                                             round(meanRiskProb*100,2), "%, is more than the maximum acceptable risk (biopsy threshold) of ",
                                                             round(riskThreshold*100,2), "%, the last biopsy was conducted less than 1 year ago.")))
      }
    }else{
      output$biopsy_decision_yes = renderText("")
      output$biopsy_decision_no = renderText("Biopsy not recommended")
      output$biopsy_decision_reason = renderUI(HTML(paste0("<b> Reason: </b> Patient's current risk of cancer progression is ",
                                                           round(meanRiskProb*100,2), "%, which is less than the maximum acceptable risk (biopsy threshold) of ",
                                                           round(riskThreshold*100,2), "%.")))
      
    }
    
    plot = riskColumnGraph(data, curVisitTime, riskThreshold, meanRiskProb)
    
    return(plot)
  })
  
  output$nb = renderPlot({
    FONT_SIZE = 20
    
    nbData = getBoxplotStatsDf(0,10, "nb")
    offsetData = getBoxplotStatsDf(0,10-0.1, "offset")
    
    if(input$risk_choice_impact=="5_PERC"){
      nbData$fill = c("", "", 
                      "", "Currently selected (5% risk)", "")
      offsetData$fill = nbData$fill
    }else if(input$risk_choice_impact=="15_PERC"){
      nbData$fill = c("", "", 
                      "Currently selected (15% risk)", "", "")
      offsetData$fill = nbData$fill
    }else{
      nbData$fill = c("", "", 
                      "", "", 
                      "Currently selected (dynamic risk)")
      offsetData$fill = nbData$fill
    }
    
    nbPlot = ggplot(data=nbData) + 
      geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, 
                       x=methodName, fill=fill),
                   stat = "identity") + coord_flip() + 
      scale_fill_manual("", values=c("White", THEME_COLOR)) +
      theme_bw() + 
      theme(text = element_text(size=FONT_SIZE), 
            axis.line = element_line(),
            legend.background = element_blank(), legend.position = "top",
            legend.text = element_text(size=FONT_SIZE-3))  +
      scale_y_continuous(breaks = seq(1,10,2)) + 
      xlab("Method") + ylab("Number of Biopsies")
    
    output$offset = renderPlot({
      ggplot(data=offsetData) + 
        geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, 
                         x=methodName, fill=fill),
                     stat = "identity") + coord_flip() + 
        scale_fill_manual("", values=c("White", THEME_COLOR)) +
        theme_bw() + 
        theme(text = element_text(size=FONT_SIZE),
              axis.line = element_line(),
              legend.background = element_blank(), legend.position = "top",
              legend.text = element_text(size=FONT_SIZE-3))  +
        scale_y_continuous(breaks = seq(0,ceiling(max(getBoxplotStatsDf(0,10-0.1, "offset")[,-6])), 0.5)) + 
        xlab("Method") + ylab("Delay in detection of\ncancer progression (years)")
    })
    
    return(nbPlot)
  })
  
})

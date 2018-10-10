#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)
library(ggplot2)
library(plotly)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  useShinyalert(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main.css")
  ),
  
  # Application title
  tags$div(class="text-center main-header", 
           headerPanel("Biopsy Recommender for Prostate Cancer Patients in Active Surveillance ")),
  sidebarLayout(
    # Sidebar with a slider input for number of bins 
    sidebarPanel(
      
      #fileInput('RDfile', 'Load the R Workspace with the fitted joint model',
      #          accept = NULL),
      
      fileInput('patientFile', 'Load patient data (CSV file)',
                accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      radioButtons('dec', 'Decimal', c(Dot = '.', Comma = ','), 
                   selected='.'),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      width=2
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Patient Data", value="patient_data", tags$br(), 
                           verticalLayout(tableOutput("table_obs_data"), tags$hr(),
                                          splitLayout(tags$h4(class="observed_graph_title", "Observed PSA measurements"),
                                                      tags$h4(class="observed_graph_title", "Observed DRE measurements")),
                                          splitLayout(plotOutput("graph_obs_psa", height = "450px"),
                                                      plotOutput("graph_obs_dre", height = "450px"))
                           )),
                  tabPanel("Predictions", tags$br(), verticalLayout(
                    radioButtons("pred_type", "Choose the type of prediction",
                                 choices = c("Risk of cancer progression"="SUMMARY",
                                             "Prostate-specific antigen" = "PSA"),
                                 selected="SUMMARY", inline=T),
                    sliderInput("visitNumber", "Number of follow-up visits:",
                        min = 2, max = 10, value = 2, step = 1,
                        animate = animationOptions(interval = 2000, loop = F, 
                                                   playButton = NULL, pauseButton = NULL)),
                    plotOutput("graph_prediction",height="500px"),
                    tags$br(),
                    plotOutput("graph_prediction_psa_velocity")
                  )),
                  tabPanel("Biopsy recommendation", tags$br(), verticalLayout(
                    radioButtons("year_gap_biopsy",
                                 "Should there be a gap of 1 year between consecutive biopsies?",
                                 choices = c("Yes"="Yes",
                                             "No" = "No"),
                                 selected="Yes", inline=T),
                    radioButtons("risk_choice_biopsy",
                                 "Choose the maximum risk of cancer progression (biopsy threshold) you/patient are willing to accept?",
                                 choices = c("Follow-up time dependent threshold based on PRIAS dataset"="AUTO", "5%"="5_PERC", "15%" = "15_PERC"),
                                 selected="AUTO", inline=T),
                    tags$span(id="impact_check_msg","Please check the impact (graphs below) of each risk threshold, before making a selection."),
                    tags$br(),
                    wellPanel(fluidRow(column(width = 5, plotOutput("graph_risk_now", height="500px")),
                                       column(width = 7, verticalLayout(
                                         textOutput("biopsy_decision_yes"), textOutput("biopsy_decision_no"), 
                                         htmlOutput("biopsy_decision_reason"), tags$br(), 
                                         tags$span(id="biopsy_threshold_reason_title","Explanation of the biopsy threshold"),
                                         htmlOutput("biopsy_threshold_reason")))))
                  ), tags$hr(), tags$br(), titlePanel("Impact of the choice of risk thresholds"),
                  tags$div(id="info_boxplots", paste("Different risk thresholds lead to different suggestions for conducting biopsies.",
                                   "In order to evaluate the impact of each threshold, we conducted an extensive and realistic simulation study based on the PRIAS dataset.",
                                   "In the simulation study we compared existing methods for biopsies, such as using PRIAS or annual schedule of biopsies, with the personalized approaches based on risk of cancer progression.",
                                   "The comparison is done on the basis of number of biopsies conducted in order to detect cancer progression, and the delay in detection of cancer progression.",
                                   "The delay is defined as the time difference between the time of last biopsy at which cancer progression is detected and the true time of cancer progression.",
                                  "For various methods, we present boxplots for number of biopsies and delay in detection of cancer progression below.",
                                  "The boxplots are obtained after the methods are tested on simulated patients.")),
                  tags$br(), radioButtons("risk_choice_impact",
                                          "Choose a risk threshold to see its impact",
                                          choices = c("Dynamic risk based on PRIAS dataset"="AUTO", "5%"="5_PERC", "15%" = "15_PERC"),
                                          selected="AUTO", inline=T),
                  splitLayout(plotOutput("nb"), plotOutput("offset")))
                  ),
                  width = 10
      )
    )
  ))
  
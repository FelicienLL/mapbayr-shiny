library(shiny)
library(mapbayr)
library(mrgsolve)
# --- First: Drug specific elements, depending on your drug, model, protocol etc...

# - mrgsolve model

my_model <- mread("mrg_901.cpp")

# - a posteriori "Adaptation" function(s): return a dose recommendation, a comment, a specific figure...

adapt_carboplatin <- function(est, target_auc, dose_d1){
  stopifnot(inherits(est, "mapbayests"))
  CL <- get_param(est, "CL")*(1000/60)
  
  AUC_D1 <- dose_d1 / CL
  AUC_D2 <- AUC_D1
  AUC_D3 <- target_auc - (AUC_D2 + AUC_D1)
  DOSE_D3 <- CL * AUC_D3
  if(DOSE_D3<=0) DOSE_D3 <- 0
  
  paste0("Patient's clearance: ", round(CL), " mL/min.", "\n",
         "AUC Day 1: ", round(AUC_D1, 2), " mg.min/mL.", "\n",
         "Dose at day 3 to reach a target AUC of ", target_auc, " mg.min/mL over 3 days: ", round(DOSE_D3), " mg.")
}

# --- Secondly: The shiny app per se



ui <- fluidPage(
  titlePanel("High dose carboplatin TDM - mapbayr"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        h4("Administration"),
        column(width = 6,
               numericInput("amt", "Dose (mg)", 1000)
        ),
        column(width = 6,
               numericInput("dur", "Duration (h)", 1, min = .1)
        )
      ),
      fluidRow(
        h4("Observations"),
        column(width = 6,
               numericInput("time1", "Time 1 (h)", .95),
               numericInput("time2", "Time 2 (h)", 1.9),
               numericInput("time3", "Time 3 (h)", 4.7)
        ),
        column(width = 6,
               numericInput("dv1", "Conc 1 (mg/L)", 35.3),
               numericInput("dv2", "Conc 2 (mg/L)", 22.1),
               numericInput("dv3", "Conc 3 (mg/L)", 6.7)
        )
      ),
      fluidRow(
        h4("Patient/Protocol"),
        column(width = 6,
               numericInput("bsa","BSA (m2)", 1.73)),
        column(width = 6,
               numericInput("auc", "Target AUC (mg.min/mL)", 24))
      ),
      fluidRow(
        h4("Perfom estimation"),
        actionButton("GO", "Let's GO !")
      )
    ),
    mainPanel(
      fluidRow(
        h4("Data"),
        tableOutput("mapbay_tab")
      ),
      fluidRow(
        h4("Results"),
        verbatimTextOutput("results")
      ),
      fluidRow(
        h4("Figures"),
        tabsetPanel(
          tabPanel("Concentrations",  plotOutput("concvstime")),
          tabPanel("Parameters",  plotOutput("distparam"))
        )
      )
    )
  )
)

server <- function(input, output) {
  my_data <- reactive({
    my_model %>%
      adm_lines(time = 0, amt = input$amt, rate = (input$amt)/(input$dur)) %>%
      obs_lines(time = c(input$time1, input$time2, input$time3),
                DV = c(input$dv1, input$dv2, input$dv3)) %>%
      add_covariates(covariates = list(BSA = input$bsa)) %>%
      get_data()
  })
  
  my_est <- eventReactive(input$GO, {
    mapbayest(my_model, my_data(), verbose = F)
  })
  
  output$mapbay_tab <- renderTable({
    if(input$GO == 0){
      my_data()
    } else {
      as.data.frame(my_est())
    }
  })
  
  output$results <- renderText({
    adapt_carboplatin(
      est = my_est(),
      target_auc = isolate(input$auc),
      dose_d1 = isolate(input$amt))
  })
  
  output$concvstime <- renderPlot({
   shiny::req(my_est())
   plot(my_est())+
     ggplot2::scale_x_continuous(limits = c(0, 12))
  })
  
  output$distparam <- renderPlot({
    mapbayr:::hist.mapbayests(my_est())
  })
  
  
}

shinyApp(ui = ui, server = server) 

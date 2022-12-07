library(shiny)
library(mapbayr)
library(mrgsolve)
# --- First: Drug specific elements, depending on your drug, model, protocol etc...

# - mrgsolve model

my_model <- mread("ibru_911.cpp")

# - a posteriori "Adaptation" function(s): to return a dose recommendation, a comment, a specific figure...

adapt_ibru <- function(est, dose){
  stopifnot(inherits(est, "mapbayests"))
  CLF <- get_param(est, "CLF")
  AUC <- 1000*dose / CLF
  
  paste0("Oral clearance: ", round(CLF, 1), " L/h.", "\n",
         "AUC: ", round(AUC, 1), " ng.h/mL.", "\n",
         "In chronic lymphocytic leukemia patients (420mg/j), mean observed AUC was 680 ng.h/mL"
  )
}

# --- Secondly: The shiny app per se



ui <- fluidPage(
  titlePanel("Ibrutinib AUC - mapbayr"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        h4("Administration"),
        column(width = 6,
               numericInput("amt", "Dose (mg)", 420, step = 140)  
        )
      ),
      fluidRow(
        h4("Observations"),
        column(width = 4,
               numericInput("time1", "T0 (h)", 0),
               numericInput("time2", "T2 (h)", 2.2),
               numericInput("time3", "T4 (h)", 3.8)
        ),
        column(width = 4,
               numericInput("dv1", "Ibru T0 (ng/mL)", 4.8),
               numericInput("dv2", "Ibru T2 (ng/mL)", 120.4),
               numericInput("dv3", "Ibru T4 (ng/mL)", 81.1)
        ),
        column(width = 4,
               numericInput("dv1met", "DHibru T0 (ng/mL)", 11.4),
               numericInput("dv2met", "DHibru T2 (ng/mL)", 160.6),
               numericInput("dv3met", "DHibru T4 (ng/mL)", 137.1)
        )
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
      adm_lines(time = 0, amt = input$amt, ss = 1, ii = 24) %>%
      adm_lines(time = 24, amt = input$amt, ss = 0) %>%
      obs_lines(time = 24+c(input$time1, input$time2, input$time3),
                DV = c(input$dv1, input$dv2, input$dv3),
                DVmet = c(input$dv1met, input$dv2met, input$dv3met)) %>%
      add_covariates(covariates = list(BSA = input$bsa)) %>%
      get_data()
  })
  
  my_est <- eventReactive(input$GO, {
    mapbayest(my_model, my_data(), verbose = F, check = F)
  })
  
  output$mapbay_tab <- renderTable({
    if(input$GO == 0){
      my_data()
    } else {
      as.data.frame(my_est()) %>% 
        dplyr::mutate(dplyr::across(c("ID", "evid", "cmt", "mdv", "rate", "amt", "ii", "ss"), as.integer))
    }
  }, digits = 2)
  
  output$results <- renderText({
    adapt_ibru(
      est = my_est(),
      dose = isolate(input$amt)
    )
  })
  
  output$concvstime <- renderPlot({
    shiny::req(my_est())
    my_est() %>% 
      augment(delta = .05) %>% 
      plot()
  })
  
  output$distparam <- renderPlot({
    mapbayr:::hist.mapbayests(my_est())
  })
}

shinyApp(ui = ui, server = server) 

library(shiny)
library(mapbayr)
library(mrgsolve)
# --- First: Drug specific elements, depending on your drug, model, protocol etc...

# - mrgsolve model

my_model <- mread("mrg_931.cpp")

# - a posteriori "Adaptation" function(s): return a dose recommendation, a comment, a specific figure...

adapt_cabo <- function(est, ii, target, ss){
  stopifnot(inherits(est, "mapbayests"))
  
  # 1- Return the trough concentration at steady-state
  tab <- est$mapbay_tab
  cii <- tab$IPRED[(tab$mdv==1 & tab$evid==0)][1]
  
  TXT <- paste0("Cabozantinib concentration ", ii, "h after last intake: ", round(cii, 3), "mg/L")
  
  # 2- Simulate a posteriori dosing regimens 
  inputsim <- expand.grid(iII = c(24, 48, 72), iAMT = c(20, 40, 60))
  inputsim$iADDL <- ((24 / (inputsim$iII)) * ss - 1)
  inputsim$iTSIM <- (inputsim$iADDL + 1) * (inputsim$iII)
  
  iCONC <- inputsim %>% 
    purrr::pmap_dbl(
      function(iAMT, iII, iADDL, iTSIM){
        est %>% 
          use_posterior() %>% 
          adm_lines(amt = iAMT, ii = iII, addl = iADDL, realize_addl = TRUE) %>%
          obsonly() %>% 
          mrgsim(start = iTSIM, end = iTSIM) %>% 
          as.data.frame() %>% 
          dplyr::pull("DV")
      })
  
  TAB <- data.frame(
    Schedule = paste0(inputsim$iAMT, "mg/", inputsim$iII, "h"),
    Concentration = iCONC
  ) %>% 
    dplyr::arrange(Concentration)
  
  TAB <- list(
    TXT = TXT, 
    TAB = TAB)
  
}

# --- Secondly: The shiny app per se

ui <- fluidPage(
  titlePanel("Cabozantinib - mapbayr"),
  sidebarLayout(
    sidebarPanel(
      h4("Administration"),
      numericInput("amt", "Dose (mg)", 60),
      numericInput("ii", "Interdose Interval (h)", 24),
      numericInput("addl", "Number of doses taken (n)", 30),
      fluidRow(
        h4("Observations"),
        column(width = 6, numericInput("time1", "Time after last dose (h)", 12.5)),
        column(width = 6, numericInput("dv1", "Concentration (mg/L)", .65))
      ),
      fluidRow(
        h4("Covariates"),
        column(width = 6, 
               numericInput("age", "Age (years)", 60), 
               numericInput("wt", "Body weight (kg)", 80), 
               radioButtons("sex", "Sex", list(Male = 0, Female = 1))
        ),
        column(width = 6, 
               radioButtons("cancer", "Cancer Type", c("RCC", "CRPC", "MTC", "GB", "HCC", "OTHER"))
        )
      ),
      fluidRow(
        h4("Target"),
        column(width = 6, numericInput("target", "Target conc (mg/L)", 1.1)),
        column(width = 6, numericInput("ss", "Steady-state (days)", 30))
      ),
      h4("Perfom estimation"),
      actionButton("GO", "Let's GO !")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data",  tableOutput("mapbay_tab")),
        tabPanel("Concentrations",  plotOutput("concvstime")),
        tabPanel("Parameters",  plotOutput("distparam")), 
        tabPanel("Results", 
                 verbatimTextOutput("results_txt"), 
                 h3("Individual predicted concentration at steady-state (mg/L)"),
                 tableOutput("results_tab"))
      )
    )
  )
)

server <- function(input, output) {
  
  my_data <- reactive({
    
    vec_cov <- c("RCC", "CRPC", "MTC", "GB", "HCC", "OTHER")
    list_cov <- c(list(AGE = input$age, 
                       SEX = as.integer(input$sex),
                       WT = input$wt), 
                  as.list(setNames(as.integer(vec_cov%in%input$cancer), vec_cov))
                  )
    
    tldos <- (input$ii)*(input$addl-1)
    
    my_model %>%
      adm_lines(time = 0, amt = input$amt, ii = input$ii, addl = (input$addl)-1) %>%
      obs_lines(time = (input$time1 + tldos), DV = input$dv1) %>%
      obs_lines(time = (input$ii + tldos), DV = NA_real_, mdv = 1) %>% 
      add_covariates(list_cov) %>%
      get_data()
  })
  
  my_est <- eventReactive(input$GO, {
    mapbayest(my_model, my_data(), verbose = F)
  })
  
  my_adapt <- reactive({
    adapt_cabo(
      est = my_est(), 
      ii = isolate(input$ii), 
      target = input$target, 
      ss = input$ss
    )
  })
  
  output$mapbay_tab <- renderTable({
    if(input$GO == 0){
      my_data()
    } else {
      as.data.frame(my_est())
    }
  })
  output$concvstime <- renderPlot({
    mapbayr:::plot.mapbayests(my_est())+ggplot2::geom_hline(yintercept = input$target, linetype = 2)
    #We can also use :
    # shiny::req(my_est())
    # plot(my_est())
  })
  
  output$distparam <- renderPlot({
    mapbayr:::hist.mapbayests(my_est())
  })
  
  output$results_txt <- renderText({
    my_adapt()[["TXT"]]
  })
  
  output$results_tab <- renderTable({
    my_adapt()[["TAB"]]
  }, digits = 3)
  
}

shinyApp(ui = ui, server = server) 

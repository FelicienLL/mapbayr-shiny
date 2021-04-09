library(shiny)
library(mapbayr)
library(mrgsolve)
# --- First: Drug specific elements, depending on your drug, model, protocol etc...

# - mrgsolve model

my_model <- mread("mrg_951.cpp")

# - a posteriori "Adaptation" function(s): return a dose recommendation, a comment, a specific figure...

adapt_mtx <- function(est, dur, target){
  stopifnot(inherits(est, "mbrests"))
  
  est2 <- augment(est, end = 300, delta = .1)

  ttarg <- dplyr::filter(est2$aug_tab, time > dur, type == "IPRED", value < target)[[1, "time"]]
  
  fig <- plot(est2)+
    ggplot2::geom_hline(yintercept = target, linetype = 2)+
    ggplot2::scale_y_log10(limits = c(target/10, NA))+
    ggplot2::scale_x_continuous(limits = c(NA, ttarg+10))
  
  list(TXT = paste0("Predicted time for MTX concentration under ", target, "\u00B5mol/L: T", round(ttarg), "h."), 
       FIG = fig
  )
}

# --- Secondly: The shiny app per se



ui <- fluidPage(
  titlePanel("Methotrexate - mapbayr"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        h4("Administration"),
        column(width = 6,
               numericInput("amt", "Dose (mg)", 5000)
        ),
        column(width = 6,
               numericInput("dur", "Duration (h)", 6, min = .1)
        )
      ),
      fluidRow(
        h4("Observations"),
        column(width = 6,
               numericInput("time1", "Time 1 (h)", 24),
               numericInput("time2", "Time 2 (h)", NA),
               numericInput("time3", "Time 3 (h)", NA),
               numericInput("time4", "Time 4 (h)", NA)
        ),
        column(width = 6,
               numericInput("dv1", "Conc 1 (\u00B5mol/L)", 3.2),
               numericInput("dv2", "Conc 2 (\u00B5mol/L)", NA),
               numericInput("dv3", "Conc 3 (\u00B5mol/L)", NA),
               numericInput("dv4", "Conc 4 (\u00B5mol/L)", NA)
        )
      ),
      fluidRow(
        h4("Covariates"), 
        column(width = 6, 
               numericInput("age", "Age (years)", 62)), 
        column(width = 6, 
               numericInput("scr", "Creatinine clearance (mL/min/1.73m2)", 67))
      ),
      h4("Patient/Protocol"),
      numericInput("target","Target Concentration (\u00B5mol/L)", 0.2),
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
      adm_lines(time = 0, amt = (input$amt), rate = (input$amt)/(input$dur)) %>%
      obs_lines(time = c(input$time1,input$time2,input$time3,input$time4),
                DV = c(input$dv1,input$dv2,input$dv3,input$dv4)) %>%
      add_covariates(list(AGE = input$age, SCR = input$scr)) %>%
      see_data() %>% 
      dplyr::filter(!((mdv==0)&(is.na(time)|is.na(DV))))
  })
  
  my_est <- eventReactive(input$GO, {
    mbrest(my_model, my_data(), verbose = F)
  })
  
  my_adapt <- reactive({
    adapt_mtx(
      est = my_est(),
      dur = isolate(input$dur),
      target = isolate(input$target)
    )
  })
  
  
  output$mapbay_tab <- renderTable({
    if(input$GO == 0){
      my_data()
    } else {
      as.data.frame(my_est())
    }
  })
  
  output$results <- renderText({
    my_adapt()[["TXT"]]
  })
  
  output$concvstime <- renderPlot({
    my_adapt()[["FIG"]]
  })
  
  output$distparam <- renderPlot({
    mapbayr:::hist.mbrests(my_est())
  })
  
  
}

shinyApp(ui = ui, server = server) 
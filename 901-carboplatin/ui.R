library(shiny)
library(mapbayr)
library(mrgsolve)

fluidPage(
  titlePanel("High dose carboplatin TDM - mapbayr"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        h4("Administration"),
        column(width = 6,
               numericInput("amt", "Dose (mg)", 500)
        ),
        column(width = 6,
               numericInput("dur", "Duration (h)", 1, min = .1)
        )
      ),
      fluidRow(
        h4("Observations"),
        column(width = 6,
               numericInput("time1", "Time 1 (h)", 1.05),
               numericInput("time2", "Time 2 (h)", 1.9),
               numericInput("time3", "Time 3 (h)", 4.7)
        ),
        column(width = 6,
               numericInput("dv1", "Conc 1 (mg/L)", 23.1),
               numericInput("dv2", "Conc 2 (mg/L)", 14.6),
               numericInput("dv3", "Conc 3 (mg/L)", 4.8)
        )
      ),
      fluidRow(
        h4("Patient/Protocol"),
        column(width = 6,
               numericInput("bsa","BSA (m2)", 1.73)),
        column(width = 6,
               numericInput("auc", "Target AUC (mg.h/mL)", 24))
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
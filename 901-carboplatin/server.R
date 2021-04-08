function(input, output) {
  # --- First: Drug specific elements, depending on your drug, model, protocol etc...
  
  # - mrgsolve model
  
  my_model <- mread("mrg_901.cpp")
  
  # - a posteriori "Adaptation" function(s): return a dose recommendation, a comment, a specific figure...
  
  adapt_carboplatin <- function(est, target_auc, dose_d1){
    stopifnot(inherits(est, "mbrests"))
    CL <- (1000/60) * as.data.frame(est)[1, "CL"]
    AUC_D1 <- dose_d1 / CL
    AUC_D2 <- AUC_D1
    AUC_D3 <- target_auc - (AUC_D2 + AUC_D1)
    DOSE_D3 <- CL * AUC_D3
    if(DOSE_D3<=0) DOSE_D3 <- 0
    
    paste0("Patient's clearance: ", round(CL), " mL/min.", "\n",
           "AUC Day 1: ", round(AUC_D1, 2), " mg.h/mL.", "\n",
           "Dose at day 3 to reach a target AUC of ", target_auc, " mg.h/mL over 3 days: ", round(DOSE_D3), " mg.")
  }
  
  # --- Secondly: The shiny app per se
  
  
  
  my_data <- reactive({
    my_model %>%
      adm_lines(time = 0, amt = input$amt, rate = (input$amt)/(input$dur)) %>%
      obs_lines(time = c(input$time1, input$time2, input$time3),
                DV = c(input$dv1, input$dv2, input$dv3)) %>%
      add_covariates(list(BSA = input$bsa)) %>%
      see_data()
  })
  
  my_est <- eventReactive(input$GO, {
    mbrest(my_model, my_data(), verbose = F)
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
    mapbayr:::plot.mbrests(my_est())
    #We can also use :
    # shiny::req(my_est())
    # plot(my_est())
  })
  
  output$distparam <- renderPlot({
    mapbayr:::hist.mbrests(my_est())
  })
  
  
}
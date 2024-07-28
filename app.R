library(shiny)
library(bslib)
library(shinyBS)
library(DT)
library(ggplot2)

# Define the parameter look-up dictionary
if (!file.exists("parameter_lookup.csv")) {
  pt <- data.frame(
    id = c(1:25),
    name = c("CMA_cost", "GP_cost", "ES_cost", "CMA_post_n", "CMA_post_p", "GP_post_n", "GP_post_p", "ES_post_n", "ES_post_p",
             "penalty", "CMA_yield", "GP_yield", "ES_yield", "expert_CMA_tn", "expert_CMA_tp", "expert_GP_tn", "expert_GP_tp",
             "expert_ES_tn", "expert_ES_tp","ai_CMA_tp" ,"ai_CMA_tn", "ai_GP_tn", "ai_GP_tp", "ai_ES_tn", "ai_ES_tp"),
    display_name = c("CMA test cost", "Gene-panel test cost", "Exome sequencing test cost",
                     "CMA post-test cost(positive)", "CMA post-test cost (negative)",
                     "Gene-panel post-test cost(negative)", "Gene-panel post-test cost(positive)",
                     "Exome Sequencing post-test cost(negative)", "Exome Sequencing post-tets cost(positive)",
                     "FN penalty cost","CMA diagnosis yield", "Gene-panel diagnosis yield", "Exome Sequencing diagnosis yield", 
                     "Expert TN rate (CMA)", "Expert TP rate (CMA)", "Expert TN rate (GP)", "Expert TP rate (GP)",
                     "Expert TN rate (ES)", "Expert TP rate (ES)", "AI TN rate(CMA)", "AI TP rate (CMA)", "AI TN rate (GP)", 
                     "AI TP rate (GP)","AI TN rate (ES)", "AI TP rate (ES)"
                     ),
    default_value = rep(0,25),
    description = c("The cost of CMA test", "The cost of gene-panel test", "The cost of exome sequencing test",
                    "Cost of post-test procedure for negative result of CMA", "Cost of post-test procedure for positive result of CMA", 
                    "Cost of post-test procedure for negative result of gene-panel", "Cost of post-test procedure for positive result of gene-panel",
                    "Cost of post-test procedure for negative result of exome sequencing", "Cost of post-test procedure for positive result of exome sequencing",
                    "Penalty cost for exiting the test procedure", 
                    "Diagnostic yield of CMA", "Diagnostic yield of gene-panel","Diagnostic yield of exome sequencing", 
                    "Expert's evaluation of CMA test result (true negative)", "Expert's evaluation of CMA test result (true positive)",
                    "Expert's evaluation of gene-panel test result (true negative)", "Expert's evaluation of gene-panel test result (true negative)",
                    "Expert's evaluation of exome sequencing test result (true negative)", "Expert's evaluation of exome sequencing test result (true negative)",
                    "AI's evaluation of CMA test result (true negative)", "AI's evaluation of CMA test result (true negative)",
                    "AI's evaluation of gene-panel test result (true negative)", "AI's evaluation of gene-panel test result (true negative)",
                    "AI's evaluation of exome sequencing test result (true negative)", "AI's evaluation of exome sequencing test result (true negative)"
                     ),
    group = c(rep("Cost", 10), rep("Diagnostic Yield", 3), rep("Expert Performance", 6), rep("AI Performance", 6)),
    stringsAsFactors = FALSE
  )
  write.csv2(pt, "parameter_lookup.csv", row.names = FALSE)
} else {
  pt <- read.csv2("parameter_lookup.csv", stringsAsFactors = FALSE)
}

# Define function logic. Don't include it in the server
cost_equation <- function(input){
  
  CMA_cost <- input$CMA_cost
  GP_cost <- input$GP_cost
  ES_cost <- input$ES_cost
  CMA_post_n <-input$CMA_post_n 
  CMA_post_p <-input$CMA_post_p
  GP_post_n <-input$GP_post_n
  GP_post_p <-input$GP_post_p
  ES_post_n <-input$ES_post_n
  ES_post_p <-input$ES_post_p
  penalty <- input$penalty
  
  CMA_yield <- input$CMA_yield
  GP_yield <- input$GP_yield
  ES_yield <- input$ES_yield
  
  expert_CMA_tn <- input$expert_CMA_tn
  expert_CMA_tp <- input$expert_CMA_tp
  expert_GP_tn <- input$expert_GP_tn
  expert_GP_tp <- input$expert_GP_tp
  expert_ES_tn <- input$expert_ES_tn
  expert_ES_tp <- input$expert_ES_tp
  ai_CMA_tn <- input$ai_CMA_tn
  ai_CMA_tp <- input$ai_CMA_tp
  ai_GP_tn <- input$ai_GP_tn
  ai_GP_tp <- input$ai_GP_tp
  ai_ES_tn <- input$ai_ES_tn
  ai_ES_tp <- input$ai_ES_tp
  
  expert_CMA_fn <- 1- expert_CMA_tp
  expert_CMA_fp <- 1- expert_CMA_tn
  ai_CMA_fn <- 1- ai_CMA_tp
  ai_CMA_fp <- 1- ai_CMA_tn
  expert_GP_fn <- 1- expert_GP_tp
  expert_GP_fp <- 1- expert_GP_tn
  expert_ES_fn <- 1- expert_ES_tp
  expert_ES_fp <- 1- expert_ES_tn
  ai_GP_fn <- 1- ai_GP_tp
  ai_GP_fp <- 1- ai_GP_tn
  ai_ES_fn <- 1- ai_ES_tp
  ai_ES_fp <- 1- ai_ES_tn
  
  expert_alone_cost <- 0
  delegation_cost <- 0
  
  if (input$xaxis == "Branch 1"){
    
    p_t12_gp <- GP_yield * (expert_GP_tp - expert_GP_fp) + expert_GP_fp
    p_t12_es <- 1-p_t12_gp
    p_t13_exit <- ES_yield * (expert_ES_fn - expert_ES_tn) + expert_ES_tn
    p_t13_es <- 1-p_t13_exit
    
    p_r1_gt_r1_star <- GP_yield*(ai_GP_tp-ai_GP_fp)+ai_GP_fp
    p_r1_sm_r1_star <- 1-p_r1_gt_r1_star
    p_r2_gt_r2_star <- ES_yield*(ai_ES_tp-ai_ES_fp)+ai_ES_fp
    p_r2_sm_r2_star <- 1-p_r2_gt_r2_star
    
    if (input$exhaust_or_not == "No"){
    
      expert_alone_cost <- CMA_cost + p_t12_gp*(GP_cost + GP_yield*GP_post_p + (1-GP_yield)*p_t13_exit*(GP_post_n + penalty)) +
        (p_t12_gp*(1-GP_yield)*p_t13_es + p_t12_es)*(ES_cost + ES_yield*ES_post_p + (1-ES_yield)*ES_post_n)
      
    } else{
      
      expert_alone_cost <- CMA_cost + p_t12_gp*(GP_cost + GP_yield*GP_post_p) +
        (p_t12_gp*(1-GP_yield) + p_t12_es)*(ES_cost + ES_yield*ES_post_p + (1-ES_yield)*ES_post_n)
    }
    
    delegation_cost <- CMA_cost + p_r1_gt_r1_star*(GP_cost + GP_yield*GP_post_p + (1-GP_yield)*GP_post_n) +
    (p_r1_gt_r1_star*(1-GP_yield) + p_r1_sm_r1_star)*(p_r2_gt_r2_star*(ES_cost + ES_yield*ES_post_p + 
    (1-ES_yield)*ES_post_n) + p_r2_sm_r2_star*expert_alone_cost)

  }
  
  else if (input$xaxis == "Branch 2"){ #Only one decision node
 
    p_t22_exit <- CMA_yield * (expert_CMA_fn - expert_CMA_tn) + expert_CMA_tn
    p_t22_cma <- 1 - p_t22_exit
    p_r0_gt_r0_star <- CMA_yield*(ai_CMA_tp - ai_CMA_fp) + ai_CMA_fp
    p_r0_sm_r0_star <- 1-p_r0_gt_r0_star
    
    if (input$exhaust_or_not == "No"){
     
      expert_alone_cost <- ES_cost + ES_yield*ES_post_p + (1-ES_yield)*p_t22_exit*(penalty+ES_post_n) + 
                          (1-ES_yield)*p_t22_cma*(CMA_yield*CMA_post_p + (1-CMA_yield)*CMA_post_n)
    } else {
      
      expert_alone_cost <- ES_cost + ES_yield*ES_post_p + (1-ES_yield)*(CMA_yield*CMA_post_p + (1-CMA_yield)*CMA_post_n)
        
    }
    
    delegation_cost <- ES_cost + ES_yield*ES_post_p + (1-ES_yield)*(p_r0_sm_r0_star)*expert_alone_cost +
                       (1-ES_yield)*p_r0_gt_r0_star*(CMA_cost + CMA_yield*CMA_post_p + (1-CMA_yield)*CMA_post_n)
  }
  
  else{
    expert_alone_cost <- ES_cost + CMA_cost
    delegation_cost <- ES_cost + CMA_cost
  }
  
  return (c(expert_alone_cost, delegation_cost))
}

# Output a Plot
plotBarPlot <-function(input){
  # Create a plot
  cost_data <- data.frame(
    Mode = c("Expert-alone mode", "Delegation mode"),
    Cost = cost_equation(input)
  )
  
  #plot(1:10, 1:10, type = "n", xlab = input$xaxis, ylab = input$yaxis)
  
  ggplot(cost_data, aes(x = Mode, y = Cost, fill = Mode)) +
    geom_bar(stat = "identity") + geom_text(aes(label = round(Cost, 2)), vjust = -0.5, size = 4) +
    labs(x = input$xaxis, y = input$yaxis) +
    theme_minimal() + theme(legend.position = "top")
}

# Output a Table
generateParameterTable <-function(input){
  working_pt = pt 
  working_pt$default_value = NULL
  value = c(input$CMA_cost, input$GP_cost, input$ES_cost, input$CMA_post_n, input$CMA_post_p, 
            input$GP_post_n, input$GP_post_p, input$ES_post_n, input$ES_post_p, input$penalty, 
            input$CMA_yield, input$GP_yield, input$ES_yield, input$expert_CMA_tn, input$expert_CMA_tp,  
            input$expert_GP_tn, input$expert_GP_tp, input$expert_ES_tn, input$expert_ES_tp, 
            input$ai_CMA_tn, input$ai_CMA_tp, input$ai_GP_tn, input$ai_GP_tp, input$ai_ES_tn, input$ai_ES_tp)
  working_pt$value = value
  working_pt
}


ui <- fluidPage(
  titlePanel("Cost Function Visualization"),
  # Add some customized CSS styles
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(12,
               div(
                 #h4("X and Y axis"),
                 # Drop Down Menu for selecting y-axis. only one option can be selected for now
                 selectInput("yaxis", "Y-axis",
                             choices = c("Expected Cost"),
                             selected = "Expected Cost"),
                 selectInput("xaxis", "X-axis",
                             choices = c("Branch 1", "Branch 2", "Branch 3"),
                             selected = "Branch 1"),
                 bsPopover(id = "xaxis", title = NULL, content = "Branch 1 starts with CMA; Branch 2 starts with ES; Branch 3 conducts ES and CMA concurrently", 
                            placement = "bottom", trigger = "hover"),
                 selectInput("exhaust_or_not", "Exhaust all possible tests or not?",
                             choices = c("Yes", "No"),
                             selected = "No"),
                 bsPopover(id = "exhaust_or_not", title = NULL, content = "If No is selected, then the patient will have the option 
                           to exit the diagnosis test procedure; otherwise all possible tests will be conducted", 
                           placement = "bottom", trigger = "hover")
               )
        )
      ),
      fluidRow(
        column(12,
               div(
                 selectInput("parameter_group", "Change Detailed Parameters",
                             choices = c("Cost", "Diagnostic Yield", 'Expert Performance', "AI Performance"),
                             selected = "Cost"),
               )
        )
      ),
      actionButton("reset", "Reset", class = "btn-warning"),
      
      conditionalPanel(
        condition = "input.parameter_group == 'Cost'",
        fluidRow(
          column(12,
                 div(
                   h4("Cost Parameters"),
                   numericInput("CMA_cost", pt$display_name[pt$name == "CMA_cost"], value = pt$default_value[pt$name == "CMA_cost"], min = 0),
                   numericInput("GP_cost", pt$display_name[pt$name == "GP_cost"], value = pt$default_value[pt$name == "GP_cost"], min = 0),
                   numericInput("ES_cost", pt$display_name[pt$name == "ES_cost"], value = pt$default_value[pt$name == "ES_cost"], min = 0),
                   numericInput("CMA_post_n", pt$display_name[pt$name == "CMA_post_n"], value = pt$default_value[pt$name == "CMA_post_n"], min = 0),
                   numericInput("CMA_post_p", pt$display_name[pt$name == "CMA_post_p"], value = pt$default_value[pt$name == "CMA_post_p"], min = 0),
                   numericInput("GP_post_n", pt$display_name[pt$name == "GP_post_n"], value = pt$default_value[pt$name == "GP_post_n"], min = 0),
                   numericInput("GP_post_p", pt$display_name[pt$name == "GP_post_p"], value = pt$default_value[pt$name == "GP_post_p"], min = 0),
                   numericInput("ES_post_n", pt$display_name[pt$name == "ES_post_n"], value = pt$default_value[pt$name == "ES_post_n"], min = 0),
                   numericInput("ES_post_p", pt$display_name[pt$name == "ES_post_p"], value = pt$default_value[pt$name == "ES_post_p"], min = 0),
                   numericInput("penalty", pt$display_name[pt$name == "penalty"], value = pt$default_value[pt$name == "penalty"], min = 0)
                 )
          )
        )
      ),
      
      conditionalPanel(
          
        condition = "input.parameter_group == 'Diagnostic Yield'",
        fluidRow(
          column(12,
                 div(
                   h4("Diagnostic Yields"),
                   numericInput("CMA_yield", pt$display_name[pt$name == "CMA_yield"], value = pt$default_value[pt$name == "CMA_yield"], min = 0),
                   numericInput("GP_yield", pt$display_name[pt$name == "GP_yield"], value = pt$default_value[pt$name == "GP_yield"], min = 0),
                   numericInput("ES_yield", pt$display_name[pt$name == "ES_yield"], value = pt$default_value[pt$name == "ES_yield"], min = 0)
                 )
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.parameter_group == 'Expert Performance'",
        fluidRow(
          column(12,
                 div(
                   h4("Expert Performance"),
                   numericInput("expert_CMA_tn", pt$display_name[pt$name == "expert_CMA_tn"], value = pt$default_value[pt$name == "expert_CMA_tn"], min = 0),
                   numericInput("expert_CMA_tp", pt$display_name[pt$name == "expert_CMA_tp"], value = pt$default_value[pt$name == "expert_CMA_tp"], min = 0),
                   numericInput("expert_GP_tn", pt$display_name[pt$name == "expert_GP_tn"], value = pt$default_value[pt$name == "expert_GP_tn"], min = 0),
                   numericInput("expert_GP_tp", pt$display_name[pt$name == "expert_GP_tp"], value = pt$default_value[pt$name == "expert_GP_tp"], min = 0),
                   numericInput("expert_ES_tn", pt$display_name[pt$name == "expert_ES_tn"], value = pt$default_value[pt$name == "expert_ES_tn"], min = 0),
                   numericInput("expert_ES_tp", pt$display_name[pt$name == "expert_ES_tp"], value = pt$default_value[pt$name == "expert_ES_tp"], min = 0),
                  )
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.parameter_group == 'AI Performance'",
        fluidRow(
          column(12,
                 div(
                   h4("AI Performance"),
                   numericInput("ai_CMA_tn", pt$display_name[pt$name == "ai_CMA_tn"], value = pt$default_value[pt$name == "ai_CMA_tn"], min = 0),
                   numericInput("ai_CMA_tp", pt$display_name[pt$name == "ai_CMA_tp"], value = pt$default_value[pt$name == "ai_CMA_tp"], min = 0),
                   numericInput("ai_GP_tn", pt$display_name[pt$name == "ai_GP_tn"], value = pt$default_value[pt$name == "ai_GP_tn"], min = 0),
                   numericInput("ai_GP_tp", pt$display_name[pt$name == "ai_GP_tp"], value = pt$default_value[pt$name == "ai_GP_tp"], min = 0),
                   numericInput("ai_ES_tn", pt$display_name[pt$name == "ai_ES_tn"], value = pt$default_value[pt$name == "ai_ES_tn"], min = 0),
                   numericInput("ai_ES_tp", pt$display_name[pt$name == "ai_ES_tp"], value = pt$default_value[pt$name == "ai_ES_tp"], min = 0),
                 )
          )
        )
      ),
    ),
    mainPanel(
      h3("Introduction"),
      conditionalPanel(
        condition = "input.xaxis == 'Branch 1'",
        img(src = "branch1.png", width = "60%", height = "60%")
      ),
      conditionalPanel(
        condition = "input.xaxis == 'Branch 2'",
        img(src = "branch2.png", width = "60%", height = "60%")
      ),
      conditionalPanel(
        condition = "input.xaxis == 'Branch 3'",
        img(src = "branch3.png", width = "60%", height = "60%")
      ),
      plotOutput("expected_cost_plot"),
      DTOutput("Parameters")
    )
  )
)

server <- function(input, output, session) {
  # Reactive value to store the updated parameter lookup
  param_data <- reactiveVal(pt)
  
  # Reset all parameters to their default values when the reset button is clicked
  observeEvent(input$reset, {
    param_data(pt)
    # Reset input fields to their default values
    for (i in 1:nrow(pt)) {
      updateNumericInput(session, pt$name[i], value = pt$default_value[i], min=0)
    }
  })
  
  # Define a reactive expression for the plot
  output$expected_cost_plot <- renderPlot({
    plotBarPlot(input)
  })
  output$Parameters <- renderDT({
    generateParameterTable(input)
  })
}

shinyApp(ui = ui, server = server)

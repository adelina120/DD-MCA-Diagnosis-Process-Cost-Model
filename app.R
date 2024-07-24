library(shiny)
library(bslib)
library(shinyBS)
library(DT)

# Define the parameter lookup dictionary
if (!file.exists("parameter_lookup.csv")) {
  pt <- data.frame(
    id = c(1:19),
    name = c("CMA_cost", "GP_cost", "ES_cost", "GP_post_n", "GP_post_p", "ES_post_n", "ES_post_p",
             "penalty", "CMA_yield", "GP_yield", "ES_yield", "expert_GP_tn", "expert_GP_tp",
             "expert_ES_tn", "expert_ES_tp", "ai_GP_tn", "ai_GP_tp", "ai_ES_tn", "ai_ES_tp"),
    display_name = c("CMA test cost", "Gene-panel test cost", "Exome sequencing test cost",
                     "Gene-panel post-test cost(negative)", "Gene-panel(positive)",
                     "Exome Sequencing(negative)", "Exome Sequencing(positive)", "FN penalty","CMA yield",
                     "Gene-panel yield", "Exome Sequencing yield", "Expert TN rate (GP)", "Expert TP rate (GP)",
                     "Expert TN rate (ES)", "Expert TP rate (ES)", "AI TN rate (GP)", "AI TP rate (GP)",
                     "AI TN rate (ES)", "AI TP rate (ES)"),
    default_value = c(500, 1500, 3000, 500, 1000, 1500, 1500, 2000, 0.05, 0.2, 0.3, 0.6, 0.65, 0.7, 0.75, 0.6, 0.65, 0.7, 0.75),
    description = c("Description of param1", "Description of param2", "Description of param3",
                    "Description of param4", "Description of param5", "Description of param6",
                    "Description of param7", "Description of param8", "Description of param9",
                    "Description of param10", "Description of param11", "Description of param12",
                    "Description of param13", "Description of param14", "Description of param15",
                    "Description of param16", "Description of param17", "Description of param18",
                    "Description of param19"),
    
    group = c(rep("Cost", 8), rep("Diagnostic Yield", 3), rep("Expert Performance", 4), rep("AI Performance", 4)),
    stringsAsFactors = FALSE
  )
  write.csv2(pt, "parameter_lookup.csv")
} else {
  pt <- read.csv2("parameter_lookup.csv", stringsAsFactors = FALSE)
}

# Define function logic. Don't include it in the server
cost <- function(CMA, GP, ES, GP_post_n, GP_post_p, ES_post_n, ES_post_p, penalty,
                 GP_yield, ES_yield, GP_tn, ES_tn, ES_tp){
  
  GP_fn <- 1-GP_tp
  GP_fp <- 1-GP_tn
  ES_fn <- 1-ES_tp
  ES_fp <- 1-ES_tn
  
  p_t12_gp <- GP_yield * (GP_tp - GP_fp) + GP_fp
  p_t12_es <- 1-p_t12_gp
  p_t13_exit <- ES_yield * (ES_fn - ES_tn) + ES_tn
  p_t13_es <- 1-p_t13_exit
  
  expected_cost <- CMA + p_t12_gp*(GP + GP_yield*GP_post_p + (1-GP_yield)*p_t13_exit*(GP_post_n + penalty)) +
    (p_t12_es*(1-GP_yield)*p_t13_es + p_t12_es)*(ES + ES_yield*ES_post_p + (1-ES_yield)*ES_post_n)
}

plotBarPlot <-function(input){
  # Create a plot
  # implement business logic here
  plot(1:10, 1:10, type = "n", xlab = input$xaxis, ylab = input$yaxis)
}

generateParameterTable <-function(input){
  working_pt = pt 
  working_pt$default_value = NULL
  value = c(input$CMA_cost, input$GP_cost, input$ES_cost, input$GP_post_n, input$GP_post_p, input$ES_post_n, input$ES_post_p,
            input$penalty, input$CMA_yield, input$GP_yield, input$ES_yield, input$expert_GP_tn, input$expert_GP_tp,
            input$expert_ES_tn, input$expert_ES_tp, input$ai_GP_tn, input$ai_GP_tp, input$ai_ES_tn, input$ai_ES_tp)
  working_pt$value = value
  working_pt
}


ui <- fluidPage(
  titlePanel("Cost Function Visualization"),
  # Add some custermized CSS styles
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(12,
               div(
                 h4("X and Y axis"),
                 # Drop Down Menu for selecting y-axis. only one option can be selected for now
                 selectInput("yaxis", "Y-axis",
                             choices = c("Expected Cost"),
                             selected = "Expected Cost"),
                 selectInput("xaxis", "X-axis",
                             choices = c("Branch 1", "Branch 2", "Branch 3"),
                             selected = "Branch 1")
                 # add some annotation. The annotation will be shown when the mouse hover over the question mark
                 # bsPopover(id = "xaxis", title = NULL, content = "Branch 1: start with ES; Branch 2: start with ES", placement = "bottom", trigger = "hover")
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
      # Add the branch figure from a file
      img(src = "placeholder.png", align = "center"),
      plotOutput("expected_cost"),
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
  output$expected_cost <- renderPlot({
    plotBarPlot(input)
  })
  output$Parameters <- renderDT({
    generateParameterTable(input)
  })
}

shinyApp(ui = ui, server = server)

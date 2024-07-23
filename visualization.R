library(shiny)
library(bslib)

ui <- fluidPage(
  titlePanel("Cost Function Visualization"),
  sidebarLayout(
    sidebarPanel(
      numericInput("CMA", "CMA test cost", value = 1000),
      numericInput("GP", "Gene-panel test cost", value = 1500, min = 0),
      numericInput("ES", "Exome sequencing test cost", value = 3000, min = 0),
      numericInput("GP_post_n", "Gene-panel post-test cost(negative)",value = 500, min = 0),
      numericInput("GP_post_p", "Gene-panel(positive)", value = 1000, min = 0),
      numericInput("ES_post_n", "Exome Sequencing(negative)", value = 1500, min = 0),
      numericInput("ES_post_p", "Exome Sequencing(positive)", value = 1500, min = 0),
      numericInput("penalty", "FN penalty", value = 2000, min = 0),
      numericInput("GP_yield", "Gene-panel", value = 0.3, min = 0),
      numericInput("ES_yield", "Exome Sequencing", value = 0.45, min = 0),
      numericInput("GP_tn", "Expert TN rate (GP)", value = 0.6, min = 0),
      numericInput("GP_tp", "Expert TP rate (GP)", value = 0.65, min = 0),
      numericInput("ES_tn", "Expert TN rate (ES)", value = 0.7, min = 0),
      numericInput("ES_tp", "Expert TP rate (ES)", value = 0.75, min = 0),
      actionButton("update", "Update")
    ),
    mainPanel(
      plotOutput("Expected Cost"),
      tableOutput("Parameters")
    )
  )
)

server <- function(input, output) {
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
}

shinyApp(ui = ui, server = server)

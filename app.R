source("00_init.R")
source("01_simulation_functions.R")
source("01_missingness_functions.R")
source("01_profile_plot_functions.R")
#reticulate::source_python("indices_output.py")
reticulate::source_python("indices.py")


ui <- fluidPage(
  titlePanel("Bandit Algorithms with Missing Outcomes"),
  sidebarLayout(
    sidebarPanel(
    
           selectInput("method_name", "Choose your bandit algorithm", Method),
           selectInput("beta_name", "Choose your outcome model", beta_names),
           selectInput("alpha_name", "Choose your missing data mechanism", alpha_names),
           selectInput("missing_name", "Choose your missing data method", missing_names),
           sliderInput("sample_size", label = "Sample Size", value=20, min = 1, max = 100), 
           
  ),
  
   mainPanel(
     tableOutput("summary"),
     plotOutput("plot"),
    )
)
)

server <- function(input, output, session){
  
  

  output$summary <- renderTable({
    beta <- beta_list[input$beta_name][[1]]
    alpha <- alpha_list[input$alpha_name][[1]]
    
      setting_summary(alpha=alpha, beta=beta)
  
   
    }, rownames=TRUE)
  
  
  output$plot <- renderPlot({
    beta <- beta_list[input$beta_name][[1]]
    alpha <- alpha_list[input$alpha_name][[1]]
    
    create_plot(method_name=input$method_name, alpha=alpha, beta=beta, missing_approach=input$missing_name, Nt=input$sample_size)
  })
}

shinyApp(ui, server)

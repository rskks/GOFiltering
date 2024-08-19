library(shiny)
library(ggplot2)
library(DT)
library(shinythemes)  # For themes
library(scales)

# Sample data for illustration (replace with your actual data)

rnaseq_data <- data.frame(
  Gene = "Gene1",
  Condition = "Control",
  Expression = rnorm(1)
)

lipidomics_data <- data.frame(
  Lipid = "Lipid1",
  Condition = "Control",
  Abundance = rnorm(1)
)

#Define UI for application
ui <- fluidPage(
  theme = shinytheme("cerulean"),  # Add a theme for better styling
  titlePanel("Nanoparticle Dataset Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3, # Make the sidebar take up 1/3 of the width
      selectInput("dataset", "Choose a Dataset:",
                  choices = c("Proteomics", "RNA-seq", "Lipidomics")),
      uiOutput("dynamicUI"),
      style = "padding: 20px;"  # Add padding for better layout
    ),
    mainPanel(
      fluidRow(
        column(12, plotOutput("plot", height = "400px"))
      ),
      fluidRow(
        column(12, DT::DTOutput("datatable"))
      ),
      style = "padding: 20px;"  # Add padding for better layout
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  output$dynamicUI <- renderUI({
    switch(input$dataset,
           "Proteomics" = selectizeInput("protein", "Select Protein:", choices = rownames(protein_data)),
           "RNA-seq" = selectizeInput("gene", "Select Gene:", choices = unique(rnaseq_data$Gene)),
           "Lipidomics" = selectizeInput("lipid", "Select Lipid:", choices = unique(lipidomics_data$Lipid))
    )
  })
  
  output$plot <- renderPlot({
    data <- switch(input$dataset,
                   "Proteomics" = protein_data,
                   "RNA-seq" = rnaseq_data,
                   "Lipidomics" = lipidomics_data)
    
    # Ensure the dataset is filtered correctly
    if (input$dataset == "Proteomics") {
      if (input$protein %in% rownames(data)) {
        data_to_plot <- as.data.frame(t(data[rownames(data) == input$protein, ]))
        colnames(data_to_plot) <- "Expression"
        data_to_plot$Condition <- rownames(data_to_plot)
        
        ggplot(data_to_plot, aes(x = Condition, y = Expression)) +
          geom_point(color = "steelblue", size = 4) +
          theme_minimal(base_size = 15) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95"),
            title = element_text(size = 18, face = "bold")
          ) +
          labs(title = paste("Expression of", input$protein)) +
          scale_y_continuous(labels = scales::number_format(accuracy = 0.1))  # Ensure correct number formatting
        
      } else {
        ggplot() + 
          labs(title = "No data available for selected protein")
      }
      
    } else if (input$dataset == "RNA-seq") {
      plot_data <- data[data$Gene == input$gene, ]
      if (nrow(plot_data) > 0) {
        ggplot(plot_data, aes(x = Condition, y = Expression)) +
          geom_point(color = "lightcoral", size = 4) +
          theme_minimal(base_size = 15) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95"),
            title = element_text(size = 18, face = "bold")
          ) +
          labs(title = paste("Expression of", input$gene)) +
          scale_y_continuous(labels = scales::number_format(accuracy = 0.1))  # Ensure correct number formatting
      } else {
        ggplot() + 
          labs(title = "No data available for selected gene")
      }
      
    } else {
      plot_data <- data[data$Lipid == input$lipid, ]
      if (nrow(plot_data) > 0) {
        ggplot(plot_data, aes(x = Condition, y = Abundance)) +
          geom_point(color = "mediumseagreen", size = 4) +
          theme_minimal(base_size = 15) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95"),
            title = element_text(size = 18, face = "bold")
          ) +
          labs(title = paste("Abundance of", input$lipid)) +
          scale_y_continuous(labels = scales::number_format(accuracy = 0.1))  # Ensure correct number formatting
      } else {
        ggplot() + 
          labs(title = "No data available for selected lipid")
      }
    }
  })
  
  output$datatable <- DT::renderDT({
    data <- switch(input$dataset,
                   "Proteomics" = protein_data,
                   "RNA-seq" = rnaseq_data,
                   "Lipidomics" = lipidomics_data)
    DT::datatable(data, options = list(
      pageLength = 10,
      lengthMenu = c(5, 10, 15, 20),
      initComplete = JS(
        "function(settings, json) {",
        "  $(this.api().table().header()).css({'background-color': '#f5f5f5', 'color': '#333'});",
        "  $(this.api().table().body()).css({'background-color': '#fff'});",
        "}"
      )
    ), style = "bootstrap4")  # Use Bootstrap 4 styling for better appearance
  }, server = FALSE)  # Ensure client-side processing for better performance
}

# Run the application 
shinyApp(ui = ui, server = server)
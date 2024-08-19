library(shiny)
library(ggplot2)
library(DT)
library(shinythemes)  # For themes
library(scales)

# Define UI for application
ui <- fluidPage(
  theme = shinytheme("cerulean"),  # Add a theme for better styling
  titlePanel("Nanoparticle Dataset Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3, # Make the sidebar take up 1/3 of the width
      selectInput("dataset", "Choose a Dataset:",
                  choices = c("Proteomics", "RNA-seq", "Lipidomics")),
      uiOutput("dynamicUI"),
      selectInput("plotType", "Select Plot Type:",
                  choices = c("Individual", "Grouped")),
      style = "padding: 20px;"  # Add padding for better layout
    ),
    mainPanel(
      fluidRow(
        column(9, plotOutput("plot", height = "400px"))  # Adjusted to fit the desired layout
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
    
    plot_type <- input$plotType
    
    if (input$dataset == "Proteomics") {
      if (plot_type == "Individual") {
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
            scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
          
        } else {
          ggplot() + 
            labs(title = "No data available for selected protein")
        }
        
      } # Check the plot type and process accordingly
      # Check the plot type and process accordingly
      # Check the plot type and process accordingly
      else if (plot_type == "Grouped") {
        group1 <- c("3DEV", "NCBEV")
        group2 <- c("2DSuper", "3DSuper", "FPLCSuper", "NCBSuper")
        group3 <- c("2DExomere", "3DExomere", "FPLCExomere", "NCBExomere")
        group4 <- c("2DEVP", "3DEVP")
        
        # Define the group assignments
        groups <- list(EVs = group1, Super = group2, Exomere = group3, EVp = group4)
        
        # Subset data to keep only the selected protein
        protein_data_subset <- data[input$protein, , drop = FALSE]
        
        # Initialize vectors to hold data for plotting
        group_names <- c()
        expression_values <- c()
        
        # Collect data for each group
        for (group_name in names(groups)) {
          group_cols <- groups[[group_name]]
          
          # Ensure the columns exist in the subset
          valid_cols <- intersect(group_cols, colnames(protein_data_subset))
          
          if (length(valid_cols) > 0) {
            # Append group names and expression values
            group_names <- c(group_names, rep(group_name, length(valid_cols)))
            expression_values <- c(expression_values, as.vector(t(protein_data_subset[, valid_cols])))
          }
        }
        
        # Create a data frame for plotting
        grouped_data <- data.frame(
          Group = factor(group_names, levels = c("EVs", "Super", "Exomere", "EVp")),
          Expression = expression_values
        )
        
        # Print grouped_data for debugging
        print(head(grouped_data))
        print(str(grouped_data))
        
        # Plot
        ggplot(grouped_data, aes(x = Group, y = Expression, fill = Group)) +
          geom_boxplot() +
          theme_minimal(base_size = 15) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95"),
            title = element_text(size = 18, face = "bold")
          ) +
          labs(title = paste("Grouped Protein Expression for", input$protein)) +
          scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
          scale_fill_manual(values = c("EVs" = "steelblue", "Super" = "lightcoral", "Exomere" = "mediumseagreen", "EVp" = "goldenrod"))
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
          scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
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
          scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
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
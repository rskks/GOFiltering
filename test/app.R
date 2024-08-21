library(shiny)
library(ggplot2)
library(DT)
library(shinythemes)  # For themes
library(scales)
library(readxl)
library(reactable)

# Define UI for application
ui <- fluidPage(
  theme = shinytheme("cerulean"),  # Add a theme for better styling
  titlePanel("EV and NVEP Dataset Explorer"),
  sidebarLayout(
    fluid = TRUE,
    sidebarPanel(
      width = 3, # Make the sidebar take up 1/3 of the width
      fluid = TRUE,
      
      selectInput("dataset", "Choose a Dataset:",
                  choices = c("Proteomics", "RNA-seq", "Lipidomics")),
      uiOutput("dynamicUI"),
      selectInput("plotType", "Select Plot Type:",
                  choices = c("Individual", "Grouped")),
      style = "padding: 20px;"  # Add padding for better layout
    ),
    mainPanel(
      fluidRow(
        column(12, plotOutput("plot", height = "400px"))  # Adjusted to fit the desired layout
      ),
      fluidRow(
        column(12, reactableOutput("datatable"))
      ),
      style = "padding: 20px;"  # Add padding for better layout
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Load the data
  protein_data <- reactive({
    file_path <- "protein_data_avgN.xlsx"
    cat("Loading data from:", file_path, "\n")  # Debugging output
    if (file.exists(file_path)) {
      cat("File exists. Reading the data...\n")
      data <- read_xlsx(file_path)
      cat("Data loaded. Converting to data frame...\n")
      data <- as.data.frame(data, check.names = FALSE)
      rownames(data) <- data$protein
      data <- data[,-1]  # Remove the $protein column
      data <- round(data, 2)  # Round data to 2 decimal places
      cat("Data processing complete. Returning data...\n")
      return(data)
    } else {
      stop("File not found: ", file_path)
    }
  })
  
  output$dynamicUI <- renderUI({
    switch(input$dataset,
           "Proteomics" = selectizeInput("protein", "Select Protein:", 
                                         choices = rownames(protein_data())),
           "RNA-seq" = selectizeInput("gene", "Select Gene:", 
                                      choices = unique(rnaseq_data$Gene)),
           "Lipidomics" = selectizeInput("lipid", "Select Lipid:", 
                                         choices = unique(lipidomics_data$Lipid))
    )
  })
  
  output$plot <- renderPlot({
    data <- switch(input$dataset,
                   "Proteomics" = protein_data(),
                   "RNA-seq" = rnaseq_data,
                   "Lipidomics" = lipidomics_data)
    
    plot_type <- input$plotType
    
    if (input$dataset == "Proteomics") {
      cat("Selected dataset is Proteomics.\n")
      cat("Checking if protein is selected...\n")
      if (!is.null(input$protein) && input$protein %in% rownames(data)) {
        cat("Protein is selected and valid.\n")
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
    } else {
      cat("Protein is NULL or not in the dataset.\n")
      ggplot() + labs(title = "No data available for selected protein")
    }
  }
  )
  
  data <- reactive({
    switch(input$dataset,
           "Proteomics" = protein_data(),
           "RNA-seq" = rnaseq_data,
           "Lipidomics" = lipidomics_data)
  })
  
  # Render the reactable table
  output$datatable <- renderReactable({
    df <- data()
    
    if (is.data.frame(df)) {
      cat("Data is a data frame.\n")
    } else {
      cat("Data is not a data frame!\n")
      str(df)  # Print the structure of the data to see what's wrong
    }
    
    reactable(
      df,
      searchable = TRUE,
      pagination = TRUE,  # Enables pagination
      defaultPageSize = 10,  # Default page length
      pageSizeOptions = c(5, 10, 20),  # Length menu options
      theme = reactable::reactableTheme(
        headerStyle = list(
          backgroundColor = '#f5f5f5',
          color = '#333'
        ),
        cellStyle = list(
          backgroundColor = '#fff'
        )
      )
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
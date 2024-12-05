library(shiny)
library(ggplot2)
library(DT)
library(shinythemes)  # For themes
library(scales)
library(reactable)
library(yarrr)
library(ggpirate)

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
      selectInput("grouping", "Grouping", choices = c("Individual", "Grouped")),
      checkboxInput("facet_isolation", "Isolation Method", value = FALSE),
      checkboxInput("facet_growth", "Growth Conditions", value = FALSE),
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
    file_path <- "protein_data.csv"
    cat("Loading data from:", file_path, "\n")  # Debugging output
    if (file.exists(file_path)) {
      cat("File exists. Reading the data...\n")
      data <- read.csv(file_path)
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
  protein_meta <- reactive({
    file_path <- "protein_meta2_fixed.csv"
    if (file.exists(file_path)) {
      data <- read.csv(file_path, row.names = 1)
      #data <- t(data)
      data <- as.data.frame(data, check.names = FALSE)
      return(data)
    } else {
      stop("File not found: ", file_path)
    }
  })
  
  # RNA data (unchanged, as it works perfectly)
  rna_data <- reactive({
    file_path <- "mirna_data_avgN.csv"
    cat("Loading data from:", file_path, "\n")  # Debugging output
    if (file.exists(file_path)) {
      cat("File exists. Reading the data...\n")
      data <- as.data.frame(read.csv(file_path), check.names = FALSE)
      rownames(data) <- data$rna
      data <- data[,-1]  # Remove the $rna column
      cat("Data processing complete. Returning data...\n")
      return(data)
    } else {
      stop("File not found: ", file_path)
    }
  })
  
  # Lipid data (updated for consistency)
  lipid_data <- reactive({
    req(input$lipidomics_type) # Ensure input exists
    
    file_path <- switch(input$lipidomics_type,
                        "Category" = "lipid_data_cat.csv",
                        "Class" = "lipid_data_cl.csv",
                        "Subclass" = "lipid_data_subcl.csv",
                        "Species" = "lipid_data_spe.csv")
    
    if (file.exists(file_path)) {
      data <- as.data.frame(read.csv(file_path, stringsAsFactors = FALSE))
      
      # Debugging output
      cat("Loaded data from:", file_path, "\n")
      cat("Column names:\n", colnames(data), "\n")
      cat("First few rows:\n")
      print(head(data))
      
      if ("lipid" %in% colnames(data)) {
        rownames(data) <- data$lipid
        data <- data[,-1]  # Remove the 'lipid' column
        cat("Row names set successfully.\n")
      } else {
        stop("Column 'lipid' not found in the data.")
      }
      
      return(data)
    } else {
      showNotification("Selected lipidomics file not found.", type = "error")
      # Return an empty data frame to avoid issues downstream
      return(data.frame())
    }
  })
  

  # Dynamic UI rendering for Proteomics, RNA-seq, or Lipidomics (with sub-types for Lipidomics)
  output$dynamicUI <- renderUI({
    if (input$dataset == "Lipidomics") {
      # Lipidomics-specific UI
      tagList(
        selectInput("lipidomics_type", "Choose Lipidomics Type:",
                    choices = c("Category", "Class", "Subclass", "Species"),
                    selected = "Category"),
        selectizeInput("lipid", "Select Lipid:", choices = NULL) # Initially empty
      )
    } else if (input$dataset == "Proteomics") {
      selectizeInput("protein", "Select Protein:", choices = NULL) # Placeholder
    } else if (input$dataset == "RNA-seq") {
      selectizeInput("rna", "Select miRNA:", choices = NULL) # Placeholder
    }
  })
  
  # Update lipidomics dropdown when lipidomics_type changes
  observeEvent(input$lipidomics_type, {
    req(input$lipidomics_type) # Ensure lipidomics type is selected
    data <- lipid_data()
    req(data) # Ensure lipid_data is loaded
    
    updateSelectizeInput(session, "lipid", choices = rownames(data), server = TRUE)
  })
  
  # Update proteomics dropdown when protein_data changes
  observe({
    data <- protein_data()
    req(data) # Ensure protein_data is loaded
    
    updateSelectizeInput(session, "protein", choices = rownames(data), server = TRUE)
  })
  
  # Update RNA-seq dropdown when rna_data changes
  observe({
    data <- rna_data()
    req(data) # Ensure rna_data is loaded
    
    updateSelectizeInput(session, "rna", choices = rownames(data), server = TRUE)
  })
  
  # Observe event to update inputs on dataset switching
  observeEvent(input$dataset, {
    # Reset the dropdown menu for proteins
    if (input$dataset == "Proteomics") {
      updateSelectInput(session, "protein", 
                        choices = rownames(protein_data()),
                        selected = NULL)
    }
    
    # Reset the dropdown menu for RNA
    if (input$dataset == "RNA-seq") {
      updateSelectInput(session, "rna", 
                        choices = rownames(rna_data()),
                        selected = NULL)
    }
    
    # Reset the dropdown menu for lipids
    if (input$dataset == "Lipidomics") {
      lipid_choices <- rownames(lipid_data())
      updateSelectInput(session, "lipid", 
                        choices = lipid_choices,
                        selected = NULL)
    }
  })
  
  
  output$plot <- renderPlot({
    data <- switch(input$dataset,
                   "Proteomics" = protein_data(),
                   "RNA-seq" = rna_data(),
                   "Lipidomics" = lipid_data())
    
    if (input$dataset == "Proteomics") {
      if (!is.null(input$protein) && input$protein %in% rownames(data)) {
        selected_protein_data <- data[input$protein, , drop = FALSE]
        grouping <- input$grouping
        plot_type <- ifelse(grouping == "Individual", "dotplot", "pirateplot")
        
        # Get the meta information for faceting
        meta_data <- protein_meta()
        
        individual_data <- as.vector(unlist(meta_data["Individual", colnames(selected_protein_data)]))
        partile_data <- as.vector(unlist(meta_data["Particle", colnames(selected_protein_data)]))
        isolation_data <- as.vector(unlist(meta_data["Isolation", colnames(selected_protein_data)]))
        growth_data <- as.vector(unlist(meta_data["Growth", colnames(selected_protein_data)]))
        
        total_data <- data.frame(
          Condition = colnames(selected_protein_data),
          Value = as.numeric(selected_protein_data),
          Individual = factor(individual_data),
          Particle = factor(partile_data),
          Isolation = factor(isolation_data),
          Growth = factor(growth_data)
        )
        
        print(head(total_data))
        
        if (input$grouping == "Individual") {
          # Create the plot
          p <- ggplot(total_data, aes(x = Individual, y = Value, color = Particle)) +
            geom_point(position = position_jitter(width = 0.1), size = 3) +
            theme_minimal(base_size = 15) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(y = "Normalized counts",
                 x = NULL,  # Remove the default x-axis label
                 title = paste("Expression of", input$protein)) +
            scale_x_discrete(labels = function(x) {
              # Create custom labels based on the levels of 'Individual'
              sapply(x, function(i) {
                paste(
                  unique(total_data$Particle[total_data$Individual == i]), 
                  "", unique(total_data$Isolation[total_data$Individual == i]), 
                  "", unique(total_data$Growth[total_data$Individual == i]), 
                  "", i
                )
              })
            }) #+
          #scale_color_manual(values = c("EVs" = "steelblue", "Group2" = "darkorange"))
          
          
          # Apply faceting based on the user's choices
          if (input$facet_isolation & input$facet_growth) {
            p <- p + facet_grid(rows = vars(Isolation), cols = vars(Growth))
          } else if (input$facet_isolation) {
            p <- p + facet_grid(cols = vars(Isolation))
          } else if (input$facet_growth) {
            p <- p + facet_grid(cols = vars(Growth))
          }
          
        } else if (input$grouping == "Grouped") {
          p <- ggplot(total_data, aes(x = Particle, y = Value)) +
            geom_pirate(aes(colour = Particle), bars = FALSE,
                        points_params = list(shape = 19, alpha = 0.2),
                        lines_params = list(size = 0.8)) +
            labs(title = paste("Grouped Protein Expression for", input$protein))
          
          # Apply faceting based on the user's choices
          if (input$facet_isolation & input$facet_growth) {
            p <- p + facet_grid(rows = vars(Isolation), cols = vars(Growth))
          } else if (input$facet_isolation) {
            p <- p + facet_grid(cols = vars(Isolation))
          } else if (input$facet_growth) {
            p <- p + facet_grid(cols = vars(Growth))
          }
        }
        
        print(p)
      } else {
        ggplot() + labs(title = "No data available for selected protein")
      }
      
    } else if (input$dataset == "RNA-seq") {
      if (!is.null(input$rna) && input$rna %in% rownames(data)) {
        if (input$grouping == "Individual") {
          data_to_plot <- as.data.frame(t(data[rownames(data) == input$rna, ]))
          colnames(data_to_plot) <- "Expression"
          data_to_plot$Condition <- rownames(data_to_plot)
          
          print(head(data_to_plot))
          
          ggplot(data_to_plot, aes(x = Condition, y = Expression)) +
            geom_point(color = "steelblue", size = 4) +
            theme_minimal(base_size = 15) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = paste("Expression of", input$rna)) +
            scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
        } else if (input$grouping == "Grouped") {
          group1 <- c("EV3D")
          group2 <- c("Super2D", "Super3D", "SuperFPLC")
          group3 <- c("Exomere2D", "Exomere3D", "ExomereFPLC")
          group4 <- c("EVP2D", "EVP3D")
          
          groups <- list(EVs = group1, Super = group2, Exomere = group3, EVp = group4)
          rna_data_subset <- data[input$rna, , drop = FALSE]
          
          group_names <- c()
          expression_values <- c()
          
          for (group_name in names(groups)) {
            group_cols <- groups[[group_name]]
            valid_cols <- intersect(group_cols, colnames(rna_data_subset))
            if (length(valid_cols) > 0) {
              group_names <- c(group_names, rep(group_name, length(valid_cols)))
              expression_values <- c(expression_values, as.vector(t(rna_data_subset[, valid_cols])))
            }
          }
          
          grouped_data <- data.frame(
            Group = factor(group_names, levels = c("EVs", "Super", "Exomere", "EVp")),
            Expression = expression_values
          )
          
          pirateplot(formula = Expression ~ Group, data = grouped_data,
                     theme = "white", 
                     pal = c("EVs" = "steelblue", "Super" = "lightcoral", "Exomere" = "mediumseagreen", "EVp" = "goldenrod"),
                     ylab = "Reads per million total reads", 
                     xlab = "Group", 
                     main = paste("Grouped miRNA Expression for", input$rna))
        }
      } else {
        ggplot() + labs(title = "No data available for selected miRNA")
      }
      
    } else if (input$dataset == "Lipidomics") {
      if (!is.null(input$lipid) && input$lipid %in% rownames(data)) {
        selected_protein_data <- data[input$lipid, , drop = FALSE]
        grouping <- input$grouping
        plot_type <- ifelse(grouping == "Individual", "dotplot", "pirateplot")
        
        # Get the meta information for faceting
        meta_data <- protein_meta()
        
        individual_data <- as.vector(unlist(meta_data["Individual", colnames(selected_protein_data)]))
        partile_data <- as.vector(unlist(meta_data["Particle", colnames(selected_protein_data)]))
        isolation_data <- as.vector(unlist(meta_data["Isolation", colnames(selected_protein_data)]))
        growth_data <- as.vector(unlist(meta_data["Growth", colnames(selected_protein_data)]))
        
        total_data <- data.frame(
          Condition = colnames(selected_protein_data),
          Value = as.numeric(selected_protein_data),
          Individual = factor(individual_data),
          Particle = factor(partile_data),
          Isolation = factor(isolation_data),
          Growth = factor(growth_data)
        )
        
        print(head(total_data))
        
        if (input$grouping == "Individual") {
          # Create the plot
          p <- ggplot(total_data, aes(x = Individual, y = Value, color = Particle)) +
            geom_point(position = position_jitter(width = 0.1), size = 3) +
            theme_minimal(base_size = 15) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(y = "Normalized counts",
                 x = NULL,  # Remove the default x-axis label
                 title = paste("Expression of", input$lipid)) +
            scale_x_discrete(labels = function(x) {
              # Create custom labels based on the levels of 'Individual'
              sapply(x, function(i) {
                paste(
                  unique(total_data$Particle[total_data$Individual == i]), 
                  "", unique(total_data$Isolation[total_data$Individual == i]), 
                  "", unique(total_data$Growth[total_data$Individual == i]), 
                  "", i
                )
              })
            }) #+
          #scale_color_manual(values = c("EVs" = "steelblue", "Group2" = "darkorange"))
          
          
          # Apply faceting based on the user's choices
          if (input$facet_isolation & input$facet_growth) {
            p <- p + facet_grid(rows = vars(Isolation), cols = vars(Growth))
          } else if (input$facet_isolation) {
            p <- p + facet_grid(cols = vars(Isolation))
          } else if (input$facet_growth) {
            p <- p + facet_grid(cols = vars(Growth))
          }
          
        } else if (input$grouping == "Grouped") {
          p <- ggplot(total_data, aes(x = Particle, y = Value)) +
            geom_pirate(aes(colour = Particle), bars = FALSE,
                        points_params = list(shape = 19, alpha = 0.2),
                        lines_params = list(size = 0.8)) +
            labs(title = paste("Grouped Lipid Expression for", input$lipid))
          
          # Apply faceting based on the user's choices
          if (input$facet_isolation & input$facet_growth) {
            p <- p + facet_grid(rows = vars(Isolation), cols = vars(Growth))
          } else if (input$facet_isolation) {
            p <- p + facet_grid(cols = vars(Isolation))
          } else if (input$facet_growth) {
            p <- p + facet_grid(cols = vars(Growth))
          }
        }
        
        print(p)
      } else {
        ggplot() + labs(title = "No data available")
      }
    }
  })
  
  data <- reactive({
    switch(input$dataset,
           "Proteomics" = protein_data(),
           "RNA-seq" = rna_data(),
           "Lipidomics" = lipid_data())
  })
  
  # Render the reactable table
  output$datatable <- renderReactable({
    df <- data()
    
    if (nrow(df) == 0 || ncol(df) == 0) {
      cat("Data is empty. Returning a placeholder table.\n")
      return(reactable(data.frame(Message = "No data available."), 
                       pagination = FALSE, searchable = FALSE))
    }
    
    cat("Rendering table with valid data.\n")
    reactable(
      df,
      searchable = TRUE,
      pagination = TRUE,
      defaultPageSize = 10,
      pageSizeOptions = c(5, 10, 20),
      theme = reactable::reactableTheme(
        headerStyle = list(backgroundColor = '#f5f5f5', color = '#333'),
        cellStyle = list(backgroundColor = '#fff')
      )
    )
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

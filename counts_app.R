## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/13/22

##################
# file description: testing the counts tab of the final project R Shiny app
##################

# import libraries
library(tidyverse)
library(BiocManager)
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(RColorBrewer)

# increase file upload size to accomodate the large files
options(shiny.maxRequestSize=30*1024^2)

# design the gui
# design the gui
ui <- fluidPage(
  # sample tab title and directions
  titlePanel("Investigation of Counts"),
  markdown(paste0("To use this tab of the application, download the CSV `GSE64810_DESeq2_norm_counts.csv` from the GitHub.")),
  
  # sample tab sidebar
  sidebarLayout(
    sidebarPanel(
      fileInput("count_file",
                label = "Load normalized counts matrix",
                accept = c(".csv", ".tsv"),
                placeholder = 'GSE64810_DESeq2_norm_counts.csv'),
      sliderInput("variance",
                  "Select the minimum variance for each gene",
                  min = 0,
                  max = 511800,
                  value = 2000,
                  step = 100),
      sliderInput("nonzeros",
                  "Select the minimum number of samples that are non-zero for each gene",
                  min = 0,
                  max = 69,
                  value = 0,
                  step = 1),
      numericInput("pc1",
                  "Select the number of the first principle component to examine",
                  1,
                  min = 1,
                  max = 69),
      numericInput("pc2",
                   "Select the number of the second principle component to examine",
                   2,
                   min = 1,
                   max = 69),
      submitButton(text = "Explore Data",
                   width = '100%')
    ),
    # show the data in the main panel
    mainPanel(tabsetPanel(
      tabPanel("Summary", {tableOutput("filter_summary")}),
      tabPanel("Diagnostics", {plotOutput("counts_diagnostic1")}, {plotOutput("counts_diagnostic2")}),
      tabPanel("Heatmap", {plotOutput("counts_heatmap")}),
      tabPanel("PCA", {plotOutput("counts_PCA")})
    )
    )
  )
)

##############
# define the server logic
server <- function(input, output, session){
  
  #' Load the data, but check if it is a tsv or csv and return an error if not
  #' 
  #' @param file A text string with the full file path to the data file
  #' 
  #' @return A dataframe containing the data within the file
  #' 
  #' @details
  #' 
  #' @example diffexp <- file_load('data/GSE64810_DESeq2_diffexp.csv')  
  file_load <- reactive({
    if (grep('.csv', input$count_file$datapath) == TRUE){
      data <- read.csv(input$count_file$datapath)
    } else if (grep('.tsv', input$count_file$datapath) == TRUE){
      data <- read.tsv(input$count_file$datapath)
    } else {
      stop('File is neither a csv nor tsv. Please upload the correct file format.')
    }
    return(data)
  })
  
  #' Filter the counts matrix based on non-zero sample counts and variance
  #' 
  #' @param gene_info A dataframe with variance, count of zeros, and median count
  #' @param nonzero_slider An input from a slider that gives the minimum number of nonzeros we want
  #' @param variance_slider An input from a slider that gives the minimum variance we want
  #'
  filter_counts <- function(gene_info, nonzero_slider, variance_slider){
    # filter for the variance first
    filtered_df <- dplyr::filter(gene_info, variance > variance_slider) %>%
      # filter for the nonzeros second
      dplyr::filter(num_zero <= (69 - nonzero_slider))
    return(filtered_df)
  }
  
  # make the output for the summary of the counts matrix after filtering
  output$filter_summary <- renderTable({
    # require a file input
    req(input$count_file)
    # load in the sample data
    counts <- file_load()
    # format as necessary
    counts_rownames <- dplyr::select(counts, -X)
    rownames(counts_rownames) <- counts$X
    # find the variance, median count, and zero count for each gene
    gene_info <- tibble(gene = counts$X,
                        variance  = apply(counts_rownames, 1, var),
                        median = apply(counts_rownames, 1, median),
                        num_zero = apply(counts_rownames, 1, function(x) sum(x==0)))
    # filter the genes
    filtered_genes <- filter_counts(gene_info, input$nonzeros, input$variance)
    # find number of samples
    num_samples <- ncol(counts_rownames)
    # find the number of genes
    num_genes <- nrow(counts_rownames)
    # find how many genes pass the filter
    num_pass <- nrow(filtered_genes)
    # find the percent of genes that pass the filter
    per_pass <- (num_pass/num_genes) * 100
    # make a tibble with all the summary statistics of the filter effects
    filter_effects <- tibble(stat = c("Number of Samples", "Total Genes", "Genes Passing", "Genes Not Passing"),
                             count = c(num_samples, num_genes, num_pass, num_genes - num_pass),
                             percent = c('NA', 'NA', round(per_pass, digits = 2), round(100 - per_pass, digits = 2)))
    return(filter_effects)
  }, striped = T)
  
  # make the first output for the diagnostics tab (plot of median vs variance)
  output$counts_diagnostic1 <- renderPlot({
    # require a file input
    req(input$count_file)
    # load in the sample data
    counts <- file_load()
    # format as necessary
    counts_rownames <- dplyr::select(counts, -X)
    rownames(counts_rownames) <- counts$X
    # find the variance, median count, and zero count for each gene
    gene_info <- tibble(gene = counts$X,
                        variance  = apply(counts_rownames, 1, var),
                        median = apply(counts_rownames, 1, median),
                        num_zero = apply(counts_rownames, 1, function(x) sum(x==0)))
    # plot the median vs the variance, where the ones above a certain variance are colored red
    medVvar_plot <- ggplot(gene_info, aes(x = log10(median), y = log10(variance))) +
      geom_point(aes(color = variance > 2000), alpha = 0.5) +
      xlab('Log10(Median)') +
      ylab('Log10(Variance)') +
      ggtitle('Log10(Variance between Counts) versus Log10(Median of Counts)') +
      scale_color_manual('Outside Filter', values = c('black', 'red')) +
      theme(legend.position = 'bottom')
    return(medVvar_plot)
  })
  
  # make the second output for the diagnostics tab (plot of median and number of zeros)
  output$counts_diagnostic2 <- renderPlot({
    # require a file input
    req(input$count_file)
    # load in the sample data
    counts <- file_load()
    # format as necessary
    counts_rownames <- dplyr::select(counts, -X)
    rownames(counts_rownames) <- counts$X
    # find the variance, median count, and zero count for each gene
    gene_info <- tibble(gene = counts$X,
                        variance  = apply(counts_rownames, 1, var),
                        median = apply(counts_rownames, 1, median),
                        num_zero = apply(counts_rownames, 1, function(x) sum(x==0)))
    #make a tibble with the count of zeros with the median count
    medVzero_plot <- ggplot(gene_info, aes(x = log10(median), y = num_zero)) +
      geom_point(aes(color = num_zero > 5), alpha = 0.5) +
      xlab('Log10(Median)') +
      ylab('Count of Zeros') +
      ggtitle('Count of Zeros versus Log10(Median of Counts)') +
      scale_color_manual('Outside Filter', values = c('black', 'red')) +
      theme(legend.position = 'bottom')
    return(medVzero_plot)
  })
  
  # make output for the heatmap tab
  output$counts_heatmap <- renderPlot({
    # require a file input
    req(input$count_file)
    # load in the sample data
    counts <- file_load()
    # format as necessary
    counts_rownames <- dplyr::select(counts, -X)
    rownames(counts_rownames) <- counts$X
    # find the variance, median count, and zero count for each gene
    gene_info <- tibble(gene = counts$X,
                        variance  = apply(counts_rownames, 1, var),
                        median = apply(counts_rownames, 1, median),
                        num_zero = apply(counts_rownames, 1, function(x) sum(x==0)))
    # filter the genes
    filtered_genes <- filter_counts(gene_info, input$nonzeros, input$variance)
    # find the counts of filtered genes
    filtered_counts <- dplyr::filter(counts, X%in%filtered_genes$gene)
    # format the data so that the heatmap function will accept it
    rownames(filtered_counts) <- filtered_counts$X
    filtered_counts <- dplyr::select(filtered_counts, -X) %>%
      as.matrix()
    # make a heatmap of the filtered counts
    hm <- heatmap(filtered_counts,
                  col= colorRampPalette(brewer.pal(9, "BuPu"))(25))
    # make a legend for the heatmap
    legend <- legend(x = 'bottomright', legend=c('low','high'),
                     fill=colorRampPalette(brewer.pal(11, "BuPu"))(2))
    return(hm)
  })
  
  # make the output for the PCA plot of two principle components
  output$counts_PCA <- renderPlot({
    # require a file input
    req(input$count_file)
    # load in the sample data
    counts <- file_load()
    # format as necessary
    counts_rownames <- dplyr::select(counts, -X)
    rownames(counts_rownames) <- counts$X
    # do pca
    pca_results <- prcomp(scale(t(counts_rownames)), center = FALSE, scale = FALSE)
    # pull out info related to all pcs
    pca_tbl <- as_tibble(pca_results$x)
    # calculate variance explained
    v <- pca_results$sdev^2
    ve <- round((v / sum(v)) * 100, digits = 2)
    # make variables for the principle component inputs
    pc1 <- paste("PC", as.character(input$pc1), sep="")
    pc2 <- paste("PC", as.character(input$pc2), sep="")
    # plot the two pcs 
    biplot <- pca_tbl %>%
      ggplot() +
      geom_point(aes(x=!!sym(pc1), y=!!sym(pc2))) +
      xlab(paste(pc1, ' (', as.character(ve[as.integer(gsub("^.*?PC","", pc1))]), "% VE)")) +
      ylab(paste(pc2, ' (', as.character(ve[as.integer(gsub("^.*?PC","", pc2))]), "% VE)")) +
      ggtitle('Percent Variance Explained by the Chosen Principle Components')
    return(biplot)
  })
}

##############
# run the application
shinyApp(ui = ui, server = server)



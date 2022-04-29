## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/13/22

##################
# file description: all apps for the four tabs of the main app compiled into one
#                   R Shiny app for the final project
##################

# import libraries
library(tidyverse)
library(BiocManager)
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(RColorBrewer)
library(biomaRt)
library(ggbeeswarm)

# increase file upload size to allow large files
options(shiny.maxRequestSize=30*1024^2)

# put choice options for any of the tabs here
sample_choices <- c("sample_info", "sample_name", "geo_accession", "diagnosis",
                    "pmi", "age_of_death", "rin", "id_ref") #used in the samples information tab
deseq_choices <- c("X", "symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj") #used in the deseq tab
ige_choices <- c("bar plot", "boxplot", "violin plot", "beeswarm plot") #used in the individual gene expression tab
sort_choices <- c('ascending', 'descending') #used in both samples and deseq tabs

################
# design the ui
ui <- shinyUI(
  
  # overall title
  navbarPage("BF591 Project: Adalee Koshiol",
  
  # start of samples information tab
  tabPanel("Samples", value="samples",
    # sample tab title and directions
    titlePanel("Investigation of Samples"),
      markdown(paste0("To use this tab of the application, download the CSV `GSE64810_processed_sample_info` from the GitHub.",
                      " Allow for ample loading time.")),
      markdown(paste0("This tab can be used to look at distinct values and distributions of the samples. ",
                      "It loads a matrix of sample information and gives the following information in separate tabs: ")),
      markdown(paste0("1) Summary of Samples: a summary table that includes data types, means of values, and example values, ")),
      markdown(paste0("2) Sample Information: a table of the samples that can be sorted (which can be sorted using the 'choose how to sort the data' selections), ")),
      markdown(paste0("3) Plots of Samples: a histogram of the number of neurologically normal samples versus ones with Huntington's.")),
    # sidebar for samples information tab
    sidebarLayout(
      sidebarPanel(
        # need csv file for sample information
        fileInput("sample_file",
                  label = "Load processed sample info",
                  accept = c(".csv", ".tsv"),
                  placeholder = 'GSE64810_processed_sample_info'),
        # have a button that allows the user to sort the data
        radioButtons(inputId = 'sort_sample',
                     label = 'Choose how to sort the data',
                     choices = sample_choices,
                     selected = 'sample_name'),
        # have a button that allows the user to sort the data in ascending or descending order
        radioButtons(inputId = 'sample_up_down',
                     label = 'Choose whether the data should be shown in ascending or descending order',
                     choices = sort_choices,
                     selected = 'ascending'),
        # have a submit button
        submitButton(text = "Explore Data",
                     width = '100%')
      ), #close for the 'sidebarpanel' function
      # main panel for the samples information tab
      mainPanel(tabsetPanel(
        tabPanel("Summary of Samples", {tableOutput("summary_info")}),
        tabPanel("Sample Information", {tableOutput("sample_data")}),
        tabPanel("Plots of Samples", plotOutput("sample_histogram"))
      )) #close of samples information mainpanel function
    ) #close of samples information sidebarlayout function 
  ), #close of samples information tab
  
  # start of counts tab
  tabPanel("Counts", value = "counts",
           # counts tab title and directions
           titlePanel("Investigation of Count Data"),
           markdown(paste0("To use this tab of the application, download the CSV `GSE64810_DESeq2_norm_counts.csv` from the GitHub.",
                           " Allow for ample loading time.")),
           markdown(paste0("This tab can be used to look at the count matrix to understand the count structure and aid in gene filtering. ",
                           "It loads a matrix of count information, as well as user input for percent variance and number of allowed non-zero samples,",
                           " and gives the following information in separate tabs: ")),
           markdown(paste0("1) Summary: a summary table that shows the effect of the user's filtering input, ")),
           markdown(paste0("2) Diagnostics: diagnositc scatter plots of median count vs variance and median count vs non-zeros (where genes that pass filters are shown in black), ")),
           markdown(paste0("3) Heatmap: a clustered heatmap of counts that remain after filtering, ")),
           markdown(paste0("4) PCA: PCA projections based on user input for which principle components to examine.")),
           # sidebar for the counts tab
           sidebarLayout(
             sidebarPanel(
               # need a file input for the normalized counts matrix
               fileInput("count_file",
                         label = "Load normalized counts matrix",
                         accept = c(".csv", ".tsv"),
                         placeholder = 'GSE64810_DESeq2_norm_counts.csv'),
               # allow the user to select the minimum variance
               sliderInput("variance",
                           "Select the minimum variance for each gene",
                           min = 0,
                           max = 511800,
                           value = 2000,
                           step = 100),
               # allow user to select the minimum number of non-zero samples for each gene
               sliderInput("nonzeros",
                           "Select the minimum number of samples that are non-zero for each gene",
                           min = 0,
                           max = 69,
                           value = 34,
                           step = 1),
               # allow the user to select a principle component for plotting
               numericInput("pc1",
                            "Select the number of the first principle component to examine",
                            1,
                            min = 1,
                            max = 69),
               # allow the user to select a second principle component to plot against the other
               numericInput("pc2",
                            "Select the number of the second principle component to examine",
                            2,
                            min = 1,
                            max = 69),
               # allow the user to select the number of principle components to plot in a beeswarm
               numericInput("beeswarm_pca",
                            "Select the number of top principle components to put in a beeswarm plot",
                            10,
                            min = 1,
                            max = 69),
               # have a submit button
               submitButton(text = "Explore Data",
                            width = '100%')
             ), #close of sidebarpanel function
             # main panel for the counts tab
             mainPanel(tabsetPanel(
               tabPanel("Summary", {tableOutput("filter_summary")}),
               tabPanel("Diagnostics", {plotOutput("counts_diagnostic1")}, {plotOutput("counts_diagnostic2")}),
               tabPanel("Heatmap", {plotOutput("counts_heatmap")}),
               tabPanel("PCA", {plotOutput("counts_PCA")}, {plotOutput("counts_beeswarm")})
               )) #close of the counts mainpanel function
          ) #close of counts sidebarlayout function
        ), #close of the counts tab
  
  # start of deseq tab
  tabPanel("DESeq", value="deseq",
           # deseq title and instructions
           titlePanel("Investigation of DESeq Data"),
           markdown(paste0("To use this tab of the application, download the CSV `GSE64810_DESeq2_diffexp.csv` from the GitHub.",
                           " Allow for ample loading time.")),
           markdown(paste0("This tab can be used to look at the results of differential expression using DESeq. ",
                           "It loads a matrix of deseq results and user input for a desired p-value",
                           " and gives the following information in separate tabs: ")),
           markdown(paste0("1) Summary: a table with the differential expression results, ")),
           markdown(paste0("2) Volcano Plot: a volcano plot of fold change values, ")),
           # sidebar for deseq tab
           sidebarLayout(
             sidebarPanel(
               fileInput("deseq_file",
                         label = "Load DESeq2 data",
                         accept = c(".csv", ".tsv"),
                         placeholder = 'GSE64810_DESeq2_diffexp.csv'),
               radioButtons("deseq_sort",
                            "Choose how to sort the data",
                            choices = deseq_choices,
                            selected = "padj"),
               sliderInput("padj_slider",
                           "Select the magnitude of the p adjusted coloring:",
                           min = -20,
                           max = 0,
                           value = -5,
                           step = 1),
               submitButton(text = "Explore Data",
                            width = '100%')
             ), #close of deseq sidebarpanel function
             # show the data in the main panel
             mainPanel(tabsetPanel(
               tabPanel("Data Table", {tableOutput("deseq_table")}),
               tabPanel("Volcano Plot", {plotOutput("volcano_plot")})
             )) #close of deseq mainPanel function
           ) #close of deseq sidebarlayout function
        ), #close of deseq tab
  
  # start of individual gene expression (IGE) tab
  tabPanel("Individual Gene Expression", value="ige",
           # ige title and instructions
           titlePanel("Invesitgation of Individual Gene Expressions from Normalized Counts Data"),
           markdown(paste0("To use this tab of the application, download the CSV `GSE64810_DESeq2_norm_counts.csv` and `GSE64810_processed_sample_info` from the GitHub.",
                           " Allow for ample loading time.")),
           markdown(paste0("This tab can be used to look at the expression of individual genes using counts. ",
                           "It loads a matrix of count information, information on samples, ",
                           "user input for categorical fields to use, ",
                           "user input for a gene name, ",
                           "and user input for the type of plot to make.")),
           # sidebar for ige tab
           sidebarLayout(
             sidebarPanel(
               # need a file input for the normalized counts matrix
               fileInput("count_file2",
                         label = "Load normalized counts matrix",
                         accept = c(".csv", ".tsv"),
                         placeholder = 'GSE64810_DESeq2_norm_counts.csv'),
               # need csv file for sample information
               fileInput("sample_file2",
                         label = "Load processed sample info",
                         accept = c(".csv", ".tsv"),
                         placeholder = 'GSE64810_processed_sample_info'),
               # have a button that allows the user to choose what to plot
               radioButtons(inputId = 'sample_categorical',
                            label = 'Choose which categorical variable to plot by',
                            choices = sample_choices,
                            selected = 'diagnosis'),
               # allow for user input for a gene name to search
               textInput(inputId = 'gene_search',
                         label = 'Choose a gene to explore. Example: ENSG00000000003.10',
                         placeholder = "ENSG00000000003.10"),
               # have a button that allows the user to choose a type of plot
               radioButtons(inputId = 'plot_type',
                            label = 'Choose what the type of plot',
                            choices = ige_choices,
                            selected = 'boxplot'),
               submitButton(text = "Explore Data",
                            width = '100%')
             ), #close of ige sidebarpanel function
           mainPanel(tabsetPanel(
             tabPanel("Individual Gene Expression Plot", {plotOutput("ige_plot")})
           )) #close of ige mainpanel function
      ) #close of ige sidebarlayout function
  ) #close of IGE tab 
  
  ) #close of navbarpage function
) #close of shiny ui function

##############
# define the server logic
server <- function(input, output, session){
  
  ############# START OF OUTPUT FOR SAMPLES INFORMATION TAB
  #' Load the data, but check if it is a tsv or csv and return an error if not
  #' 
  #' @param file A text string with the full file path to the data file
  #' 
  #' @return A dataframe containing the data within the file
  sample_load <- reactive({
    if (grep('.csv', input$sample_file$datapath) == TRUE){
      data <- read.csv(input$sample_file$datapath)
    } else if (grep('.tsv', input$sample_file$datapath) == TRUE){
      data <- read.tsv(input$sample_file$datapath)
    } else {
      stop('File is neither a csv nor tsv. Please upload the correct file format.')
    }
    return(data)
  })
  
  # create the output for the summary information table tab
  output$summary_info <- renderTable({
    # require a file input
    req(input$sample_file)
    # load in the sample data
    samples <- sample_load()
    # process the sample information
    sample_info <- data.frame(t(samples))
    sample_info$sample_name <- rownames(sample_info)
    sample_info <- rename(sample_info, geo_accession = X1)
    sample_info <- rename(sample_info, diagnosis = X2)
    sample_info <- rename(sample_info, pmi = X3)
    sample_info$pmi <- as.integer(sample_info$pmi)
    sample_info <- rename(sample_info, age_of_death = X4)
    sample_info$age_of_death <- as.integer(sample_info$age_of_death)
    sample_info <- rename(sample_info, rin = X5)
    sample_info$rin <- as.numeric(sample_info$rin)
    sample_info <- rename(sample_info, id_ref = X6)
    sample_info <- sample_info[2:70,]
    sample_info <- dplyr::select(sample_info, sample_name, geo_accession, diagnosis,
                                 pmi, age_of_death, rin, id_ref)
    # make the summary table of the sample information
    column_names <- colnames(sample_info)
    type <- c(typeof(sample_info$sample_name), typeof(sample_info$geo_accession),
              typeof(sample_info$diagnosis), typeof(sample_info$pmi),
              typeof(sample_info$age_of_death), typeof(sample_info$rin),
              typeof(sample_info$id_ref))
    mean <- c('NA', 'NA', 'NA',
              round(mean(sample_info$pmi, na.rm = TRUE), digits = 2),
              round(mean(sample_info$age_of_death), digits = 2),
              round(mean(sample_info$rin), digits = 2),
              'NA')
    example <- c(sample_info$sample_name[1],
                 sample_info$geo_accession[1],
                 sample_info$diagnosis[1],
                 sample_info$pmi[1],
                 sample_info$age_of_death[1],
                 sample_info$rin[1],
                 sample_info$id_ref[1])
    # put all the summary information into a tibble
    summary_info <- tibble(column_names, type, mean, example)
    return(summary_info)
  }, striped = T)
  
  # create the output for the whole data table of samples, sorted how the user instructs
  output$sample_data <- renderTable({
    # require an input file
    req(input$sample_file)
    # load in the data
    samples <- sample_load()
    # process the sample information
    sample_info <- data.frame(t(samples))
    sample_info$sample_name <- rownames(sample_info)
    sample_info <- rename(sample_info, geo_accession = X1)
    sample_info <- rename(sample_info, diagnosis = X2)
    sample_info <- rename(sample_info, pmi = X3)
    sample_info$pmi <- as.integer(sample_info$pmi)
    sample_info <- rename(sample_info, age_of_death = X4)
    sample_info$age_of_death <- as.integer(sample_info$age_of_death)
    sample_info <- rename(sample_info, rin = X5)
    sample_info$rin <- as.numeric(sample_info$rin)
    sample_info <- rename(sample_info, id_ref = X6)
    sample_info <- sample_info[2:70,]
    sample_info <- dplyr::select(sample_info, sample_name, geo_accession, diagnosis,
                                 pmi, age_of_death, rin, id_ref)
    # look at how the user wants the data sorted (asc or desc, and which column)
    if (input$sample_up_down == 'descending'){
      sample_arranged <- dplyr::arrange(sample_info, desc(!!sym(input$sort_sample)))
    } else{
      sample_arranged <- dplyr::arrange(sample_info, !!sym(input$sort_sample))
    }
    return(sample_arranged)
  }, striped = T)
  
  # create the example histogram for the sample data
  output$sample_histogram <- renderPlot({
    # require an input file
    req(input$sample_file)
    # load in the data
    samples <- sample_load()
    # process the sample information
    sample_info <- data.frame(t(samples))
    sample_info$sample_name <- rownames(sample_info)
    sample_info <- rename(sample_info, geo_accession = X1)
    sample_info <- rename(sample_info, diagnosis = X2)
    sample_info <- rename(sample_info, pmi = X3)
    sample_info$pmi <- as.integer(sample_info$pmi)
    sample_info <- rename(sample_info, age_of_death = X4)
    sample_info$age_of_death <- as.integer(sample_info$age_of_death)
    sample_info <- rename(sample_info, rin = X5)
    sample_info$rin <- as.numeric(sample_info$rin)
    sample_info <- rename(sample_info, id_ref = X6)
    sample_info <- sample_info[2:70,]
    sample_info <- dplyr::select(sample_info, sample_name, geo_accession, diagnosis,
                                 pmi, age_of_death, rin, id_ref)
    # create a histogram of the sample data
    h <- ggplot(sample_info, aes(diagnosis)) +
      geom_histogram(stat = 'count',
                     color = 'black',
                     fill = 'white') +
      xlab('Diagnosis') +
      ylab('Count') +
      ggtitle("Count of Normal vs Huntington's Disease")
    return(h)
  })
  ############# END OF SAMPLES INFORMATION TAB OUTPUTS
  
  ############# START OF OUTPUT FOR COUNTS TAB 
  
  #' Load the data, but check if it is a tsv or csv and return an error if not
  #' 
  #' @param file A text string with the full file path to the data file
  #' 
  #' @return A dataframe containing the data within the file
  #' 
  #' @details
  #' 
  #' @example diffexp <- file_load('data/GSE64810_DESeq2_diffexp.csv')  
  counts_load <- reactive({
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
    counts <- counts_load()
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
    counts <- counts_load()
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
    counts <- counts_load()
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
    counts <- counts_load()
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
    # load in the counts data
    counts <- counts_load()
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
    
    # make the output for the top N components in a beeswarm plot
    output$counts_beeswarm <- renderPlot({
      # require a file input
      req(input$count_file)
      # load in the counts data
      counts <- counts_load()
      # format as necessary
      counts_rownames <- dplyr::select(counts, -X)
      rownames(counts_rownames) <- counts$X
      # do pca
      pca_results <- prcomp(scale(t(counts_rownames)), center = FALSE, scale = FALSE)
      # pull out info related to all pcs
      pca_tbl <- as_tibble(pca_results$x)
      # take only the top N components
      bee_pca_data <- pca_tbl[,1:input$beeswarm_pca]
      # pivot the data for plotting
      bee_pca_data <- bee_pca_data %>%
        tidyr::pivot_longer(colnames(bee_pca_data), names_to="PCA_num", values_to="values")
      # plot the data in a beeswarm
      bee_pca <- bee_pca_data %>%
        ggplot() +
        geom_beeswarm(aes(x = PCA_num, y = values, color = PCA_num)) +
        ggtitle(paste0('Top ', as.character(input$beeswarm_pca), ' Principle Components')) +
        xlab('Principle Component') +
        ylab('') +
        labs(color = "Principle Component")
      return(bee_pca)
    })
  
  ############# END OF COUNTS OUTPUTS
  
  ############# START OF OUTPUT FOR DESEQ TAB
  #' Load the data, but check if it is a tsv or csv and return an error if not
  #' 
  #' @param file A text string with the full file path to the data file
  #' 
  #' @return A dataframe containing the data within the file
  #' 
  #' @details
  #' 
  #' @example diffexp <- file_load('data/GSE64810_DESeq2_diffexp.csv')  
  deseq_load <- reactive({
    if (grep('.csv', input$deseq_file$datapath) == TRUE){
      data <- read.csv(input$deseq_file$datapath)
    } else if (grep('.tsv', input$deseq_file$datapath) == TRUE){
      data <- read.tsv(input$deseq_file$datapath)
    } else {
      stop('File is neither a csv nor tsv. Please upload the correct file format.')
    }
    return(data)
  })
  
  # make the table for the table tab and have sort it by input
  output$deseq_table <- renderTable({
    # require a file input
    req(input$deseq_file)
    # load the data
    deseq_data <- deseq_load()
    # sort the data based on radio button input
    deseq_sorted <- dplyr::arrange(deseq_data, input$deseq_sort)
    return(deseq_sorted)
  }, striped= T)
  
  # make the volcano plot for the plot tab
  output$volcano_plot <- renderPlot({
    # require a file input
    req(input$deseq_file)
    # load the data
    deseq_data <- deseq_load()
    # plot the data with the padj coloring based on the slider
    vp <- deseq_data %>%
      ggplot(aes(x=log2FoldChange, y= -log10(padj))) +
      geom_point(aes(color = padj < 1 * 10 ^ (as.numeric(input$padj_slider)))) +
      theme_bw() +
      scale_color_manual(values = c("black", "red")) +
      theme(legend.position = "bottom") +
      xlab("log2FoldChange") +
      ylab("-log10(padj)") +
      ggtitle("Adjusted P-value versus log2FoldChange")
    return(vp)
  }, height = 700)
  ############# END OF DESEQ OUTPUTS
  
  ############# START OF OUTPUT FOR INDIVIDUAL GENE EXPRESSION TAB
  #' Load the data, but check if it is a tsv or csv and return an error if not
  #' 
  #' @param file A text string with the full file path to the data file
  #' 
  #' @return A dataframe containing the data within the file
  sample_load2 <- reactive({
    if (grep('.csv', input$sample_file2$datapath) == TRUE){
      data <- read.csv(input$sample_file2$datapath)
    } else if (grep('.tsv', input$sample_file2$datapath) == TRUE){
      data <- read.tsv(input$sample_file2$datapath)
    } else {
      stop('File is neither a csv nor tsv. Please upload the correct file format.')
    }
    return(data)
  })
  
  #' Load the data, but check if it is a tsv or csv and return an error if not
  #' 
  #' @param file A text string with the full file path to the data file
  #' 
  #' @return A dataframe containing the data within the file
  #' 
  #' @details
  #' 
  #' @example diffexp <- file_load('data/GSE64810_DESeq2_diffexp.csv')  
  counts_load2 <- reactive({
    if (grep('.csv', input$count_file2$datapath) == TRUE){
      data <- read.csv(input$count_file2$datapath)
    } else if (grep('.tsv', input$count_file2$datapath) == TRUE){
      data <- read.tsv(input$count_file2$datapath)
    } else {
      stop('File is neither a csv nor tsv. Please upload the correct file format.')
    }
    return(data)
  })
  
  output$ige_plot <- renderPlot({
    # require the files for the plotting
    req(input$count_file2)
    req(input$sample_file2)
    # load the count data
    counts2 <- counts_load2()
    # load in the sample information data
    samples2 <- sample_load2()
    # process the sample information
    sample_info <- data.frame(t(samples2))
    sample_info$sample_name <- rownames(sample_info)
    sample_info <- rename(sample_info, geo_accession = X1)
    sample_info <- rename(sample_info, diagnosis = X2)
    sample_info <- rename(sample_info, pmi = X3)
    sample_info$pmi <- as.integer(sample_info$pmi)
    sample_info <- rename(sample_info, age_of_death = X4)
    sample_info$age_of_death <- as.integer(sample_info$age_of_death)
    sample_info <- rename(sample_info, rin = X5)
    sample_info$rin <- as.numeric(sample_info$rin)
    sample_info <- rename(sample_info, id_ref = X6)
    sample_info <- sample_info[2:70,]
    sample_info <- dplyr::select(sample_info, sample_name, geo_accession, diagnosis,
                                 pmi, age_of_death, rin, id_ref)
    # filter for data only related to the chosen gene
    counts_filtered <- dplyr::filter(counts2, input$gene_search==X)
    # transpose the data
    counts_transposed <- t(counts_filtered)
    # get rid of the first row (the gene name)
    counts_transposed <- counts_transposed[2:nrow(counts_transposed),]
    # add the counts to the sample information dataframe
    sample_info$counts <- as.numeric(counts_transposed)
    # look at what the choice is and plot the corresponding plot type
    if (input$plot_type == 'bar plot'){
      # make a bar plot with the data
      ige_plot <- ggplot(sample_info, aes(x=!!sym(input$sample_categorical), y=counts)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
        xlab(input$sample_categorical) +
        ylab("counts") +
        ggtitle(paste("Counts vs.", input$sample_categorical, "for gene", input$gene_search))
    } else if (input$plot_type == 'boxplot'){
      # make a box plot with the data
      ige_plot <- ggplot(sample_info) +
        geom_boxplot(aes(x=!!sym(input$sample_categorical), y=counts)) +
        xlab(input$sample_categorical) +
        ylab("counts") +
        ggtitle(paste("Counts vs.", input$sample_categorical, "for gene", input$gene_search))
    } else if (input$plot_type == 'violin plot'){
      # make a violin plot with the data
      ige_plot <- ggplot(sample_info) +
        geom_violin(aes(x=!!sym(input$sample_categorical), y=counts), fill = "steelblue") +
        xlab(input$sample_categorical) +
        ylab("counts") +
        ggtitle(paste("Counts vs.", input$sample_categorical, "for gene", input$gene_search))
    } else if (input$plot_type == 'beeswarm plot'){
      # make a beeswarm plot with the data
      ige_plot <- ggplot(sample_info) +
        geom_beeswarm(aes(x=!!sym(input$sample_categorical), y=counts)) +
        xlab(input$sample_categorical) +
        ylab("counts") +
        ggtitle(paste("Counts vs.", input$sample_categorical, "for gene", input$gene_search))
    }
    return (ige_plot)
  })
  ############# END OF INDIVIDUAL GENE EXPRESSION OUTPUTS
  
} #close of the server function

##############
# run the application
shinyApp(ui = ui, server = server)


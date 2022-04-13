## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/12/22

##################
# file description: testing the samples tab of the final project R Shiny app
##################

# import libraries
library(tidyverse)
library(BiocManager)
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)

sample_choices <- colnames(sample_info)
sort_choices <- c('ascending', 'descending')

# design the gui
ui <- fluidPage(
  # sample tab title and directions
  titlePanel("Investigation of Samples"),
  markdown(paste0("To use this tab of the application, download the CSV `GSE64810_processed_sample_info` from the GitHub.")),
  
  # sample tab sidebar
  sidebarLayout(
    sidebarPanel(
      fileInput("sample_file",
                label = "Load processed sample info",
                accept = c(".csv", ".tsv"),
                placeholder = 'GSE64810_processed_sample_info'),
      radioButtons(inputId = 'sort_sample',
                   label = 'Choose how to sort the data',
                   choices = sample_choices,
                   selected = 'sample_name'),
      radioButtons(inputId = 'sample_up_down',
                   label = 'Choose whether the data should be shown in ascending or descending order',
                   choices = sort_choices,
                   selected = 'ascending'),
      submitButton(text = "Explore Data",
                   width = '100%')
    ),
    # show the data in the main panel
    mainPanel(tabsetPanel(
      tabPanel("Summary", {tableOutput("summary_info")}),
      tabPanel("Sample Information", {tableOutput("sample_data")}),
      tabPanel("Plots", plotOutput("sample_histogram"))
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
  #' @example diffexp <- file_load('/BF591-project/GSE64810_DESeq2_diffexp.csv')  
  file_load <- reactive({
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
    samples <- file_load()
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
    samples <- file_load()
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
    samples <- file_load()
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
}

##############
# run the application
shinyApp(ui = ui, server = server)




## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/14/22

##################
# file description: a testing file for uploading and exploring DESeq data;
#                   code that works here will be modified for the final 
#                   R Shiny app
##################

# set working directory
setwd("~/Desktop/BU/BF591/BF591-project")

# import libraries
library(tidyverse)
library(BiocManager)
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)

# increase file upload size to accomodate the large files
options(shiny.maxRequestSize=30*1024^2)

# define choices for deseq sorting
deseq_choices = c("X", "symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# Define UI for application that draws a histogram
ui <- fluidPage(
  # deseq title and instructions
  titlePanel("Investigation of DESeq Data"),
  markdown(paste0("To use this tab of the application, download the CSV `GSE64810_DESeq2_diffexp.csv` from the GitHub.")),
  
  # deseq tab sidebar
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
    ),
    # show the data in the main panel
    mainPanel(tabsetPanel(
      tabPanel("Data Table", {tableOutput("deseq_table")}),
      tabPanel("Volcano Plot", {plotOutput("volcano_plot")})
    )
    )
  )

)

# define server logic
server <- function(input, output, session) {
  
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
    deseq_data <- file_load()
    # sort the data based on radio button input
    deseq_sorted <- dplyr::arrange(deseq_data, input$deseq_sort)
    return(deseq_sorted)
  }, striped= T)
  
  # make the volcano plot for the plot tab
  output$volcano_plot <- renderPlot({
    # require a file input
    req(input$deseq_file)
    # load the data
    deseq_data <- file_load()
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
}

# Run the application
shinyApp(ui = ui, server = server)




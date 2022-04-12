## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/12/22

##################
# file description: the final project R Shiny app that explores the Post-mortem
#                   Huntington's Disease data set (GSE64810 on NCBI)
##################

# load necessary libraries
library(shiny)
library(ggplot2)
library(colourpicker)
library(tidyverse)

# design the gui
ui <- fluidPage()

# define the server logic
server <- function(input, output, session){
  
}

# run the application
shinyApp(ui = ui, server = server)
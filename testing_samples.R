# import necessary libraries
library(tidyverse)
library(BiocManager)

# set working directory
setwd("~/Desktop/BU/BF591/BF591-project")

#' Load the data, but check if it is a tsv or csv and return an error if not
#' 
#' @param file A text string with the full file path to the data file
#' 
#' @return A dataframe containing the data within the file
#' 
#' @details
#' 
#' @example diffexp <- file_load('/BF591-project/GSE64810_DESeq2_diffexp.csv')  
file_load <- function(file){
  if (grep('.csv', file) == TRUE){
    data <- read.csv(file)
  } else if (grep('.tsv', file) == TRUE){
    data <- read.tsv(file)
  } else {
    stop('File is neither a csv nor tsv. Please upload the correct file format.')
  }
  return(data)
}

# load in both the counts and differential expression data files
diffexp <- file_load('data/GSE64810_DESeq2_diffexp.csv')
counts <- file_load('data/GSE64810_DESeq2_norm_counts.csv')




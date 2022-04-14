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

#######
# for the table tab
#######

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


# load in the deseq data
deseq_data <- file_load('data/GSE64810_DESeq2_diffexp.csv') %>%
  # only keep the important information in the table
  dplyr::select(X, symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# enable the sorting
deseq_sorted <- dplyr::arrange(deseq_data, pvalue)

#########
# for the plot tab
#########

# make a volcano plot
volcano_plot <- ggplot(deseq_data, aes(x=log2FoldChange, y= -log10(padj))) +
  geom_point(aes(color = padj < 1 * 10 ^ (-20))) +
  theme_bw() +
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "bottom") +
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  ggtitle("Adjusted P-value versus log2FoldChange")

# note: make the slider go from -20 to 0; anything below -20 doesn't really show up




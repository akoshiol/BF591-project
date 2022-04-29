## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/17/22

##################
# file description: a testing file for uploading and exploring
#                   individual gene expression
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
#' @example diffexp <- file_load('data/GSE64810_DESeq2_diffexp.csv')  
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

# load the count data
counts <- file_load('data/GSE64810_DESeq2_norm_counts.csv')
counts_rownames <- dplyr::select(counts, -X)
rownames(counts_rownames) <- counts$X

# load in both the sample files
samples <- file_load('data/GSE64810_processed_sample_info.csv')
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

# have an example gene for testing
ex_gene <- "ENSG00000000003.10"
# have an example categorical variable for testing
ex_variable <- 'diagnosis'

# filter for data only related to that gene
counts_filtered <- dplyr::filter(counts, ex_gene==X)
# transpose the data
counts_transposed <- t(counts_filtered)
# get rid of the first row (the gene name)
counts_transposed <- counts_transposed[2:nrow(counts_transposed),]
# add the counts to the sample information dataframe
sample_info$counts <- as.numeric(counts_transposed)

choice <- 'beeswarm plot'

if (choice == 'bar plot'){
    # make a bar plot with the data
    ige_plot <- ggplot(sample_info, aes(x=!!sym(ex_variable), y=counts)) +
      geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
      xlab(ex_variable) +
      ylab("counts") +
      ggtitle(paste("Counts vs.", ex_variable, "for gene", ex_gene))
  } else if (choice == 'boxplot'){
    # make a box plot with the data
    ige_plot <- ggplot(sample_info) +
      geom_boxplot(aes(x=!!sym(ex_variable), y=counts)) +
      xlab(ex_variable) +
      ylab("counts") +
      ggtitle(paste("Counts vs.", ex_variable, "for gene", ex_gene))
  } else if (choice == 'violin plot'){
    # make a violin plot with the data
    ige_plot <- ggplot(sample_info) +
      geom_violin(aes(x=!!sym(ex_variable), y=counts), fill = "steelblue") +
      xlab(ex_variable) +
      ylab("counts") +
      ggtitle(paste("Counts vs.", ex_variable, "for gene", ex_gene))
  } else if (choice == 'beeswarm plot'){
    # make a beeswarm plot with the data
    ige_plot <- ggplot(sample_info) +
      geom_beeswarm(aes(x=!!sym(ex_variable), y=counts)) +
      xlab(ex_variable) +
      ylab("counts") +
      ggtitle(paste("Counts vs.", ex_variable, "for gene", ex_gene))
  }







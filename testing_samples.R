## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/12/22

##################
# file description: a testing file for uploading and exploring DE data;
#                   code that works here will be modified for the final 
#                   R Shiny app
##################

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
             sample_info$id_ref[1]
             )

summary_info <- tibble(column_names, type, mean, example)

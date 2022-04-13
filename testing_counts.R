## Author: Adalee Koshiol
## Class: BF591
## Assignment: Final Project
## Start Date: 4/13/22

##################
# file description: a testing file for uploading and exploring count data;
#                   code that works here will be modified for the final 
#                   R Shiny app
##################

# import libraries
library(tidyverse)
library(BiocManager)
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(RColorBrewer)

# set working directory
setwd("~/Desktop/BU/BF591/BF591-project")

########
# for inputs
########

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

#' Filter the counts matrix based on non-zero sample counts and variance
#' 
#' @param gene_info A dataframe containing variance, median count, and zero count for each gene
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



# take in the counts matrix
counts <- file_load('data/GSE64810_DESeq2_norm_counts.csv')
counts_rownames <- dplyr::select(counts, -X)
rownames(counts_rownames) <- counts$X

# find the variance, median count, and zero count for each gene
gene_info <- tibble(gene = counts$X,
                    variance  = apply(counts_rownames, 1, var),
                    median = apply(counts_rownames, 1, median),
                    num_zero = apply(counts_rownames, 1, function(x) sum(x==0)))

###########
# for tab with table summarizing effect of filtering
###########
num_samples <- ncol(counts_rownames)
num_genes <- nrow(counts_rownames)
num_pass <- nrow(filtered_genes)
per_pass <- (num_pass/num_genes) * 100

# make a tibble with all the summary statistics of the filter effects
filter_effects <- tibble(stat = c("Number of Samples", "Total Genes", "Genes Passing", "Genes Not Passing"),
                         count = c(num_samples, num_genes, num_pass, num_genes - num_pass),
                         percent = c('NA', 'NA', round(per_pass, digits = 2), round(100 - per_pass, digits = 2)))


##########
# for tab with diagnostic scatter plots
##########

# plot the median vs the variance, where the ones above a certain variance are colored red
medVvar_plot <- ggplot(gene_info, aes(x = log10(median), y = log10(variance))) +
  geom_point(aes(color = variance > 2000), alpha = 0.5) +
  xlab('Log10(Median)') +
  ylab('Log10(Variance)') +
  ggtitle('Log10(Variance between Counts) versus Log10(Median of Counts)') +
  scale_color_manual('Outside Filter', values = c('black', 'red')) +
  theme(legend.position = 'bottom')

#make a tibble with the count of zeros with the median count
medVzero_plot <- ggplot(gene_info, aes(x = log10(median), y = num_zero)) +
  geom_point(aes(color = num_zero > 5), alpha = 0.5) +
  xlab('Log10(Median)') +
  ylab('Count of Zeros') +
  ggtitle('Count of Zeros versus Log10(Median of Counts)') +
  scale_color_manual('Outside Filter', values = c('black', 'red')) +
  theme(legend.position = 'bottom')

###########
# for tab with heatmap
###########

# find the counts of filtered genes
filtered_counts <- dplyr::filter(counts, X%in%filtered_genes$gene)
# format the data so that the heatmap function will accept it
rownames(filtered_counts) <- filtered_counts$X
filtered_counts <- dplyr::select(filtered_counts, -X) %>%
  as.matrix()

# make a heatmap of the filtered counts
hm <- heatmap(filtered_counts,
              col= colorRampPalette(brewer.pal(9, "BuPu"))(25))

legend <- legend(x = 'bottomright', legend=c('low','high'),
                 fill=colorRampPalette(brewer.pal(11, "BuPu"))(2))


############
# for tab with pca plots
############


# do pca
pca_results <- prcomp(scale(t(counts_rownames)), center = FALSE, scale = FALSE)
# pull out info related to the pcs
pca_tbl <- as_tibble(pca_results$x)
# calculate variance explained
v <- pca_results$sdev^2
ve <- round((v / sum(v)) * 100, digits = 2)

# plot the two pcs 
biplot <- pca_tbl %>%
  ggplot() +
  geom_point(aes(x=PC2, y=PC3)) +
  xlab(paste('PC2 (', as.character(ve[as.integer(gsub("^.*?PC","","PC2"))]), "% VE)")) +
  ylab(paste('PC3 (', as.character(ve[as.integer(gsub("^.*?PC","","PC3"))]), "% VE)")) +
  ggtitle('Percent Variance Explained by the Chosen Principle Components')





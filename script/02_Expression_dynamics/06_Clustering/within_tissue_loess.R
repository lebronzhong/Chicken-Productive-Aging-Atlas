#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Within tissue loess.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript SCRIPT MATRIX_TSV SAMPLE_MANIFEST_TSV TISSUE OUTPUT_TSV")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "read_age_samples.R"))
aligned <- read_age_samples(args[1], args[2], args[3])
output_file <- args[4]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
data_long <- as.data.frame(aligned$matrix, check.names = FALSE) %>%
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "value", -gene)
data_long$week <- aligned$samples$age[match(data_long$sample, aligned$samples$sample)]
data_long <- na.omit(data_long)
models <- data_long %>%
  group_by(gene) %>%
  do(model = loess(value ~ week, data = .))
#########################################
weeks <- c(20,50,80,100)
#########################################
predictions <- models %>%
  group_by(gene) %>%
  do(data.frame(week = weeks, 
                prediction = predict(.$model[[1]], newdata = data.frame(week = weeks))))

library(tidyverse)
predictions <- predictions %>%
  ungroup() %>%
  spread(key = week, value = prediction)
prediction_matrix <- predictions[,-1]
rownames(prediction_matrix) <- predictions$gene
#library(proxy)
#dist_matrix <- proxy::dist(prediction_matrix, method = "Euclidean")
dist_matrix <- dist(prediction_matrix, diag = TRUE, upper = TRUE)
hc <- hclust(dist_matrix, method = "complete")
clusters <- cutree(hc, k = 15)
predictions$cluster <- clusters


write.table(predictions,output_file,quote = F,col.names = T,row.names = F,sep = "\t")

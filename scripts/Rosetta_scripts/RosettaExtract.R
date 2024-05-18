# ROSETTA SCORE EXTRACTION

# Load the required packages
library(tidyr)
library(data.table)
library(lubridate)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggalt)
library(ggrepel)

# Set the parent directory where the subfolders are located
parent_directory <- "~/Documents/Github/master_thesis/rosetta_scores/"

# Get a list of all subfolders in the parent directory
subfolders <- list.dirs(parent_directory, recursive = FALSE)

# Initialize an empty list to store the data
data_list <- list()

# Loop through each subfolder
for (subfolder in subfolders) {
  sc_files <- list.files(path = subfolder, pattern = "scores_", full.names = TRUE)
 
  if (length(sc_files) > 0) {
    data <- rbindlist(lapply(sc_files, function(file) {
      if (file.exists(file) && !file.size(file) == 0) {
        read.table(file, fill = TRUE, header = TRUE)
      } else {
        message("Skipping empty file or file not found: ", file)
        NULL
      }
    }), fill = TRUE)
    data_list[[subfolder]] <- data
  }
}

# Combine all the data.tables into a single data.table
combined_data <- rbindlist(data_list, fill = T)

piv_combined_data <- combined_data %>% pivot_longer(cols = c(dG_AB, dG_AC, dG_AD, dG_AE, dG_AF, dG_AG, dG_AH, dG_AI, dG_AJ), names_to = "dG_IA", values_to = "dG") %>% 
  filter(!is.na(dG))

piv_combined_data <- piv_combined_data %>% 
  mutate(FILE = unlist(lapply(strsplit(piv_combined_data$description, "_model", fixed=TRUE), function(x) return(x[1])))) %>%
  unique()

# Write the filtered data to a summary file without recycles
write.csv(x = piv_combined_data, file = paste0(parent_directory, "/", "summary.csv"))


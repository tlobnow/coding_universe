library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr)

# Specify the file path
main_path <- "/Users/u_lobnow/Documents/etc/Desktop_Cleanup/BDLD_57_WP_086756477.1__Nostoc_sp_106C_x10_model_2_renamed_chains/"
file_path <- file.path(main_path,"default.sc")

# Read in the data using read.table
data <- read.table(file_path, header = TRUE, sep = "", fill = T) %>% distinct(description, .keep_all = T)

# select + filter data
fil_data <- data %>% 
  select(-SCORE.) %>%
  filter(grepl("model", x = description))

# filter for best models (here: lowest score)
best_scoring_models <- fil_data %>% 
  arrange(desc(score)) %>% tail(n = 10) %>%
  mutate(file_name = gsub("_0001$", ".pdb", description))

# copy best models into extra folder
for (i in best_scoring_models$file_name) {
  
  # provide file names with path
  source_path = file.path(main_path, i)
  
  # provide a destination folder (will be created if new)
  dest_folder = file.path(main_path, "best_scoring_models")
  if (!dir.exists(dest_folder)) {dir.create(dest_folder)}
  dest_path = file.path(main_path, "best_scoring_models", i)
  
  # copy files into destination folder
  file.copy(source_path, dest_path, overwrite = TRUE)
}

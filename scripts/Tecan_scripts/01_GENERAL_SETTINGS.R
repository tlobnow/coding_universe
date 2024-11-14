########################################################
### THIS IS RUN AS PART OF THE TECAN ANALYSIS SCRIPT ###
###    See Tecan_Analysis.Rmd for the full script.   ###
########################################################

READ_DATA = TRUE
PREP_KEY  = TRUE
JOIN_DATA = TRUE
PLOT_DATA = TRUE

SAVE      = TRUE

# Set knitr options
knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo = FALSE)
knitr::opts_knit$set(root.dir = MAIN)

# Set global options
options(stringsAsFactors = FALSE)

# Set ggplot2 theme
theme_set(cowplot::theme_cowplot())
theme_update(panel.grid.minor = element_blank())

# source the functions and read the CL key
source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")
cl_key <- fread(file.path(MAIN, "TECAN_CL_KEY.csv"), na.strings = c("NA", "NaN", ""), header = T)

# Set the input and output directories
INPUT  <- file.path(MAIN, FOLDER)
OUTPUT <- file.path(INPUT, "output") ; if (!dir.exists(OUTPUT)) {dir.create(OUTPUT)}

# Function to find matching sheets in an Excel file and read them into data frames
read_matching_sheets <- function(file, patterns) {
  sheets <- excel_sheets(file)
  matching_sheets <- sheets[str_detect(sheets, regex(paste(patterns, collapse="|"), ignore_case = TRUE))]
  
  # Read each matching sheet into a data frame
  sheet_dfs <- lapply(matching_sheets, function(sheet) read_excel(file, sheet = sheet, col_names = F))
  
  # Return a list of data frames
  names(sheet_dfs) <- matching_sheets
  return(sheet_dfs)
}


# Define column categories
basic_columns <- c("Well", "Cycle", "Time_s", "Time_min", "Plate", "Date", 
                   # any columns that follow the pattern "Replicate_"
                   paste0("Replicate_", 1:3),
                   "CELL_LINE", "ARABINOSE", "T6_PHAGE", "REPLICATE", "CONDITION",
                   "mean_OD600_normalized", "OD600_normalized", "OD600", "OD600_blank")

colors <- c("black", "red", "#E69F00","#56B4E9")
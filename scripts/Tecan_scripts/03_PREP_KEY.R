########################################################
### THIS IS RUN AS PART OF THE TECAN ANALYSIS SCRIPT ###
###    See Tecan_Analysis.Rmd for the full script.   ###
########################################################

allplates <- data.frame()

# Define patterns to match the relevant sheets
patterns <- list(
  CELL_LINE  = c("CELL", "LINE", "COHORT"),
  ARABINOSE   = c("ARA", "ARABINOSE"),
  T6_PHAGE    = c("T6", "PHAGE"),
  MOI         = c("MOI"),
  REPLICATE   = c("REPLICATE")
)

# Define the 96-well layout as well IDs in the correct order
well_ids <- outer(LETTERS[1:8], 1:12, paste0) %>% as.vector()

# Function to read and assign well IDs to each pattern
read_and_label_sheet <- function(key_file, pattern, well_ids) {
  # Read data matching the pattern
  sheet_data <- read_matching_sheets(key_file, patterns = pattern)
  flat_data <- unlist(sheet_data) # Flatten any nested list structure
  
  # Pad with "EMPTY" if sheet_data is shorter than 96
  if (length(flat_data) < 96) {
    flat_data <- c(flat_data, rep("EMPTY", 96 - length(flat_data)))
  } else if (length(flat_data) > 96) {
    stop("Sheet data has more than 96 entries; please check the input file.")
  }
  
  # Ensure the length is exactly 96, then assign well IDs
  data.frame(Well = well_ids, DATA = flat_data)
}

# Initialize the main PLATE data frame with Well IDs
PLATE <- data.frame(Well = well_ids)

# Loop through each pattern and add data for each category
for (name in names(patterns)) {
  # Read and label the data
  labeled_sheet <- read_and_label_sheet(key_file, patterns[[name]], well_ids)
  PLATE[[name]] <- labeled_sheet$DATA
}

# Filter out rows with "EMPTY" for CELL_LINE (optional, if needed)
PLATE <- PLATE %>% filter(CELL_LINE != "EMPTY")

# Combine all plates into the main allplates data frame
allplates <- rbind(allplates, PLATE)

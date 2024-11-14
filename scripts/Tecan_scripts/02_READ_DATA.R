########################################################
### THIS IS RUN AS PART OF THE TECAN ANALYSIS SCRIPT ###
###    See Tecan_Analysis.Rmd for the full script.   ###
########################################################

# Define the regex pattern for the plate data files
pattern <- "^\\d{8}_Plate_\\d{1,2}\\.xlsx$"

# List files that match the pattern
files <- list.files(INPUT, pattern = pattern, full.names = TRUE)

# Read the files, skip the first 57 lines, and assign the plate number to the data
# Define the function to process each file
data <- lapply(files, function(file) {
  
  # Debugging 
  # file = files[1]
  
  # Extract the plate number as a number without underscores
  plate_number <- str_extract(file, "(?<=_Plate_)\\d{1,2}")
  date_number  <- str_extract(file, "\\d{8}")
  
  # Read the file, skipping the first 57 lines, and add the Plate column
  read_xlsx(file, skip = SKIP) %>%
    mutate(Plate = as.numeric(plate_number),
           Date  = as.Date(date_number, format = "%Y%m%d"))
}) %>% bind_rows() %>%
  # reorder columns
  select(`Cycles / Well`, "Plate", "Date", 
         colnames(.)[grepl(x = colnames(.), pattern = "\\.\\.\\.")])

# determine the number of cycles based on the column names
max_cycle <- length(colnames(data)[grepl(x = colnames(data), pattern = "\\.\\.\\.")])

# rename the columns to match the expected format
colnames(data) <- c("X", "Plate", "Date", as.character(1:max_cycle))

# list all files in the input directory that end with "_plate_setup.xlsx" and read the first one (give warning if multiple files were found)
key_file <- list.files(INPUT, pattern = "_plate_setup.xlsx", full.names = T) %>% 
  {if (length(.) > 1) {warning("Multiple key files found; reading the first one.")} else 
    if (length(.) == 0) {stop("No key file found.")} else {.}} %>% .[1]

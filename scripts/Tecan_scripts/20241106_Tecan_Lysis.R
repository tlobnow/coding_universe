### plot TECAN data ### 

GENERAL = TRUE

MAIN   <- "/Users/u_lobnow/Documents/etc/PhD/Phage_Phase"

list.dirs(MAIN, full.names = F, recursive = F)
# FOLDER <- "20241106"
FOLDER <- "20241108"
# FOLDER <- "20241109"

# Latest Phage Titer Details
n_plaques           = 280
dilution_step       = -4
phage_volume_plated = 0.1

if (GENERAL) {
  
  READ_DATA = TRUE
  PREP_KEY  = TRUE
  JOIN_DATA = TRUE
  PLOT_DATA = TRUE
  
  SAVE      = TRUE
  
  # Load required packages
  pacman::p_load(tidyverse, data.table, stringr, glue, readxl, cowplot, lubridate,
                 purrr, matrixStats)
  
  # source the functions and read the CL key
  source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")
  cl_key <- fread(file.path(MAIN, "TECAN_CL_KEY.csv"), na.strings = c("NA", "NaN", ""), header = T)
  
  # Set the working directory
  setwd(MAIN)
  
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
  
}

if (READ_DATA) {
  
  # Define the regex pattern for the plate data files
  pattern <- "^\\d{8}_Plate_\\d{1,2}\\.xlsx$"
  
  # List files that match the pattern
  files <- list.files(INPUT, pattern = pattern, full.names = TRUE)
  
  # Read the files, skip the first 57 lines, and assign the plate number to the data
  data <- lapply(files, function(file) {
    
    # Debugging 
    # file = files[1]
    
    # Extract the plate number as a number without underscores
    plate_number <- str_extract(file, "(?<=_Plate_)\\d{1,2}")
    date_number  <- str_extract(file, "\\d{8}")
    
    # Read the file, skipping the first 57 lines, and add the Plate column
    read_xlsx(file, skip = 57) %>%
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
} ; glimpse(data)

if (PREP_KEY) {
  
  allplates <- data.frame()
  
  # Define patterns to match the relevant sheets
  patterns <- list(
    CELL_LINE  = c("CELL", "LINE", "COHORT"),
    ARABINOSE   = c("ARA", "ARABINOSE"),
    T6_PHAGE    = c("T6", "PHAGE"),
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
  
} ; glimpse(allplates)

if (JOIN_DATA) {
  
  # Filter for relevant rows based on pattern matching in the first column
  data_filtered <- data %>%
    # Create a helper column 'Well' to identify well ID rows (e.g., A1, B2, etc.)
    mutate(Well = if_else(str_detect(X, "^[A-H][1-9][0-2]?$"), X, NA_character_)) %>%
    # Use fill() to propagate well IDs to the subsequent rows (Mean and StDev)
    fill(Well, .direction = "down") %>%
    # Filter rows that are either well IDs, "Time [s]", "Mean", or "StDev"
    # filter(str_detect(X, "^[A-H][1-9][0-2]?$") | X == "Time [s]" | X == "Mean" | X == "StDev") %>%
    filter(!is.na(X), X != "Cycles / Well") %>%
    # Keep only the rows where Well is not NA to ensure proper pairing
    filter(!is.na(Well))
  
  data_filtered$`1`<- as.numeric(data_filtered$`1`)
  
  pos_cols <- unique(data_filtered$X[grepl(x = data_filtered$X, pattern = "^\\d+;\\d+$")])
  
  data_transformed <- data_filtered %>%
    filter(X %in% c("Mean", "StDev", pos_cols)) %>%
    pivot_longer(
      cols = as.character(1:max_cycle),
      names_to = "Cycle",
      values_to = "Value"
    ) %>%
    # Add columns to identify data type for easier selection later
    mutate(
      Measurement = case_when(
        X == "Mean" ~ "Mean",
        X == "StDev" ~ "StDev",
        TRUE ~ "Raw"
      ),
      Position = ifelse(Measurement == "Raw", X, NA)
    ) %>%
    # Convert Cycle to two-digit format for consistency
    mutate(Cycle = str_pad(Cycle, 2, pad = "0")) %>%
    # Spread out Mean, StDev, and Raw values
    pivot_wider(
      names_from = Measurement,
      values_from = Value,
      values_fill = NA#, names_repair = "unique"
    ) %>%
    group_by(Date, Plate, Well, Cycle) %>%
    fill(Mean, StDev, .direction = "downup") %>%
    # Remove rows without positions to avoid issues in the next steps
    filter(!is.na(Raw)) %>%
    select(-X) %>%
    unique()
  
  # Determine the grid size for subgrid calculations
  grid_size <- unique(data_transformed$Position) %>% strsplit(";") %>% unlist() %>% unique() %>% max() %>% as.integer() + 1
  
  # Apply subgrid calculations to the data
  mean_df <- data_transformed %>%
    group_by(Well, Plate, Date, Cycle) %>%
    nest() %>%
    mutate(original_mean      = map_dbl(data, ~ mean(.x$Mean, na.rm = TRUE)),  # Original Mean
           original_sd        = map_dbl(data, ~ mean(.x$StDev, na.rm = TRUE)),   # Original StDev
           subgrid_stats      = map(data, ~ calculate_subgrid_stats(.x, grid_size = grid_size)),  # Sub-grid calculation
           refined_OD600_mean = map_dbl(subgrid_stats, ~ .x$mean)) %>% 
    ungroup() %>%
    select(Well, Plate, Date, Cycle, original_mean, original_sd, refined_OD600_mean)
  
  # extract the rows that match the "Time [s]" pattern
  time_df <- data_filtered %>%
    # Filter for rows with "Time [s]" in the first column
    filter(X == "Time [s]") %>% 
    select(-Well) %>% unique() %>%
    # pivot_longer() to match the structure of mean_df
    pivot_longer(names_to = "Cycle", values_to = "Time_s", cols = as.character(1:max_cycle)) %>%
    mutate(Cycle = str_pad(Cycle, 2, pad = "0")) %>%
    select(-X)
  
  # extract the blanks to normalize the data
  blanks <- allplates %>% left_join(mean_df) %>% 
    filter(CELL_LINE == "BLANK") %>% 
    group_by(Cycle) %>%
    #summarise(OD600_blank = mean(OD600, na.rm = T)) %>%
    summarise(OD600_blank = mean(refined_OD600_mean, na.rm = T)) %>%
    select(Cycle, OD600_blank)
  
  # Join the mean and time data frames
  dp <- mean_df %>%
    left_join(time_df) %>%
    left_join(allplates) %>%
    left_join(blanks) %>%
    filter(CELL_LINE != "EMPTY") %>%
    group_by(Plate, Date, Cycle) %>%
    mutate(
      T6_PHAGE = case_when(
        T6_PHAGE == 0 ~ NA,
        T6_PHAGE > 0  ~ T6_PHAGE,
        TRUE ~ "Unknown"),
      CONDITION = case_when(
        T6_PHAGE == 0 ~ paste0(CELL_LINE, "_", ARABINOSE, "Ara_-T6"),
        TRUE ~ paste0(CELL_LINE, "_", ARABINOSE, "Ara_+10e-", T6_PHAGE)),
      CONDITION_2 = case_when(
        ARABINOSE == 0 ~ paste0(CELL_LINE, "-Ara"),
        TRUE ~ paste0(CELL_LINE, "+Ara")),
      OD600_normalized = as.numeric(case_when(
        # OD600 - OD600_blank > 0 ~ OD600 - OD600_blank,
        refined_OD600_mean - OD600_blank > 0 ~ refined_OD600_mean - OD600_blank,
        TRUE ~ 0)),
      Time_s = round(Time_s),
      Time_min = round(Time_s / 60),
      # Rename replicate numbers to `replicate_1`, `replicate_2`, etc.
      REPLICATE = paste0("replicate_", REPLICATE)
    ) %>%
    ungroup()
  
  # Calculate OD600 mean and SD for Cycle 01
  OD600_C01_df <- dp %>%
    filter(Cycle == "01" & !grepl("BLANK", CONDITION_2)) %>%
    group_by(CONDITION_2) %>%
    mutate(
      OD600 = unique(refined_OD600_mean),
      mean = round(mean(refined_OD600_mean, na.rm = TRUE), 2),
      stdev = round(sd(refined_OD600_mean, na.rm = TRUE), 2)
    ) %>%
    ungroup()
  
  # Ensure T6_PHAGE is treated as numeric
  OD600_C01_df <- OD600_C01_df %>% mutate(T6_PHAGE = as.numeric(T6_PHAGE)) ; unique(OD600_C01_df$T6_PHAGE)
  
  # Calculate PFU based on latest Phage Titer

  # Calculate MOI for each condition based on the mean OD600 and T6_PHAGE dilution steps
  OD600_C01_df <- OD600_C01_df %>%
    rowwise() %>%
    mutate(
      MOI = calculate_MOI_lysis_assay(
        OD600 = 0.2, conversion_factor = 8E+08,
        T6_PHAGE = 2, 
        PFU = calculate_PFU(n_plaques = n_plaques, dilution_step = dilution_step, phage_volume_plated = phage_volume_plated),
        bacterial_volume = 0.18, phage_volume_added = 0.02
      )
    ) %>%
    ungroup()

  
  # View the final data frame with OD600 values, mean, standard deviation, and MOI
  OD600_C01_df

  dp_replicates <- dp %>%dp_reCONDITION_2plicates <- dp %>%
    filter(CELL_LINE != "BLANK") %>%
    # select(Date, Plate, Cycle, REPLICATE, CELL_LINE, OD600, OD600_normalized, T6_PHAGE, ARABINOSE) %>%
    # pivot_wider(names_from = REPLICATE, values_from = OD600_normalized) %>%
    # mutate(mean_OD600_normalized = rowMeans(select(., starts_with("replicate_")), na.rm = TRUE)) %>%
    # select(-starts_with("replicate_"))
    select(Date, Plate, Cycle, REPLICATE, CELL_LINE, refined_OD600_mean, T6_PHAGE, ARABINOSE) %>%
    pivot_wider(names_from = REPLICATE, values_from = refined_OD600_mean) %>%
    mutate(mean_OD600 = rowMeans(select(., starts_with("replicate_")), na.rm = TRUE)) %>%
    select(-starts_with("replicate_"))
  
  dp_final <- dp %>%
    select(-REPLICATE) %>% unique() %>%
    # left_join(dp_replicates) %>%
    left_join(dp) %>%
    distinct(CONDITION, Cycle, Date, .keep_all = T)
  
}

if (PLOT_DATA) {
  
  data_for_plotting <- dp_final %>%
    filter(CELL_LINE != "BLANK") %>%
    filter(T6_PHAGE %in% c(1,3,5,7, NA)) %>%
    mutate(T6_PHAGE = ifelse(is.na(T6_PHAGE), "Ctrl", T6_PHAGE),
           ) %>%
    left_join(cl_key, relationship = "many-to-many")
  
  # data_for_plotting %>% 
  #   arrange(Cycle) %>%
  #   # ggplot(aes(x = Time_min, y = mean_OD600_normalized, color = ARABINOSE, group = CONDITION, alpha = -(as.numeric(T6_PHAGE)))) +
  #   ggplot(aes(x = Time_min, y = mean_OD600, color = ARABINOSE, group = CONDITION, alpha = -(as.numeric(T6_PHAGE)))) +
  #   # geom_line() +
  #   geom_smooth() +
  #   facet_wrap(~CELL_LINE+ARABINOSE, nrow = 1) +
  #   theme_cowplot() +
  #   labs(title = "OD600 over time", x = "time [min]", y = "OD600")
  # 
  # data_for_plotting %>% 
  #   arrange(Cycle) %>%
  #   # ggplot(aes(x = Time_min, y = mean_OD600_normalized, color = ARABINOSE, group = CONDITION, linetype = paste("10e-", T6_PHAGE))) +
  #   ggplot(aes(x = Time_min, y = mean_OD600, color = ARABINOSE, group = CONDITION, linetype = paste("10e-", T6_PHAGE))) +
  #   # geom_line() +
  #   geom_smooth() +
  #   facet_wrap(~CELL_LINE, nrow = 2) +
  #   theme_cowplot() +
  #   labs(title = "OD600 over time", x = "time [min]", y = "OD600")
  
  plot <- data_for_plotting %>% 
    arrange(Cycle) %>%
    # ggplot(aes(x = Time_min, y = mean_OD600_normalized, color = CONDITION_2, group = CONDITION)) +
    # ggplot(aes(x = Time_min, y = mean_OD600, color = CONDITION_2, group = CONDITION)) +
    ggplot(aes(x = Time_min, y = original_mean, color = CONDITION_2, group = CONDITION)) +
    geom_smooth() +
    facet_wrap(~T6_PHAGE, nrow = 1) +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    labs(title = "OD600 over time", x = "time [min]", y = "OD600") +
    scale_color_manual(values = c("WT -Ara" = "gray50",
                                  "WT +Ara" = "gray30",
                                  "bGSDM -Ara" = "red",
                                  "bGSDM +Ara" = "#F99999"))
  
  plot
}



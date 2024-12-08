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

# Join the mean and time data frames, checking for the presence of MOI or T6_PHAGE
# Check if MOI or T6_PHAGE columns contain data other than "EMPTY"
use_moi <- any(allplates$MOI != "EMPTY", na.rm = TRUE)
use_t6_phage <- any(allplates$T6_PHAGE != "EMPTY", na.rm = TRUE)

# Join the mean, time, and other data frames, without setting CONDITION initially
dp <- mean_df %>%
  left_join(time_df) %>%
  left_join(allplates) %>%
  left_join(blanks) %>%
  filter(CELL_LINE != "EMPTY") %>%
  group_by(Plate, Date, Cycle) %>%
  mutate(
    # Set CONDITION based on whether we have MOI or T6_PHAGE data, using only relevant column
    CONDITION = if (use_moi) {
      case_when(
        MOI == 0 ~ paste0(CELL_LINE, "_", ARABINOSE, "Ara_-T6"),
        TRUE ~ paste0(CELL_LINE, "_", ARABINOSE, "Ara_+MOI_", MOI)
      )
    } else {
      case_when(
        T6_PHAGE == 0 ~ paste0(CELL_LINE, "_", ARABINOSE, "Ara_-T6"),
        TRUE ~ paste0(CELL_LINE, "_", ARABINOSE, "Ara_+10e-", T6_PHAGE)
      )
    },
    CONDITION_2 = case_when(
      ARABINOSE == 0 ~ paste0(CELL_LINE, "-Ara"),
      TRUE ~ paste0(CELL_LINE, "_", ARABINOSE, "%Ara")
    ),
    # Calculate OD600 normalized with positive values
    OD600_normalized = as.numeric(case_when(
      refined_OD600_mean - OD600_blank > 0 ~ refined_OD600_mean - OD600_blank,
      TRUE ~ 0
    )),
    Time_s = round(Time_s),
    Time_min = round(Time_s / 60),
    # Rename replicate numbers to `replicate_1`, `replicate_2`, etc.
    REPLICATE = paste0("replicate_", REPLICATE)
  ) %>%
  ungroup()


# Calculate PFU based on latest Phage Titer
current_PFU <- calculate_PFU(n_plaques           = n_plaques, 
                             dilution_step       = dilution_step, 
                             phage_volume_plated = phage_volume_plated) 
cat("Current PFU:", current_PFU, "/ mL\n")

# Calculate OD600 mean and SD for Cycle 01, accounting for either MOI or T6_PHAGE
OD600_C01_df <- dp %>%
  filter(Cycle == "01" & !grepl("BLANK", CONDITION_2)) %>%
  group_by(CONDITION_2) %>%
  mutate(
    OD600 = refined_OD600_mean,
    mean  = round(mean(refined_OD600_mean, na.rm = TRUE), 2),
    stdev = round(sd(refined_OD600_mean, na.rm = TRUE), 2)
  ) %>%
  ungroup() %>%
  mutate(T6_PHAGE = as.numeric(T6_PHAGE)) %>%
  group_by(Plate, Date, T6_PHAGE) %>%
  mutate(mean_OD600_per_T6 = mean(OD600, na.rm = TRUE)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    # Calculate MOI only if T6_PHAGE is used and MOI was not predefined
    MOI = if_else(use_moi, as.numeric(MOI),
                  calculate_MOI_lysis_assay(
                    OD600              = mean_OD600_per_T6, # Use mean OD600
                    conversion_factor  = 8E+08,             # E. coli conversion factor
                    T6_PHAGE           = T6_PHAGE,
                    PFU                = current_PFU,
                    bacterial_volume   = bacterial_volume,
                    phage_volume_added = phage_volume_added
                  ))
  ) %>%
  ungroup()

# Check for duplicates
dp_duplicates <- dp %>%
  filter(CELL_LINE != "BLANK") %>%
  group_by(Date, Plate, Cycle, CELL_LINE, T6_PHAGE, MOI, ARABINOSE, REPLICATE) %>%
  summarise(n = n()) %>%
  filter(n > 1)

# If duplicates need to be handled, resolve them with pivot_wider and values_fn
dp_replicates <- dp %>%
  filter(CELL_LINE != "BLANK") %>%
  select(Date, Plate, Cycle, REPLICATE, CELL_LINE, refined_OD600_mean, T6_PHAGE, ARABINOSE, MOI) %>%
  pivot_wider(
    names_from = REPLICATE, 
    values_from = refined_OD600_mean,
    values_fn = mean # Take the mean if duplicates exist
  ) %>%
  mutate(
    # Calculate mean OD600 across replicates
    mean_OD600 = rowMeans(select(., starts_with("replicate_")), na.rm = TRUE)
  ) %>%
  select(-starts_with("replicate_"))

# Final data preparation
dp_final <- dp %>%
  select(-REPLICATE) %>%
  unique() %>%
  distinct(CONDITION, Cycle, Date, .keep_all = TRUE) %>%
  filter(CELL_LINE != "BLANK") %>%
  mutate(T6_PHAGE = as.numeric(T6_PHAGE))

# Data for plotting, adjusting fields and filling missing values

if (use_moi) {
  data_for_plotting <- dp_final %>%
    mutate(MOI = as.numeric(MOI))
} else {
  data_for_plotting <- dp_final %>%
    select(-MOI)
}

data_for_plotting <- data_for_plotting %>%
  left_join(cl_key, relationship = "many-to-many") %>%
  left_join(OD600_C01_df %>% select(-REPLICATE)) %>%
  group_by(Date, Plate) %>%
  fill(OD600, mean, stdev, PURPOSE, Promoter, MOI, mean_OD600_per_T6, .direction = "downup") %>%
  ungroup() %>%
  unique()

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

fill_columns_by_id <- function(data, id, target_cols) {
  data <- data %>%
    group_by(!!sym(id)) %>%
    # Apply fill to each target column individually
    purrr::reduce(target_cols, ~ fill(.x, !!sym(.y), .direction = "downup"), .init = .) %>%
    ungroup()
  return(data)
}

fill_subset <- function(data, id, cols_to_join, new_col_name) {
  temp_data <- data %>%
    select(!!sym(id), all_of(cols_to_join)) %>%
    pivot_longer(names_to = "Temp", values_to = new_col_name, cols = all_of(cols_to_join)) %>%
    arrange(!!sym(id)) %>%
    group_by(!!sym(id)) %>%
    fill(everything(), .direction = "downup") %>%
    ungroup() %>%
    select(!!sym(id), !!sym(new_col_name)) %>%
    distinct(!!sym(id), .keep_all = TRUE)
  
  return(temp_data)
}

join_subset_to_main <- function(data, temp_data, id, cols_to_join) {
  data <- data %>%
    full_join(temp_data, by = id) %>%
    select(-all_of(cols_to_join)) %>%
    arrange(!!sym(id)) %>%
    group_by(!!sym(id)) %>%
    fill(everything(), .direction = "downup") %>%
    ungroup() %>%
    distinct(!!sym(id), .keep_all = TRUE)
  
  return(data)
}

join_columns <- function(data, id, cols_to_join, new_col_name) {
  temp_data <- fill_subset(data, id, cols_to_join, new_col_name)
  data <- join_subset_to_main(data, temp_data, id, cols_to_join)
  
  return(data)
}

ELISA_Fx <- function(Input_Directory, Output_Directory = Input_Directory) {
  # Initialize an empty data frame to store all plates' data
  All_plates_data = data.frame()
  
  # Get a list of Excel files in the main directory matching the pattern
  excel_files <- list.files(Input_Directory, recursive = T, full.names = TRUE, pattern = "\\d{8}_Plate_\\d+\\.xlsx$")
  
  # Check if there are plates found
  if (length(excel_files) > 0) {
    print("Plates exist!")
    # Iterate through each Excel file
    # file = "SOTA/01_raw_data/ELISA_PLATES/20220623_Plate_1.xlsx"
    for (file in excel_files) {
      
      print(paste("Processing", file))
      
      # Define patterns for reading sheets in files
      MEASUREMENTS_PATTERN <- c("MEASURE", "VALUE")
      CELL_LINES_PATTERN   <- c("CELL", "LINE", "COHORT")
      CONDITIONS_PATTERN   <- c("COND")
      DILUTIONS_PATTERN    <- c("DIL")
      STIM_DAYS_PATTERN    <- c("DAY")
      STIM_TIMES_PATTERN   <- c("TIME")
      PATHWAYS_PATTERN     <- c("PATHWAY")
      STIMULANTS_PATTERN   <- c("STIMULANT")
      STIM_CONCENTRATIONS_PATTERN <- c("CONC")
      
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
      
      # Initialize lists for storing results
      MEASUREMENTS <- list()
      CELL_LINES   <- list()
      CONDITIONS   <- list()
      DILUTIONS    <- list()
      STIM_DAYS    <- list()
      STIM_TIMES   <- list()
      PATHWAYS     <- list()
      STIMULANTS   <- list()
      STIM_CONCENTRATIONS <- list()
      
      # Extract sheets based on patterns
      MEASUREMENTS <- c(MEASUREMENTS, read_matching_sheets(file, patterns = MEASUREMENTS_PATTERN))
      CELL_LINES   <- c(CELL_LINES,   read_matching_sheets(file, patterns = CELL_LINES_PATTERN))
      CONDITIONS   <- c(CONDITIONS,   read_matching_sheets(file, patterns = CONDITIONS_PATTERN))
      DILUTIONS    <- c(DILUTIONS,    read_matching_sheets(file, patterns = DILUTIONS_PATTERN))
      STIM_DAYS    <- c(STIM_DAYS,    read_matching_sheets(file, patterns = STIM_DAYS_PATTERN))
      STIM_TIMES   <- c(STIM_TIMES,   read_matching_sheets(file, patterns = STIM_TIMES_PATTERN))
      PATHWAYS     <- c(PATHWAYS,     read_matching_sheets(file, patterns = PATHWAYS_PATTERN))
      STIMULANTS   <- c(STIMULANTS,   read_matching_sheets(file, patterns = STIMULANTS_PATTERN))
      STIM_CONCENTRATIONS <- c(STIM_CONCENTRATIONS, read_matching_sheets(file, patterns = STIM_CONCENTRATIONS_PATTERN))
      
      # Setting defaults
      default_dilution <- 5           # 1:5 dilution
      default_stim_time <- 24         # 24 hours
      default_stim_concentration <- 5 # 5ng/µL for IL-1ß stimulation
      
      DILUTIONS           <- lapply(DILUTIONS,  function(dilutions)     if (is.null(dilutions)) default_dilution else dilutions)
      STIM_TIMES          <- lapply(STIM_TIMES, function(time)          if (is.null(time)) default_stim_time else time)
      STIM_CONCENTRATIONS <- lapply(STIM_CONCENTRATIONS, function(conc) if (is.null(conc)) default_stim_concentration else conc)
      
      # Create Plate data.table
      Plate <- data.table(
        MEASUREMENT = unlist(MEASUREMENTS),
        CELL_LINE   = unlist(CELL_LINES),
        CONDITION   = unlist(CONDITIONS),
        DILUTION    = as.numeric(unlist(DILUTIONS)),
        STIM_DAY    = as.numeric(unlist(STIM_DAYS)),
        STIM_TIME   = as.numeric(unlist(STIM_TIMES)),
        PATHWAY     = unlist(PATHWAYS),
        STIMULANT   = unlist(STIMULANTS),
        STIM_CONCENTRATION = as.numeric(unlist(STIM_CONCENTRATIONS))
      )
      
      # Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINE != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate[Plate$CONDITION == "CALIBRATION"]
      Plate_Standards$CELL_LINE <- as.numeric(Plate_Standards$CELL_LINE)
      
      Plate_Standards <- Plate_Standards %>%
        group_by(CELL_LINE) %>%
        summarise(MEASUREMENT_mean = mean(MEASUREMENT)) %>%
        mutate(CELL_LINE = as.numeric(CELL_LINE),
               Date  = as_date(str_extract(basename(file), "\\d{8}"))) %>%
        arrange(CELL_LINE)
      
      # Fit the linear model
      Fit <- lm(CELL_LINE ~ MEASUREMENT_mean - 1, data = Plate_Standards[Plate_Standards$MEASUREMENT_mean <= 1.1, ])
      
      R       <- summary(Fit)$r.squared
      Rsquare <- signif(R, digits = 4)
      
      print(paste0("Secretion = slope*Intensity"))
      print(paste0("Secretion = ", Fit$coefficients[1],"*Intensity"))
      Plate_Standards <- Plate_Standards %>% mutate(Fit_Test = (Fit$coefficients[1]*MEASUREMENT_mean))
      
      # Plotting Standard Curve
      p <- ggplot(data = Plate_Standards) +
        geom_point(aes(x = MEASUREMENT_mean, y = CELL_LINE, col = MEASUREMENT_mean >  1.1), size = 5) +
        geom_line(aes(x = MEASUREMENT_mean, y = Fit_Test), linetype = "dashed") +
        annotate('text', x = 0.15, y = 700, label = paste0("R^2 = ", Rsquare), size = 10) +
        annotate('text',
                 x = max(Plate_Standards$MEASUREMENT_mean) - (0.25 * max(Plate_Standards$MEASUREMENT_mean)),
                 y = 150, label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        labs(x = "Measured Values",
             y = "IL-Concentration (pg/mL)") +
        ggtitle(label = paste0(basename(file)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        scale_color_manual(values = c("#79d2a3", "salmon"), guide = FALSE) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text  = element_text(size = 20)) +
        theme(legend.position = "none")
      
      # Saving the plot
      Save_Name <- file.path(Output_Directory, paste0(basename(file), "_Standard_Curve.pdf"))
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      # Fitting Data To Standard Curve ----------------------------------------
      # tmp_Plate <- Plate
      # Plate <- tmp_Plate
      Plate <- Plate %>%
        filter(CONDITION != "CALIBRATION") %>%
        mutate(
          Plate = as.numeric(unlist(lapply(strsplit(gsub(".xlsx", "", x = basename(file)),   "_", fixed=TRUE), function(x) return(x[3])))),
          Date  = as_date(str_extract(basename(file), "\\d{8}")),
          MEASUREMENT = as.numeric(MEASUREMENT),
          
          # Adjust measurements and concentrations
          Concentration = (Fit$coefficients[1] * (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT))),
          Concentration = Concentration * DILUTION,
          
          Is_Dose_Response = ifelse(str_detect(basename(file), "^DR_"), TRUE, FALSE)
        )
      
      All_plates_data <- rbind(All_plates_data, Plate)
    }
  } else {
    print("No plates found!")
  }
  return(All_plates_data)
}

process_data <- function(plate_data, pattern, negative_ctrl, positive_ctrl) {
  filtered_dates <- plate_data %>%
    filter(grepl(x = CL_NAME_ON_PLOT, pattern)) %>%
    pull(Date) %>%
    unique()
  
  filtered_data <- plate_data %>%
    filter(grepl(x = CL_NAME_ON_PLOT, pattern) | (PURPOSE == "CONTROL" & Date %in% filtered_dates))
  
  processed_data <- process_ELISA_data(DF = filtered_data, NEGATIVE_CTRL = negative_ctrl, POSITIVE_CTRL = positive_ctrl)
  
  return(processed_data)
}

process_plate_data <- function(data, date, plate_condition, id, df_name, negative_ctrl, positive_ctrl) {
  filtered_data <- data %>% filter(Date == date & {{ plate_condition }})
  processed_data <- process_data(filtered_data, id, negative_ctrl, positive_ctrl) %>% 
    mutate(ID = id, DF = df_name, NEGATIVE_CTRL = negative_ctrl, POSITIVE_CTRL = positive_ctrl) %>% 
    filter(PURPOSE == "SAMPLE")
  rm(filtered_data)
  return(processed_data)
}

get_species <- function(taxid) {
  classification_data <- classification(taxid, db = "ncbi")
  
  # Check if classification data is retrieved successfully
  if (is.null(classification_data[[1]])) {
    return(NA)
  }
  
  # Extract the species information
  species_info <- classification_data[[1]]
  
  # Find the species level
  species <- species_info[species_info$rank == "species", "name"]
  
  # Return species name (if exists) or NA
  if (length(species) == 1) {
    return(species)
  } else {
    return(NA)
  }
}

# Define the updated function to retrieve full taxonomy info
get_taxonomy_info <- function(taxid) {
  classification_data <- classification(taxid, db = "ncbi")
  
  # Check if classification data is retrieved successfully
  if (is.null(classification_data[[1]])) {
    return(list(species = NA, genus = NA, family = NA, order = NA, 
                class = NA, phylum = NA, clade = NA, superkingdom = NA))
  }
  
  # Extract the species and other taxonomic information
  tax_info <- classification_data[[1]]
  
  # Extract different ranks
  species <- tax_info[tax_info$rank == "species", "name"]
  genus <- tax_info[tax_info$rank == "genus", "name"]
  family <- tax_info[tax_info$rank == "family", "name"]
  order <- tax_info[tax_info$rank == "order", "name"]
  class <- tax_info[tax_info$rank == "class", "name"]
  phylum <- tax_info[tax_info$rank == "phylum", "name"]
  clade <- tax_info[tax_info$rank == "clade", "name"]
  superkingdom <- tax_info[tax_info$rank == "superkingdom", "name"]
  
  # Create a named list to return
  result <- list(
    species = ifelse(length(species) == 1, species, NA),
    genus = ifelse(length(genus) == 1, genus, NA),
    family = ifelse(length(family) == 1, family, NA),
    order = ifelse(length(order) == 1, order, NA),
    class = ifelse(length(class) == 1, class, NA),
    phylum = ifelse(length(phylum) == 1, phylum, NA),
    clade = ifelse(length(clade) == 1, clade, NA),
    superkingdom = ifelse(length(superkingdom) == 1, superkingdom, NA)
  )
  
  return(result)
}

# Function to get taxonomy and convert to data frame
get_taxonomy_df <- function(taxid) {
  tax_info <- get_taxonomy_info(taxid)
  return(as.data.frame(t(tax_info), stringsAsFactors = FALSE))
}

safe_get_taxonomy_df <- purrr::safely(get_taxonomy_df)

get_name_from_ncbi <- function(accession_number) {
  # Use the rentrez function to search NCBI for the accession number
  search_result <- tryCatch({
    entrez_search(db = "protein", term = accession_number, retmode = "xml", use_history = TRUE)
  }, error = function(e) {
    return(NA)
  })
  
  # Check if the search_result is empty or NULL
  if (is.null(search_result) || length(search_result$ids) == 0) {
    return(NA)
  }
  
  # Fetch the summary of the first search result
  summary_result <- tryCatch({
    entrez_summary(db = "protein", id = search_result$ids[1])
  }, error = function(e) {
    return(NA)
  })
  
  # Check if the summary_result is NULL
  if (is.null(summary_result)) {
    return(NA)
  }
  
  # Extract the gene name from the summary result
  gene_name <- summary_result$title
  
  return(gene_name)
}

# Define the function to write a DataFrame to a FASTA file
write_fasta <- function(df, id_col, seq_col, file_path) {
  # Ensure the data frame contains the necessary columns
  if (!all(c(id_col, seq_col) %in% colnames(df))) {
    stop("The DataFrame must contain the specified ID and sequence columns.")
  }
  
  # Open a connection to the file
  file_conn <- file(file_path, "w")
  
  # Iterate over each row in the DataFrame
  for (i in 1:nrow(df)) {
    # Write the ID as the FASTA header
    writeLines(paste0(">", df[[id_col]][i]), file_conn)
    # Write the sequence
    writeLines(df[[seq_col]][i], file_conn)
  }
  
  # Close the file connection
  close(file_conn)
}

################################################################################
################################################################################
################################################################################

sem <- function(x) sd(x)/sqrt(length(x))

################################################################################
################################################################################
################################################################################

# DF = CHARMS
# NEGATIVE_CTRL = "3xKO"
# POSITIVE_CTRL = "Wild Type"

process_ELISA_data <- function(DF, NEGATIVE_CTRL, POSITIVE_CTRL) {
  
  group_vars <- c("STIM_DAY", "Date")
  
  # Get baseline
  get_baseline <- function(DF, NEGATIVE_CTRL) {
    if (any((DF$CELL_LINE %in% NEGATIVE_CTRL & DF$CONDITION == "UNSTIM") | unique(DF$CL_NAME_ON_PLOT %in% NEGATIVE_CTRL & DF$CONDITION == "UNSTIM"))) {
      baseline <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        filter((CELL_LINE %in% NEGATIVE_CTRL & CONDITION == "UNSTIM") | (CL_NAME_ON_PLOT %in% NEGATIVE_CTRL & CONDITION == "UNSTIM")) %>%
        summarise(baseline_control_value = mean(Concentration))
    } else if (any(DF$CELL_LINE %in% NEGATIVE_CTRL)) {
      baseline <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(baseline_control_value = min(Concentration[CONDITION == "UNSTIM"]))
    } else {
      baseline <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(baseline_control_value = min(Concentration))
    }
    
    # Join the calculated values with the dataset
    DF_baseline_adj <- left_join(DF, baseline) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value, TRUE ~ Concentration))
    
    return(DF_baseline_adj)
  }
  
  DF_baseline_adj <- get_baseline(DF = DF, NEGATIVE_CTRL = NEGATIVE_CTRL)
  
  get_normalization_value <- function(DF, POSITIVE_CTRL) {
    
    if (any(DF$CELL_LINE %in% POSITIVE_CTRL & DF$CONDITION %in% "STIM" | DF$CL_NAME_ON_PLOT %in% POSITIVE_CTRL & DF$CONDITION %in% "STIM")) {
      normalization_control_value <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        filter((CELL_LINE %in% POSITIVE_CTRL & CONDITION == "STIM") | (CL_NAME_ON_PLOT %in% POSITIVE_CTRL & CONDITION == "STIM")) %>%
        summarise(normalization_control_value = case_when(mean(Concentration_REDUCED) > 0 ~ mean(Concentration_REDUCED), TRUE ~ -Inf))
    } else {
      normalization_control_value <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(normalization_control_value = max(Concentration))
    }
    
    # Join the calculated control means
    DF_normalization_adj <- left_join(DF, normalization_control_value)
    
    return(DF_normalization_adj)
  }
  
  DF_normalization_adj <- get_normalization_value(DF = DF_baseline_adj, POSITIVE_CTRL = POSITIVE_CTRL)
  
  
  # Normalize ELISA data
  normalize_ELISA <- function(DF) {
    
    DATA_NORMALIZED <- DF %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / normalization_control_value < 0 ~ 0, TRUE ~ Concentration_REDUCED / normalization_control_value),
             triplicate_mean_per_day  = mean(Concentration_NORMALIZED)) %>%
      ungroup()
    return(DATA_NORMALIZED)
  }
  
  DATA_NORMALIZED <- DF_normalization_adj %>%
    group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
    mutate(Concentration_NORMALIZED = ifelse(Concentration_REDUCED / normalization_control_value < 0, 0, Concentration_REDUCED / normalization_control_value),
           triplicate_mean_per_day  = mean(Concentration_NORMALIZED),
           POSITIVE_CTRL = POSITIVE_CTRL,
           NEGATIVE_CTRL = NEGATIVE_CTRL) %>%
    ungroup()
  
  return(DATA_NORMALIZED)
}

################################################################################
################################################################################
################################################################################

perform_statistical_analysis <- function(DATA, GROUP_BY_COLUMN, TESTING_COLUMN) {
  # Internal function for pairwise t-test
  unpaired_ttest <- function(DATA, return_annotation = FALSE) {
    
    if (nrow(DATA) < 3) {
      if (return_annotation) {
        return("")
      } else {
        return("")
      }
    }
    
    # unpaired t-test
    # p_values <- t.test(DATA[[TESTING_COLUMN]]~ DATA$CONDITION)$p.value
    formula  <- as.formula(paste(TESTING_COLUMN, "~ CONDITION"))
    p_values <- ggpubr::compare_means(formula,
                                      data = DATA,
                                      method = "t.test", 
                                      paired = FALSE)$p.adj
    
    if (return_annotation) {
      p_annotation <- case_when(
        p_values < 1e-4 ~ '****',
        p_values < 1e-3 ~ '***',
        p_values < 1e-2 ~ '**',
        p_values < 0.05 ~ '*',
        TRUE ~ 'ns'
      )
      return(p_annotation)
    } else {
      return(formatC(p_values, format = "e", digits = 3))
    }
  }
  
  # Calculate statistical significance using a t-test for each group
  annotations <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), unpaired_ttest)
  p_values    <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), function(DATA) unpaired_ttest(DATA, return_annotation = TRUE))
  
  filtered_annotations <- unlist(lapply(annotations, function(annotation) annotation[annotation != ""]))
  filtered_p_values    <- unlist(lapply(p_values, function(p_value) p_value[p_value != ""]))
  
  # return(list(annotations = annotations, p_values = p_values))
  return(list(annotations = filtered_annotations, p_values = filtered_p_values))
}

################################################################################
################################################################################
################################################################################

process_data_for_plot <- function(data, change_unstim_plt_col = T, unstim_plt_col = "#BEBDBD", unstim_plt_col_lightest = F) {
  # Reorder cell lines for plotting
  data$CL_NAME_ON_PLOT <- reorder(data$CL_NAME_ON_PLOT, -data$ORDER_NO)
  
  # Reformat condition for legend text
  data$CONDITION <- factor(data$CONDITION, levels = c("UNSTIM", "STIM"))
  if (unstim_plt_col_lightest) {
    data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- data$PLT_LIGHTEST[data$CONDITION == "UNSTIM"] 
  } else if (change_unstim_plt_col) {
    data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- unstim_plt_col 
  }
  
  return(data)
}

################################################################################
################################################################################
################################################################################

prepare_plotting_means <- function(data, group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "POSITIVE_CTRL", "NEGATIVE_CTRL")) {
  # Group and summarize data for plotting_means
  plotting_means <- data %>%
    group_by(!!!syms(group_var)) %>%
    summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
              Relative_Intensity_mean = mean(Concentration_NORMALIZED)) %>%
    as.data.table()
  
  # Round Relative_Intensity_mean
  plotting_means$Relative_Intensity_mean <- round(plotting_means$Relative_Intensity_mean, 3)
  
  # Reorder CL_NAME_ON_PLOT
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  
  return(plotting_means)
}

################################################################################
################################################################################
################################################################################

process_statistical_analysis <- function(data, group_var, value_var) {
  # Perform statistical analysis
  statistical_significance <- perform_statistical_analysis(data, group_var, value_var)
  
  # Turn statistical_significance list into data.table
  stat_significance_dt <- data.table(
    CL_NAME_ON_PLOT = names(statistical_significance$annotations),
    p_value = statistical_significance$annotations,
    significance = statistical_significance$p_values
  )
  return(stat_significance_dt)
}

################################################################################
################################################################################
################################################################################

# data = plotting_data_main
# group_var = c("CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO")
# mean_var = "triplicate_mean_per_day"
# change_unstim_plt_col = T
# unstim_plt_col = "#BEBDBD"
# unstim_plt_col_lightest = F
# fold_change_option = F

prepare_plotting_stats <- function(data, stat_significance_dt, 
                                   group_var = c("CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", 
                                                 "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", 
                                                 "PLOTTING_COLOR", "ORDER_NO",
                                                 "POSITIVE_CTRL", "NEGATIVE_CTRL"),
                                   mean_var = "triplicate_mean_per_day",
                                   change_unstim_plt_col = T,
                                   unstim_plt_col = "#BEBDBD",
                                   unstim_plt_col_lightest = F,
                                   fold_change_option = F) {
  
  # Determine grouping variables based on the condition
  if (unstim_plt_col_lightest) {
    group_vars <- c(group_var, "PLT_LIGHTEST")
    # group_vars <- group_var
  } else {
    group_vars <- group_var
  }
  
  # Group and summarize data for plotting_stats
  plotting_stats_main <- data %>%
    group_by(!!!syms(group_vars)) %>%
    summarise(
      IL2_concentration_Dilution_Factor_mean = mean(Concentration),
      IL2_concentration_Dilution_Factor_sem = sem(Concentration),
      Relative_Intensity_mean = mean(!!!syms(mean_var)),
      Relative_Intensity_sem = sem(!!!syms(mean_var))
    ) %>%
    as.data.table() %>%
    left_join(stat_significance_dt)
  
  plotting_stats_main$CONDITION <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
  plotting_stats_main$CL_NAME_ON_PLOT <- reorder(plotting_stats_main$CL_NAME_ON_PLOT, -plotting_stats_main$ORDER_NO)
  
  if (unstim_plt_col_lightest) {
    plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- plotting_stats_main$PLT_LIGHTEST[plotting_stats_main$CONDITION == "UNSTIM"]
    # plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- unstim_plt_col
  } else if (change_unstim_plt_col) {
    plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- unstim_plt_col
  }
  
  return(plotting_stats_main)
}

################################################################################
################################################################################
################################################################################

process_alignments <- function(file_path, 
                               output_dir, 
                               key, 
                               min.fraction=0.5, 
                               min.block.width=2, 
                               OVERWRITE=TRUE,
                               base = paste0(temp_dir, "/SEQS_AND_ALNS/", basename(file_path) %>% str_replace_all("\\.fasta", ""))){
  ## load the multiple alignment
  origMAlign <- Biostrings::readAAMultipleAlignment(file_path)
  
  ## list the names of the sequences in the alignment
  rownames(origMAlign)
  
  ## rename the sequences to the species names
  species_nms <- c(key[match(rownames(origMAlign), key$ID), "species"])
  rownames(origMAlign) <- species_nms[["species"]]
  
  # ## automatically masks (removes or marks as unreliable) columns in the alignment that have a high proportion of gaps.
  # ## Improves the reliability of subsequent analyses by excluding poorly aligned regions.
  autoMasked <- maskGaps(origMAlign, min.fraction = min.fraction, min.block.width = min.block.width)
  colmask(autoMasked)
  
  ## Calculates the frequency of each amino acid in the masked alignment.
  ## Provides insights into the composition of the sequences
  alph_freq <- Biostrings::alphabetFrequency(autoMasked)
  
  ## Provides insights into the consensus sequence across the alignment.
  cons_mat  <- Biostrings::consensusMatrix(autoMasked, baseOnly=TRUE)[, 1:nrow(autoMasked)]
  
  ## Generates a consensus sequence from the alignment
  ## "#" represents a masked column
  ## "-" represents a gap
  ## "?" represents an ambiguous position
  cons_string <- consensusString(autoMasked)
  
  # Summarizes the most common amino acids at each position, highlighting conserved regions.
  cons_view   <- consensusViews(autoMasked)
  
  # Cluster the sequences based on their similarity using hierarchical clustering.
  # Identifies groups of similar sequences, which can reveal evolutionary relationships or functional similarities.
  sdist      <- stringDist(as(autoMasked,"AAStringSet"), method="hamming")
  clust      <- hclust(sdist, method = "single")
  
  # cut a hierarchical clustering tree (dendrogram) into a specified number of groups (clusters)
  fourgroups <- cutree(clust, k = 4)
  
  if (OVERWRITE) {
    ## write out the results to CSV files
    # base <- paste0(temp_dir, "/SEQS_AND_ALNS/", basename(file_path) %>% str_replace_all("\\.fasta", ""))
    
    fwrite(as.data.frame(alph_freq),   paste0(base, "_alph_freq.csv"))
    fwrite(as.data.frame(cons_mat),    paste0(base, "_cons_mat.csv"))
    fwrite(as.data.frame(cons_string), paste0(base, "_cons_string.csv"))
    fwrite(as.data.frame(cons_view),   paste0(base, "_cons_view.csv"))
    fwrite(as.matrix(sdist),           paste0(base, "_sdist.csv"))
    png(filename = paste0(base, "_clust.png"), width = 1600, height = 1200); plot(clust <- hclust(sdist, method = "single")); dev.off()
    
  }
  ## write out the alignement object (with current masks) to Phylip format
  # write.phylip(x = autoMasked, filepath = file.path(output_dir, "DLD_sequences_aln.phy"))
}

################################################################################

iterative_msa_masking <- function(file_path, 
                                  output_dir, 
                                  key, 
                                  min_fraction_start_end, 
                                  min_block_width_start_end, 
                                  OVERWRITE = TRUE) {
  
  for (min_fraction in min_fraction_start_end) {
    for (min_block_width in min_block_width_start_end) {
      # Construct the base variable inside the loops
      base <- paste0(output_dir, "/SEQS_AND_ALNS/iterative_msa_masking/", 
                     basename(file_path) %>% str_replace_all("\\.fasta", ""), 
                     "_min_fraction_", min_fraction, 
                     "_min_block_width_", min_block_width)
      
      process_alignments(file_path, 
                         output_dir, 
                         key, 
                         min.fraction = min_fraction, 
                         min.block.width = min_block_width,
                         OVERWRITE,
                         base)
    }
  }
}

################################################################################

consensus_char <- function(column, min_pct, max_pct) {
  # Remove "-" and count frequencies
  column <- column[column != ""]
  if(length(column) == 0) {
    return("-")
  }
  freq_table <- table(column[column != "-"])
  total_non_gap <- sum(freq_table)
  
  # If the column has fewer non-gap characters than the min percentage, return "-"
  if (total_non_gap < length(column) * min_pct) {
    return("-")
  }
  
  # Find the most common character
  if (length(freq_table) == 0) {
    return("-")
  } else {
    consensus <- names(freq_table)[which.max(freq_table)]
  }
  
  # If the most common character appears in more than max_pct of rows, return it, otherwise return "-"
  if (freq_table[consensus] >= length(column) * max_pct) {
    return(consensus)
  } else {
    return("-")
  }
}

################################################################################
# plotting_data = plotting_data; name_key = NAME_KEY
prepare_fold_change_plots <- function(plotting_data, name_key) {
  
  # Select Data of interest
  cohort_data <- plotting_data
  
  # Set negative measured values to zero
  cohort_data$MEASUREMENT <- ifelse(cohort_data$MEASUREMENT < 0, 0, cohort_data$MEASUREMENT)
  
  # Normalization for fold change from UNSTIM to STIM
  normalized_to_control <- cohort_data %>%
    group_by(Date, STIM_DAY, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PLOTTING_COLOR, CL_NUMBER) %>%
    summarise(mean_per_day = mean(MEASUREMENT),
              Concentration = mean(Concentration),
              .groups = 'drop')
  
  # First, calculate the mean for the UNSTIM condition separately.
  unstim_means <- cohort_data %>%
    filter(CONDITION == "UNSTIM") %>%
    group_by(CELL_LINE, CL_NUMBER) %>%
    summarise(mean_unstim = ifelse(mean(MEASUREMENT) == 0, 1, mean(MEASUREMENT)),
              .groups = 'drop')
  
  # Now, join this back to the main dataset.
  cohort_subset <- cohort_data %>% left_join(unstim_means, by = c("CELL_LINE", "CL_NUMBER"), relationship = "many-to-many")
  
  # Compute the fold change.
  normalized_to_control <- cohort_subset %>%
    mutate(fold_change = case_when(CONDITION == "STIM" ~ MEASUREMENT / mean_unstim, 
                                   TRUE ~ NA_real_)) %>%
    group_by(CELL_LINE, CONDITION, PLOTTING_COLOR) %>%
    mutate(trip_mean = mean(fold_change, na.rm = TRUE),
           fold_change_sd = sd(fold_change, na.rm = TRUE)) %>%
    unique() %>%
    ungroup()
  
  # Helper function to perform a t-test and return p-value and annotation
  # perform_ttest <- function(data) {
  #   ttest_result <- t.test(data$MEASUREMENT[data$CONDITION == "STIM"],
  #                          data$MEASUREMENT[data$CONDITION == "UNSTIM"],
  #                          paired = FALSE)
  #   p_value <- ttest_result$p.value
  #   annotation <- ifelse(p_value < 0.0001, '****', 
  #                        ifelse(p_value < 0.001, '***', 
  #                               ifelse(p_value < 0.01, '**',
  #                                      ifelse(p_value < 0.05, '*', 'ns'))))
  #   return(list(p_value = p_value, annotation = annotation))
  # }
  
  # Helper function to perform a t-test and return p-value and annotation
  perform_ttest <- function(data) {
    ttest_result <- t.test(,
                           paired = FALSE)
    p_value <- ttest_result$p.value
    annotation <- ifelse(p_value < 0.0001, '****',
                         ifelse(p_value < 0.001, '***',
                                ifelse(p_value < 0.01, '**',
                                       ifelse(p_value < 0.05, '*', 'ns'))))
    return(list(p_value = p_value, annotation = annotation))
  }
  
  control_subset_fc <- cohort_subset %>% 
    group_by(CL_NAME_ON_PLOT, Date, CL_NUMBER) %>%
    filter(CELL_LINE == "cl204") %>% 
    summarise(fold_change_neg_ctrl = round(mean(MEASUREMENT[CONDITION == "UNSTIM"])))

  cohort_subset  <- cohort_subset %>% filter(CELL_LINE != "cl204") %>% 
    left_join(control_subset_fc)
  
  results <- cohort_subset %>% group_by(CL_NAME_ON_PLOT, Date) %>%
    do({
      data <- .
      fold_change_data <- data %>%
        filter(MEASUREMENT != 0) %>%
        summarise(fold_change = round(mean(MEASUREMENT[CONDITION == "STIM"]) / mean(MEASUREMENT[CONDITION == "UNSTIM"]), digits = 2), .groups = 'drop')
      #ttest_results <- perform_ttest(data)
      #cbind(fold_change_data, ttest_results)
    }) %>%
    ungroup() %>%
    left_join(name_key[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR", "CL_NUMBER")], 
              by = c("CL_NAME_ON_PLOT" = "CL_NAME_ON_PLOT"), relationship = "many-to-many") %>%
    unique()
  
  # Extract annotations for plotting
  # annotations <- results$annotation
  # names(annotations) <- results$CL_NAME_ON_PLOT
  
  # Reorder by your preferred visualization order
  results$CL_NAME_ON_PLOT <- reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  
  return(list(results = results#,
              #annotations = annotations
              ))
}

################################################################################

# # extract_af3_json_files <- function(af_json_files) {
# #   # Initialize a list to store the large data frames for each af_json_file
# #   large_af3_dfs <- list()
# #   
# #   for (af_json_file in af_json_files) {
# #     
# #     # af_json_file = af_json_files[1] # for debugging
# #     json <- fromJSON(af_json_file)
# #     
# #     af3_name     <- json[[1]][[1]]
# #     af3_seed     <- as.numeric(json[[2]][[1]])
# #     af3_seq      <- json[[3]][[1]]$proteinChain$sequence
# #     af3_n_chains <- json[[3]][[1]]$proteinChain$count
# #     
# #     summary_confidences_files <- list.files(af_dir, pattern = "summary_confidences", full.names = TRUE) %>%
# #       str_subset(paste0("(?i)", gsub("-", "_", af3_name), "_")) %>%
# #       str_subset("\\.json$")
# #     
# #     # Initialize a list to store the extracted information from each summary_confidences file
# #     sum_conf_extract_list <- list()
# #     
# #     for (summary_confidences_file in summary_confidences_files) {
# #       
# #       # summary_confidences_file <- summary_confidences_files[1] # for debugging
# #       af3_model               <- basename(summary_confidences_file)
# #       
# #       json <- fromJSON(summary_confidences_file)
# #       
# #       af3_chain_iptm          <- json[["chain_iptm"]]
# #       af3_chain_pair_iptm     <- json[["chain_pair_iptm"]]
# #       af3_chain_pair_pae_min  <- json[["chain_pair_pae_min"]]
# #       af3_chain_ptm           <- json[["chain_ptm"]]
# #       af3_fraction_disordered <- json[["fraction_disordered"]]
# #       af3_has_clash           <- json[["has_clash"]]
# #       
# #       ##########################################################################
# #       
# #       af3_iptm                <- json[["iptm"]]
# #       af3_n_recycles          <- json[["num_recycles"]]
# #       af3_ptm                 <- json[["ptm"]]
# #       af3_ranking_score       <- json[["ranking_score"]]
# #       
# #       sum_conf_extract <- cbind(af3_name,
# #                                 af3_seed,
# #                                 af3_seq,
# #                                 af3_n_chains,
# #                                 af3_model,
# #                                 af3_chain_iptm,
# #                                 af3_chain_pair_pae_min,
# #                                 af3_chain_ptm,
# #                                 af3_fraction_disordered,
# #                                 af3_has_clash,
# #                                 af3_iptm,
# #                                 af3_n_recycles,
# #                                 af3_ptm,
# #                                 af3_ranking_score
# #       ) %>%  as.data.frame()
# #       
# #       # Add the data frame to the list
# #       sum_conf_extract_list <- append(sum_conf_extract_list, list(sum_conf_extract))
# #       
# #       # Combine the af_json_extract data frame with all the summary_confidences data frames
# #       large_af3_df <- bind_rows(sum_conf_extract_list)
# #       
# #       # Add the large data frame to the list
# #       large_af3_dfs <- append(large_af3_dfs, list(large_af3_df))
# #     }
# #     
# #     # Combine all large data frames into one
# #     mega_af3_df <- bind_rows(large_af3_dfs) %>%
# #       mutate(ID = str_extract(af3_name, "(?i)^(?:[^_]*_)?([A-Z]+[12]?DLD\\d?)_\\d{2}"),
# #              ID = str_replace(ID, "(DLD|BDLD[12])_(\\d{2})", "\\1_\\2"))
# #     
# #     # Subset the data to save a mini version with unique rows based on specific columns
# #     
# #     ### TODO: find a way to pivot the remaining columns into a wide format!
# #     mini_af3_df <- mega_af3_df %>%
# #       select(
# #         ID,
# #         af3_name,
# #         af3_seed,
# #         af3_seq,
# #         af3_n_chains,
# #         af3_model,
# #         af3_iptm,
# #         af3_n_recycles,
# #         af3_ptm,
# #         af3_ranking_score
# #       ) %>%
# #       distinct()
# #     
# #     
# #     AF3_SUMMARY <- mini_af3_df %>%
# #       group_by(af3_name) %>%
# #       filter(grepl("summary_confidences_0", af3_model))
# #   }
# #   return(list(mega_af3_df, mini_af3_df, AF3_SUMMARY))
# # }
# 
# 
# # af_json_files <- af_json_files[1:2] # for debugging
# # af_json_file <- af_json_files[1]    # for debugging
# # af_dir <- "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/01_raw_data/AF3_PREDICTIONS/OPERONS"
# # af_json_files <- list.files(af_dir, full.names = TRUE, recursive = TRUE, pattern = "job_request") %>% str_subset("(?i)dld") %>% str_subset("\\.json$")
# 
# # Function to extract information from job_request JSON file
# 
# extract_job_request_info <- function(af_json_file) {
#   
#   json <- fromJSON(af_json_file)
#   
#   af3_name  <- json[[1]][[1]]
#   af3_seed  <- as.numeric(json[[2]][[1]])
#   sequences <- json[[3]]
#   
#   # seq <- sequences[[1]]
#   chain_info <- lapply(sequences, function(seq) {
#     list(
#       sequence = seq$proteinChain[[1]],
#       count = seq$proteinChain$count
#     )
#   })
#   
#   chain_sequences <- sapply(chain_info, function(x) x$sequence)
#   af3_n_chains <- sapply(chain_info, function(x) x$count)
#   
#   return(list(
#     af3_name = af3_name,
#     af3_seed = af3_seed,
#     chain_sequences = chain_sequences,
#     # af3_n_chains = af3_n_chains
#     af3_n_chains = sum(af3_n_chains)  # Sum to get the total number of chains
#   ))
# }
# 
# # Function to process chain pair data
# process_chain_pair_data <- function(af3_n_chains) {
#   chain_combinations <- expand.grid(1:af3_n_chains, 1:af3_n_chains)
#   colnames(chain_combinations) <- c("Chain_i", "Chain_j")
#   return(chain_combinations)
# }
# 
# # summary_confidences_file <- summary_confidences_files[1] # for debugging
# # Function to extract information from summary_confidences JSON files
# extract_summary_confidences_info <- function(summary_confidences_file, af3_n_chains) {
#   
#   json <- fromJSON(summary_confidences_file)
#   
#   chain_pair_iptm <- json[["chain_pair_iptm"]]
#   chain_pair_pae_min <- json[["chain_pair_pae_min"]]
#   
#   chain_combinations <- process_chain_pair_data(af3_n_chains)
#   
#   chain_combinations <- chain_combinations %>%
#     mutate(
#       iptm = mapply(function(i, j) chain_pair_iptm[i, j], Chain_i, Chain_j),
#       pae_min = mapply(function(i, j) chain_pair_pae_min[i, j], Chain_i, Chain_j)
#     )
#   
#   return(chain_combinations)
# }
# 
# 
# job_request_info <- extract_job_request_info(af_json_file)
# af3_n_chains     <- job_request_info$af3_n_chains
# 
# summary_confidences_files <- list.files(af_dir, pattern = "summary_confidences", full.names = TRUE, recursive = T) %>%
#   str_subset(paste0("(?i)", gsub("-", "_", af3_name), "_")) %>%
#   str_subset("\\.json$")
# 
# # Read and combine all summary_confidences data
# all_chain_pair_data <- lapply(summary_confidences_files, function(summary_confidences_file) {
#   extract_summary_confidences_info(summary_confidences_file, af3_n_chains)
# })
# 
# # Combine all chain pair data into one data frame
# combined_chain_pair_data <- bind_rows(all_chain_pair_data)
# 
# # pivot the information from combined_chain_pair_data into a wide format
# wide_chain_pair_data <- combined_chain_pair_data %>%
#   pivot_wider(id_cols = c(Chain_i, Chain_j, iptm, pae_min), names_from = c(Chain_i, Chain_j, iptm, pae_min), values_from = c(iptm, pae_min))


################################################################################

# operons <- operon_df$operon[1:10]
# ids     <- operon_df$ID[1:10]
# operon = operons[2]               # for debugging
# ID     = ids[2]                   # for debugging

# Function to interpret directions in an operon string
interpret_operon_direction <- function(ID, operon) {

  operon <- gsub("<-", "←", operon)

  # Extract components with direction symbols using stringr::str_split
  # components <- unlist(str_split(operon, "(?<=[←→<\\|])")) ; components
  components <- unlist(str_split(operon, "(?<=→|←|\\|\\|)")) ; components

  # Initialize vectors for parts and directions
  parts <- character(length(components))
  directions <- character(length(components))

  # Temporary storage for the last direction sign encountered
  last_direction <- NA

  # Loop through components to separate ?|| into individual parts
  for (i in seq_along(components)) {
    if (components[i] == "?") {
      parts[i] <- "?"
      directions[i] <- last_direction  # Assign the last stored direction
    } else if (components[i] == "||") {
      parts[i] <- "||"
      directions[i] <- NA  # No direction associated with ||
    } else if (str_detect(components[i], "\\?\\|\\|")) {
      parts[i] <- "?"
      directions[i] <- last_direction  # Assign the last stored direction
    } else {
      parts[i] <- components[i]
      directions[i] <- str_extract(components[i], "←|→")
      if (!is.na(directions[i])) {
        last_direction <- directions[i]  # Update last stored direction
      }
    }
  }

  # Determine the direction for the last gene correctly
  if (!is.na(directions[length(directions)])) {
    if (directions[length(directions)] == "←" || str_ends(parts[length(parts)], "REV")) {
      directions[length(directions)] <- "REV"
    } else {
      directions[length(directions)] <- "FWD"
    }
  } else {
    directions[length(directions)] <- "REV"  # Default direction as "<-" for the last gene if its direction is NA
  }

  # Remove empty parts
  operon_data <- data.frame(ID = ID, part = parts, direction = directions, stringsAsFactors = FALSE)
  operon_data <- operon_data[!is.na(operon_data$part) & operon_data$part != "", ]

  # Replace direction symbols with "FWD" and "REV"
  operon_data$direction[operon_data$direction == "→"] <- "FWD"
  operon_data$direction[operon_data$direction == "←"] <- "REV"

  operon_data$part <- str_replace_all(operon_data$part, "→|←|\\|", "")

  # Return the interpreted operon data
  return(operon_data)
}

################################################################################

# Function to search for matches in search_col1 and search_col2 and return the target_col
search_target_by_pattern_col <- function(df, pattern, search_col1, search_col2, target_col) {
  # Search in search_col1 column
  match_col1 <- df %>%
    filter(grepl(pattern, .[[search_col1]])) %>%
    select(all_of(target_col))
  
  # If match found in search_col1, return the target_col name(s)
  if (nrow(match_col1) > 0) {
    return(match_col1[[target_col]])
  }
  
  # If no match found in search_col1, search in search_col2 column
  match_assembly_2 <- df %>%
    filter(grepl(pattern, .[[search_col2]])) %>%
    select(all_of(target_col))
  
  # Return the target_col name(s) found in search_col2 or NULL if not found
  if (nrow(match_assembly_2) > 0) {
    return(match_assembly_2[[target_col]])
  } else {
    return(NULL)  # Return NULL if no match is found in either column
  }
}

################################################################################

# gene = "WP_044290893.1"
# Function to split the gene name into character and numeric parts
split_gene_name <- function(gene) {
  parts <- regmatches(gene, regexec("^([A-Za-z_]+)([0-9]+)\\.([0-9]+)$", gene))
  if (length(parts[[1]]) > 0) {
    return(list(prefix = parts[[1]][2], number = parts[[1]][3], suffix = parts[[1]][4]))
  } else {
    return(NULL)
  }
}

# Function to get upstream genes
get_upstream_genes <- function(gene, n) {
  parts <- split_gene_name(gene)
  if (!is.null(parts)) {
    num_length <- nchar(parts$number)
    upstream_genes <- sapply((as.numeric(parts$number) - n):(as.numeric(parts$number) - 1), function(num) {
      sprintf("%s%0*d.%s", parts$prefix, num_length, num, parts$suffix)
    })
    return(upstream_genes)
  } else {
    return(NULL)
  }
}

# Function to get downstream genes
get_downstream_genes <- function(gene, n) {
  parts <- split_gene_name(gene)
  if (!is.null(parts)) {
    num_length <- nchar(parts$number)
    downstream_genes <- sapply((as.numeric(parts$number) + 1):(as.numeric(parts$number) + n), function(num) {
      sprintf("%s%0*d.%s", parts$prefix, num_length, num, parts$suffix)
    })
    return(downstream_genes)
  } else {
    return(NULL)
  }
}

################################################################################

extract_operon <- function(genome, start, end, strand) {
  # Calculate cumulative lengths of contigs
  cumulative_lengths <- cumsum(sapply(genome, length))
  
  # Find the contig that contains the start position
  contig_start <- which(start <= cumulative_lengths & start > cumulative_lengths - sapply(genome, length))
  # Find the contig that contains the end position
  contig_end <- which(end <= cumulative_lengths & end > cumulative_lengths - sapply(genome, length))
  
  if (length(contig_start) == 0 || length(contig_end) == 0) {
    stop("Start or end position is out of bounds.")
  }
  
  # If both positions are within the same contig
  if (contig_start == contig_end) {
    local_start <- start - (cumulative_lengths[contig_start] - length(genome[[contig_start]]))
    local_end   <- end   - (cumulative_lengths[contig_end]   - length(genome[[contig_end]]))
    operon_sequence <- subseq(genome[[contig_start]], start=local_start, end=local_end)
  } else {
    # If the operon spans multiple contigs
    operon_sequence <- DNAStringSet()
    
    # Extract the sequence from the start position to the end of the start contig
    local_start <- start - (cumulative_lengths[contig_start] - length(genome[[contig_start]]))
    operon_sequence <- c(operon_sequence, DNAStringSet(subseq(genome[[contig_start]], start=local_start, end=length(genome[[contig_start]]))))
    
    # Extract the sequences from any intermediate contigs
    if (contig_start + 1 < contig_end) {
      for (i in (contig_start + 1):(contig_end - 1)) {
        operon_sequence <- c(operon_sequence, DNAStringSet(genome[[i]]))
      }
    }
    
    # Extract the sequence from the beginning of the end contig to the end position
    local_end <- end - (cumulative_lengths[contig_end] - length(genome[[contig_end]]))
    operon_sequence <- c(operon_sequence, DNAStringSet(subseq(genome[[contig_end]], start=1, end=local_end)))
    
    # Concatenate all the sequences into a single DNAString
    operon_sequence <- unlist(operon_sequence)
  }
  
  # Reverse complement if the strand is negative
  if (strand == "-") {
    operon_sequence <- reverseComplement(operon_sequence)
  }
  
  return(operon_sequence)
}

################################################################################

# Function to find all genes in the operon starting from a given gene
find_operon_genes <- function(start_gene, operonic_gp_list) {
  operon_genes <- c(start_gene)  # Initialize with the start gene
  new_genes <- c(start_gene)     # Genes to be checked
  
  while (length(new_genes) > 0) {
    current_genes <- new_genes
    new_genes <- c()
    
    # Find all pairs where current_genes are involved
    for (gene in current_genes) {
      new_pairs <- operonic_gp_list %>%
        filter((IdGene %in% gene | IdGene2 %in% gene) & CLASS == "Operon" & Probability >= 0.5)
      
      # Add new genes to operon_genes
      new_genes <- c(new_genes, new_pairs$IdGene, new_pairs$IdGene2)
    }
    
    # Keep only unique genes that are not already in operon_genes
    new_genes <- setdiff(unique(new_genes), operon_genes)
    operon_genes <- unique(c(operon_genes, new_genes))
  }
  
  return(operon_genes)
}

################################################################################

# Define a function to convert DNAStringSet to a character vector
convert_to_character <- function(dna_string_set) {
  if (length(dna_string_set) > 0) {
    as.character(dna_string_set)
  } else {
    NULL
  }
}

################################################################################

# Function to read files and add ID and db columns
read_and_add_ID_and_db <- function(folder_path, file_names, db_names) {
  ID <- unlist(lapply(strsplit(basename(folder_path), "_mmseqs", fixed=TRUE), function(x) return(x[1])))
  file_paths <- file.path(folder_path, file_names)
  data_list <- lapply(seq_along(file_paths), function(i) {
    data <- fread(file_paths[i], sep = "\t")
    data$ID <- ID
    data$db <- db_names[i]
    return(data)
  })
  return(data_list)
}

################################################################################

extract_operon_components_from_list <- function(main_df = sota, genome, list, output_path = file.path(output_dir, "OPERONS", genome_name)) {
  output_path <- output_path
  if (!dir.exists(output_path)) {dir.create(output_path, showWarnings = FALSE)}
  
  genome_path           <- genome_fasta_files[sapply(genome_fasta_files, function(x) grepl(genome_name, x))]
  gff_path              <- genome_gff_files[sapply(genome_gff_files, function(x) grepl(genome_name, x))]
  operon_list_path      <- operon_lists[sapply(operon_lists, function(x) grepl(genome_name, x))]
  operonic_gp_list_path <- operonic_gene_pairs[sapply(operonic_gene_pairs, function(x) grepl(genome_name, x))]
  protein_seq_path      <- predicted_protein_sequences[sapply(predicted_protein_sequences, function(x) grepl(genome_name, x))]
  ORFs_coordinates_path <- fread(ORFs_coordinates_main[sapply(ORFs_coordinates_main, function(x) grepl(genome_name, x))])
  
  genome            <- readDNAStringSet(genome_path)
  gff               <- ape::read.gff(gff_path)
  operon_list       <- fread(operon_list_path)
  operonic_gp_list  <- fread(operonic_gp_list_path)
  protein_sequences <- readAAStringSet(protein_seq_path)
  ORFs_coordinates  <- ORFs_coordinates_path
  
  # Retrieve the bDD AA sequences and convert to AAStringSet
  cat("*** Retrieving bDD AA Sequences *************************************** \n")
  print(aa_seq_patterns <- main_df %>% filter(genbank_gca == genome_name | genbank_gca_2 == genome_name | refseq_gcf == genome_name | refseq_gcf_2 == genome_name) %>% pull(seq_DLD_aa) %>% unique())
  print(IDs             <- main_df %>% filter(genbank_gca == genome_name | genbank_gca_2 == genome_name | refseq_gcf == genome_name | refseq_gcf_2 == genome_name) %>% select(ID, seq_DLD_aa) %>% unique())
  
  # Initialize a vector to store names of matching sequences
  matching_names <- character()
  
  # Loop through each pattern and find matches
  pattern <- aa_seq_patterns[1]
  for (pattern in aa_seq_patterns) {
    ID <- IDs %>% filter(seq_DLD_aa == pattern) %>% pull(ID) %>% unique() %>% paste(collapse = "_")
    aa_pattern <- AAString(pattern)
    matches <- vmatchPattern(aa_pattern, protein_sequences, max.mismatch = 1, fixed = TRUE)
    cat(matching_names <- unique(c(matching_names, names(protein_sequences)[sapply(matches, function(x) length(x) > 0)])), "\n")
    
    # save the matching sequences to a file - the file name should have the genome and matching name
    if (length(matching_names) > 0) {
      cat("Match found! \n")
      for (matching_name in matching_names) {
        # save the matching sequences to a file - the file name should have the genome and matching name
        matching_seq        <- protein_sequences[[matching_name]]
        print(matching_seq        <- AAStringSet(matching_seq))
        
        names(matching_seq) <- paste0(genome_name, "_", matching_name, "_", ID)
        writeXStringSet(matching_seq, file.path(output_path, paste0(names(matching_seq), "_matching_seq.faa")))
        cat("Overwriting", file.path(output_path, paste0(names(matching_seq), "_matching_seqs.faa")), "\n")
      }
    }
    
    for (protein in gene_list) {
      matching_seq        <- protein_sequences[[paste0("cds-", protein)]]
      matching_seq        <- AAStringSet(matching_seq)
      if (length(matching_seq) > 0 && OVERWRITE) {
        names(matching_seq) <- paste0(genome_name, "_", protein, "_", ID)
        writeXStringSet(matching_seq, file.path(output_path, paste0(names(matching_seq), ".fna")))
      }
    }
  }
}

################################################################################


extract_specific_operons <- function(genome_name, list, output_path = file.path(output_dir, "OPERONS", genome_name), OVERWRITE = TRUE) {
  
  # Initialize a list to store results for each genome
  genome_fasta_files <- list.files(file.path(input_dir, "GENOMES"), pattern = "\\.fna$", full.names = TRUE, recursive = TRUE)
  genome_path <- genome_fasta_files[sapply(genome_fasta_files, function(x) grepl(genome_name, x))]
  genome      <- readDNAStringSet(genome_path) 
  gff_path    <- genome_gff_files[sapply(genome_gff_files, function(x) grepl(genome_name, x))]    ; gff    <- ape::read.gff(gff_path)
  gff_data    <- read_tsv(gff_path, comment = "#", col_names = c("seqid","source","type","start","end","score" ,"strand","phase","attributes"), col_types = cols()) %>% filter(type == "CDS")
  
  # Create a data frame to store the operon components
  df <- data.frame(gene = list,
                   start = NA,
                   end = NA,
                   seqid = NA,
                   type = NA,
                   strand = NA,
                   stringsAsFactors = FALSE)
  
  
  # Iterate over df and match the 'gene' column with the 'attributes' column in gff_data
  for (i in 1:nrow(df)) {
    gene <- df$gene[i]
    matched_row <- gff_data %>% filter(grepl(gene, attributes))
    if (nrow(matched_row) > 0) {
      df$start[i]  <- matched_row$start[1]
      df$end[i]    <- matched_row$end[1]
      df$seqid[i]  <- matched_row$seqid[1]
      df$type[i]   <- matched_row$type[1]
      df$strand[i] <- matched_row$strand[1]
    }
  }
  
  # Display the updated df
  print(df)
  
  # extract the operon from the genome based on the start and end positions
  print(special_operon <- extract_operon(genome = genome,
                                         start  = min(df$start), 
                                         end    = max(df$end),
                                         strand = unique(df$strand)))
  special_operon <- DNAStringSet(special_operon)
  names(special_operon) <- paste0(genome_name, "_", paste0(list, collapse = "_"), "_operon")
  
  df$operon <- convert_to_character(special_operon)
  if (OVERWRITE) {
    # create ouutput dir if it does not exist
    if (!dir.exists(output_path)) {dir.create(output_path, showWarnings = FALSE)}
    fwrite(df, file.path(output_path, paste0(genome_name, "_", paste0(list, collapse = "_"), "_operon.csv")), append = FALSE)
    writeXStringSet(DNAStringSet(special_operon), file.path(output_path, paste0(genome_name, "_", paste0(list, collapse = "_"), "_operon.fna")))
    cat("Files written to", file.path(getwd(), output_path), "\n")
  }
  
  # add the fasta files of the individual genes based on each start and end position
  for (i in 1:nrow(df)) {
    gene <- df$gene[i]
    start <- df$start[i]
    end <- df$end[i]
    strand <- df$strand[i]
    gene_seq <- extract_operon(genome = genome, start = start, end = end, strand = strand)
    gene_seq <- DNAStringSet(gene_seq)
    names(gene_seq) <- paste0(genome_name, "_", gene)
    if (OVERWRITE) {
      writeXStringSet(gene_seq, file.path(output_path, paste0(genome_name, "_", gene, ".fna")))
    }
  }
  return(df)
}

##############################################################################
# Function to plot data
# *** DEBUGGING *** 
# plotting_means = plotting_means; plotting_stats = plotting_stats_relative; x_mean    = "Relative_Intensity_mean"; x_sem     = "Relative_Intensity_sem"; x_label   = "Relative IL-2 secretion"; cl_label = "CL_NAME_ON_PLOT"
prepare_and_plot <- function(plotting_means, plotting_stats, x_mean, x_sem, x_label, cl_label = "CL_NAME_ON_PLOT") {
  
  # print if stats has significance column
  if("significance" %in% colnames(plotting_stats)){ 
    cat("Yep, stats were calculated for ya!\n") 
  } else { 
    cat("Sorry, couldn't calculate stats for ya!\n") 
  }
  
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats$CL_NAME_ON_PLOT <- reorder(plotting_stats$CL_NAME_ON_PLOT, -plotting_stats$ORDER_NO)
  
  # !!sym(x_mean) dynamically references the column name stored in x_mean, allowing flexible use of different columns in plots or calculations.
  ELISA_PLOT <- ggplot(data = plotting_stats, aes(x = !!sym(x_mean), y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
    # geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_point(data = plotting_means, aes(x = !!sym(x_mean), y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
               col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = cl_label, values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    scale_color_manual(name = cl_label, values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    labs(x = x_label, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  if ("significance" %in% colnames(plotting_stats)) {
    ELISA_PLOT <- ELISA_PLOT +
      geom_errorbar(aes(y = !!sym(cl_label),
                        xmin = !!sym(x_mean) - !!sym(x_sem),
                        xmax = !!sym(x_mean) + !!sym(x_sem)),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_text(data = plotting_stats, aes(x = 1.2 * max(!!sym(x_mean)), y = CL_NAME_ON_PLOT, label = significance), 
                hjust = .5, vjust = 1, size = TEXT, 
                #angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)
      )
  }
  
  if (length(unique(plotting_means$POSITIVE_CTRL)) > 1 & x_mean == "Relative_Intensity_mean") {
    ELISA_PLOT <- ELISA_PLOT +
      facet_wrap(~POSITIVE_CTRL, scales = "free", ncol = 1)
  }
  
  print(ELISA_PLOT)
  
}

##############################################################################
  # Function to generate statistics from means
  run_statistics <- function(plotting_means, x_mean) {
    
    tryCatch({
      stat_significance_dt <- process_statistical_analysis(data = plotting_means, group_var = "CL_NAME_ON_PLOT", value_var = x_mean)

    }, error = function(e) {
      print("Statistical analysis failed.")
      print("The plot will be generated without statistical analysis.")
    })
    
    
    if (exists("stat_significance_dt")) {
      plotting_stats <- prepare_plotting_stats(plotting_data_main, stat_significance_dt)
      
      if (sum(plotting_stats$Relative_Intensity_sem) == 0) {
        group_vars_1 <- c("CL_NUMBER", "CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", 
                          "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", 
                          "ORDER_NO", "Plate", "POSITIVE_CTRL", "NEGATIVE_CTRL")
        
        plotting_stats_main <- plotting_data_main %>%
          group_by(!!!syms(group_vars_1)) %>%
          summarise(
            IL2_concentration_Dilution_Factor_mean = mean(Concentration),
            IL2_concentration_Dilution_Factor_sem  = sem(Concentration),
            Relative_Intensity_mean                = mean(Concentration_NORMALIZED),
            Relative_Intensity_sem                 = sem(Concentration_NORMALIZED)) %>%
          as.data.table() %>%
          left_join(stat_significance_dt) %>%
          distinct(CL_NUMBER, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, .keep_all = T)
        
        plotting_stats_main$CONDITION       <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
        plotting_stats_main$CL_NAME_ON_PLOT <- reorder(plotting_stats_main$CL_NAME_ON_PLOT, -plotting_stats_main$ORDER_NO)
        plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- "#BEBDBD"
        
        plotting_stats <- plotting_stats_main
      }
    } else {
      plotting_stats <- plotting_means
    }
    return(plotting_stats)
  }



################################################################################
### 20241031 FX TO EXTRACT DOMAINS FROM PROTEIN SEQUENCES ######################
################################################################################

# Function to find extended boundaries with adjustments for N and C terminal conditions
# extend_boundaries <- function(df) {
#   df %>%
#     group_by(accession) %>%
#     arrange(start) %>%
#     mutate(
      # extended_start = case_when(
      #   domain_position == "N" ~ start,
      #   domain_position == "C" ~ lag(end, default = NA, order_by = start) + 1,
      #   TRUE ~ lag(start, default = NA, order_by = start)
      # ),
      # extended_end = case_when(
      #   domain_position == "N" ~ lead(start, default = NA, order_by = start) - 1,
      #   domain_position == "C" ~ length,
      #   TRUE ~ lead(start, default = NA, order_by = start) - 1
      # ),
#       
      # # Fill NA defaults with original boundaries if needed
      # extended_start = ifelse(is.na(extended_start), start, extended_start),
      # extended_end = ifelse(is.na(extended_end), end, extended_end),
      # 
      # # Extract the extended sequence based on adjusted boundaries
      # fasta_extended = substr(fasta_FL, extended_start, extended_end),
#       
#       # Calculate medium boundaries with ±10 amino acids, respecting sequence and domain limits
#       medium_start = pmax(start - 10, extended_start),
#       medium_end = pmin(end + 10, extended_end),
#       
#       # Extract the medium sequence based on medium boundaries
#       medium_fasta = substr(fasta_FL, medium_start, medium_end)
#     ) %>%
#     ungroup()
# }

# 20241104 optimized to include slightly extended bDDs
extend_boundaries <- function(df) {
  df %>%
    group_by(accession) %>%
    arrange(start) %>%
    mutate(
      # Calculate extended boundaries based on position within the sequence
      extended_start = case_when(
        domain_position == "N" ~ start,
        domain_position == "C" ~ lag(end, default = NA, order_by = start) + 1,
        TRUE ~ lag(start, default = NA, order_by = start)
      ),
      extended_end = case_when(
        domain_position == "N" ~ lead(start, default = NA, order_by = start) - 1,
        domain_position == "C" ~ length,
        TRUE ~ lead(start, default = NA, order_by = start) - 1
      ),
      
      # Fill NA defaults with original boundaries if needed
      extended_start = ifelse(is.na(extended_start), start, extended_start),
      extended_end = ifelse(is.na(extended_end), end, extended_end),
      
      # Extract the extended sequence based on adjusted boundaries
      fasta_extended = substr(fasta_FL, extended_start, extended_end),
      
      # Calculate medium boundaries with ±10 amino acids, constrained by protein limits
      medium_start = pmax(1,      start - 10),
      medium_end   = pmin(length, end   + 10),
      
      # Ensure medium boundaries do not exceed extended boundaries
      medium_start = pmax(extended_start, medium_start),
      medium_end = pmin(extended_end, medium_end),
      
      # Extract the medium sequence based on medium boundaries
      fasta_medium = substr(fasta_FL, medium_start, medium_end)
    ) %>%
    ungroup()
}

# Prepare headers and sequences for strict and extended FASTA files
prepare_fasta_data <- function(df, domain_name = "bDLD3") {
  # Strict domain FASTA
  df_strict <- df %>%
    filter(pfam_name == domain_name) %>%
    mutate(
      header = glue(">{accession} | {domain_name} aa({start}-{end} | {tax_name} | TaxID={tax_id} | Mode=strict"),
      sequence = fasta_domain
    )
  
  # Extended domain FASTA
  df_extended <- df %>%
    filter(pfam_name == domain_name) %>%
    mutate(
      header = glue(">{accession} | {domain_name} aa{extended_start}-{extended_end} | {tax_name} | TaxID={tax_id} | Mode=extended"),
      sequence = fasta_extended
    )
  
  # Combine strict and extended data into a list for easy writing
  list(strict = df_strict, extended = df_extended)
}

# Write domain sequences to FASTA files
write_domain_fasta <- function(df_list, domain_name = "bDLD3", strict_filename, extended_filename) {
  # Check if the sequences are DNA or protein
  if (all(grepl("^[ACGTacgt]*$", df_list$strict$sequence))) {
    strict_fasta <- DNAStringSet(df_list$strict$sequence)
  } else {
    strict_fasta <- AAStringSet(df_list$strict$sequence)
  }
  names(strict_fasta) <- df_list$strict$header
  writeXStringSet(strict_fasta, filepath = strict_filename)
  
  if (all(grepl("^[ACGTacgt]*$", df_list$extended$sequence))) {
    extended_fasta <- DNAStringSet(df_list$extended$sequence)
  } else {
    extended_fasta <- AAStringSet(df_list$extended$sequence)
  }
  names(extended_fasta) <- df_list$extended$header
  writeXStringSet(extended_fasta, filepath = extended_filename)
}

################################################################################
### 20241101 FX TO BACKTRANSLATE SEQUENCES BASED ON CODON USAGE ################
################################################################################

# Define E. coli codon usage table for each amino acid
ecoli_codon_table <- list(
  A = "GCT", C = "TGT", D = "GAT", E = "GAA", F = "TTT",
  G = "GGT", H = "CAT", I = "ATT", K = "AAA", L = "CTG",
  M = "ATG", N = "AAT", P = "CCT", Q = "CAA", R = "CGT",
  S = "TCT", T = "ACT", V = "GTT", W = "TGG", Y = "TAT",
  # Stop codons: Let's use TAA as the default
  "*" = "TAA"
)  

# Define Mus musculus codon usage table for each amino acid
mouse_codon_table <- list(
  A = "GCC", C = "TGC", D = "GAC", E = "GAG", F = "TTC",
  G = "GGC", H = "CAC", I = "ATC", K = "AAG", L = "CTG",
  M = "ATG", N = "AAC", P = "CCC", Q = "CAG", R = "CGC",
  S = "AGC", T = "ACC", V = "GTC", W = "TGG", Y = "TAC",
  # Stop codons: Let's use TAA as the default
  "*" = "TAA"
)

# Function to backtranslate an amino acid sequence to nucleotide sequence
backtranslate <- function(aa_seq, codon_table) {
  # Split sequence into individual amino acids and map to codons
  codons <- sapply(strsplit(aa_seq, NULL)[[1]], function(aa) codon_table[[aa]])
  # Concatenate codons to form the nucleotide sequence
  paste(codons, collapse = "")
}

################################################################################
################################################################################
################################################################################
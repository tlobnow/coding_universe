################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

ELISA_Fx <- function(Input_Directory, Output_Directory) {
  # Initialize as an empty data frame
  All_plates_data = data.frame()
  # Input_Directory = "~/Desktop/ALL"
  # Output_Directory = "~/Desktop/ALL"
  
  # Get the list of subdirectories matching the pattern "Plate_"
  subdirs <- list.files(Input_Directory, recursive = FALSE, full.names = TRUE, pattern = "Plate_\\d+_\\d{8}$")
  # input_plate_dir = "~/Desktop/ALL/Plate_1_20220609"
  
  if (length(subdirs) > 0) {
    print("Plates exist!")
    for (input_plate_dir in subdirs) {
      
      Input_plate <- input_plate_dir
      
      #Reading Plate Treatment 
      MEASUREMENTS <- fread(paste0(Input_plate, "/MEASUREMENTS.csv"), header = F)
      CELL_LINES   <- fread(paste0(Input_plate, "/CELL_LINES.csv"),   header = F)
      CONDITIONS   <- fread(paste0(Input_plate, "/CONDITIONS.csv"),   header = F)
      STIM_DAYS    <- fread(paste0(Input_plate, "/STIM_DAYS.csv"),    header = F)
      DILUTIONS    <- fread(paste0(Input_plate, "/DILUTIONS.csv"),    header = F)
      
      #Converting tables into vector for to make a single table
      MEASUREMENT <- as.vector(as.matrix(MEASUREMENTS))
      CELL_LINE   <- as.vector(as.matrix(CELL_LINES))
      CONDITION   <- as.vector(as.matrix(CONDITIONS))
      STIM_DAY    <- as.vector(as.matrix(STIM_DAYS))
      DILUTION    <- as.vector(as.matrix(DILUTIONS))
      
      #Creating Table containing all plate Information
      Plate <- NULL
      Plate$MEASUREMENT <- MEASUREMENT
      Plate$CELL_LINE   <- CELL_LINE
      Plate$CONDITION   <- CONDITION
      Plate$STIM_DAY    <- STIM_DAY
      
      if (exists("DILUTION")) {
        Plate$DILUTION <- DILUTION
      } else {
        Plate$DILUTION <- ifelse(Plate$CONDITION != "CALIBRATION", 5, NA_real_)
      }
      Plate$DILUTION <- as.numeric(Plate$DILUTION)
      
      rm(MEASUREMENT, CELL_LINE, CONDITION, STIM_DAY, DILUTION)
      
      Plate <- Plate %>% as.data.table()
      
      #Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINE != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate[Plate$CONDITION == "CALIBRATION"]
      Plate_Standards$CELL_LINE <- as.numeric(Plate_Standards$CELL_LINE)
      
      Plate_Standards <- Plate_Standards %>%
        group_by(CELL_LINE) %>% 
        summarise(MEASUREMENT_mean = mean(MEASUREMENT)) %>%
        mutate(CELL_LINE = as.numeric(CELL_LINE),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}"))) %>% 
        arrange(CELL_LINE)
      
      # We will only use standard curve values of 1 and below (machine is optimized to measure absorptions between 0 and 1.1)
      Fit <- lm(CELL_LINE ~ MEASUREMENT_mean - 1, data = Plate_Standards[Plate_Standards$MEASUREMENT_mean <= 1.1, ])
      
      R       <- summary(Fit)$r.squared
      Rsquare <- signif(R, digits = 4)
      
      print(paste0("Secretion = slope*Intensity"))
      print(paste0("Secretion = ", Fit$coefficients[1],"*Intensity"))
      Plate_Standards <- Plate_Standards %>% mutate(Fit_Test = (Fit$coefficients[1]*MEASUREMENT_mean))
      
      # Plotting Standard Curve
      p <- ggplot(data = Plate_Standards) +
        geom_point(aes(x = MEASUREMENT_mean, y = CELL_LINE, col = MEASUREMENT_mean <= 1.5), size = 5) +
        geom_line(aes(x = MEASUREMENT_mean, y = Fit_Test), linetype = "dashed") +
        annotate('text', x = 0.15, y = 700, label = paste0("R^2 = ", Rsquare), size = 10) +
        annotate('text', 
                 x = max(Plate_Standards$MEASUREMENT_mean) - (0.25 * max(Plate_Standards$MEASUREMENT_mean)),
                 y = 150, label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        labs(x = "Measured Values",
             y = "IL-Concentration (pg/mL)") +
        ggtitle(label = paste0(basename(Input_plate)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text  = element_text(size = 20)) +
        theme(legend.position = "none")
      
      # Saving the plot
      Save_Name <- paste0(basename(Input_plate), "_Standard_Curve.pdf")
      Save_Name <- file.path(Output_Directory, Save_Name)
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      
      # Fitting Data To Standard Curve ----------------------------------------
      Plate <- Plate %>% 
        filter(CONDITION != "CALIBRATION") %>% 
        mutate(Plate = as.numeric(gsub("(DR_)?Plate_(\\d+)_\\d{8}$", "\\2", basename(input_plate_dir))),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}")),
               MEASUREMENT = as.numeric(MEASUREMENT),
               
               # METHOD I: Correct negative values to zero, multiply measurements by dilution factor and THEN extrapolate values based on SC
               # ########  20231206 @Fakun: We should adjust for the dilution factor BEFORE extrapolating values based on the Standard Curve..
               # ########
               # MEASUREMENT_DIL_ADJ = (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT)*DILUTION),
               # Concentration = (Fit$coefficients[1]*MEASUREMENT_DIL_ADJ),
               
               # METHOD II: Extrapolate values based on SC and THEN multiply measurements by dilution factor
               # ########  20231206 @Finn: set machine measurement errors to zero (raw values below zero should be set to zero before extrapolating etc)
               # ########  20231206 @Finn:  We remove the column name *Concentration_DILUTION_FACTOR* from downstream analysis and simply stick to *Concentration*!
               # ########
               # Concentration = (Fit$coefficients[1]*MEASUREMENT),
               # Concentration_DILUTION_FACTOR = Concentration*DILUTION,
               Concentration = (Fit$coefficients[1] * (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT))),
               Concentration = Concentration * DILUTION,
               
               Is_Dose_Response = ifelse(str_detect(basename(input_plate_dir), "^DR_"), TRUE, FALSE)
        )
      
      All_plates_data <- rbind(All_plates_data, Plate)
    }
  } else {
    print("No plates found!")
  }
  return(All_plates_data)
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Define the function to save plots with specified dimensions
save_plots <- function(folder_name, plots, path = save_to, width = 12, height = 8) {
  # Create the full path for the folder
  folder_path <- file.path(save_to, folder_name)
  
  # Check if the folder exists, and if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # Loop through the plots and save them
  for (i in seq_along(plots)) {
    file_base <- file.path(folder_path, paste0(folder_name, "_plot_", i))
    
    # Save as SVG
    ggsave(paste0(file_base, ".svg"), plot = plots[[i]], device = "svg", width = width, height = height)
    
    # Save as PNG
    ggsave(paste0(file_base, ".png"), plot = plots[[i]], device = "png", width = width, height = height)
  }
}


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Vectorized function for pattern matching
matches_any_pattern_vec <- Vectorize(function(string, patterns) {
  any(sapply(patterns, function(pattern) grepl(pattern, string)))
}, vectorize.args = "string")

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# FILTER_VALUES = c("MyD88-GFP-synTRAF6-BD-1x", "MyD88-GFP-synTRAF6-BD-3x", "MyD88-GFP-synTRAF6-BD-5x")
# FILTER_TYPE = "COHORT"        # Should be "COHORT" or "DAY"
# COLOR = "salmon"
# SEED = 600
# plot_results = TRUE
# plot_pval = TRUE
# run_anova = FALSE
# plot_faceted_by_date = FALSE
# save_to = "~/Desktop"
# POSITIVE_CTRL = c("WT")
# NEGATIVE_CTRL = 0
# DATA = All_plates_data

filter_data <- function(DATA, FILTER_VALUES, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL) {
  filter_pattern_func <- if (FILTER_TYPE %in% "COHORT") {
    function(value, df) matches_any_pattern_vec(df$CELL_LINE, value) | matches_any_pattern_vec(df$CL_NAME_ON_PLOT, value)
  } else if (FILTER_TYPE %in% "DAY") {
    function(value, df) matches_any_pattern_vec(df$Date, value)
  } else {
    stop("Invalid FILTER_TYPE. Must be either 'COHORT' or 'DAY'.")
  }
  
  # Generate lists for plates, dates, and stim_days based on the filter type
  plates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Plate) %>%
      pull(Plate)
  })
  
  dates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Date) %>%
      pull(Date)
  })
  
  stim_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(STIM_DAY) %>%
      pull(STIM_DAY)
  })
  
  condition_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(CONDITION) %>%
      pull(CONDITION)
  })
  
  # Adding names to the list elements
  names(plates_list)    <- FILTER_VALUES
  names(dates_list)     <- FILTER_VALUES
  names(stim_list)      <- FILTER_VALUES
  names(condition_list) <- FILTER_VALUES
  
  # Create subset based on the lists
  subset_list <- lapply(names(plates_list), function(value) {
    plates <- plates_list[[value]]
    dates  <- dates_list[[value]]
    stim   <- stim_list[[value]]
    cond   <- condition_list[[value]]
    
    DATA %>%
      filter(Date %in% dates,
             Plate %in% plates,
             STIM_DAY %in% stim,
             CONDITION %in% cond,
             filter_pattern_func(c(POSITIVE_CTRL, NEGATIVE_CTRL, value), DATA))
  })
  
  # Combine the subsets and return
  FILTERED_SUBSET <- bind_rows(subset_list)
  return(FILTERED_SUBSET)
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# DATA = FILTERED_SUBSET
calculate_baseline_and_control <- function(DATA, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL) {
  # Determine the group_by and filter parameters based on FILTER_TYPE
  if (FILTER_TYPE == "COHORT") {
    group_vars <- c("Date", "STIM_DAY")
    filter_var <- "CELL_LINE"
  } else if (FILTER_TYPE == "DAY") {
    group_vars <- c("Date", "STIM_DAY", "Plate")
    filter_var <- "CELL_LINE"  # Assuming the correct variable for DAY filter
  } else {
    stop("Invalid FILTER_TYPE. Must be either 'cohort' or 'day'.")
  }
  
  # Debugging information
  # print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  # print(paste("Filtering using variable:", filter_var))
  
  # Calculate the baseline control value
  if (is.character(NEGATIVE_CTRL)) {
    baseline <- DATA %>%
      group_by(!!!syms(group_vars)) %>%
      filter(matches_any_pattern_vec(!!sym(filter_var), NEGATIVE_CTRL), CONDITION == "UNSTIM") %>% 
      summarise(baseline_control_value = mean(Concentration))
    
    # Join the calculated values with the dataset
    data <- left_join(DATA, baseline, by = group_vars) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value, TRUE ~ Concentration))
  } else {
    # baseline <- DATA %>%
    #   group_by(!!!syms(group_vars)) %>%
    #   summarise(baseline_control_value = min(Concentration))
    baseline <- data %>%
      group_by(Date) %>%
      filter(Date == FILTER_VALUE) %>%
      summarise(baseline_control_value = min(Concentration))
    
    # Join the calculated values with the dataset
    data <- left_join(DATA, baseline) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value,
                                               TRUE ~ Concentration))
  }
  
  
  if (is.character(POSITIVE_CTRL)) {
    control_mean_per_day <- data %>%
      filter(matches_any_pattern_vec(!!sym(filter_var), POSITIVE_CTRL), CONDITION == "STIM") %>% 
      group_by(!!!syms(group_vars)) %>%
      summarise(control_mean_MEASUREMENT = case_when(mean(Concentration_REDUCED) > 0 ~ mean(Concentration_REDUCED), TRUE ~ -Inf))
    
    # Join the calculated control means
    data <- left_join(data, control_mean_per_day, by = group_vars)
    
    # Perform normalization
    DATA_NORMALIZED <- data %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_mean_MEASUREMENT < 0 ~ 0,
                                                  TRUE ~ Concentration_REDUCED / control_mean_MEASUREMENT),
             triplicate_mean_per_day = mean(Concentration_NORMALIZED)) %>%
      ungroup()
    
  } else {
    # Calculate control mean using an alternative method, e.g., mean value per group
    # control_mean_per_day <- data %>%
    #   group_by(!!!syms(group_vars)) %>%
    #   summarise(control_mean_MEASUREMENT = mean(Concentration_REDUCED))
    control_mean_per_day <- data %>%
      group_by(STIM_DAY) %>%
      reframe(Concentration_REDUCED = Concentration_REDUCED) %>%
      ungroup() %>%
      group_by(STIM_DAY) %>%
      reframe(control_MEASUREMENT = max(Concentration_REDUCED),
              STIM_DAY = STIM_DAY) %>%
      unique()
    
    # Join the calculated control means
    data <- left_join(data, control_mean_per_day)
    
    # Perform normalization
    DATA_NORMALIZED <- data %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_MEASUREMENT < 0 ~ 0,
                                                  TRUE ~ Concentration_REDUCED / control_MEASUREMENT),
             triplicate_mean_per_day = mean(Concentration_NORMALIZED)) %>%
      ungroup()
  }
  
  return(DATA_NORMALIZED)
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# DATA = DATA_NORMALIZED
perform_statistical_analysis <- function(DATA, GROUP_BY_COLUMN) {
  # Internal function for pairwise t-test
  pairwise_ttest <- function(DATA, return_annotation = FALSE) {
    p_values <- pairwise.t.test(DATA$triplicate_mean_per_day, DATA$CONDITION, p.adjust.method = "none")$p.value
    
    if (return_annotation) {
      p_annotation <- ifelse(p_values < 0.001, '***', 
                             ifelse(p_values < 0.01, '**', 
                                    ifelse(p_values < 0.05, '*', 
                                           'ns')))
      return(p_annotation)
    } else {
      return(formatC(p_values, format = "e", digits = 3))
    }
  }
  
  # Calculate statistical significance using a t-test for each group
  annotations <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), pairwise_ttest)
  p_values    <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), function(DATA) pairwise_ttest(DATA, return_annotation = TRUE))
  
  return(list(annotations = annotations, p_values = p_values))
}
# Example usage
# results <- perform_statistical_analysis(MEANS, "CL_NAME_ON_PLOT")

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

create_plot <- function(FILTER_VALUES, FILTER_TYPE, MEANS, MOM_SUBSET, STATISTICAL_RESULTS, COLOR = "salmon", 
                        x_label = xlabel, y_label = y_label, plot_title = plot_title,
                        subtitle = paste0("Assay performed for ", FILTER_VALUES, ", based on filter type ", FILTER_TYPE),
                        SEED = 600, plot_pval = T, plot_faceted_by_date = F) {
  plot <- ggplot(MEANS, aes(x = CL_NAME_ON_PLOT)) +
    geom_col(     data = MOM_SUBSET,  aes(y = triplicate_mean_per_day, fill = CONDITION), position = position_dodge(width = 1), alpha = 0.5) +
    geom_point(   data = MEANS,       aes(y = triplicate_mean_per_day, group = CONDITION, shape = STIM_DAY),               position = position_jitterdodge(jitter.height = 0, jitter.width = 1.2, seed = SEED), col = "white", size = 4) +
    geom_point(   data = MEANS,       aes(y = triplicate_mean_per_day, group = CONDITION, shape = STIM_DAY), position = position_jitterdodge(jitter.height = 0, jitter.width = 1.2, seed = SEED), 
                  col = "black",
                  size = 3) +
    geom_errorbar(data = MOM_SUBSET, aes(ymin = triplicate_mean_per_day - triplicate_sd_per_day, 
                                         ymax = triplicate_mean_per_day + triplicate_sd_per_day, group = CONDITION), width = 0.25, position = position_dodge(width = 1)) +
    labs(x = x_label,
         y = y_label) +
    scale_fill_manual(values = c("UNSTIM" = "gray50", "STIM" = COLOR)) +
    ggtitle(plot_title) +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (FILTER_TYPE == "DAY") {
    plot_faceted_by_date = T
    plot_pval = F
  }
  
  if (plot_pval == T & plot_faceted_by_date == F) {
    plot <- plot +
      geom_text(data = data.frame(CL_NAME_ON_PLOT = names(STATISTICAL_RESULTS[["annotations"]]), y = max(MEANS$triplicate_mean_per_day) + 0.1, label = STATISTICAL_RESULTS[["annotations"]]), aes(y = y, label = label), col = "gray40", size = 4) +
      geom_text(data = data.frame(CL_NAME_ON_PLOT = names(STATISTICAL_RESULTS[["p_values"]]), y = max(MEANS$triplicate_mean_per_day) + 0.2, label = STATISTICAL_RESULTS[["p_values"]]), aes(y = y, label = label), size = 7)
  } else if (plot_pval == F & plot_faceted_by_date == T) {
    plot <- plot  + 
      facet_wrap(~Date, scales = "free_x")
  } else {
    plot_pval = F
    plot <- plot  + 
      facet_wrap(~Date, scales = "free_x")
    print("Overriding p-value plotting. Option not available for date-faceted plots.")
  }
  
  return(plot)
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# I am providing an example cohort of the All_plates_data data frame to debug the functions:
# FILTER_VALUES = c("BDLD_57", "BDLD_6-")  # Example cohort names
# FILTER_TYPE = "COHORT"        # Should be "COHORT" or "DAY"
# save_to = "~/Desktop"
# POSITIVE_CTRL = c("3E10_GFP", "MyD88_GFP")
# NEGATIVE_CTRL = c("204_TRIPLE_KO", "tKO_EL4")

plot_ELISA <- function(FILTER_VALUES, FILTER_TYPE, POSITIVE_CTRL = 0, NEGATIVE_CTRL = 0,
                       COLOR = COLOR,#"salmon", 
                       SEED = 600, 
                       plot_pval = TRUE,  plot_results = TRUE, 
                       run_anova = FALSE, plot_faceted_by_date = FALSE, 
                       x_label = "", y_label = "relative IL-2 conc.", plot_title = "IL-2 ELISA",
                       subtitle = paste0("Assay performed for ", FILTER_VALUES, ", based on filter type ", FILTER_TYPE)) {
  
  ALL_COHORT_DATA <- lapply(FILTER_VALUES, function(FILTER_VALUE) {
    # Step 1: Filter & Subset the main data into lists
    COHORT_DATA <- filter_data(All_plates_data, FILTER_VALUE, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL)
    # Step 2 and Step 3
    NORMALIZED_DATA <- calculate_baseline_and_control(COHORT_DATA, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL)
    return(NORMALIZED_DATA)
  })
  
  # Step 4: Combine the subset lists into one data frame
  COMBINED_DATA <- bind_rows(ALL_COHORT_DATA)
  # print(COMBINED_DATA)  # Debugging line
  
  # Step 5: Calculate the means per cohort
  MEANS <- COMBINED_DATA %>%
    group_by(CELL_LINE, CONDITION, Date) %>%
    distinct(triplicate_mean_per_day, STIM_DAY, .keep_all = TRUE) %>%
    ungroup()
  # print(MEANS)  # Debugging line
  
  # Step 6: Calculate the mean of means (MOM) per cohort
  MOM_SUBSET <- MEANS %>%
    group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION, Date) %>%
    summarise(triplicate_sd_per_day = sd(triplicate_mean_per_day),
              triplicate_mean_per_day = mean(triplicate_mean_per_day)) %>%
    distinct(.keep_all = TRUE)
  # print(MOM_SUBSET)  # Debugging line
  
  # Step 7: Perform statistics per cohort
  STATISTICAL_RESULTS <- perform_statistical_analysis(MEANS, "CL_NAME_ON_PLOT")
  # print(STATISTICAL_RESULTS)  # Debugging line
  
  # Step 8: Plot the results
  if (plot_results) {
    if (all(POSITIVE_CTRL == 0) & all(NEGATIVE_CTRL == 0)) {
      PLOT <- ggplot(MEANS, aes(x = CL_NAME_ON_PLOT)) +
        geom_col(     data = MOM_SUBSET,  aes(y = triplicate_mean_per_day, fill = CONDITION, group = CONDITION), position = position_dodge(width = 1), alpha = 0.5) +
        labs(x = x_label,
             y = y_label) +
        scale_fill_manual(values = c("UNSTIM" = "gray50", "STIM" = COLOR)) +
        labs(x = x_label, y = y_label) +
        ggtitle(label = plot_title, subtitle = subtitle) +
        theme_cowplot() +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~STIM_DAY)
      
      PLOT
      
    } else {
      PLOT <- create_plot(FILTER_VALUES, FILTER_TYPE, MEANS, MOM_SUBSET, STATISTICAL_RESULTS, COLOR = "salmon", 
                          SEED = 600, x_label = x_label, y_label = y_label, plot_title = plot_title, subtitle = subtitle)
    }
    list(PLOT, COMBINED_DATA, MEANS, MOM_SUBSET, MEANS)
  }
}


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

plot_dose_response_ELISA <- function(DATA, FILTER_VALUES, 
                                     COLOR_PALETTE_KEY = color_palette_key,
                                     SEED = 600, 
                                     plot_pval = TRUE, plot_results = TRUE, 
                                     run_anova = FALSE, plot_faceted_by_date = FALSE, 
                                     x_label = "Stimulant Conc.", 
                                     y_label = "Relative Response", 
                                     plot_title = "Dose Response IL-2 ELISA",
                                     subtitle = paste0("Assay performed on ", FILTER_VALUES),
                                     without_scientific_notation = FALSE) {
  
  group_vars <- c("CELL_LINE", "CL_NAME_ON_PLOT", "CONDITION")
  extra_var <- "STIM_DAY" 
  cat("group_vars used:\n")
  cat(group_vars, "\n")
  cat(paste("extra_var used:", extra_var, "\n"))
  
  ALL_COHORT_DATA <- lapply(FILTER_VALUES, function(FILTER_VALUE) {
    
    cat("Step 1: Filter & Subset the main data into lists\n")
    COHORT_DATA <- All_plates_data %>%
      group_by(CELL_LINE, STIM_DAY, CONDITION, Date, Plate) %>%
      filter(Date == FILTER_VALUE)
    
    cat("Step 2: Calculate baseline values (lowest value per plate)\n")
    baseline <- COHORT_DATA %>%
      group_by(Plate) %>%
      filter(Date == FILTER_VALUE) %>%
      summarise(baseline_control_value = min(Concentration))
    
    cat("Step 3: Subtract the lowest value from all measurements to normalize per plate (Machine Normalization).\n")
    COHORT_DATA <- left_join(COHORT_DATA, baseline) %>%
      mutate(Concentration_REDUCED    = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value,
                                                  TRUE ~ Concentration))
    cat("Step 4: Baseline-corrected per cell line, per condition, per stimulation and retrieve the max values to normalize for day-to-day differences.\n")
    control_mean_per_day <- COHORT_DATA %>%
      group_by(CELL_LINE, CONDITION, Date, STIM_DAY, Plate) %>%
      reframe(Concentration_REDUCED = Concentration_REDUCED) %>%
      ungroup() %>%
      group_by(STIM_DAY) %>%
      reframe(control_MEASUREMENT = max(Concentration_REDUCED))
    
    cat("Step 5: Normalize the baseline-corrected values by dividing each value by the max value per stimulation day.\n")
    COHORT_DATA <- left_join(COHORT_DATA, control_mean_per_day) %>%
      group_by(CELL_LINE, CONDITION, Date, STIM_DAY) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_MEASUREMENT < 0 ~ 0, TRUE ~ Concentration_REDUCED / control_MEASUREMENT),
             triplicate_mean_per_day = mean(Concentration_NORMALIZED)) %>%
      ungroup()
    
    # Make sure that the conditions are now read as numbers for plotting!
    COHORT_DATA$CONDITION <- as.numeric(as.character(COHORT_DATA$CONDITION))
    
    # Additional processing steps can be added here
    return(COHORT_DATA)
  })
  
  COHORT_DATA <- ALL_COHORT_DATA[[1]]
  
  group_and_summarize <- function(COHORT_DATA, group_vars = c("CELL_LINE", "CL_NAME_ON_PLOT", "CONDITION"), extra_var = "STIM_DAY") {
    
    cat("Step 6: Calculate the normalized means per grouped variables (group_vars + extra_vars)\n")
    MEANS <- COHORT_DATA %>%
      group_by_at(c(group_vars, extra_var)) %>%
      distinct(triplicate_mean_per_day, STIM_DAY, .keep_all = TRUE) %>%
      ungroup()
    # print(MEANS)  # Debugging line
    
    cat("Step 7: Calculate the mean of means (MOM) per cell line and per condition.\n")
    MOM_SUBSET <- COHORT_DATA %>%
      group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION) %>%
      summarise(triplicate_sd_per_day = sd(triplicate_mean_per_day),
                triplicate_mean_per_day = mean(triplicate_mean_per_day)) %>%
      distinct(.keep_all = TRUE)
    # print(MOM_SUBSET)  # Debugging line
    
    return(list(MEANS = MEANS, MOM_SUBSET = MOM_SUBSET))
  }
  
  # Group, Summarize, and Extract Results from group_and_summarize
  result_list <- group_and_summarize(COHORT_DATA, group_vars, extra_var)
  MEANS       <- result_list$MEANS
  MOM_SUBSET  <- result_list$MOM_SUBSET
  
  cat("Step 8: Plot the results.\n")
  
  if (plot_results) {
    create_dose_response_plot <- function(MEANS = MEANS, MOM_SUBSET = MOM_SUBSET, plot_title = plot_title, subtitle = subtitle, x_label = x_label, y_label = y_label, COLOR_PALETTE_KEY = COLOR_PALETTE_KEY) {
      
      # Deriving breaks from unique CONDITION values
      condition_breaks <- sort(unique(MOM_SUBSET$CONDITION))
      
      PLOT <- ggplot(data = MOM_SUBSET, aes(x = CONDITION, y = triplicate_mean_per_day, group = CELL_LINE, color = CELL_LINE)) +
        geom_path(size = 1.5) +
        color_palette(palette = COLOR_PALETTE_KEY) +
        fill_palette(palette  = COLOR_PALETTE_KEY) +
        geom_point(data = MEANS, aes(x = CONDITION, y = triplicate_mean_per_day, fill = CELL_LINE), size = 3, color = "black", shape = 21) +
        geom_errorbar(data = MOM_SUBSET, col = "black", linewidth = .4, width = 0.25,
                      aes(x = CONDITION, y = triplicate_mean_per_day, 
                          ymin = triplicate_mean_per_day - triplicate_sd_per_day, 
                          ymax = triplicate_mean_per_day + triplicate_sd_per_day)) +
        labs(x = x_label, y = y_label) +
        ggtitle(label = plot_title, subtitle = subtitle) +
        theme_cowplot() +
        theme(legend.position = "bottom", axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
      
      if (without_scientific_notation) {
        PLOT <- PLOT +
          scale_x_continuous(trans = "log10", breaks = condition_breaks,
                             # removes the scientific notation
                             labels = scales::label_comma())
      } else {
        PLOT <- PLOT +
          scale_x_continuous(trans = "log10", breaks = condition_breaks)
      }
      
      return(PLOT)
    }
    
    DR_PLOT <- create_dose_response_plot(MEANS = MEANS, MOM_SUBSET = MOM_SUBSET, plot_title = plot_title, subtitle = subtitle, x_label = x_label, y_label = y_label, COLOR_PALETTE_KEY = COLOR_PALETTE_KEY)
    list(DR_PLOT, MOM_SUBSET, MEANS, COHORT_DATA)
  }
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Plot for a specific day
# FILTER_VALUES = c("2023-11-21")
# FILTER_VALUES = c("2023-11-16")
# FILTER_TYPE = "DAY"
# POSITIVE_CTRL = 0
# NEGATIVE_CTRL = 0
# DATA = All_plates_data

plot_dose_response_ELISA_2 <- function(FILTER_VALUES, FILTER_TYPE = "DAY", DATA = All_plates_data, POSITIVE_CTRL = 0, NEGATIVE_CTRL = 0,
                                       x_label = "X-Axis", y_label = "Y-Axis", 
                                       plot_title = "Dose Response Plot", 
                                       subtitle = paste0("Assay performed on ", FILTER_VALUES),
                                       COLOR_PALETTE_KEY = color_palette_key) {
  cat("This function is optimized for FILTER_TYPE = 'DAY'. Using this as default for subsetting data.\n")
  
  filter_pattern_func <- if (FILTER_TYPE %in% "COHORT") {
    # function(value, df) matches_any_pattern_vec(df$CELL_LINE, value) | matches_any_pattern_vec(df$CL_NAME_ON_PLOT, value)
    print("plot_dose_response_ELISA_2() is currently optizimed for 'DAY' inputs only. Please adjust your input accordingly to use this function.")
  } else if (FILTER_TYPE %in% "DAY") {
    function(value, df) matches_any_pattern_vec(df$Date, value)
    # cat(paste0("Filter Type is 'DAY'. Subsetting for the Date supplied in FILTER_VALUES: ", FILTER_VALUES, "\n"))
  } else {
    stop("Invalid FILTER_TYPE. Must be either 'COHORT' or 'DAY'.")
  }
  
  # Generate lists for plates, dates, and stim_days based on the filter type
  plates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Plate) %>%
      pull(Plate)
  })
  
  dates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Date) %>%
      pull(Date)
  })
  
  stim_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(STIM_DAY) %>%
      pull(STIM_DAY)
  })
  
  condition_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(CONDITION) %>%
      pull(CONDITION)
  })
  
  # Adding names to the list elements
  names(plates_list)    <- FILTER_VALUES
  names(dates_list)     <- FILTER_VALUES
  names(stim_list)      <- FILTER_VALUES
  names(condition_list) <- FILTER_VALUES
  
  # Create subset based on the lists
  subset_list <- lapply(names(plates_list), function(value) {
    plates <- plates_list[[value]]
    dates  <- dates_list[[value]]
    stim   <- stim_list[[value]]
    cond   <- condition_list[[value]]
    
    DATA %>%
      filter(Date %in% dates,
             Plate %in% plates,
             STIM_DAY %in% stim,
             CONDITION %in% cond,
             filter_pattern_func(c(POSITIVE_CTRL, NEGATIVE_CTRL, value), DATA))
  })
  
  # Combine the subsets and return
  FILTERED_SUBSET <- bind_rows(subset_list)
  # return(FILTERED_SUBSET)
  cat("Filtered Data.\n")
  
  # DATA = FILTERED_SUBSET
  
  group_vars <- c("Date", "STIM_DAY", "Plate")
  filter_var <- "CELL_LINE"  # Assuming the correct variable for DAY filter
  
  # Debugging information
  # print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  # print(paste("Filtering using variable:", filter_var))
  
  # DATA_NORMALIZED <- DATA %>% # Debugging line
  DATA_NORMALIZED <- FILTERED_SUBSET %>%
    group_by(Date, Plate) %>%
    filter(Date == FILTER_VALUES) %>%
    mutate(baseline_control_value =  min(Concentration),
           Concentration_REDUCED = Concentration - baseline_control_value,
           Concentration_REDUCED = case_when(Concentration_REDUCED < 0 ~ 0, TRUE ~ Concentration_REDUCED)) %>%
    ungroup() %>%
    group_by(Plate) %>%
    mutate(max_control_value = max(Concentration),
           Concentration_REDUCED_NORMALIZED = Concentration_REDUCED/max_control_value)
  cat("Extracted min value per plate to perform baseline-normalization.\n")
  cat("Extracted max value per plate to perform upper limit normalization for relative secretion.\n")
  
  DATA_NORMALIZED <- DATA_NORMALIZED %>%
    group_by(CELL_LINE, CONDITION, STIM_DAY) %>%
    mutate(triplicate_mean_per_day = mean(Concentration_REDUCED_NORMALIZED),
           triplicate_sd_per_day = sd(Concentration_REDUCED_NORMALIZED))
  # return(DATA_NORMALIZED)
  cat("Calculated means per cell line, condition, and stimulation.\n")
  
  # Set the dodge width
  dodge_width  <- 0.9
  jitter_width <- 0.3
  
  cat("Plotting results:\n")
  DR2_PLOT_1 <- ggplot(DATA_NORMALIZED, aes(x = CONDITION, y = triplicate_mean_per_day, fill = CELL_LINE, group = CELL_LINE)) +
    geom_col(position = position_dodge(width = dodge_width)) +
    geom_point(aes(y = Concentration_REDUCED_NORMALIZED, shape = as.factor(DILUTION)), 
               position = position_jitterdodge(dodge.width = dodge_width, jitter.width = jitter_width), 
               size = 3,
               show.legend = T) +
    geom_errorbar(aes(ymin = triplicate_mean_per_day - triplicate_sd_per_day,
                      ymax = triplicate_mean_per_day + triplicate_sd_per_day),
                  position = position_dodge(width = dodge_width)) +
    facet_wrap(~paste0(STIM_DAY, "h"), nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  DR2_PLOT_2 <- ggplot(DATA_NORMALIZED, aes(x = CONDITION, y = Concentration_REDUCED_NORMALIZED, fill = CELL_LINE, group = CELL_LINE)) +
    geom_col(position = position_dodge(width = dodge_width)) +
    facet_wrap(~paste0(STIM_DAY, "h")+DILUTION, nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  DR2_PLOT_3 <- ggplot(DATA_NORMALIZED, aes(x = CONDITION, y = triplicate_mean_per_day, fill = STIM_DAY, group = STIM_DAY)) +
    geom_col(position = position_dodge(width = 1)) +
    facet_wrap(~CELL_LINE, nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  DR2_PLOT_4 <- ggplot(DATA_NORMALIZED, aes(x = STIM_DAY, y = triplicate_mean_per_day, fill = CONDITION, group = CONDITION)) +
    geom_col(position = position_dodge(width = 1)) +
    facet_wrap(~CELL_LINE+CONDITION, nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  return(list(FILTERED_SUBSET = FILTERED_SUBSET, 
              DATA_NORMALIZED = DATA_NORMALIZED, 
              DR2_PLOT_1 = DR2_PLOT_1, 
              DR2_PLOT_2 = DR2_PLOT_2, 
              DR2_PLOT_3 = DR2_PLOT_3, 
              DR2_PLOT_4 = DR2_PLOT_4))
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

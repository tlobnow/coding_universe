library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, ggforce, ggbreak, patchwork, lemon)

################################################################################

run_settings_and_prep     <- TRUE
run_processing_and_subset <- TRUE

if (run_settings_and_prep) {
  
  # GENERAL SETTINGS
  figure <- "Figure Amyloids"
  SAVE   <- TRUE
  
  RELATIVE_SECRETION      <- TRUE
  RELATIVE_SECRETION_ZOOM <- F
  REAL_SECRETION          <- F
  FOLD_CHANGE             <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 6
  TEXT   <- 8
  
  
  # GATHER DATA & FUNCTIONS
  Input_Directory <- ifelse(dir.exists(file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)),
                            file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure), 
                            file.path("~/Desktop/Source files/", figure))
  
  # NAME_KEY  <- fread("https://raw.githubusercontent.com/tlobnow/coding_universe/main/ELISA_CL_KEY.csv", header = T) 
  NAME_KEY  <- fread("~/Documents/Github/coding_universe/ELISA_CL_KEY.csv", header = T) 
  
  source("~/Documents/Github/Analysis_Synthethic_Myddosome/FigureScripts/functions.R")
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
}

################################################################################
################################################################################
################################################################################

if (run_processing_and_subset) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  #### EXTRA FILTER ####
  # plate_data <- plate_data %>% filter(STIM_DAY != 5)
  #### EXTRA FILTER ####
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  AMYLOIDS            <- plate_data
  AMYLOIDS_normalized <- process_ELISA_data(DF = AMYLOIDS, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069")
  
  # check completeness
  table(AMYLOIDS_normalized$CELL_LINE, AMYLOIDS_normalized$STIM_DAY)

  plotting_data        <- AMYLOIDS_normalized
  
  # Manipulate data for plotting
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main)
  
  plotting_means$CL_NAME_ON_PLOT = reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
}

if (RELATIVE_SECRETION) {
  
  # statistical analysis based on relative secretion values
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  plotting_stats$CL_NAME_ON_PLOT = reorder(plotting_stats$CL_NAME_ON_PLOT, -plotting_stats$ORDER_NO)
  
  # plotting relative values
  figure_X_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                      xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, 
                                          shape = as.factor(STIM_DAY), color = PLOTTING_COLOR), 
               # col = "black",
               # shape = 21, 
               size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), 
               # show.legend = FALSE
               ) +
    geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.15, y = CL_NAME_ON_PLOT, label = significance), 
              hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) + 
    scale_x_continuous(breaks = seq(from = 0, to = 3.2, by = 0.5), position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom", #c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  print(figure_X_ELISA)
}

if (RELATIVE_SECRETION_ZOOM) {
  
  # statistical analysis based on relative secretion values
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting relative values
  figure_X_ELISA_zoom <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                      xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
               shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.05, y = CL_NAME_ON_PLOT, label = significance), 
              hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
    scale_x_continuous(breaks = seq(from = 0, to = 3.2, by = 0.5), position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom", #c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, 0.3), zoom.data = ifelse(a <= 0.3, NA, FALSE))
  
  print(figure_X_ELISA_zoom)
}


if (SAVE) {
  
  save_to_data_tay <- file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)
  save_to_local    <- file.path("~/Desktop/Source files/", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  
  if (RELATIVE_SECRETION) {
    ggsave(file.path(save_to, "figure_X_ELISA.svg"), plot = figure_X_ELISA, device = "svg", width = 12, height = 7)
  }
  
  if (RELATIVE_SECRETION_ZOOM) {
    ggsave(file.path(save_to, "figure_X_ELISA_zoom.svg"), plot = figure_X_ELISA_zoom, device = "svg", width = 14, height = 21)
  }

  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
  fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt.csv"))
}

################################################################################
################################################################################
################################################################################

subset_Plate_1_20240227 = F

if (subset_Plate_1_20240227) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  Plate_1_20240227    <- plate_data %>% filter(Date == "2024-02-27")
  Plate_1_20240227_normalized <- process_ELISA_data(DF = Plate_1_20240227, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069")
  
  plotting_data        <- Plate_1_20240227_normalized
  
  # Manipulate data for plotting
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main)
  
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "IL2_concentration_Dilution_Factor_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting real values
  figure_X_ELISA_real_subset_Plate_1_20240227 <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                      xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
               shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(IL2_concentration_Dilution_Factor_mean) + 1000, y = CL_NAME_ON_PLOT, label = significance), 
              hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
    scale_x_continuous(position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = "Real IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, 500), zoom.data = ifelse(a <= 500, NA, FALSE))
  
  print(figure_X_ELISA_real_subset_Plate_1_20240227)
  
}

if (FOLD_CHANGE) {
  
  # Select Data
  COHORT_DATA <- plotting_data
  
  # Set negative values to 1
  COHORT_DATA$MEASUREMENT <- ifelse(COHORT_DATA$MEASUREMENT < 0, yes = 0, COHORT_DATA$MEASUREMENT)
  
  # Normalization for fold change from UNSTIM to STIM
  NORMALIZED_TO_CONTROL <- COHORT_DATA %>%
    group_by(Date, STIM_DAY, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PLOTTING_COLOR) %>%
    summarise(mean_per_day = mean(MEASUREMENT),
              Concentration = mean(Concentration))
  
  # First, calculate the mean for the UNSTIM condition separately.
  unstim_means <- COHORT_DATA %>%
    filter(CONDITION == "UNSTIM") %>%
    group_by(CELL_LINE) %>%
    summarise(mean_unstim = ifelse(mean(MEASUREMENT) == 0, 1, mean(MEASUREMENT)))
  
  # Now, join this back to the main dataset.
  COHORT_SUBSET <- COHORT_DATA %>%
    left_join(unstim_means, by = "CELL_LINE")
  
  # Compute the fold change. 
  # The `case_when` logic ensures that the fold change is calculated only for the STIM condition.
  NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
    mutate(fold_change = case_when(CONDITION == "STIM" ~ MEASUREMENT / mean_unstim, TRUE ~ NA_real_ )) %>% # set NA for non-STIM conditions
    ungroup() %>%
    group_by(CELL_LINE, CONDITION, PLOTTING_COLOR) %>%
    mutate(trip_mean = mean(fold_change),
           fold_change_sd = sd(fold_change)) %>%
    unique()
  
  # Helper function to perform a t-test and return p-value and annotation
  perform_ttest <- function(data) {
    ttest_result <- t.test(data$MEASUREMENT[data$CONDITION == "STIM"],
                           data$MEASUREMENT[data$CONDITION == "UNSTIM"],
                           paired = F)
    p_value <- ttest_result$p.value
    annotation <- ifelse(p_value < 0.0001, '****', 
                         ifelse(p_value < 0.001, '***', 
                                ifelse(p_value < 0.01, '**',
                                       ifelse(p_value < 0.05, '*', 'ns'))))
    return(list(p_value = p_value, annotation = annotation))
  }
  
  
  results <- COHORT_SUBSET %>% group_by(CL_NAME_ON_PLOT, Date) %>%
    do({
      data <- .
      fold_change_data <- data %>%
        summarise(fold_change = round(mean(MEASUREMENT[CONDITION == "STIM"]) / mean(MEASUREMENT[CONDITION == "UNSTIM"]), digits = 2))
      ttest_results <- perform_ttest(data)
      cbind(fold_change_data, ttest_results)
    }) %>%
    ungroup() %>%
    left_join(NAME_KEY[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")], relationship = "many-to-many") %>%
    unique()
  
  # Extract annotations for plotting
  annotations        <- results$annotation
  names(annotations) <- results$CL_NAME_ON_PLOT
  
  # Reorder
  # results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$fold_change)
  results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  
  # Plotting the fold changes
  figure_X_fold_change_plot_subset_Plate_1_20240227 <- ggplot(results, aes(x = CL_NAME_ON_PLOT, y = fold_change)) +
    geom_col(aes(col = PLOTTING_COLOR, fill = PLOTTING_COLOR), position = position_dodge(width = 0.8), alpha = 0.3) +
    geom_text(aes(label = annotations), size = TEXT, hjust = -0.3) +
    geom_hline(yintercept = 2, linetype = "longdash") +
    scale_color_identity() +
    scale_fill_identity() +
    labs(x = "", y = "fold change") +
    coord_flip() +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  figure_X_fold_change_plot_subset_Plate_1_20240227
}

subset_Plate_1_20240405 = T

if (subset_Plate_1_20240405) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  Plate_1_20240405    <- plate_data %>% filter(Date == "2024-04-05")
  Plate_1_20240405_normalized <- process_ELISA_data(DF = Plate_1_20240405, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069")
  
  plotting_data        <- Plate_1_20240405_normalized
  
  # Manipulate data for plotting
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main)
  
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "IL2_concentration_Dilution_Factor_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting real values
  figure_X_ELISA_real_subset_Plate_1_20240405 <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                      xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
               shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(IL2_concentration_Dilution_Factor_mean) + 1000, y = CL_NAME_ON_PLOT, label = significance), 
              hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
    scale_x_continuous(position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = "Real IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, 500), zoom.data = ifelse(a <= 500, NA, FALSE))
  
  print(figure_X_ELISA_real_subset_Plate_1_20240405)
  
}

if (FOLD_CHANGE) {
  
  # Select Data
  COHORT_DATA <- plotting_data
  
  # Set negative values to 1
  COHORT_DATA$MEASUREMENT <- ifelse(COHORT_DATA$MEASUREMENT < 0, yes = 0, COHORT_DATA$MEASUREMENT)
  
  # Normalization for fold change from UNSTIM to STIM
  NORMALIZED_TO_CONTROL <- COHORT_DATA %>%
    group_by(Date, STIM_DAY, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PLOTTING_COLOR) %>%
    summarise(mean_per_day = mean(MEASUREMENT),
              Concentration = mean(Concentration))
  
  # First, calculate the mean for the UNSTIM condition separately.
  unstim_means <- COHORT_DATA %>%
    filter(CONDITION == "UNSTIM") %>%
    group_by(CELL_LINE) %>%
    summarise(mean_unstim = ifelse(mean(MEASUREMENT) == 0, 1, mean(MEASUREMENT)))
  
  # Now, join this back to the main dataset.
  COHORT_SUBSET <- COHORT_DATA %>%
    left_join(unstim_means, by = "CELL_LINE")
  
  # Compute the fold change. 
  # The `case_when` logic ensures that the fold change is calculated only for the STIM condition.
  NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
    mutate(fold_change = case_when(CONDITION == "STIM" ~ MEASUREMENT / mean_unstim, TRUE ~ NA_real_ )) %>% # set NA for non-STIM conditions
    ungroup() %>%
    group_by(CELL_LINE, CONDITION, PLOTTING_COLOR) %>%
    mutate(trip_mean = mean(fold_change),
           fold_change_sd = sd(fold_change)) %>%
    unique()
  
  # Helper function to perform a t-test and return p-value and annotation
  perform_ttest <- function(data) {
    ttest_result <- t.test(data$MEASUREMENT[data$CONDITION == "STIM"],
                           data$MEASUREMENT[data$CONDITION == "UNSTIM"],
                           paired = F)
    p_value <- ttest_result$p.value
    annotation <- ifelse(p_value < 0.0001, '****', 
                         ifelse(p_value < 0.001, '***', 
                                ifelse(p_value < 0.01, '**',
                                       ifelse(p_value < 0.05, '*', 'ns'))))
    return(list(p_value = p_value, annotation = annotation))
  }
  
  
  results <- COHORT_SUBSET %>% group_by(CL_NAME_ON_PLOT, Date) %>%
    do({
      data <- .
      fold_change_data <- data %>%
        summarise(fold_change = round(mean(MEASUREMENT[CONDITION == "STIM"]) / mean(MEASUREMENT[CONDITION == "UNSTIM"]), digits = 2))
      ttest_results <- perform_ttest(data)
      cbind(fold_change_data, ttest_results)
    }) %>%
    ungroup() %>%
    left_join(NAME_KEY[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")], relationship = "many-to-many") %>%
    unique()
  
  # Extract annotations for plotting
  annotations        <- results$annotation
  names(annotations) <- results$CL_NAME_ON_PLOT
  
  # Reorder
  # results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$fold_change)
  results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  
  # Plotting the fold changes
  figure_X_fold_change_plot_subset_Plate_1_20240405 <- ggplot(results, aes(x = CL_NAME_ON_PLOT, y = fold_change)) +
    geom_col(aes(col = PLOTTING_COLOR, fill = PLOTTING_COLOR), position = position_dodge(width = 0.8), alpha = 0.3) +
    geom_text(aes(label = annotations), size = TEXT, hjust = -0.3) +
    geom_hline(yintercept = 2, linetype = "longdash") +
    scale_color_identity() +
    scale_fill_identity() +
    labs(x = "", y = "fold change") +
    coord_flip() +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  figure_X_fold_change_plot_subset_Plate_1_20240405
}

subset_Plate_1_20240518 = T

if (subset_Plate_1_20240518) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  Plate_1_20240518    <- plate_data %>% filter(Date == "2024-05-18")
  Plate_1_20240518_normalized <- process_ELISA_data(DF = Plate_1_20240518, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069")
  
  plotting_data        <- Plate_1_20240518_normalized
  
  # Manipulate data for plotting
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main)
  
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "IL2_concentration_Dilution_Factor_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting real values
  figure_X_ELISA_real_subset_Plate_1_20240518 <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                      xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
               shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(IL2_concentration_Dilution_Factor_mean) + 1000, y = CL_NAME_ON_PLOT, label = significance), 
              hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
    scale_x_continuous(position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = "Real IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, 500), zoom.data = ifelse(a <= 500, NA, FALSE))
  
  print(figure_X_ELISA_real_subset_Plate_1_20240518)
  
}

if (FOLD_CHANGE) {
  
  # Select Data
  COHORT_DATA <- plotting_data
  
  # Set negative values to 1
  COHORT_DATA$MEASUREMENT <- ifelse(COHORT_DATA$MEASUREMENT < 0, yes = 0, COHORT_DATA$MEASUREMENT)
  
  # Normalization for fold change from UNSTIM to STIM
  NORMALIZED_TO_CONTROL <- COHORT_DATA %>%
    group_by(Date, STIM_DAY, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PLOTTING_COLOR) %>%
    summarise(mean_per_day = mean(MEASUREMENT),
              Concentration = mean(Concentration))
  
  # First, calculate the mean for the UNSTIM condition separately.
  unstim_means <- COHORT_DATA %>%
    filter(CONDITION == "UNSTIM") %>%
    group_by(CELL_LINE) %>%
    summarise(mean_unstim = ifelse(mean(MEASUREMENT) == 0, 1, mean(MEASUREMENT)))
  
  # Now, join this back to the main dataset.
  COHORT_SUBSET <- COHORT_DATA %>%
    left_join(unstim_means, by = "CELL_LINE")
  
  # Compute the fold change. 
  # The `case_when` logic ensures that the fold change is calculated only for the STIM condition.
  NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
    mutate(fold_change = case_when(CONDITION == "STIM" ~ MEASUREMENT / mean_unstim, TRUE ~ NA_real_ )) %>% # set NA for non-STIM conditions
    ungroup() %>%
    group_by(CELL_LINE, CONDITION, PLOTTING_COLOR) %>%
    mutate(trip_mean = mean(fold_change),
           fold_change_sd = sd(fold_change)) %>%
    unique()
  
  # Helper function to perform a t-test and return p-value and annotation
  perform_ttest <- function(data) {
    ttest_result <- t.test(data$MEASUREMENT[data$CONDITION == "STIM"],
                           data$MEASUREMENT[data$CONDITION == "UNSTIM"],
                           paired = F)
    p_value <- ttest_result$p.value
    annotation <- ifelse(p_value < 0.0001, '****', 
                         ifelse(p_value < 0.001, '***', 
                                ifelse(p_value < 0.01, '**',
                                       ifelse(p_value < 0.05, '*', 'ns'))))
    return(list(p_value = p_value, annotation = annotation))
  }
  
  
  results <- COHORT_SUBSET %>% group_by(CL_NAME_ON_PLOT, Date) %>%
    do({
      data <- .
      fold_change_data <- data %>%
        summarise(fold_change = round(mean(MEASUREMENT[CONDITION == "STIM"]) / mean(MEASUREMENT[CONDITION == "UNSTIM"]), digits = 2))
      ttest_results <- perform_ttest(data)
      cbind(fold_change_data, ttest_results)
    }) %>%
    ungroup() %>%
    left_join(NAME_KEY[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")], relationship = "many-to-many") %>%
    unique()
  
  # Extract annotations for plotting
  annotations        <- results$annotation
  names(annotations) <- results$CL_NAME_ON_PLOT
  
  # Reorder
  # results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$fold_change)
  results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  
  # Plotting the fold changes
  figure_X_fold_change_plot_subset_Plate_1_20240518 <- ggplot(results, aes(x = CL_NAME_ON_PLOT, y = fold_change)) +
    geom_col(aes(col = PLOTTING_COLOR, fill = PLOTTING_COLOR), position = position_dodge(width = 0.8), alpha = 0.3) +
    geom_text(aes(label = annotations), size = TEXT, hjust = -0.3) +
    geom_hline(yintercept = 2, linetype = "longdash") +
    scale_color_identity() +
    scale_fill_identity() +
    labs(x = "", y = "fold change") +
    coord_flip() +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  figure_X_fold_change_plot_subset_Plate_1_20240518
}


if (SAVE) {
  
  save_to_data_tay <- file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)
  save_to_local    <- file.path("~/Desktop/Source files/", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  if (subset_Plate_1_20240227) {
    ggsave(file.path(save_to, "figure_X_ELISA_real_subset_Plate_1_20240227.svg"), plot = figure_X_ELISA_real_subset_Plate_1_20240227, device = "svg", width = 12, height = 7)
    ggsave(file.path(save_to, "figure_X_fold_change_plot_subset_Plate_1_20240227.svg"), plot = figure_X_fold_change_plot_subset_Plate_1_20240227, device = "svg", width = 12, height = 7)
  }
  
  if (subset_Plate_1_20240405) {
    ggsave(file.path(save_to, "figure_X_ELISA_real_subset_Plate_1_20240405.svg"), plot = figure_X_ELISA_real_subset_Plate_1_20240405, device = "svg", width = 12, height = 7)
    ggsave(file.path(save_to, "figure_X_fold_change_plot_subset_Plate_1_20240405.svg"), plot = figure_X_fold_change_plot_subset_Plate_1_20240405, device = "svg", width = 12, height = 7)
  }
  if (subset_Plate_1_20240518) {
    ggsave(file.path(save_to, "figure_X_ELISA_real_subset_Plate_1_20240518.svg"), plot = figure_X_ELISA_real_subset_Plate_1_20240518, device = "svg", width = 12, height = 7)
    ggsave(file.path(save_to, "figure_X_fold_change_plot_subset_Plate_1_20240518.svg"), plot = figure_X_fold_change_plot_subset_Plate_1_20240518, device = "svg", width = 12, height = 7)
  }
}

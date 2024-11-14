################################################################################
### ELISA ANALYSIS SCRIPT ######################################################
################################################################################

library(pacman) ; pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, 
                                 dplyr, cowplot, readxl, scales, knitr, tidyr, ggforce, 
                                 ggbreak, patchwork, lemon, openxlsx, here)

################################################################################

GENERAL_SETTINGS          <- TRUE
RUN_SETTINGS_AND_PREP     <- TRUE
RUN_PROCESSING_AND_SUBSET <- TRUE
RELATIVE_SECRETION        <- TRUE
REAL_SECRETION            <- TRUE
FOLD_CHANGE               <- TRUE

if (GENERAL_SETTINGS) {
  
  # GENERAL SETTINGS
  figure <- "Figure_Amyloids"
  SAVE   <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 6
  TEXT   <- 8
  
  # SET WD
  setwd(paste0("~/Desktop/", figure))
  # if (exists(paste0("/Users/u_lobnow/Desktop/", figure))) {
  #   setwd(paste0("~/Desktop/", figure))
  #   } else {setwd(file.path(here()))}
  print(paste0("Working directory set to ", getwd()))
  
  if (file.exists("x/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")) {
    source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")
  } else {
    source("https://raw.githubusercontent.com/tlobnow/coding_universe/refs/heads/main/scripts/SOTA_FX.R")
  }

  # NAME_KEY
  if (file.exists("x~/Documents/Github/coding_universe/ELISA_CL_KEY.csv")) {
    NAME_KEY  <- fread("~/Documents/Github/coding_universe/ELISA_CL_KEY.csv")
  } else {
    NAME_KEY  <- fread("https://raw.githubusercontent.com/tlobnow/coding_universe/refs/heads/main/ELISA_CL_KEY.csv")
  }

  # GATHER DATA & FUNCTIONS
  Input_Directory  <- file.path(getwd())
  Output_Directory <- file.path(Input_Directory, as.character(Sys.Date()))
  
  if (dir.exists(Input_Directory)) {print(paste0("Data will be loaded from ", Input_Directory))} else {dir.create(Input_Directory, recursive = TRUE)}
  
}

if (RUN_SETTINGS_AND_PREP) {
  
  print("Running settings and preparation.")
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory)
  
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA")
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
}

if (RUN_PROCESSING_AND_SUBSET) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA")
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  ##############################################################################
  # subset and normalize data
  
  AMYLOIDS_20240531 <- plate_data %>% filter(Date == "2024-05-31")
  AMYLOIDS_20240822 <- plate_data %>% filter(Date == "2024-08-22", Plate != 1)
  AMYLOIDS_20241004 <- plate_data %>% filter(Date == "2024-10-04", STIM_DAY > 2)
  AMYLOIDS_20241018 <- plate_data %>% filter(Date == "2024-10-18")
  AMYLOIDS_20241101 <- plate_data %>% filter(Date == "2024-11-01")
  AMYLOIDS_20241107 <- plate_data %>% filter(Date == "2024-11-07")
  AMYLOIDS_20241114 <- plate_data %>% filter(Date == "2024-11-14")
  
  
  AMYLOIDS_20240531_normalized  <- process_ELISA_data(DF = AMYLOIDS_20240531, 
                                                      NEGATIVE_CTRL = "cl204", 
                                                      POSITIVE_CTRL = "cl069")
  
  AMYLOIDS_20240822_normalized  <- process_ELISA_data(DF = AMYLOIDS_20240822, 
                                                      NEGATIVE_CTRL = "cl204", 
                                                      POSITIVE_CTRL = "cl69")
  
  AMYLOIDS_20241004_normalized  <- process_ELISA_data(DF = AMYLOIDS_20241004, 
                                                      NEGATIVE_CTRL = "cl204", 
                                                      POSITIVE_CTRL = "cl069")
  
  AMYLOIDS_20241018_normalized  <- process_ELISA_data(DF = AMYLOIDS_20241018, 
                                                      NEGATIVE_CTRL = "cl204", 
                                                      POSITIVE_CTRL = "cl069")
  
  AMYLOIDS_20241101_normalized  <- process_ELISA_data(DF = AMYLOIDS_20241101, 
                                                      NEGATIVE_CTRL = "cl204", 
                                                      POSITIVE_CTRL = "cl069")
  
  AMYLOIDS_20241107_normalized  <- process_ELISA_data(DF = AMYLOIDS_20241107, 
                                                      NEGATIVE_CTRL = "cl204", 
                                                      POSITIVE_CTRL = "cl069")
  
  AMYLOIDS_20241114_normalized  <- process_ELISA_data(DF = AMYLOIDS_20241114,
                                                      NEGATIVE_CTRL = "cl204",
                                                      POSITIVE_CTRL = "cl069")
  
  ##############################################################################
  # join data frames of interest
  plotting_data <- rbind(#AMYLOIDS_20240822_normalized,
                         # AMYLOIDS_20241004_normalized#,
                         # AMYLOIDS_20241018_normalized,
                         # AMYLOIDS_20241101_normalized,
                         # AMYLOIDS_20241107_normalized
                         AMYLOIDS_20241114_normalized)
  
  # check completeness
  table(plotting_data$CELL_LINE, plotting_data$STIM_DAY)
  
  ##############################################################################
  ##############################################################################
  ##############################################################################

  plotting_data2 <- plotting_data
  plotting_data2$CL_NAME_ON_PLOT <- paste(plotting_data$CL_NAME_ON_PLOT, plotting_data$Date, plotting_data$EXPRESSION_LVL, sep = "_")
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  
  
  # Manipulate plotting_data2 for plotting
  plotting_data_main <- process_data_for_plot(plotting_data2, change_unstim_plt_col = T)
  plotting_means     <- prepare_plotting_means(data = plotting_data_main)
  plotting_means$CL_NAME_ON_PLOT = reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  
  ##############################################################################
  # Snippet to streamline stats for saving one large data frame
  plotting_stats_relative <- run_statistics(plotting_means = plotting_means, x_mean = "Relative_Intensity_mean")
  if (!exists("plotting_stats_relative")) {plotting_stats_relative <- plotting_means}
  plotting_stats_real     <- run_statistics(plotting_means = plotting_means, x_mean = "IL2_concentration_Dilution_Factor_mean")
  if (!exists("plotting_stats_real")) {plotting_stats_real <- plotting_means}
  colnames_relative       <- colnames(plotting_stats_relative) %>% str_remove_all("significance|p_value")
  colnames_relative       <- colnames_relative[colnames_relative != ""]
  plotting_stats          <- left_join(plotting_stats_relative, plotting_stats_real, by = colnames_relative, suffix = c("_relative", "_real"))
}

if (RELATIVE_SECRETION) {
  
  ELISA_PLOT_RELATIVE <- 
    prepare_and_plot(
      plotting_means = plotting_means,
      plotting_stats = plotting_stats_relative,
      x_mean    = "Relative_Intensity_mean", 
      x_sem     = "Relative_Intensity_sem", 
      x_label   = "Relative IL-2 secretion"
    )
  
  # ELISA_PLOT_RELATIVE +
  #   facet_zoom(xlim = c(0, 0.2), zoom.data = ifelse(a <= 0.2, NA, FALSE))
  
}

if (REAL_SECRETION) {
  
  ELISA_PLOT_REAL <- 
    prepare_and_plot(
      plotting_means = plotting_means,
      plotting_stats = plotting_stats_real,
      x_mean    = "IL2_concentration_Dilution_Factor_mean", 
      x_sem     = "IL2_concentration_Dilution_Factor_sem", 
      x_label   = "Real IL-2 secretion"
    )
}

if (SAVE) {
  
  if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
  
  # save figure
  if (RELATIVE_SECRETION)      {ggsave(file.path(Output_Directory, "ELISA_rltv.svg"), plot = ELISA_PLOT_RELATIVE, device = "svg", width = 12, height = 12)}
  if (REAL_SECRETION)          {ggsave(file.path(Output_Directory, "ELISA_real.svg"), plot = ELISA_PLOT_REAL,     device = "svg", width = 12, height = 12)}
  
  # save tables
  # fwrite(plate_data,           file.path(Output_Directory, "ELISA_plate_data.csv"))
  
  plate_data$Date    <- as.character(plate_data$Date)
  plotting_data$Date <- as.character(plotting_data$Date)
  
  excel_data <- list('plate_data'     = plate_data, 
                     'plotting_data'  = plotting_data, 
                     'plotting_means' = plotting_means,
                     'plotting_stats' = plotting_stats)
  
  write.xlsx(excel_data, file = paste0(Output_Directory, "/ELISA_ANALYSIS.xlsx"))
  
}

if (FOLD_CHANGE) {
  
  plotting_data <- read.xlsx(paste0(Output_Directory,"/ELISA_ANALYSIS.xlsx"), sheet = "plotting_data")
  
  NEGATIVE_CONTROL <- "cl204"
  
  ### **************************************************************************
  ### **************************************************************************
  ### **************************************************************************
  
  # Select Data of interest & set negative measured values to zero
  cohort_data             <- plotting_data 
  cohort_data$MEASUREMENT <- ifelse(cohort_data$MEASUREMENT < 0, 0, cohort_data$MEASUREMENT)
  
  # Calculate the mean for the UNSTIM condition separately.
  unstim_means <- cohort_data %>%
    filter(CONDITION == "UNSTIM") %>% group_by(CELL_LINE, CL_NAME_ON_PLOT, STIM_DAY) %>%
    summarise(mean_unstim = ifelse(mean(MEASUREMENT) == 0, 1, mean(MEASUREMENT)), .groups = 'drop')
  
  # Join this back to the main dataset.
  cohort_subset <- cohort_data %>% 
    left_join(unstim_means, relationship = "many-to-many") %>%
    filter(MEASUREMENT != 0)
  # add individual fold change??
  
  prepare_fold_change_plots <- function(plotting_data, name_key) {
    
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
      group_by(CL_NAME_ON_PLOT, Date, PLOT_ID) %>%
      filter(CELL_LINE == "cl204") %>%
      summarise(fold_change_neg_ctrl = round(mean(MEASUREMENT[CONDITION == "UNSTIM"])))
    
    cohort_subset  <- cohort_subset %>% 
      left_join(control_subset_fc)
    
    results <- cohort_subset %>% group_by(CL_NAME_ON_PLOT, Date, PLOT_ID) %>%
      do({
        data <- .
        fold_change_data <- data %>%
          filter(MEASUREMENT != 0) %>%
          summarise(fold_change = round(mean(MEASUREMENT[CONDITION == "STIM"]) / mean(MEASUREMENT[CONDITION == "UNSTIM"]), digits = 2), .groups = 'drop')
        #ttest_results <- perform_ttest(data)
        #cbind(fold_change_data, ttest_results)
      }) %>%
      ungroup() %>%
      left_join(name_key[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")],
               by = c("CL_NAME_ON_PLOT" = "CL_NAME_ON_PLOT"), relationship = "many-to-many") %>%
      distinct(CL_NAME_ON_PLOT, fold_change, .keep_all = T)
    
    # Extract annotations for plotting
    # annotations <- results$annotation
    # names(annotations) <- results$CL_NAME_ON_PLOT
    
    # Reorder by your preferred visualization order
    results$CL_NAME_ON_PLOT <- reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
    
    return(list(results = results,
                cohort_subset = cohort_subset#,
                #annotations = annotations
    ))
  }
  
  ### **************************************************************************
  ### **************************************************************************
  ### **************************************************************************
  
  fold_change_data <- prepare_fold_change_plots(plotting_data = plotting_data, 
                                                name_key = NAME_KEY)
  
  
  
  results     <- fold_change_data$results %>% distinct(CL_NAME_ON_PLOT, fold_change, .keep_all = T)
  # annotations <- fold_change_data$annotations
  
  # Reorder by your preferred visualization order
  # results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$fold_change)
  results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  results$PLOT_ID = reorder(results$PLOT_ID, -results$ORDER_NO)
  
  # Plotting the fold changes
  fold_change <- ggplot(results, aes(x = PLOT_ID, y = fold_change)) +
    geom_col(aes(col = PLOTTING_COLOR, fill = PLOTTING_COLOR), position = position_dodge(width = 0.8), alpha = 0.3) +
    # geom_text(aes(label = annotations), size = TEXT, hjust = -0.3) +
    geom_hline(yintercept = 2, linetype = "longdash") +
    scale_color_identity() +
    scale_fill_identity() +
    labs(x = "", y = "fold change") +
    coord_flip() +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position.inside = c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  fold_change
  
  # add a plot that is faceted by assay date
  fc_date <- fold_change  +
    facet_wrap(~Date, scales = "free_y", ncol = 1)
  
  fc_date
  
  if (SAVE) {
    
    if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
    
    results$Date    <- as.character(results$Date)
    
    excel_data <- list('plate_data'     = plate_data, 
                       'plotting_data'  = plotting_data, 
                       'plotting_means' = plotting_means,
                       'plotting_stats' = plotting_stats,
                       'fold_change'    = results)
    
    write.xlsx(excel_data, file = paste0(Input_Directory, "/ELISA_ANALYSIS.xlsx"))
    
    # save figures
    ggsave(file.path(Output_Directory, "fold_change.svg"),      plot = fold_change, device = "svg", width = 12, height = 7)
    ggsave(file.path(Output_Directory, "fold_change_date.svg"), plot = fc_date, device = "svg", width = 12, height = 15)
    # ggsave(file.path(Output_Directory, "fold_change_date_facet.svg"), plot = fc_date,     device = "svg", width = 12, height = 7)
  }
  
  fc_date
}



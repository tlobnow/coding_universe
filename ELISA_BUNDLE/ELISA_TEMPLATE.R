################################################################################
### ONE ELISA SCRIPT TO RULE THEM ALL ##########################################
################################################################################
### This script is designed to be a one-stop-shop for all ELISA data analysis.

### The script is divided into sections, each of which can be run independently.
### It is designed to be run from top to bottom, but you can run any
### section independently by setting the corresponding flag to TRUE.

### We will look at data sets from two different days that used different positive controls
### The negative control is the same for all plates (cl204).

### Dataset 1: 2024-03-04, plates 1 and 2, positive control cl028
### Dataset 2: 2024-03-04, plate 3, positive control cl069
### Dataset 3: 2023-07-04, plate 3, positive control cl028

# load necessary libraries
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
  
  # determine the path of this script (you should not need to change this)
  MAIN   <- dirname(rstudioapi::getSourceEditorContext()$path)
  SAVE   <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 6
  TEXT   <- 8
  
  # Set the working directory to the location of the script or provide your own location of the input data
  if (dir.exists(MAIN)) {setwd(MAIN) ; print(paste0("Working directory set to ", getwd()))} else {print("No working directory set.")}
  # setwd("your/path/to/input/data")
  
  # Load the SOTA_FX.R script from the provided link or provide the path to your own script
  try(source("https://raw.githubusercontent.com/tlobnow/coding_universe/refs/heads/main/scripts/SOTA_FX.R"))
  # try(source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R"))
  
  # Use the provided link to read the cell line key or supply the name and path of your own key file
  try(NAME_KEY  <- fread("https://raw.githubusercontent.com/tlobnow/coding_universe/refs/heads/main/ELISA_CL_KEY.csv"))
  # NAME_KEY <- fread("/Users/u_lobnow/Documents/Github/coding_universe/ELISA_CL_KEY.csv")
  # NAME_KEY <- ("your/path/to/ELISA_CL_KEY.csv")
  
  # GATHER DATA & FUNCTIONS
  Input_Directory  <- file.path(getwd(), "example_data")
  # Input_Directory  <- file.path(getwd(), "input")
  Output_Directory <- file.path(MAIN, paste0("output_", as.character(Sys.Date())))
  
  if (dir.exists(Input_Directory))  {print(paste0("Data will be loaded from ", Input_Directory))}
  if (dir.exists(Output_Directory)) {print(paste0("Data will be saved to ",    Output_Directory))} else {dir.create(Output_Directory, recursive = TRUE)}  
}

if (RUN_SETTINGS_AND_PREP) {
  
  print("Running settings and preparation.")
  
  # process raw data into a usable format and save the standard curves in PDF format
  plate_data_raw <- ELISA_Fx(Input_Directory)
  
  # join the raw data with the cell line key
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many") %>% filter(CELL_LINE != "NA")
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
}

if (RUN_PROCESSING_AND_SUBSET) {
  
  ##############################################################################
  ## 1 ADD YOUR DATA & NORMALIZE IT ############################################
  ##############################################################################
  # subset and normalize data to your needs
  
  # This is a way to subset the data by date and conditions
  # You can use this to compare different plates or conditions
  
  data_20240304_p1_p2 <- plate_data %>% filter(Date == "2024-03-04",
                                               ### only include plates 1 and 2 (plate 3 was from a different experiment)
                                               Plate %in% c(1, 2))
  
  data_20240304_p3 <- plate_data %>% filter(Date == "2024-03-04", 
                                             ### only include plates 1 and 2 (plate 3 was from a different experiment)
                                             Plate %in% c(3),
                                             ### exclude outliers
                                             MEASUREMENT != 0.1087, 
                                             ### subset specific cell lines
                                             CL_NUMBER %in% c("cl321", "cl204", "cl069"))
  
  data_20230704 <- plate_data %>% filter(Date == "2023-07-04")
  
  # Normalize the data based on the provided negative and positive controls
  # Check out the cell line key to understand the naming convention
      # controls as specified in the CL_NUMBER column,       e.g., cl204 and cl028
      # controls as specified in the CL_NAME_ON_PLOT column, e.g., tKO_EL4 and 3E10
      # controls as specified in the raw CELL_LINE column,   e.g. cl204 and 3E_10 (or whatever you specified in the input excel file CELL_LINES sheet)
  
  normalized_data_20240304_p1_p2 <- process_ELISA_data(DF = data_20240304_p1_p2, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  normalized_data_20240304_p3    <- process_ELISA_data(DF = data_20240304_p3,    NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl069")
  normalized_data_20230704       <- process_ELISA_data(DF = data_20230704,       NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")

  ##############################################################################
  ## 2 BIND THE DATA OF INTEREST ###############################################
  ##############################################################################
  # join data frames of interest 
  # I would always include all subsets and hash out (#) what you don't need at the moment
  # this way you can always go back and check the data
  # Go ahead and play around to see how the data changes when you include or exclude certain subsets
  # !!! If you hash out some lines, make sure that there are no empty arguments (don't forget to remove unnecessary commas after the last argument in the function call) !!!
  
  plotting_data <- rbind(
    normalized_data_20240304_p1_p2,
    # normalized_data_20240304_p3
    normalized_data_20230704
    )
  
  # check completeness per cell line and stimulation day (ideally 3 replicates per condition)
  table(plotting_data$CELL_LINE, plotting_data$STIM_DAY, plotting_data$CONDITION)
  
  # this pastes the cell line name with expression levels and date for better plotting!
  plotting_data2 <- plotting_data
  plotting_data2$CL_NAME_ON_PLOT <- paste(plotting_data$CL_NAME_ON_PLOT, plotting_data$Date, plotting_data$EXPRESSION_LVL, sep = "_")
  
  # Manipulate plotting_data2 for plotting
  plotting_data_main <- process_data_for_plot(plotting_data2, change_unstim_plt_col = T)
  
  # calculate means for plotting and reorder the cell lines based on the order in the key file 
  # >> higher order number ~ higher on the plot
  plotting_means <- prepare_plotting_means(plotting_data_main) ; glimpse(plotting_means)
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  
  ##############################################################################
  # Snippet to streamline stats for saving one large data frame
  try(rm(plotting_stats_relative))
  try(plotting_stats_relative <- run_statistics(plotting_means = plotting_means, x_mean = "Relative_Intensity_mean"))
  if (!exists("plotting_stats_relative")) {plotting_stats_relative <- plotting_means}
  try(rm(plotting_stats_real))
  try(plotting_stats_real     <- run_statistics(plotting_means = plotting_means, x_mean = "IL2_concentration_Dilution_Factor_mean"))
  if (!exists("plotting_stats_real")) {plotting_stats_real <- plotting_means}
  colnames_relative       <- colnames(plotting_stats_relative) %>% str_remove_all("significance|p_value")
  colnames_relative       <- colnames_relative[colnames_relative != ""]
  plotting_stats_relative$CL_NAME_ON_PLOT <- reorder(plotting_stats_relative$CL_NAME_ON_PLOT, -plotting_stats_relative$ORDER_NO)
  plotting_stats_real$CL_NAME_ON_PLOT <- reorder(plotting_stats_real$CL_NAME_ON_PLOT, -plotting_stats_real$ORDER_NO)
  plotting_stats          <- left_join(plotting_stats_relative, plotting_stats_real, by = colnames_relative, suffix = c("_relative", "_real"))
}

if (RELATIVE_SECRETION) {
  
  print("Running relative secretion.")
  
  ## easiest use case:
  
  # ELISA_PLOT_RELATIVE <- 
  #   prepare_and_plot(
  #     plotting_means = plotting_means,
  #     plotting_stats = plotting_stats_relative,
      # x_mean    = "Relative_Intensity_mean",
  #     x_sem     = "Relative_Intensity_sem", 
  #     x_label   = "Relative IL-2 secretion"
  #   )
  
  ## add some zoom to the plot if you want to focus on a specific range (a)
  # zoom_max <- 0.2 ; ELISA_PLOT_RELATIVE + facet_zoom(xlim = c(0, zoom_max), zoom.data = ifelse(a <= zoom_max, NA, FALSE))

  ## or plot the data however you like
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats_relative$CL_NAME_ON_PLOT <- reorder(plotting_stats_relative$CL_NAME_ON_PLOT, -plotting_stats_relative$ORDER_NO)

  ELISA_PLOT_RELATIVE <- ggplot(data = plotting_stats_relative, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, 
             alpha = 0.5) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR),
               col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "Relative IL-2 secretion", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    scale_color_manual(name = "Relative IL-2 secretion", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) ; ELISA_PLOT_RELATIVE

  if ("significance" %in% colnames(plotting_stats_relative)) {
    ELISA_PLOT_RELATIVE <- ELISA_PLOT_RELATIVE +
      geom_errorbar(aes(xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                        xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_text(data = plotting_stats_relative, aes(x = 1.2 * max(Relative_Intensity_mean), y = CL_NAME_ON_PLOT, label = significance),
                hjust = .5, vjust = 1, size = TEXT,
                angle = case_when(plotting_stats_relative$significance == "ns" ~ 0, T ~ 90)
      )
  } ; ELISA_PLOT_RELATIVE

  if (length(unique(plotting_means$POSITIVE_CTRL)) > 1) {
    ELISA_PLOT_RELATIVE <- ELISA_PLOT_RELATIVE +
      facet_wrap(~POSITIVE_CTRL, scales = "free", ncol = 1)
  } ; ELISA_PLOT_RELATIVE
}

if (REAL_SECRETION) {
  
  print("Running real secretion.")
  
  ELISA_PLOT_REAL <- 
    prepare_and_plot(
      plotting_means = plotting_means,
      plotting_stats = plotting_stats_real,
      x_mean    = "IL2_concentration_Dilution_Factor_mean", 
      x_sem     = "IL2_concentration_Dilution_Factor_sem", 
      x_label   = "Real IL-2 secretion"
    )
  
  # # ## add zoom to the plot if you want to focus on a specific range (a)
  # # zoom_max <- 0.2 ; ELISA_PLOT_REAL_REAL + facet_zoom(xlim = c(0, zoom_max), zoom.data = ifelse(a <= zoom_max, NA, FALSE))
  # 
  # # ## or plot the data however you like
  # plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  # plotting_stats_real$CL_NAME_ON_PLOT <- reorder(plotting_stats_real$CL_NAME_ON_PLOT, -plotting_stats_real$ORDER_NO)
  # 
  # ELISA_PLOT_REAL <- ggplot(data = plotting_stats_real, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
  #   geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
  #   # geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  #   geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
  #              col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
  #   scale_y_discrete(expand = c(0, 0)) +
  #   scale_fill_manual(name  = "Real IL-2 secretion", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
  #   scale_color_manual(name = "Real IL-2 secretion", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
  #   labs(x = "Real IL-2 secretion", y = "") +
  #   guides(color = "none", fill = guide_legend(reverse = TRUE)) +
  #   theme_cowplot(font_size = SIZE, font_family = FONT) +
  #   theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
  #         axis.title.y      = element_blank(),
  #         legend.position   = "bottom",
  #         legend.title      = element_blank(),
  #         legend.text       = element_text(size = SIZE),
  #         legend.key.size   = unit(9, "mm")) ; ELISA_PLOT_REAL
  # 
  # if ("significance" %in% colnames(plotting_stats_real)) {
  #   ELISA_PLOT_REAL <- ELISA_PLOT_REAL +
  #     geom_errorbar(aes(xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
  #                       xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
  #                   linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
  #     geom_text(data = plotting_stats_relative, aes(x = 1.2 * max(IL2_concentration_Dilution_Factor_mean), y = CL_NAME_ON_PLOT, label = significance), 
  #               hjust = .5, vjust = 1, size = TEXT, 
  #               #angle = case_when(plotting_stats_relative$significance == "ns" ~ 0, T ~ 90)
  #     )
  # } ; ELISA_PLOT_REAL
  
}

if (SAVE) {
  
  if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
  
  # save figure
  if (RELATIVE_SECRETION)      {ggsave(file.path(Output_Directory, "ELISA_rltv.svg"), plot = ELISA_PLOT_RELATIVE, device = "svg", width = 12, height = 12)}
  if (REAL_SECRETION)          {ggsave(file.path(Output_Directory, "ELISA_real.svg"), plot = ELISA_PLOT_REAL,     device = "svg", width = 12, height = 12)}
  
  # save tables
  if (RUN_SETTINGS_AND_PREP & RUN_PROCESSING_AND_SUBSET) {
    plate_data$Date    <- as.character(plate_data$Date)
    plotting_data$Date <- as.character(plotting_data$Date)
    
    excel_data <- list('plate_data'     = plate_data, 
                       'plotting_data'  = plotting_data, 
                       'plotting_means' = plotting_means,
                       'plotting_stats' = plotting_stats)
    
    write.xlsx(excel_data, file = paste0(Output_Directory, "/ELISA_ANALYSIS.xlsx"))
  }
}

if (FOLD_CHANGE) {
  
  print("Running fold change.")
  
  # if plotting_data is not found in the environment, load it from the ELISA_ANALYSIS.xlsx file
  if (!exists("plotting_data")) {plotting_data <- read.xlsx(paste0(Output_Directory,"/ELISA_ANALYSIS.xlsx"), sheet = "plotting_data")}
  
  fold_change_data <- plotting_data %>%
    group_by(CL_NUMBER, STIM_DAY, PLOTTING_COLOR, PLOT_ID) %>%
    distinct(fold_change, .keep_all = T) %>%
    ungroup() %>%
    
    ###########
    ###########

    ## There can be a special case if the controls were generated in separate experiments, on separate days, or even by separate people, perhaps with different control cell batches..
    ## Usually you can simply group by Date, but if that's insufficient because both experiments were done on the same day
    ## you can unhash out the following two lines and adjust the experimental details to your needs:
  
    # mutate(EXPERIMENT = ifelse(Plate %in% c(1, 2), "exp1", "exp2")) %>%
    # group_by(EXPERIMENT) %>%
  
    ###########
    ###########
  
    group_by(CL_NUMBER, Date) %>%
      mutate(fold_change_mean = mean(fold_change),
             fold_change_sem  = sem(fold_change))
  
  # reorder the cell lines based on the order in the key file
  fold_change_data$CL_NAME_ON_PLOT <- reorder(fold_change_data$CL_NAME_ON_PLOT, -fold_change_data$ORDER_NO)
  
  FC_ELISA_PLOT <- ggplot(data = fold_change_data, aes(x = fold_change_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR)) +
    geom_col(data = fold_change_data %>% distinct(CL_NUMBER, .keep_all = T), aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
    geom_vline(xintercept = 2, linetype = "longdash") +
    geom_point(data = fold_change_data, aes(x = fold_change, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR),
               col = "black", shape = 21, size = POINTS, alpha = 0.7,
               position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "Fold Change", values = fold_change_data$PLOTTING_COLOR, breaks = fold_change_data$PLOTTING_COLOR, labels = ifelse(fold_change_data$CONDITION == "UNSTIM", paste0("- ", fold_change_data$STIMULANT), paste0("+ ", fold_change_data$STIMULANT))) +
    scale_color_manual(name = "Fold Change", values = fold_change_data$PLOTTING_COLOR, breaks = fold_change_data$PLOTTING_COLOR, labels = ifelse(fold_change_data$CONDITION == "UNSTIM", paste0("- ", fold_change_data$STIMULANT), paste0("+ ", fold_change_data$STIMULANT))) +
    labs(x = "Fold Change", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")
          )
  
  if ("fc_significance" %in% colnames(fold_change_data)) {
    FC_ELISA_PLOT <- FC_ELISA_PLOT +
      geom_errorbar(aes(xmin = fold_change_mean - fold_change_sem,
                        xmax = fold_change_mean + fold_change_sem),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_text(data = fold_change_data, aes(x = 1.2 * max(fold_change_mean), y = CL_NAME_ON_PLOT, label = fc_significance),
                hjust = .5, vjust = 1, size = TEXT)
  } ; FC_ELISA_PLOT
  
  if (SAVE) {
    ggsave(file.path(Output_Directory, "ELISA_fold_change.svg"), plot = FC_ELISA_PLOT, device = "svg", width = 12, height = 12)
  }
  
  FC_ELISA_PLOT
}

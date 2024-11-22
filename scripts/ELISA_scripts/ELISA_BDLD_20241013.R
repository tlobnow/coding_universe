################################################################################
### ELISA bDD SCRIPTS ##########################################################
################################################################################

pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, 
               dplyr, cowplot, readxl, scales, knitr, tidyr,
               ggforce, ggbreak, patchwork, lemon, openxlsx)

################################################################################

GENERAL_SETTINGS          <- TRUE
RUN_SETTINGS_AND_PREP     <- TRUE
RUN_PROCESSING_AND_SUBSET <- TRUE
RELATIVE_SECRETION        <- TRUE
REAL_SECRETION            <- TRUE
FOLD_CHANGE               <- TRUE

if (GENERAL_SETTINGS) {
  # GENERAL SETTINGS
  figure <- "Figure bDDs"
  SAVE   <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 6
  TEXT   <- 8
  
  # SET WD
  setwd(paste0("~/Desktop/", figure)) ; getwd()
  
  # source("https://raw.githubusercontent.com/tlobnow/coding_universe/main/scripts/SOTA_FX.R")
  source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")
  
  # GATHER DATA & FUNCTIONS
  Input_Directory  <- file.path(getwd())
  Output_Directory <- file.path(Input_Directory, as.character(Sys.Date()))
  
  if (dir.exists(Input_Directory)) {print(paste0("Data will be loaded from ", Input_Directory))} else {dir.create(Input_Directory, recursive = TRUE)}
  
}

if (RUN_SETTINGS_AND_PREP) {
  
  print("Running settings and preparation.")
  
  # NAME_KEY  <- fread("ELISA_CL_KEY.csv", header = T) 
  NAME_KEY  <- fread("~/Documents/Github/coding_universe/ELISA_CL_KEY.csv")
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory)
  
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA")
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
}

if (RUN_PROCESSING_AND_SUBSET) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") #%>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  ##############################################################################
  # subset and normalize data
  
  BDLD_06 <- plate_data %>% filter(Date == "2023-07-04" & Plate == 3, CL_NUMBER %in% c("cl328", "cl204", "cl028"))
  BDLD_06_normalized <- process_ELISA_data(DF = BDLD_06, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_57            <- plate_data %>% filter(Date == "2024-03-04", MEASUREMENT != 0.1087, CL_NUMBER %in% c("cl321", "cl204", "cl069"))
  BDLD_57_normalized <- process_ELISA_data(DF = BDLD_57, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl069")
  
  BDLD_87_88            <- plate_data %>% filter(Date == "2024-10-04", STIM_DAY %in% c(1,2), CL_NUMBER %in% c("cl518", "cl519", "cl520", "cl521", "cl204", "cl069"))
  BDLD_87_88_normalized <- process_ELISA_data(DF = BDLD_87_88, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl069")
  
  BDLD_27_20230504            <- plate_data %>% filter(Date == "2023-05-04", CL_NUMBER %in% c("cl310", "cl311", "cl011", "cl028", "cl069", "cl204"))
  BDLD_27_20230504_normalized <- process_ELISA_data(DF = BDLD_27_20230504, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")

  BDLD_27_20230512            <- plate_data %>% filter(Date == "2023-05-12", CL_NUMBER %in% c("cl310", "cl311", "cl011", "cl028", "cl069", "cl204"))
  BDLD_27_20230512_normalized <- process_ELISA_data(DF = BDLD_27_20230512, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")

  BDLD_27_20230529            <- plate_data %>% filter(Date == "2023-05-29", CL_NUMBER %in% c("cl310", "cl311", "cl011", "cl028", "cl069", "cl204", "cl320"))
  BDLD_27_20230529_normalized <- process_ELISA_data(DF = BDLD_27_20230529, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_14_20230504            <- plate_data %>% filter(Date == "2023-05-04", CL_NUMBER %in% c("cl314", "cl315", "cl311", "cl011", "cl028", "cl069", "cl204"))
  BDLD_14_20230504_normalized <- process_ELISA_data(DF = BDLD_14_20230504, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_14_20230512            <- plate_data %>% filter(Date == "2023-05-12", CL_NUMBER %in% c("cl314", "cl315", "cl011", "cl028", "cl069", "cl204"))
  BDLD_14_20230512_normalized <- process_ELISA_data(DF = BDLD_14_20230512, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_14_20230529            <- plate_data %>% filter(Date == "2023-05-29", CL_NUMBER %in% c("cl314", "cl315", "cl011", "cl028", "cl069", "cl204", "cl320"))
  BDLD_14_20230529_normalized <- process_ELISA_data(DF = BDLD_14_20230529, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_50_20230504            <- plate_data %>% filter(Date == "2023-05-04", CL_NUMBER %in% c("cl312", "cl313", "cl011", "cl028", "cl069", "cl204"))
  BDLD_50_20230504_normalized <- process_ELISA_data(DF = BDLD_50_20230504, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_50_20230512            <- plate_data %>% filter(Date == "2023-05-12", CL_NUMBER %in% c("cl312", "cl313", "cl011", "cl028", "cl069", "cl204"))
  BDLD_50_20230512_normalized <- process_ELISA_data(DF = BDLD_50_20230512, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  BDLD_50_20230529            <- plate_data %>% filter(Date == "2023-05-29", CL_NUMBER %in% c("cl312", "cl313", "cl011", "cl028", "cl069", "cl204", "cl320"))
  BDLD_50_20230529_normalized <- process_ELISA_data(DF = BDLD_50_20230529, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "3E10_GFP")
  
  CHARMS              <- plate_data %>% filter(Date == "2022-06-23")
  CHARMS_normalized   <- process_ELISA_data(CHARMS, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028") %>% filter(CL_NAME_ON_PLOT %in% c("CHARMS", "3E10"))

  BDLD_11            <- plate_data %>% filter(Date == "2024-03-04" & Plate != 3, CL_NUMBER %in% c("cl332", "cl204", "cl028"))
  BDLD_11_normalized <- process_ELISA_data(BDLD_11, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  
  BDLD_58            <- plate_data %>% filter(Date == "2024-03-04" & Plate != 3, CL_NUMBER %in% c("cl333", "cl204", "cl028"))
  BDLD_58_normalized <- process_ELISA_data(BDLD_58, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  
  BDLD_38            <- plate_data %>% filter(Date == "2024-03-04" & Plate != 3, CL_NUMBER %in% c("cl338", "cl204", "cl028"))
  BDLD_38_normalized <- process_ELISA_data(BDLD_38, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  
  BDLD_69            <- plate_data %>% filter(Date == "2023-07-04" & Plate != 3, CL_NUMBER %in% c("cl326", "cl204", "cl028"))
  BDLD_69_normalized <- process_ELISA_data(BDLD_69, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  
  BDLD_62            <- plate_data %>% filter(Date == "2023-06-21", CL_NUMBER %in% c("cl318", "cl204", "cl028"))
  BDLD_62_normalized <- process_ELISA_data(BDLD_62, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  
  BDLD_67            <- plate_data %>% filter(Date == "2023-07-04" & Plate == 3, CL_NUMBER %in% c("cl327", "cl204", "cl028"))
  BDLD_67_normalized <-process_ELISA_data(BDLD_67, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  
  # join data frames of interest
  plotting_data        <- rbind(
                                BDLD_87_88_normalized,
                                BDLD_57_normalized,
                                BDLD_06_normalized,
                                BDLD_27_20230504_normalized,
                                BDLD_27_20230512_normalized,
                                BDLD_27_20230529_normalized,
                                BDLD_14_20230504_normalized,
                                BDLD_14_20230512_normalized,
                                BDLD_14_20230529_normalized,
                                BDLD_50_20230504_normalized,
                                BDLD_50_20230512_normalized,
                                BDLD_50_20230529_normalized,
                                CHARMS_normalized,
                                BDLD_11_normalized,
                                BDLD_58_normalized,
                                BDLD_38_normalized,
                                BDLD_69_normalized,
                                BDLD_62_normalized,
                                BDLD_67_normalized
                                )

  # *** PROBLEM: Stats calculation based on two stimulations.. ***
  #     Adding a repeat of day 1 or 2, significance is assessed more correctly.
  #     The n is simply too small for significance calculations!
  #     At least 9 values per cell line and stimulation condition required!
  
  # BDLD_87_88_normalized_1 <- BDLD_87_88_normalized %>%
  #   filter(STIM_DAY == 1) %>%
  #   mutate(STIM_DAY = 3)
  # 
  # BDLD_87_88_normalized_2 <- BDLD_87_88_normalized %>%
  #   filter(STIM_DAY == 2) %>%
  #   mutate(STIM_DAY = 3)
  # 
  # BDLD_87_88_normalized <- rbind(BDLD_87_88_normalized, 
  #                                BDLD_87_88_normalized_1
  #                                # BDLD_87_88_normalized_2
  #                                )
  
  # check completeness
  table(plotting_data$CL_NAME_ON_PLOT, plotting_data$STIM_DAY)
  
  ##############################################################################
  
  plotting_data2 <- plotting_data
  plotting_data2$CL_NAME_ON_PLOT <- paste(plotting_data$CL_NAME_ON_PLOT, plotting_data$Date, plotting_data$EXPRESSION_LVL, sep = "_")
  
  ##############################################################################
  
  # Manipulate plotting_data for plotting
  plotting_data_main <- process_data_for_plot(plotting_data, change_unstim_plt_col = T)
  # plotting_data_main <- process_data_for_plot(plotting_data2, change_unstim_plt_col = T)
  plotting_means     <- prepare_plotting_means(data = plotting_data_main, group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "POSITIVE_CTRL", "NEGATIVE_CTRL", "Date"))
  plotting_means$CL_NAME_ON_PLOT = reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  
  ##############################################################################
  # Snippet to streamline stats for saving one large data frame
  plotting_stats_relative <- run_statistics(plotting_means = plotting_means, x_mean = "Relative_Intensity_mean")
  colnames_relative       <- colnames(plotting_stats_relative) %>% str_remove_all("significance|p_value")
  colnames_relative       <- colnames_relative[colnames_relative != ""]
  
  plotting_stats_real     <- run_statistics(plotting_means = plotting_means, x_mean = "IL2_concentration_Dilution_Factor_mean")
  colnames_real           <- colnames(plotting_stats_real) %>% str_remove_all("significance|p_value")
  colnames_real           <- colnames_real[plotting_stats_real != ""]
  
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
  
  # ELISA_PLOT_REAL + 
  #   xlim(c(0, 100))
    # facet_zoom(xlim = c(0, 100), zoom.data = ifelse(a <= 100, NA, FALSE))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  
  ### 
  # 20241115
  # Snippet to streamline stats for saving the ELISA plot for the Lunch Seminar slides
  ###
  
  # plotting_means_MAIN <- plotting_means
  # plotting_stats_real_MAIN <- plotting_stats_real
  # 
  # x_mean    = "IL2_concentration_Dilution_Factor_mean"; 
  # x_sem     = "IL2_concentration_Dilution_Factor_sem"; 
  # x_label   = "Real IL-2 secretion"
  # cl_label = "CL_NAME_ON_PLOT"
  # 
  # 
  # # print if stats has significance column
  # if("significance" %in% colnames(plotting_stats_real_MAIN)){ 
  #   cat("Yep, stats were calculated for ya!\n") 
  # } else { 
  #   cat("Sorry, couldn't calculate stats for ya!\n") 
  # }
  # 
  # # Define the custom order for CELL_LINE
  # custom_order <- c("590_L", 
  #                   "BDLD_6", 
  #                   "BDLD_11", 
  #                   "BDLD_58", 
  #                   "BDLD_57", 
  #                   "BDLD_38", 
  #                   "BDLD_69", 
  #                   "BDLD_62", 
  #                   "BDLD_67",
  #                   "cl518",
  #                   "cl519")
  # 
  # # Filter and set factor levels based on the custom order
  # plotting_means <- plotting_means_MAIN %>% 
  #   filter(CELL_LINE %in% custom_order | CL_NAME_ON_PLOT %in% custom_order) %>%
  #   mutate(CELL_LINE = factor(CELL_LINE, levels = custom_order))
  # 
  # plotting_stats <- plotting_stats_real_MAIN %>% 
  #   filter(CELL_LINE %in% custom_order | CL_NAME_ON_PLOT %in% custom_order) %>%
  #   mutate(CELL_LINE = factor(CELL_LINE, levels = custom_order))
  # 
  # # Your plotting code
  # ELISA_PLOT <- ggplot(data = plotting_stats, aes(x = !!sym(x_mean), y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
  #   geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
  #   scale_y_discrete(expand = c(0, 0)) +
  #   scale_fill_manual(name  = cl_label, values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
  #   scale_color_manual(name = cl_label, values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
  #   labs(x = x_label, y = "") +
  #   guides(color = "none", fill = guide_legend(reverse = TRUE)) +
  #   theme_cowplot(font_size = SIZE, font_family = FONT) +
  #   theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
  #         axis.title.y      = element_blank(),
  #         legend.position   = "bottom",
  #         legend.title      = element_blank(),
  #         legend.text       = element_text(size = SIZE),
  #         legend.key.size   = unit(9, "mm"))
  # 
  # # Add error bars and text if significance column is present
  # if ("significance" %in% colnames(plotting_stats)) {
  #   ELISA_PLOT <- ELISA_PLOT +
  #     geom_errorbar(aes(y = CL_NAME_ON_PLOT,
  #                       xmin = !!sym(x_mean) - !!sym(x_sem),
  #                       xmax = !!sym(x_mean) + !!sym(x_sem)),
  #                   linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
  #     geom_text(data = plotting_stats, aes(x = 1.2 * max(!!sym(x_mean)), y = CL_NAME_ON_PLOT, label = significance), 
  #               hjust = .5, vjust = 1, size = TEXT)
  # }
  # 
  # # Add facet if needed
  # if (length(unique(plotting_means$POSITIVE_CTRL)) > 1 & x_mean == "Relative_Intensity_mean") {
  #   ELISA_PLOT <- ELISA_PLOT +
  #     facet_wrap(~POSITIVE_CTRL, scales = "free", ncol = 1)
  # }
  # 
  # # Print the plot
  # ELISA_PLOT
  # 
  # ggsave("ELISA_PLOT.svg", plot = ELISA_PLOT, width = 8, height = 6, dpi = 300)
  
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
  
  results_plus <- results %>%
    left_join(NAME_KEY[NAME_KEY$CELL_LINE %in% unique(plotting_data$CELL_LINE)]) %>%
    distinct(PLOT_ID, fold_change, .keep_all = T)
  
  results_plus$PLOT_ID = reorder(results_plus$PLOT_ID, -results_plus$ORDER_NO)
  
  # Plotting the fold changes
  fold_change <- ggplot(results_plus, aes(x = PLOT_ID, y = fold_change)) +
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
  
  fold_change_plus <- results_plus %>%
    filter(!CL_NUMBER %in% c("cl204", "cl011", "cl028"),
           ifelse(CL_NUMBER != "cl069", CL_NAME_ON_PLOT %in% unique(results_plus$CL_NAME_ON_PLOT[results_plus$EXPRESSION_LVL == "HIGH"]), TRUE)) %>%
    arrange(#desc(fold_change), 
            # ORIGIN
            ) %>%  # Arrange by ORIGIN, then by fold_change in descending order
    mutate(PLOT_ID = factor(PLOT_ID, levels = unique(PLOT_ID)),
           ORIGIN  = factor(ORIGIN, levels = c("Mus musculus", "Richelia sp. SL_2_1", "Thiotrichaceae bacterium", "Nostoc sp. 106C")),
           STATUS  = case_when(CL_NUMBER %in% c("cl311", "cl320", "cl321") ~ paste0(ORIGIN, "_strict"),
                               CL_NUMBER %in% c("cl518", "cl519", "cl520", "cl521") ~ paste0(ORIGIN, "_extended"),
                               TRUE ~ paste0(ORIGIN, "_", CL_NAME_ON_PLOT))) %>%
    ggplot(aes(x = STATUS, y = fold_change)) +
    geom_col(aes(col = ORIGIN, fill = ORIGIN), position = position_dodge(width = 0.8), alpha = 0.3) +
    geom_hline(yintercept = 2, linetype = "longdash") +
    labs(x = "", y = "fold change") +
    coord_flip() +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          legend.position   = "bottom") +
    scale_color_manual(values = c("#E69F00", "dodgerblue", "lightblue", "lightgreen")) +
    scale_fill_manual(values  = c("#E69F00", "dodgerblue", "lightblue", "lightgreen")) +
    # add axis ticks every 5 units
    scale_y_continuous(breaks = c(2, seq(0, 25, by = 5)))

  fold_change_plus
  
  # add a plot that is faceted by assay date
  fc_date <- fold_change  +
    facet_wrap(~Date, scales = "free_y", ncol = 1)
  
  fc_date
  
  if (SAVE) {
    
    if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
    
    results$Date    <- as.character(results$Date)
    
    excel_data <- list('plate_data'       = plate_data, 
                       'plotting_data'    = plotting_data, 
                       'plotting_means'   = plotting_means,
                       'plotting_stats'   = plotting_stats,
                       'fold_change'      = results,
                       "fold_change_plus" = results_plus)
    
    write.xlsx(excel_data, file = paste0(Input_Directory, "/ELISA_ANALYSIS.xlsx"))
    
    # save figures
    ggsave(file.path(Output_Directory, "fold_change.svg"),      plot = fold_change, device = "svg", width = 12, height = 7)
    ggsave(file.path(Output_Directory, "fold_change_date.svg"), plot = fc_date, device = "svg", width = 12, height = 15)
    # ggsave(file.path(Output_Directory, "fold_change_date_facet.svg"), plot = fc_date,     device = "svg", width = 12, height = 7)
    ggsave(file.path(Output_Directory, "fold_change_plus.svg"), plot = fold_change_plus, device = "svg", width = 12, height = 7)
  }
  
  fc_date
}

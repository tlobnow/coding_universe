################################################################################
### ELISA ANALYSIS SCRIPT FOR AMYLOIDS #########################################
################################################################################

library(pacman) ; pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, 
                                 dplyr, cowplot, readxl, scales, knitr, tidyr, 
                                 ggforce, ggbreak, patchwork, lemon)

################################################################################

run_settings_and_prep     <- TRUE
run_processing_and_subset <- TRUE

if (run_settings_and_prep) {
  
  # GENERAL SETTINGS
  figure <- "Figure Amyloids"
  SAVE   <- TRUE
  
  RELATIVE_SECRETION      <- TRUE
  RELATIVE_SECRETION_ZOOM <- FALSE
  REAL_SECRETION          <- TRUE
  FOLD_CHANGE             <- TRUE
  FOLD_CHANGE_SGL         <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 6
  TEXT   <- 8
  
  # SET WD
  setwd("/Users/u_lobnow/Documents/Github/coding_universe")
  source("scripts/SOTA_FX.R")
  
  NAME_KEY  <- fread("ELISA_CL_KEY.csv", header = T) 

  # GATHER DATA & FUNCTIONS
  Input_Directory <- ifelse(dir.exists(file.path("/Volumes/TAYLOR-LAB/Finn_v2/Other_Project/AMYLOIDS/", figure)),
                            file.path("/Volumes/TAYLOR-LAB/Finn_v2/Other_Project/AMYLOIDS/", figure), 
                            file.path(figure))
  
  Output_Directory <- Input_Directory
  
  if (dir.exists(Input_Directory)) {print(paste0("Data will be loaded from ", Input_Directory))} else {dir.create(Input_Directory, recursive = TRUE)}
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory)
}

################################################################################

if (run_processing_and_subset) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  #### EXTRA FILTER ####
  plate_data <- plate_data %>% filter(Date != "2024-05-18")
  #### EXTRA FILTER ####
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  AMYLOIDS            <- plate_data
  AMYLOIDS_normalized <- process_ELISA_data(DF = AMYLOIDS, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl069")
  
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
  
  if (nrow(stat_significance_dt) > 0) {
    plotting_stats                 <- prepare_plotting_stats(plotting_data, stat_significance_dt)
    plotting_stats$CL_NAME_ON_PLOT <- reorder(plotting_stats$CL_NAME_ON_PLOT, -plotting_stats$ORDER_NO)
    
    # plotting relative values
    Amyloids_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
      geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
      geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                        xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                        xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
                 col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
      geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.15, y = CL_NAME_ON_PLOT, label = significance), 
                hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) + 
      scale_y_discrete(expand = c(0, 0)) +
      scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
      scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
      labs(x = "Relative IL-2 secretion", y = "") +
      guides(color = "none", fill = guide_legend(reverse = TRUE)) +
      theme_cowplot(font_size = SIZE, font_family = FONT) +
      theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
            axis.title.y      = element_blank(),
            legend.position   = "bottom", #c(x = 0.4, y = 0.73),
            legend.title      = element_blank(),
            legend.text       = element_text(size = SIZE),
            legend.key.size   = unit(9, "mm"))
    
    print(Amyloids_ELISA)
    
  } else if (nrow(stat_significance_dt) == 0) {
    print("stat_significance_dt has no row content.")
    print("the plot will be generated without statisticial analysis.")

    Amyloids_ELISA <- ggplot(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
      geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
      geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
                 col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
      scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
      scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
      labs(x = "Relative IL-2 secretion", y = "") +
      guides(color = "none", fill = guide_legend(reverse = TRUE)) +
      theme_cowplot(font_size = SIZE, font_family = FONT) +
      theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
            axis.title.y      = element_blank(),
            legend.position   = "bottom", #c(x = 0.4, y = 0.73),
            legend.title      = element_blank(),
            legend.text       = element_text(size = SIZE),
            legend.key.size   = unit(9, "mm"))
    
    print(Amyloids_ELISA)
  }
}

if (REAL_SECRETION) {
  
  # statistical analysis based on relative secretion values
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "IL2_concentration_Dilution_Factor_mean")
  
  if (nrow(stat_significance_dt) > 0) {
    plotting_stats                 <- prepare_plotting_stats(plotting_data, stat_significance_dt)
    plotting_stats$CL_NAME_ON_PLOT <- reorder(plotting_stats$CL_NAME_ON_PLOT, -plotting_stats$ORDER_NO)
    
    # plotting relative values
    Amyloids_ELISA_real <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
      geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
      geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                        xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                        xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
                 col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
      geom_text(data = plotting_stats, aes(x = max(IL2_concentration_Dilution_Factor_mean) + 0.15, y = CL_NAME_ON_PLOT, label = significance), 
                hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) + 
      scale_y_discrete(expand = c(0, 0)) +
      scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
      scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
      labs(x = "Real IL-2 secretion", y = "") +
      guides(color = "none", fill = guide_legend(reverse = TRUE)) +
      theme_cowplot(font_size = SIZE, font_family = FONT) +
      theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
            axis.title.y      = element_blank(),
            legend.position   = "bottom", #c(x = 0.4, y = 0.73),
            legend.title      = element_blank(),
            legend.text       = element_text(size = SIZE),
            legend.key.size   = unit(9, "mm"))
    
    print(Amyloids_ELISA_real)
    
  } else if (nrow(stat_significance_dt) == 0) {
    print("stat_significance_dt has no row content.")
    print("the plot will be generated without statisticial analysis.")
    
    Amyloids_ELISA_real <- ggplot(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
      geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
      geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
                 col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
      scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
      scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
      labs(x = "Real IL-2 secretion", y = "") +
      guides(color = "none", fill = guide_legend(reverse = TRUE)) +
      theme_cowplot(font_size = SIZE, font_family = FONT) +
      theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
            axis.title.y      = element_blank(),
            legend.position   = "bottom", #c(x = 0.4, y = 0.73),
            legend.title      = element_blank(),
            legend.text       = element_text(size = SIZE),
            legend.key.size   = unit(9, "mm"))
    
    print(Amyloids_ELISA)
  
  }
}

if (RELATIVE_SECRETION_ZOOM) {
  
  # statistical analysis based on relative secretion values
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting relative values
  Amyloids_ELISA_zoom <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
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
  
  print(Amyloids_ELISA_zoom)
}

if (SAVE) {
  
  if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
  
  # save figure
  
  if (RELATIVE_SECRETION)      {ggsave(file.path(Output_Directory, "Amyloids_ELISA_rltv.svg"), plot = Amyloids_ELISA,      device = "svg", width = 12, height = 7)}
  if (RELATIVE_SECRETION_ZOOM) {ggsave(file.path(Output_Directory, "Amyloids_ELISA_zoom.svg"), plot = Amyloids_ELISA_zoom, device = "svg", width = 14, height = 21)}
  if (REAL_SECRETION)          {ggsave(file.path(Output_Directory, "Amyloids_ELISA_real.svg"), plot = Amyloids_ELISA_real, device = "svg", width = 12, height = 7)}

  # save tables
  fwrite(plate_data,           file.path(Output_Directory, "Amyloids_ELISA_plate_data.csv"))
  fwrite(plotting_data,        file.path(Output_Directory, "Amyloids_ELISA_plotting_data.csv"))
  fwrite(plotting_means,       file.path(Output_Directory, "Amyloids_ELISA_plotting_means.csv"))
  fwrite(plotting_stats,       file.path(Output_Directory, "Amyloids_ELISA_plotting_stats.csv"))
  fwrite(stat_significance_dt, file.path(Output_Directory, "Amyloids_ELISA_stat_significance_dt.csv"))
}

################################################################################

if (FOLD_CHANGE) {
  
  plotting_data <- fread(file.path(Output_Directory, "Amyloids_ELISA_plotting_data.csv"))

  fold_change_data <- prepare_fold_change_plots(plotting_data, NAME_KEY)
  # fold_change_data <- prepare_fold_change_plots(plotting_data %>% filter(Date == "2024-05-18"), NAME_KEY)
  # fold_change_data <- prepare_fold_change_plots(plotting_data %>% filter(Date == "2024-05-24"), NAME_KEY)
  
  results     <- fold_change_data$results
  annotations <- fold_change_data$annotations
  
  # Reorder by your preferred visualization order
  # results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$fold_change)
  results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  
  # Plotting the fold changes
  Amyloids_fold_change <- ggplot(results, aes(x = CL_NAME_ON_PLOT, y = fold_change)) +
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
  
  Amyloids_fold_change
  
  # add a plot that is faceted by assay date
  Amyloids_fc_date <- Amyloids_fold_change  +
    facet_wrap(~Date, scales = "free", ncol = 1)
  
  Amyloids_fc_date
  
  if (SAVE) {
    
    if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
    
    # save tables
    fwrite(as.data.frame(results),     file.path(Output_Directory, "Amyloids_fold_change.csv"), quote = F, row.names = F, col.names = T)

    # save figures
    ggsave(file.path(Output_Directory, "Amyloids_fold_change.svg"),            plot = Amyloids_fold_change, device = "svg", width = 12, height = 7)
    ggsave(file.path(Output_Directory, "Amyloids_fold_change_date_facet.svg"), plot = Amyloids_fc_date,     device = "svg", width = 12, height = 7)
  }
}



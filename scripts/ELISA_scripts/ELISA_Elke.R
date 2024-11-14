# DEEPIKA ELISA ANALYSIS

################################################################################
### ELISA ANALYSIS SCRIPT ######################################################
################################################################################

library(pacman) ; pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, ggforce, ggbreak, patchwork, lemon)

################################################################################

RUN_SETTINGS_AND_PREP     <- TRUE
RUN_PROCESSING_AND_SUBSET <- TRUE

if (RUN_SETTINGS_AND_PREP) {
  
  # GENERAL SETTINGS
  figure <- "ELISA_2"
  SAVE   <- TRUE
  
  REAL_SECRETION_WITHOUT_WT  <- TRUE
  REAL_SECRETION_WITH_WT     <- TRUE
  
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
  Input_Directory <- ifelse(dir.exists(file.path("~/Desktop/ELKE/", figure)),
                            file.path("~/Desktop/ELKE/", figure), 
                            file.path(figure))
  
  Output_Directory <- Input_Directory
  
  if (dir.exists(Input_Directory)) {print(paste0("Data will be loaded from ", Input_Directory))} else {dir.create(Input_Directory, recursive = TRUE)}
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory)
}

################################################################################

if (RUN_PROCESSING_AND_SUBSET) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY)
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  ################################################################################
  # subset and normalize TLR data
  
  TLR4      <- plate_data %>% filter(PATHWAY == "TLR4")
  TLR4_data <- process_ELISA_data(DF = TLR4, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  TLR7      <- plate_data %>% filter(PATHWAY == "TLR7")
  TLR7_data <- process_ELISA_data(DF = TLR7, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  TLR9      <- plate_data %>% filter(PATHWAY == "TLR9")
  TLR9_data <- process_ELISA_data(DF = TLR9, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  ################################################################################
  
  # join data frames of interest
  data <- rbind(TLR4_data, TLR7_data, TLR9_data)
  plotting_data_main <- process_data_for_plot(data, change_unstim_plt_col = F)
  
  # add column with CL_NAME_ON_PLOT and PATHWAY
  plotting_data_main$CL_NAME_ON_PLOT_PLUS_PATHWAY = paste0(plotting_data_main$CL_NAME_ON_PLOT, "_", plotting_data_main$PATHWAY)
  
  # check completeness
  check <- as.data.table(table(plotting_data_main$CELL_LINE, plotting_data_main$STIM_DAY, plotting_data_main$PATHWAY, plotting_data_main$CONDITION))

  
  # add stimulant concentration unit
  plotting_data_main <- plotting_data_main %>%
    mutate(Concentration_Unit = case_when(PATHWAY == "TLR4" & STIMULANT == "LPS"  ~ "ng/ml",
                                          PATHWAY == "TLR7" & STIMULANT == "R848" ~ "ng/ml",
                                          PATHWAY == "TLR9" & STIMULANT == "Cpg-B" ~ "µM"),
           # Some TLR7 plates have unstimulated wells that are controls for the 5ng/ml stimulated wells
           # Although we ended up using 50ng/ml for evaluation, the unstimulated controls are still valid controls
           STIM_CONCENTRATION = case_when(PATHWAY == "TLR7" & CONDITION == "UNSTIM" & STIM_CONCENTRATION == 5 ~ 50,
                                          PATHWAY == "TLR9" & CONDITION == "UNSTIM" & STIM_CONCENTRATION == 10 ~ 1,
                                          TRUE ~ STIM_CONCENTRATION),
           PLOTTING_COLOR = case_when(PATHWAY == "TLR4" & CONDITION == "STIM"    ~ "#57bec7",
                                      PATHWAY == "TLR4" & CONDITION == "UNSTIM"  ~ "#ddfcff",
                                      PATHWAY == "TLR7" & CONDITION == "STIM"    ~ "#ea6061",
                                      PATHWAY == "TLR7" & CONDITION == "UNSTIM"  ~ "#fccdcd",
                                      PATHWAY == "TLR9" & CONDITION == "STIM"    ~ "#f0bf37",
                                      PATHWAY == "TLR9" & CONDITION == "UNSTIM"  ~ "#fff5d7",
                                      TRUE ~ "#BEBDBD"
           ))
  
  
  plotting_means_main       <- prepare_plotting_means(data = plotting_data_main, group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "CL_NAME_ON_PLOT_PLUS_PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "Concentration_Unit", "PLOTTING_COLOR", "ORDER_NO"))
  statistical_significance  <- process_statistical_analysis(data = plotting_means_main, group_var = "CL_NAME_ON_PLOT_PLUS_PATHWAY", value_var = "IL2_concentration_Dilution_Factor_mean")
  
  # turn statistical_significance list into data table
  stat_significance_dt <- data.table(
    CL_NAME_ON_PLOT_PLUS_PATHWAY = statistical_significance$CL_NAME_ON_PLOT,
    p_value         = statistical_significance$p_value,
    significance    = statistical_significance$significance)
  
  plotting_stats_main <- plotting_data_main %>%
    group_by(CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PATHWAY, CL_NAME_ON_PLOT_PLUS_PATHWAY, STIMULANT, STIM_CONCENTRATION, Concentration_Unit, PLOTTING_COLOR) %>%
    summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
              IL2_concentration_Dilution_Factor_sem = sem(Concentration),
              Relative_Intensity_mean = mean(triplicate_mean_per_day),
              Relative_Intensity_sem = sem(triplicate_mean_per_day)) %>%
    as.data.table() %>%
    left_join(stat_significance_dt)
  
  plotting_stats_main$CONDITION <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
  
  plotting_stats <- plotting_stats_main
}

if (REAL_SECRETION_WITHOUT_WT) {
  plotting_data  <- plotting_data_main  %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  plotting_means <- plotting_means_main %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  plotting_stats <- plotting_stats_main %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  
  ELISA_real_no_WT <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), col = "black", width = 0.68, alpha = 0.8) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, 
               size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = F) +
    geom_text(data = plotting_stats, aes(x = case_when(PATHWAY == "TLR4" ~ max(plotting_means$IL2_concentration_Dilution_Factor_mean[plotting_means$PATHWAY == "TLR4"]) * 1.5, 
                                                       PATHWAY == "TLR7" ~ max(plotting_means$IL2_concentration_Dilution_Factor_mean[plotting_means$PATHWAY == "TLR7"]) * 1.5, 
                                                       PATHWAY == "TLR9" ~ max(plotting_means$IL2_concentration_Dilution_Factor_mean[plotting_means$PATHWAY == "TLR9"]) * 1.5), 
                                         y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    labs(x = "Real secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom", #c(0.1,0.7),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          # panel.spacing.x = unit(-5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  print(ELISA_real_no_WT)
}

if (REAL_SECRETION_WITH_WT) {
  plotting_data  <- plotting_data_main
  plotting_means <- plotting_means_main
  plotting_stats <- plotting_stats_main
  
  ELISA_real <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = case_when(PATHWAY == "TLR4" ~ max(plotting_means$IL2_concentration_Dilution_Factor_mean[plotting_means$PATHWAY == "TLR4"]) * 1.5, 
                                                       PATHWAY == "TLR7" ~ max(plotting_means$IL2_concentration_Dilution_Factor_mean[plotting_means$PATHWAY == "TLR7"]) * 1.5, 
                                                       PATHWAY == "TLR9" ~ max(plotting_means$IL2_concentration_Dilution_Factor_mean[plotting_means$PATHWAY == "TLR9"]) * 1.5), 
                                         y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    labs(x = "Real secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          # panel.spacing.x = unit(-5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  ELISA_real
}

if (SAVE) {
  
  if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
  
  # save figure
  
  if (REAL_SECRETION_WITH_WT)    {ggsave(file.path(Output_Directory, "ELISA_real.svg"),       plot = ELISA_real, device = "svg", width = 12, height = 7)}
  if (REAL_SECRETION_WITHOUT_WT) {ggsave(file.path(Output_Directory, "ELISA_real_wo_WT.svg"), plot = ELISA_real, device = "svg", width = 12, height = 7)}
  
  # save tables
  fwrite(plate_data,           file.path(Output_Directory, "ELISA_plate_data.csv"))
  fwrite(plotting_data,        file.path(Output_Directory, "ELISA_plotting_data.csv"))
  fwrite(plotting_means,       file.path(Output_Directory, "ELISA_plotting_means.csv"))
  
  if (nrow(stat_significance_dt) == 0) {
    fwrite(plotting_stats,       file.path(Output_Directory, "ELISA_plotting_stats.csv"))
    fwrite(stat_significance_dt, file.path(Output_Directory, "ELISA_stat_significance_dt.csv"))
  }
}

################################################################################

if (FOLD_CHANGE) {
  
  plotting_data <- fread(file.path(Output_Directory, "ELISA_plotting_data.csv"))
  
  fold_change_data <- prepare_fold_change_plots(plotting_data, NAME_KEY)
  # fold_change_data <- prepare_fold_change_plots(plotting_data %>% filter(Date == "2024-05-18"), NAME_KEY)
  # fold_change_data <- prepare_fold_change_plots(plotting_data %>% filter(Date == "2024-05-24"), NAME_KEY)
  
  results     <- fold_change_data$results
  annotations <- fold_change_data$annotations
  
  # Reorder by your preferred visualization order
  # results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$fold_change)
  results$CL_NAME_ON_PLOT = reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  
  # Plotting the fold changes
  fold_change <- ggplot(results, aes(x = CL_NAME_ON_PLOT, y = fold_change)) +
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
  
  fold_change
  
  # add a plot that is faceted by assay date
  fc_date <- fold_change  +
    facet_wrap(~Date, scales = "free", ncol = 1)
  
  fc_date
  
  if (SAVE) {
    
    if (dir.exists(Output_Directory)) {print(paste0("Files will be saved to ", Output_Directory))} else {dir.create(Output_Directory, recursive = T); print(paste0("Files will be saved to ", Output_Directory))}
    
    # save tables
    fwrite(as.data.frame(results),     file.path(Output_Directory, "fold_change.csv"), quote = F, row.names = F, col.names = T)
    
    # save figures
    ggsave(file.path(Output_Directory, "fold_change.svg"),            plot = fold_change, device = "svg", width = 12, height = 7)
    ggsave(file.path(Output_Directory, "fold_change_date_facet.svg"), plot = fc_date,     device = "svg", width = 12, height = 7)
  }
}



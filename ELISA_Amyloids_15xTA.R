library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, ggforce, ggbreak, patchwork, lemon)

################################################################################

# SETTINGS
save   <- TRUE
figure <- "Figure X"

################################################################################

# PREPROCESSING
Input_Directory <- file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)
NAME_KEY  <- fread("https://raw.githubusercontent.com/tlobnow/coding_universe/main/ELISA_CL_KEY.csv", header = T) 

# load functions
source("~/Documents/Github/Analysis_Synthethic_Myddosome/FigureScripts/functions.R")

# process ELISA
plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
plate_data     <- left_join(plate_data_raw, NAME_KEY)
plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()

# ensure correct column type assignment
plate_data$CONDITION <- as.factor(plate_data$CONDITION)

################################################################################

# subset data
AMYLOIDS            <- plate_data
AMYLOIDS_normalized <- process_ELISA_data(DF = AMYLOIDS, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069")

# subset data for CHARMS cell lines
AMYLOIDS_data     <- AMYLOIDS_normalized %>% filter(CL_NAME_ON_PLOT %in% c("CHARMS-10xTA", "CHARMS-12xTA", "CHARMS-15xTA", "cl069"))

# join data frames of interest
plotting_data_main <- AMYLOIDS_data

# reorder cell lines for plotting
plotting_data_main$CL_NAME_ON_PLOT <- reorder(plotting_data_main$CL_NAME_ON_PLOT, -plotting_data_main$ORDER_NO)

# reformat condition for legend text
plotting_data_main$CONDITION  <- as.factor(plotting_data_main$CONDITION)

plotting_data <- plotting_data_main
################################################################################

# calculate means concentrations for plotting
plotting_means <- plotting_data %>%
  group_by(CELL_LINE, CONDITION, STIM_DAY, CL_NAME_ON_PLOT) %>%
  summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
            Relative_Intensity_mean = mean(Concentration_NORMALIZED)) %>%
  as.data.table()

################################################################################
# calculate statistical significance

statistical_significance <- perform_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")

# turn statistical_significance list into data table 
stat_significance_dt <- data.table(CL_NAME_ON_PLOT = names(statistical_significance$annotations), 
                                   p_value         = statistical_significance$annotations,
                                   significance    = statistical_significance$p_values)

################################################################################

# calculate stats
plotting_stats <- plotting_data %>%
  group_by(CELL_LINE, CONDITION, CL_NAME_ON_PLOT) %>%
  summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
            IL2_concentration_Dilution_Factor_sd = sd(Concentration),
            Relative_Intensity_mean = mean(triplicate_mean_per_day),
            Relative_Intensity_sem = sem(triplicate_mean_per_day)) %>%
  as.data.table() %>%
  left_join(stat_significance_dt) %>%
  ungroup()

################################################################################
# plotting settings

COLOR_ELISA <- c("STIM" = "#A6529A", "UNSTIM" = "#BEBDBD")
LABELS      <- c("+ IL-1", "- IL-1")
SIZE = 25
X_AXIS <- "Relative IL-2 secretion"

################################################################################
plotting_means$STIM_DAY <- as.factor(plotting_means$STIM_DAY)

# plot
figure_X_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = CONDITION)) +
  geom_col(position = position_dodge(width = 0.7), aes(color = CONDITION), width = 0.68, alpha = 0.3) +
  geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                           xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                                           xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
  geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT), shape = 21, 
             size = 4, 
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
  geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.05, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = 5, angle = 90) +
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.5), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(values = COLOR_ELISA, labels = LABELS) +
  scale_fill_manual(values  = COLOR_ELISA, labels = LABELS) +
  labs(x = X_AXIS, y = "") +
  guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 2)) +
  theme_cowplot(font_size = SIZE) +
  theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
        axis.title.y      = element_blank(),
        # legend.position   = c(x = 0.4, y = 0.73),
        legend.title      = element_blank(),
        legend.text       = element_text(size = SIZE),
        legend.key.size   = unit(9, "mm"))

figure_X_ELISA

if (save) {
  
  save_to_data_tay <- file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)
  save_to_local    <- file.path("~/Desktop/Source files/", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(file.path(save_to, "figure_X_ELISA.svg"), plot = last_plot(), device = "svg", width = 12, height = 7)
  
  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
  fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt.csv"))
}


################################################################################
# Plotting w/o cl069

plotting_data  <- plotting_data_main
plotting_data  <- plotting_data  %>% filter(CL_NAME_ON_PLOT != "cl069")
plotting_means <- plotting_means %>% filter(CL_NAME_ON_PLOT != "cl069")
plotting_stats <- plotting_stats %>% filter(CL_NAME_ON_PLOT != "cl069")

# plot
figure_X_ELISA_no_cl069 <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = CONDITION)) +
  geom_col(position = position_dodge(width = 0.7), aes(color = CONDITION), width = 0.68, alpha = 0.3) +
  geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                           xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                                           xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
  geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT), shape = 21, 
             size = 4, 
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
  geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.05, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = 5, angle = 90) +
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.01), position = "bottom") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(values = COLOR_ELISA, labels = LABELS) +
  scale_fill_manual(values  = COLOR_ELISA, labels = LABELS) +
  labs(x = X_AXIS, y = "") +
  guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 2)) +
  theme_cowplot(font_size = SIZE) +
  theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
        axis.title.y      = element_blank(),
        # legend.position   = c(x = 0.4, y = 0.73),
        legend.title      = element_blank(),
        legend.text       = element_text(size = SIZE),
        legend.key.size   = unit(9, "mm"))

figure_X_ELISA_no_cl069

if (save) {
  
  save_to_data_tay <- file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)
  save_to_local    <- file.path("~/Desktop/Source files/", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(file.path(save_to, "figure_X_ELISA_no_cl069.svg"), plot = last_plot(), device = "svg", width = 12, height = 7)
  
  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
  fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt.csv"))
}

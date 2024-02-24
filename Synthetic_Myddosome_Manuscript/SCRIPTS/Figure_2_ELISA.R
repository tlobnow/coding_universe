library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr)

Input_Directory <- "./Synthetic_Myddosome_Manuscript/ELISA/"

NAME_KEY <- fread("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/Figure 1/1B/ELISA_CL_KEY.csv", header = T) 

# load functions
source("~/Documents/Github/coding_universe/Synthetic_Myddosome_Manuscript/SCRIPTS/functions.R")
source("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/Figure 1/1B/functions.R")

# process ELISA
All_plates_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
All_plates_data     <- left_join(All_plates_data_raw, NAME_KEY)
All_plates_data     <- All_plates_data %>% filter(CELL_LINE != "NA")

################################################################################
################################################################################
################################################################################
# FIGURE 1

save    <- T
save_to <- "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/Figure 2/2F"

All_days_data <- All_plates_data %>%
  filter(Date == "2022-06-09" | 
           Date == "2022-06-23" | 
           Date == "2022-07-01" |
           Date == "")

All_days_data$CONDITION <- as.factor(All_days_data$CONDITION)

################################################################################
# subset and normalize CHARMS data

CHARMS                   <- All_days_data %>% filter(Date == "2022-06-23" |  Date == "2022-07-01")
CHARMS_baseline_adj      <- get_baseline(CHARMS, NEGATIVE_CTRL = "3xKO")
CHARMS_normalization_adj <- get_normalization_value(CHARMS_baseline_adj, POSITIVE_CTRL = "Wild Type")
CHARMS_nomalized         <- normalize_ELISA(CHARMS_normalization_adj)

# subset data for CHARMS cell lines
CHARMS_data <- CHARMS_nomalized %>% filter(CL_NAME_ON_PLOT %in% c("CHARMS-MyD88", "CHARMS-1xTIR", "Wild Type"), Date == "2022-06-23")

# subset data for CHARMS cell line with alanine-mutated TRAF6 Binding sites
CHARMS_3xA_data <- CHARMS_nomalized %>% filter(CL_NAME_ON_PLOT %in% c("CHARMS-3xA"), Date == "2022-07-01")

################################################################################
# subset data for triple KO cell line (MyD88-/- IRAK4-/-, IRAK1-/-)

tKO_EL4               <- All_days_data %>% filter(Date == "2022-06-09")
tKO_baseline_adj      <- get_baseline(tKO_EL4, NEGATIVE_CTRL = "3xKO")
tKO_normalization_adj <- get_normalization_value(tKO_baseline_adj, POSITIVE_CTRL = "WT_EL4")
tKO_nomalized         <- normalize_ELISA(tKO_normalization_adj)

tKO_EL4_data <- tKO_nomalized %>% filter( CL_NAME_ON_PLOT == "3xKO")

################################################################################
# join data frames of interest
plotting_data <- rbind(CHARMS_data, tKO_EL4_data, CHARMS_3xA_data)

# reorder cell lines for plotting
plotting_data$CL_NAME_ON_PLOT = reorder(plotting_data$CL_NAME_ON_PLOT, -plotting_data$ORDER_NO)

################################################################################

# calculate means concentrations for plotting
plotting_means <- plotting_data %>%
  group_by(CELL_LINE, CONDITION, STIM_DAY, CL_NAME_ON_PLOT) %>%
  summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
            Relative_Intensity_mean = mean(Concentration_NORMALIZED)) %>%
  as.data.table()

################################################################################

# calculate stats
plotting_stats <- plotting_data %>%
  group_by(CELL_LINE, CONDITION, CL_NAME_ON_PLOT) %>%
  summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
            IL2_concentration_Dilution_Factor_sd = sd(Concentration),
            Relative_Intensity_mean = mean(triplicate_mean_per_day),
            Relative_Intensity_sem = sem(triplicate_mean_per_day)) %>%
  as.data.table()

################################################################################

# set plotting colors
color_elisa <- c("UNSTIM" = "gray50", "STIM" = "#A6529A")

################################################################################

SIZE = 25

# plot
figure_1_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = CONDITION)) +
  geom_col(position = position_dodge(width = 0.7), color = "black", width = 0.7, alpha = 0.7) +
  geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                           xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                                           xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
  geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT),
             shape = 21, size = 4,
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600)) +
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.5), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  scale_color_manual(labels = c("+ IL-1", "- IL-1"), values = c("#A6529A", "gray50")) +
  scale_fill_manual(labels  = c("+ IL-1", "- IL-1"), values = c("#A6529A", "gray50")) +
  labs(x = "Relative IL-2 secretion", y = "") +
  guides(color = "none",fill = guide_legend(nrow = 2)) +
  theme_cowplot(
    # font_size = SIZE
  ) +
  theme(axis.text.x = element_text(#size = SIZE,
    colour = "black", angle = 0, vjust = 0.6),
    axis.title.y = element_blank(),
    legend.position = c(0.6, 0.1),
    legend.title = element_blank(),
    legend.text = element_text(#size = SIZE,
      color = "black"),
    legend.key.size = unit(9, "mm"))

figure_1_ELISA

if (save) {
  # save figure
  ggsave(file.path(save_to, "Figure_1_ELISA.svg"), plot = figure_1_ELISA, device = "svg", width = 6, height = 10)
  # save tables
  fwrite(All_days_data,  file.path(save_to, "All_days_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv")) 
}

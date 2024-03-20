library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, 
               ggforce, ggbreak, patchwork, lemon)

################################################################################

# SETTINGS
save   <- T
prep_and_process <- TRUE
figure <- "Figure X"

# plotting settings

FONT   <- "Helvetica"
SIZE   <- 25
POINTS <- 8
TEXT   <- 8
ZOOM_LOC_X = 0.05

################################################################################

# PREPROCESSING
Input_Directory <- file.path("/Users/u_lobnow/Documents/Github/coding_universe/Synthetic_Myddosome_Manuscript/ELISA/")
NAME_KEY        <- fread("https://raw.githubusercontent.com/tlobnow/coding_universe/main/ELISA_CL_KEY.csv", header = T) 

# load functions
source("~/Documents/Github/Analysis_Synthethic_Myddosome/FigureScripts/functions.R")

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

# process ELISA
plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()

# ensure correct column type assignment
plate_data$CONDITION <- as.factor(plate_data$CONDITION)

################################################################################

# subset data

if (prep_and_process) {
  pre_CHARMS        <- plate_data %>% filter(Date == "2022-06-23" |  Date == "2022-07-01") %>% mutate(ORIGIN_2 = "CHARMS")
  CHARMS_norm       <- process_ELISA_data(DF = pre_CHARMS, NEGATIVE_CTRL = "3xKO", POSITIVE_CTRL = "3E10") %>% filter(CL_NAME_ON_PLOT %in% c("CHARMS", "3E10"), Date == "2022-06-23")
  
  pre_BDLD_6        <- plate_data %>% filter(Date == "2023-07-04" & Plate == 3)
  CHARMS_BDD_6_norm <- process_data(pre_BDLD_6, "BDLD_6$", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_6")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_11        <- plate_data %>% filter(Date == "2024-03-04" & Plate != 3)
  CHARMS_BDD_11_norm <- process_data(pre_BDLD_11, "BDLD_11", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_11")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_57        <- plate_data %>% filter(Date != "2024-03-04")
  CHARMS_BDD_57_norm <- process_data(pre_BDLD_57, "BDLD_57", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_57")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_58        <- plate_data %>% filter(Date == "2024-03-04" & Plate != 3)
  CHARMS_BDD_58_norm <- process_data(pre_BDLD_58, "BDLD_58", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_58")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_38        <- plate_data %>% filter(Date == "2024-03-04" & Plate != 3)
  CHARMS_BDD_38_norm <- process_data(pre_BDLD_38, "BDLD_38", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_38")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_69        <- plate_data %>% filter(Date == "2023-07-04" & Plate != 3)
  CHARMS_BDD_69_norm <- process_data(pre_BDLD_69, "BDLD_69", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_69")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_62        <- plate_data
  CHARMS_BDD_62_norm <- process_data(pre_BDLD_62, "BDLD_62", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_62")  %>% filter(PURPOSE == "SAMPLE")
  
  pre_BDLD_67        <- plate_data
  CHARMS_BDD_67_norm <- process_data(pre_BDLD_67, "BDLD_67", "3xKO", "3E10") %>% mutate(ORIGIN_2 = "BDLD_67")  %>% filter(PURPOSE == "SAMPLE")
  
  # join data frames of interest
  plotting_data        <- rbind(CHARMS_norm, 
                                CHARMS_BDD_6_norm, 
                                CHARMS_BDD_57_norm, 
                                CHARMS_BDD_58_norm, 
                                CHARMS_BDD_11_norm, 
                                CHARMS_BDD_38_norm, 
                                CHARMS_BDD_69_norm, 
                                CHARMS_BDD_62_norm,
                                CHARMS_BDD_67_norm)
  
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main) %>% left_join(NAME_KEY)
  
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt) %>% left_join(NAME_KEY)
  
  # reordering
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats$CL_NAME_ON_PLOT <- reorder(plotting_stats$CL_NAME_ON_PLOT, -plotting_stats$ORDER_NO)
  
}


################################################################################

# plot
figure_BDLD_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
  geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
  geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                    xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                    xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
  geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
             shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
  geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.5, y = CL_NAME_ON_PLOT, label = significance), 
            hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1.2, by = 0.5), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
  scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
  labs(x = "Relative IL-2 secretion", y = "") +
  guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
  theme_cowplot(font_size = SIZE, font_family = FONT) +
  theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
        axis.title.y      = element_blank(),
        legend.position   = "bottom",
        legend.title      = element_blank(),
        legend.text       = element_text(size = SIZE),
        legend.key.size   = unit(9, "mm"))

figure_BDLD_ELISA



if (save) {
  
  save_to_data_tay <- file.path("/Volumes/TAYLOR-LAB/Finn_v2/DLD_Project/ELISA/figures/")
  save_to_local    <- file.path("~/Desktop/Source files/")
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(file.path(save_to, "figure_BDLD_ELISA.svg"), plot = figure_BDLD_ELISA, device = "svg", width = 12, height = 24)

  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
  fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt.csv"))
}

library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr)

Input_Directory <- "./Synthetic_Myddosome_Manuscript/ELISA/"
save_to         <- Input_Directory

# adjust the names as needed!
NAME_KEY_DOC <- fread("./Synthetic_Myddosome_Manuscript/SCRIPTS/ELISA_CL_KEY_DOC.csv", header = T, sep = ",", quote = "")
head(NAME_KEY_DOC, n = 10)

NAME_KEY     <- fread("./Synthetic_Myddosome_Manuscript/SCRIPTS/ELISA_CL_KEY.csv", header = T) 
head(NAME_KEY, n = 10)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/coding_universe/main/scripts/ELISA_Fx.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/coding_universe/main/scripts/ELISA_Fx.R",
                     no  =  "./scripts/ELISA_Fx.R"))


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

All_plates_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
All_plates_data     <- left_join(All_plates_data_raw, NAME_KEY)
All_plates_data     <- All_plates_data %>% filter(CELL_LINE != "NA")
DF                  <- All_plates_data

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Figure 1

DF <- All_plates_data %>% filter(Date == "2022-06-23")

FILTER_VALUES = "CHARMS-DHF91 (H)"
FILTER_VALUE  = FILTER_VALUES
FILTER_TYPE   = "COHORT"
POSITIVE_CTRL = c("cl069")
NEGATIVE_CTRL = c("tKO_EL4", "204_TRIPLE_KO")
COLOR         = "cornflowerblue"
x_label       = ""
y_label       = "relative IL-2 conc."
plot_title    = "IL-2 ELISA"
subtitle      = paste0("Assay performed for ", FILTER_VALUES, ", based on filter type ", FILTER_TYPE)
save          = T

# Automatic settings
# DF = All_plates_data
DF$CL_NAME_ON_PLOT = reorder(DF$CL_NAME_ON_PLOT, -DF$ORDER_NO)

SEED = 600
plot_pval = T
plot_faceted_by_date = F

# Step 1: Filter & Subset the main data into lists
group_vars  <- c("Date", "STIM_DAY")

DATA <- DF

baseline <- DATA %>%
  group_by(Date) %>%
  summarise(baseline_control_value = min(Concentration))

# Join the calculated values with the dataset
data <- left_join(DATA, baseline) %>%
  mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value, TRUE ~ Concentration))

control_mean_per_day <- data %>%
  filter((CELL_LINE %in% POSITIVE_CTRL & CONDITION == "STIM") | (CL_NAME_ON_PLOT %in% POSITIVE_CTRL & CONDITION == "STIM")) %>%
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


COMBINED_DATA <- DATA_NORMALIZED

# Step 5: Calculate the means per cohort
MEANS <- COMBINED_DATA %>%
  group_by(CELL_LINE, CONDITION, Date) %>%
  distinct(triplicate_mean_per_day, STIM_DAY, .keep_all = TRUE) %>%
  ungroup()

MEANS$STIM_DAY <- as.factor(MEANS$STIM_DAY)


# Step 6: Calculate the mean of means (MOM) per cohort
MOM_SUBSET <- MEANS %>%
  group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION, Date) %>%
  summarise(triplicate_sd_per_day = sd(triplicate_mean_per_day),
            triplicate_mean_per_day = mean(triplicate_mean_per_day)) %>%
  distinct(.keep_all = TRUE)

# Step 7: Perform statistics per cohort
STATISTICAL_RESULTS <- perform_statistical_analysis(MEANS, "CL_NAME_ON_PLOT")

# Step 8: Plotting
my_plot <- ggplot(MEANS, aes(x = CL_NAME_ON_PLOT)) +
  geom_col(data = MOM_SUBSET,
           aes(y = triplicate_mean_per_day, fill = CONDITION), position = position_dodge(width = 1), alpha = 0.5) +
  geom_point(data = MEANS,
             aes(y = triplicate_mean_per_day, group = CONDITION, shape = STIM_DAY),
             position = position_jitterdodge(jitter.height = 0, jitter.width = 1.2, seed = SEED),
             col = "white", size = 4) +
  geom_point(data = MEANS, aes(y = triplicate_mean_per_day, group = CONDITION, shape = STIM_DAY),
             position = position_jitterdodge(jitter.height = 0, jitter.width = 1.2, seed = SEED),
             col = "black", size = 3) +
  geom_errorbar(data = MOM_SUBSET, aes(ymin = triplicate_mean_per_day - triplicate_sd_per_day,
                                       ymax = triplicate_mean_per_day + triplicate_sd_per_day, group = CONDITION),
                width = 0.25, position = position_dodge(width = 1)) +
  labs(x = x_label, y = y_label) +
  scale_fill_manual(values = c("UNSTIM" = "gray50", "STIM" = COLOR)) +
  ggtitle(plot_title) +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = data.frame(CL_NAME_ON_PLOT = names(STATISTICAL_RESULTS[["annotations"]]), 
                              y = max(MEANS$triplicate_mean_per_day) + 0.1, 
                              label = STATISTICAL_RESULTS[["annotations"]]), 
            aes(y = y, label = label), col = "gray40", size = 4) +
  geom_text(data = data.frame(CL_NAME_ON_PLOT = names(STATISTICAL_RESULTS[["p_values"]]), 
                              y = max(MEANS$triplicate_mean_per_day) - 0.2, 
                              label = STATISTICAL_RESULTS[["p_values"]]), 
            aes(y = y, label = label), size = 7) +
  coord_flip() +
  facet_wrap(~Date)

my_plot

















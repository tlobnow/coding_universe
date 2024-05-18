# ELISA FOR BDLD_57 & BDLD_6
#          .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^

#Load libraries
library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot)

# ****************************************************************************
# ****************************************************************************
# ****************************************************************************

MAIN <- "~/Documents/Github/master_thesis/ELISA"
Input_Directory <- file.path(MAIN, "ALL/") #Make folder with desired Plates to be analysed.  Plate name Format Plate_1_20220623
Output_Directory <- file.path(Input_Directory, "Output")
save_to = "~/Documents/Github/master_thesis/thesis_figs"

POSITIVE_CTRL = "MyD88-GFP"
NEGATIVE_CTRL = "204_TRIPLE_KO"

ifelse(test = file.exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
       yes  = source("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
       no   = source("~/Documents/Github/master_thesis/scripts/functions.R"))

CELL_LINES    <- "CELL_LINES.csv" # COHORT
STIM_DAYS     <- "STIMULATION_DAYS.csv" # SAMPLE_DAY
CONDITIONS    <- "CONDITIONS.csv" # STIMULATION_CONDITION
MEASUREMENTS  <- "MEASUREMENTS.csv" # VALUES_MEASURED

#Specify Dilution factor of supernatant
DILUTION_FACTOR_2  <- 2
DILUTION_FACTOR_5  <- 5
DILUTION_FACTOR_10 <- 10

# Run ELISA_Fx()
All_plates_data <- ELISA_Fx(Input_Directory, Output_Directory)

# Prep Analysis
All_plates_data$CELL_LINES <- gsub("_H$|_E$", "", All_plates_data$CELL_LINES)
All_plates_data$CELL_LINES <- ifelse(grepl("^BDLD_", All_plates_data$CELL_LINES), 
                                     paste0("MyD88-", All_plates_data$CELL_LINES, "-T6BM"),
                                     All_plates_data$CELL_LINES)

All_plates_data <- All_plates_data %>% 
  mutate(CELL_LINES = case_when(CELL_LINES == "MyD88-BDLD_10-T6BM"  ~ "MyD88-CNTRL-T6BM",
                                CELL_LINES == "3E10_GFP"  ~ "MyD88-GFP",
                                CELL_LINES == "69"  ~ "WT_EL4",
                                CELL_LINES == "614" ~ "DHF91_H",
                                CELL_LINES == "614_L" ~ "DHF91_L",
                                CELL_LINES == "590_L" ~ "MyD88-T6BM",
                                CELL_LINES == "591_L" ~ "MyD88-T6BM-3xAla",
                                CELL_LINES == "631" ~ "MyD88-TIR-T6BM-3xAla",
                                CELL_LINES == "632" ~ "DHF91-T6BM-3xAla",
                                CELL_LINES != "614" ~ CELL_LINES),
         # MEASUREMENTS = case_when(MEASUREMENTS < 0 ~ 0, TRUE ~ MEASUREMENTS),
         # IL2_concentration = case_when(IL2_concentration < 0 ~ 0, TRUE ~ IL2_concentration),
         # IL2_concentration_DILUTION_FACTOR = case_when(IL2_concentration_DILUTION_FACTOR < 0 ~ 0, TRUE ~ IL2_concentration_DILUTION_FACTOR),
         CONTROL_CELL_LINE = case_when(CELL_LINES  %in% c("WT_EL4", "MyD88-GFP", "204_TRIPLE_KO") ~ T,
                                       !CELL_LINES %in% c("WT_EL4", "MyD88-GFP", "204_TRIPLE_KO") ~ F))


# EXCLUDE PLATES THAT HAD TOO HIGH EXPRESSION MEASUREMENTS (Threshold will be 1.1) OF POSITIVE CONTROLS (3E10/WT_EL4)
# DropPlates <- All_plates_data %>%
#   filter(MEASUREMENTS > 1.1 & CONTROL_CELL_LINE == TRUE) %>%
#   distinct(Date)
# All_plates_data <- All_plates_data[!All_plates_data$Date %in% DropPlates, ]

# ADJUST COLUMN TYPES
All_plates_data$CONDITIONS <- as.factor(All_plates_data$CONDITIONS)
All_plates_data$STIM_DAYS  <- as.factor(All_plates_data$STIM_DAYS)

# COHORT = "BDLD_57"
COHORT = "MyD88-BDLD_57-T6BM"

plot4COHORT <- function(COHORT, COLOR = "salmon", SEED = 600) {
  if (!is.na(COHORT)) {
    # PLOTTING COLORS
    fill_colors <- c("UNSTIM" = "gray50")
    # fill_colors <- c("STIM" = "gray40")
    fill_colors["STIM"]   <- COLOR
    
    # Data Filtering
    # Step 1: Data Filtering and Control Computation
    # Unique dates corresponding to the given cohort.
    COHORT_DATES <- All_plates_data %>% filter(CELL_LINES == COHORT) %>% distinct(Date) %>% pull(Date)
    
    # Unique plate names corresponding to the given cohort.
    COHORT_PLATES <- All_plates_data %>% filter(CELL_LINES == COHORT) %>% distinct(Plate) %>% pull(Plate)
    
    # Unique stim days corresponding to the given cohort.
    STIM_DAYS_COHORT <- All_plates_data[All_plates_data$CELL_LINES %in% COHORT, ] %>% pull(STIM_DAYS) %>% unique()
    
    # Subset the dataset.
    COHORT_SUBSET <- All_plates_data %>%
      filter(Date %in% COHORT_DATES,
             Plate %in% COHORT_PLATES,
             STIM_DAYS %in% STIM_DAYS_COHORT,
             CELL_LINES %in% c(POSITIVE_CTRL,
                               NEGATIVE_CTRL,
                               COHORT))
    
    # Normalization against the negative control
      baseline <- COHORT_SUBSET %>%
        group_by(Date, STIM_DAYS) %>%
        filter(CELL_LINES == NEGATIVE_CTRL) %>% 
        summarise(baseline_control_value = mean(IL2_concentration_DILUTION_FACTOR))
      COHORT_SUBSET <- left_join(COHORT_SUBSET, baseline)
    
    # Compute Control Mean
      control_mean_per_day <- COHORT_SUBSET %>% 
        filter(CELL_LINES == POSITIVE_CTRL,
               CONDITIONS == "STIM") %>%
        group_by(Date, STIM_DAYS) %>%
        summarise(control_mean_measurements = mean(IL2_concentration_DILUTION_FACTOR))
      COHORT_SUBSET <- left_join(COHORT_SUBSET, control_mean_per_day, by = c("Date", "STIM_DAYS"))
      
      
      COHORT_SUBSET <- COHORT_SUBSET %>%
        group_by(Date, Plate, STIM_DAYS, CELL_LINES, CONDITIONS) %>%
        mutate(CNTRL_MINUS_BASELINE = control_mean_measurements-baseline_control_value,
               IL2_concentration_DILUTION_FACTOR_REDUCED = IL2_concentration_DILUTION_FACTOR - baseline_control_value,
               IL2_concentration_DILUTION_FACTOR_NORMALIZED = IL2_concentration_DILUTION_FACTOR_REDUCED/CNTRL_MINUS_BASELINE)
      
      # Mean of means (MOM) 
      MOM_by_night <- COHORT_SUBSET %>%
        group_by(CELL_LINES, CONDITIONS) %>%
        mutate(triplicate_mean_per_day = mean(IL2_concentration_DILUTION_FACTOR_NORMALIZED),
               triplicate_sd_per_day   = sd(IL2_concentration_DILUTION_FACTOR_NORMALIZED)) %>%
        unique()
      
      # Mean of means (MOM) 
      MOM_by_Day <- COHORT_SUBSET %>%
        group_by(CELL_LINES, CONDITIONS, STIM_DAYS) %>%
        summarise(triplicate_mean_per_day = mean(IL2_concentration_DILUTION_FACTOR_NORMALIZED),
                  triplicate_sd_per_day   = sd(IL2_concentration_DILUTION_FACTOR_NORMALIZED))
      
      
      Plate_Summary_Day_Rel2 <- MOM_by_night %>%
        group_by(CELL_LINES, CONDITIONS) %>%
        summarise(IL2_concentration_DILUTION_FACTOR_REDUCED = mean(IL2_concentration_DILUTION_FACTOR_REDUCED),
                  IL2_concentration_DILUTION_FACTOR_REDUCED_sd = IL2_concentration_DILUTION_FACTOR_REDUCED) %>%
        distinct(IL2_concentration_DILUTION_FACTOR_REDUCED_sd, .keep_all = TRUE)
      
      # Pairwise t-test
      pairwise_ttest <- function(data, return_annotation = FALSE) {
        p_values <- pairwise.t.test(data$IL2_concentration_DILUTION_FACTOR_REDUCED, data$CONDITIONS, p.adjust.method = "none")$p.value
        
        if (return_annotation) {
          p_annotation <- ifelse(p_values < 0.001, '***', 
                                 ifelse(p_values < 0.01, '**', 
                                        ifelse(p_values < 0.05, '*', 
                                               'ns')))
          return(p_annotation)
        } else {
          return(formatC(p_values, format = "e", digits = 3))
        }
      }
      
      # calculate statistical significance using a t-test for each cell line (df of interest is split into groups based on cell lines)
      annotations <- sapply(split(MOM_by_night, MOM_by_night$CELL_LINES), pairwise_ttest)
      p_val       <- sapply(split(MOM_by_night, MOM_by_night$CELL_LINES), function(data) pairwise_ttest(data, return_annotation = TRUE))
      
      # --- VISUALIZATION ---
      plot <- ggplot(NORMALIZED_TO_CONTROL, aes(x = CELL_LINES)) +
        geom_col(aes(y = trip_mean, fill = CONDITIONS), position = position_dodge(width = 1), alpha = 0.5) +
        geom_point(aes(y = nrmlz_fct_mean, group = CONDITIONS, shape = STIM_DAYS), 
                   position = position_jitterdodge(jitter.height = 0, 
                                                   jitter.width = 1.2, seed = SEED),
                   col = "white", size = 4) +
        geom_point(aes(y = nrmlz_fct_mean, group = CONDITIONS, shape = STIM_DAYS), 
                   position = position_jitterdodge(jitter.height = 0, 
                                                   jitter.width = 1.2, seed = SEED),
                   col = "black", size = 3) +
        geom_errorbar(data = Plate_Summary_Day_Rel2, 
                      aes(ymin = nrmlz_fct_mean - nrmlz_fct_sd, ymax = nrmlz_fct_mean + nrmlz_fct_sd, group = CONDITIONS),
                      width = 0.25, position = position_dodge(width = 1)) +
        geom_text(data = data.frame(CELL_LINES = names(annotations), y = max(NORMALIZED_TO_CONTROL$nrmlz_fct_mean) + 0.1, label = annotations),
                  aes(y = y, label = label), col = "gray40",
                  size = 4) +
        geom_text(data = data.frame(CELL_LINES = names(p_val), y = max(NORMALIZED_TO_CONTROL$nrmlz_fct_mean) + 0.2, label = p_val),
                  aes(y = y, label = label), size = 7) +
        theme_cowplot() +
        theme(legend.position = "bottom") +
        labs(x = "",
             y = "relative IL-2 conc.") +
        scale_fill_manual(values = fill_colors) +
        ggtitle(paste(COHORT, "ELISA")) +
        theme(plot.title = element_text(hjust = 0.5), 
              # legend.position = 0
        )
      
      ggsave(filename = paste0(save_to, paste0("/",COHORT,"_ELISA.svg")), 
             plot = plot, device = "svg", width = 6, height = 4)
      
      return(plot)
      
      # Normalization for fold change from UNSTIM to STIM
      NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
        group_by(Date, Plate, STIM_DAYS, CELL_LINES, CONDITIONS) %>%
        summarise(mean_per_day = mean(MEASUREMENTS),
                  IL2_concentration_mean = mean(IL2_concentration))
      
      # First, calculate the mean for the UNSTIM condition separately.
      unstim_means <- COHORT_SUBSET %>%
        filter(CONDITIONS == "UNSTIM") %>%
        group_by(CELL_LINES) %>%
        summarise(mean_unstim = mean(MEASUREMENTS))
      
      # Now, join this back to the main dataset.
      COHORT_SUBSET <- COHORT_SUBSET %>%
        left_join(unstim_means, by = "CELL_LINES")
      
      # Compute the fold change. 
      # The `case_when` logic ensures that the fold change is calculated only for the STIM condition.
      NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
        mutate(fold_change = case_when(
          # set NA for non-STIM conditions
          CONDITIONS == "STIM" ~ MEASUREMENTS / mean_unstim, TRUE ~ NA_real_ )) %>%
        ungroup() %>%
        group_by(CELL_LINES, CONDITIONS) %>%
        mutate(trip_mean = mean(fold_change),
               fold_change_sd = sd(fold_change)) %>%
        unique()
      
      # Helper function to perform a t-test and return p-value and annotation
      perform_ttest <- function(data) {
        ttest_result <- t.test(data$MEASUREMENTS[data$CONDITIONS == "STIM"], 
                               data$MEASUREMENTS[data$CONDITIONS == "UNSTIM"],
                               paired = TRUE)  # or paired = FALSE if the data are independent
        
        p_value <- ttest_result$p.value
        annotation <- ifelse(p_value < 0.001, '***', 
                             ifelse(p_value < 0.01, '**', 
                                    ifelse(p_value < 0.05, '*', 'ns')))
        return(list(p_value = p_value, annotation = annotation))
      }
      
      # Compute fold changes and perform t-tests
      results <- COHORT_SUBSET %>%
        group_by(CELL_LINES) %>%
        do({
          data <- .
          fold_change_data <- data %>%
            summarise(fold_change = mean(MEASUREMENTS[CONDITIONS == "STIM"]) / 
                        mean(MEASUREMENTS[CONDITIONS == "UNSTIM"]))
          
          ttest_results <- perform_ttest(data)
          cbind(fold_change_data, ttest_results)
        }) %>%
        ungroup()
      
      # Extract annotations for plotting
      annotations <- results$annotation
      names(annotations) <- results$CELL_LINES
      
      # Plotting the fold changes
      plot <- ggplot(results, aes(x = CELL_LINES, y = fold_change)) +
        geom_col(aes(fill = CONDITIONS), position = position_dodge(width = 0.8)) +
        # Add other layers like geom_errorbar if you have standard deviations to include
        geom_text(aes(label = annotations), vjust = -0.5) +
        theme_cowplot() +
        theme(legend.position = "bottom") +
        labs(x = "",
             y = "fold change") +
        # scale_fill_manual(values = fill_colors) +
        ggtitle(paste(COHORT, "ELISA (fold change)")) +
        theme(plot.title = element_text(hjust = 0.5), 
              # legend.position = 0
        )
      
      ggsave(filename = paste0(save_to, paste0("/",COHORT,"_ELISA_foldchange.svg")), 
             plot = plot, device = "svg", width = 6, height = 4)
      
      return(plot)
  } else print("Please supply a COHORT")
}

### Usage ----
plot4COHORT(COHORT = paste0("MyD88-","BDLD_6","-T6BM"),
            # COLOR = "#2171b5",
            COLOR = "#2d8659",
            SEED = 600, 
            normalize_to_control = F, 
            normalize_as_foldchange = T)

# plot4COHORT(paste0("MyD88-","BDLD_57","-T6BM"), "#d72323")
plot4COHORT(paste0("MyD88-","BDLD_57","-T6BM"), "#cddc39", normalize_to_control = T, normalize_as_foldchange = F)


plot4COHORT(paste0("MyD88-","BDLD_2","-T6BM"), "gray10", normalize_to_control = T, normalize_as_foldchange = F)
plot4COHORT(paste0("MyD88-","BDLD_7","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","CNTRL","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_13","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_14","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_27","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_46","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_50","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_62","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_67","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_69","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
plot4COHORT(paste0("MyD88-","BDLD_73","-T6BM"), "gray10", normalize_to_control = F, normalize_as_foldchange = T)
# plot4COHORT(paste0("DHF91_L"), "lightgreen")
# plot4COHORT(paste0("DHF91_H"), "lightgreen")
plot4COHORT(paste0("MyD88-T6BM"), "#2171b5", normalize_to_control = F, normalize_as_foldchange = T)



# plot4COHORT("BDLD_2")
# plot4COHORT("BDLD_7")
# plot4COHORT("BDLD_10")
# plot4COHORT("BDLD_13")
# plot4COHORT("BDLD_14")
# plot4COHORT("BDLD_27")
# plot4COHORT("BDLD_46")
# plot4COHORT("BDLD_62")
# plot4COHORT("BDLD_67")
# plot4COHORT("BDLD_69")
# plot4COHORT("BDLD_73")

# DAYS = All_plates_data %>% distinct(Date) %>% pull(Date)
DAY="2023-07-04"

plot4Day <- function(DAY, COLOR = "salmon", SEED = 600) {
  if (!is.na(DAY)) {
    # PLOTTING COLORS
    fill_colors <- c("UNSTIM" = "gray40")
    fill_colors["STIM"]   <- COLOR
    
    # Data Filtering
    
    # Step 1: Data Filtering and Control Computation
    # Unique dates corresponding to the given cohort.
    COHORT_DATES <- All_plates_data %>% filter(Date == DAY) %>% distinct(Date) %>% pull(Date)
    # Unique cohorts corresponding to the given date
    COHORTS <- All_plates_data %>% filter(Date == DAY) %>% distinct(CELL_LINES) %>% pull(CELL_LINES)
    # Unique plate names corresponding to the given date
    COHORT_PLATES <- All_plates_data %>% filter(Date == DAY) %>% distinct(Plate) %>% pull(Plate)
    # Unique stim days corresponding to the given date
    STIM_DAYS_COHORT <- All_plates_data[All_plates_data$Date %in% as.IDate(DAY), ] %>% pull(STIM_DAYS) %>% unique()
    
    # Subset the dataset.
    COHORT_SUBSET <- All_plates_data %>%
      filter(Date %in% COHORT_DATES,
             Plate %in% COHORT_PLATES,
             STIM_DAYS %in% STIM_DAYS_COHORT,
             CELL_LINES %in% c("MyD88-GFP",
                               "204_TRIPLE_KO",
                               COHORTS
                               # "MyD88-T6BM",
                               # "WT_EL4",
             ))
    
    # Compute Control Mean
    control_mean_per_day <- COHORT_SUBSET %>% 
      filter(CELL_LINES == "MyD88-GFP",#(if ("WT_EL4" %in% COHORT_SUBSET$CELL_LINES) "WT_EL4" else "MyD88-GFP"), 
             CONDITIONS == "STIM") %>%
      group_by(Date, STIM_DAYS) %>%
      summarise(control_mean_measurements = mean(MEASUREMENTS))
    COHORT_SUBSET <- left_join(COHORT_SUBSET, control_mean_per_day, by = c("Date", "STIM_DAYS"))
    
    
    # Normalization
    NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
      group_by(Date, Plate, STIM_DAYS, CELL_LINES, CONDITIONS) %>%
      summarise(mean_per_day = mean(MEASUREMENTS),
                control_mean_measurements = mean(control_mean_measurements),
                IL2_concentration_mean = mean(IL2_concentration)) %>%
      mutate(nrmlz_fct_mean = mean_per_day/control_mean_measurements) %>%
      ungroup() %>%
      group_by(CELL_LINES, CONDITIONS, Date) %>%
      mutate(trip_mean = mean(nrmlz_fct_mean),
             nrmlz_fct_sd = sd(nrmlz_fct_mean)) %>%
      unique()
    
    Plate_Summary_Day_Rel2 <- NORMALIZED_TO_CONTROL %>%
      group_by(CELL_LINES, CONDITIONS, Date
      ) %>%
      summarise(nrmlz_fct_mean = mean(nrmlz_fct_mean),
                nrmlz_fct_med = median(nrmlz_fct_mean),
                nrmlz_fct_sd = nrmlz_fct_sd) %>%
      distinct(nrmlz_fct_sd, .keep_all = TRUE)
    
    NORMALIZED_TO_CONTROL$CONDITIONS <- factor(NORMALIZED_TO_CONTROL$CONDITIONS, levels = c("UNSTIM", "STIM"))
    Plate_Summary_Day_Rel2$CONDITIONS <- factor(Plate_Summary_Day_Rel2$CONDITIONS, levels = c("UNSTIM", "STIM"))
    
    # Pairwise t-test
    pairwise_ttest <- function(data, return_annotation = FALSE) {
      p_values <- pairwise.t.test(data$nrmlz_fct_mean, data$CONDITIONS, p.adjust.method = "none")$p.value
      
      if (return_annotation) {
        p_annotation <- ifelse(p_values < 0.001, '***', 
                               ifelse(p_values < 0.01, '**', 
                                      ifelse(p_values < 0.05, '*', 
                                             'ns')))
        return(p_annotation)
      } else {
        return(formatC(p_values, format = "e", digits = 3))
      }
    }
    
    
    # Function to plot given data
    create_plot <- function(data_to_plot, summary_data) {
      annotations <- sapply(split(data_to_plot, data_to_plot$CELL_LINES), pairwise_ttest)
      p_val       <- sapply(split(data_to_plot, data_to_plot$CELL_LINES), function(data) pairwise_ttest(data, return_annotation = TRUE))
      
      
      # --- VISUALIZATION ---
      # add a specific color for cell lines that show significant rescue ability:
      
      # Add a marker for significant cell lines
      significant_names <- names(p_val)[!(p_val %in% c(NA, "ns"))]
      if (length(significant_names) > 0) {
        significant_lines <- data.frame(CELL_LINES = significant_names, is_significant = TRUE)
      } else {
        significant_lines <- data.frame(CELL_LINES = character(0), is_significant = logical(0))
      }
      
      # Join with the main dataset
      data_to_plot <- left_join(data_to_plot, significant_lines, by = "CELL_LINES")
      data_to_plot$is_significant[is.na(data_to_plot$is_significant)] <- FALSE
      
      # Update the color assignment
      data_to_plot <- data_to_plot %>%
        mutate(fill_color = ifelse(CONDITIONS == "UNSTIM", fill_colors["UNSTIM"], 
                                   ifelse(CONDITIONS == "STIM" & is_significant == TRUE, "turquoise", fill_colors["STIM"])))
      
      # Check if there are plates with numbers greater than 2
      plates_gt_2 <- unique(COHORT_SUBSET$Plate[COHORT_SUBSET$Plate > 2])
      
      plot <- ggplot(data_to_plot, aes(x = CELL_LINES)) +
        geom_col(aes(y = trip_mean, fill = fill_color), position = position_dodge(width = 1), alpha = 0.5) +
        geom_point(aes(y = nrmlz_fct_mean, group = CONDITIONS, shape = STIM_DAYS), 
                   position = position_jitterdodge(jitter.height = 0, 
                                                   jitter.width = 1.2, seed = SEED),
                   col = "white", size = 4) +
        geom_point(aes(y = nrmlz_fct_mean, group = CONDITIONS, shape = STIM_DAYS), 
                   position = position_jitterdodge(jitter.height = 0, 
                                                   jitter.width = 1.2, seed = SEED),
                   col = "black", size = 3) +
        geom_errorbar(data = summary_data, 
                      aes(ymin = nrmlz_fct_mean - nrmlz_fct_sd, ymax = nrmlz_fct_mean + nrmlz_fct_sd, group = CONDITIONS),
                      width = 0.25, position = position_dodge(width = 1)) +
        geom_text(data = data.frame(CELL_LINES = names(annotations), 
                                    y = ifelse(is.finite(max(data_to_plot$nrmlz_fct_mean)), max(data_to_plot$nrmlz_fct_mean) + 0.2, 1), 
                                    label = p_val),
                  aes(y = y, label = label)) +
        geom_text(data = data.frame(CELL_LINES = names(p_val), 
                                    y = ifelse(is.finite(max(data_to_plot$nrmlz_fct_mean)), max(data_to_plot$nrmlz_fct_mean) + 0.1, 1), 
                                    label = annotations),
                  aes(y = y, label = label), col = "gray40") +
        theme_cowplot() +
        facet_rep_wrap(~Date) +
        theme(legend.position = "bottom") +
        labs(x = "",
             y = "relative IL-2 conc.") +
        scale_fill_manual(values = c(fill_colors, "turquoise" = "turquoise"))
      return(plot)
    }
    
    summary_data <- Plate_Summary_Day_Rel2 %>% filter(CELL_LINES %in% unique(NORMALIZED_TO_CONTROL$CELL_LINES))
    plot <- create_plot(NORMALIZED_TO_CONTROL, summary_data)
    ggsave(filename = paste0(save_to, paste0("/", as.IDate(DAY), "_ELISA.svg")), 
           plot = plot, device = "svg", width = 12, height = 8)
    return(plot)
  } else print("Please supply a DAY")
}

DAYS = All_plates_data %>% distinct(Date) %>% pull(Date)
for (DAY in DAYS) {
  plot4Day(DAY, COLOR = "salmon", SEED = 600)
}

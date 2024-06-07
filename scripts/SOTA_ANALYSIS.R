### SOTA ANALYSIS
SETTINGS = TRUE

if (SETTINGS) {
  
  LOAD_LIBS = TRUE
  LOAD_DIRS = TRUE
  LOAD_COLS = TRUE
  
  COMPARE_MSA_ALIGNERS = TRUE
  
  GENERATE_CORR_PLOTS  = FALSE
  GENERATE_AF3_PLOTS   = TRUE
      EXTRACT_AF3_OPERON_JSON_FILES = TRUE
  
  GENERATE_LCI_PLOTS   = TRUE
      GENERATE_TRAF6_VIOLIN_PLOTS = TRUE
      GENERATE_TRAF6_LT_PLOTS     = TRUE
      
  GENERATE_OPERON_TABLE = TRUE
  
  OVERWRITE = TRUE
  
  if (LOAD_LIBS) {library(pacman)
    pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, 
                   readxl, scales, knitr, tidyr, visdat, purrr, fs, jsonlite, ggrepel,
                   scico)
    }
  
  if (LOAD_DIRS) {
    setwd("/Users/u_lobnow/Documents/Github/coding_universe/")
    main       <- "SOTA"
    input_dir  <- file.path(main, "01_raw_data")
    temp_dir   <- file.path(main, "02_temp_data")
    output_dir <- file.path(main, "03_output_data")
    
    source("scripts/SOTA_FX.R")
    source("scripts/SOTA_COLNAMES.R") ; ls(pattern = "_cols")
  }
  
  sota <- fread(file.path(output_dir, "SOTA.csv"), na.strings = c("", "NA", NA))
}

af3 <- sota %>%
  select(any_of(c(
    basics_cols,
    # taxonomy_cols,
    # accession_cols,
    basics_extended_cols,
    # environment_cols,
    # sequence_cols,
    # annotation_cols,
    alphaFold_cols,
    # gBlock_cols,
    # cl_gen_cols,
    elisa_cols#,
    # lci_cols,
    # miscellaneous,
    # ncbi_genome_cols
    )))

# elisa <- sota %>%
#   select(any_of(c(
#     basics_cols,
#     taxonomy_cols,
#     # accession_cols,
#     # basics_extended_cols,
#     # environment_cols,
#     # sequence_cols,
#     # annotation_cols,
#     # alphaFold_cols,
#     # gBlock_cols,
#     # cl_gen_cols,
#     elisa_cols#,
#     # lci_cols,
#     # miscellaneous,
#     # ncbi_genome_cols
#   )))
# 
if (GENERATE_CORR_PLOTS) {
  # start putting together final outcome df that provides a good evaluation per assay

  # vis_cor(dplyr::select_if(sota, is.numeric))
  vis_cor(dplyr::select_if(af3, is.numeric), )
  # vis_cor(dplyr::select_if(elisa, is.numeric))
}

if (GENERATE_AF3_PLOTS) {
  af3_plt_subset <- af3 %>%
    filter(af3_n_chains == 10) %>%
    distinct(ID, .keep_all = T)
  
  # grab all IDs that follow the pattern BDLD and subset af3_plt_subset
  bdld_ids <- af3_plt_subset %>%
    filter(str_detect(ID, "BDLD_"),
           str_detect(ID, "(BDLD_87|BDLD_88)") == F) %>%
    pull(ID)
  
  # subset af3_plt_subset
  af3_plt_subset <- af3_plt_subset %>%
    filter(ID %in% bdld_ids)
  
  af3_plt_label_subset <- af3_plt_subset %>% filter(af3_ranking_score < 0.5)
  af3_plt_EDLD_25      <- af3 %>% filter(ID == "EDLD_25")
  cl_generated         <- af3_plt_subset %>% filter(!is.na(cl_name_high) | !is.na(cl_name_endo))
  
  
  # plot(plt_2 <- af3_plt_subset %>% select(ID, af3_iptm, af3_ptm, af3_ranking_score, af_confidence, species, fold_change))

  plot(af3_plt <- af3_plt_subset %>%
         ggplot(aes(af3_ranking_score, af3_iptm, label = ID)) +
         geom_point(show.legend = F, size = 3, col = "#009E73", alpha = 0.6) +
         geom_point(data = af3_plt_EDLD_25, fill = "#E4AC49", size = 4, shape = 21) + 
         geom_label_repel(data = af3_plt_label_subset, force = 0.5, vjust = 0.8, col = "gray30", max.overlaps = 50) +
         geom_label_repel(data = af3_plt_EDLD_25, col = "#E4AC49", vjust = 0.8) +
         theme_bw() +
         theme(axis.text = element_text(size = 15),
               axis.title = element_text(size = 15)) +
         lims(x = c(0, 1), y = c(0, 1)))
  
    
  # af3_plt_subset %>% ggplot(aes(af3_ranking_score, fill = taxend)) + geom_histogram(alpha = 0.5, show.legend = F) + facet_wrap(~taxend, ncol = 1)
  
  plot(af3_plt_no_labels <- af3_plt_subset %>%
         ggplot(aes(af3_ranking_score, af3_iptm, label = ID)) +
         geom_point(show.legend = F, size = 3, col = "#009E73", alpha = 0.6) +
         geom_point(data = af3_plt_EDLD_25, fill = "#E4AC49", size = 4, shape = 21) + 
         scale_color_scico_d(palette = "roma", direction = -1) +
         theme_bw() +
         theme(axis.text = element_text(size = 15),
               axis.title = element_text(size = 15)) +
         lims(x = c(0, 1), y = c(0, 1)))
  
  plot(af3_plt_no_labels_cl_gen <- af3_plt_subset %>%
         ggplot(aes(af3_ranking_score, af3_iptm, label = ID)) +
         geom_point(data = cl_generated, fill = "blue", size = 5, shape = 21) +
         geom_point(show.legend = F, size = 3, col = "#009E73", alpha = 0.6) +
         geom_point(data = af3_plt_EDLD_25, fill = "#E4AC49", size = 4, shape = 21) + 
         scale_color_scico_d(palette = "roma", direction = -1) +
         theme_bw() +
         theme(axis.text = element_text(size = 15),
               axis.title = element_text(size = 15)) +
         lims(x = c(0, 1), y = c(0, 1)))
  
  
  if (OVERWRITE) {
    ggsave(file.path(temp_dir, "figures/af3_ranking_score_vs_iptm.svg"), plot = af3_plt, device = "svg", width = 6, height = 3.5, dpi = 800)
    ggsave(file.path(temp_dir, "figures/af3_ranking_score_vs_iptm_no_labels.svg"), plot = af3_plt_no_labels, device = "svg", width = 6, height = 3.5, dpi = 800)
    ggsave(file.path(temp_dir, "figures/af3_ranking_score_vs_iptm_no_labels_cl_gen.svg"), plot = af3_plt_no_labels_cl_gen, device = "svg", width = 6, height = 3.5, dpi = 800)
  }
  
  
  if (EXTRACT_AF3_OPERON_JSON_FILES) {
    
    large_af3_dfs <- list()
    
    af_dir <- "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/01_raw_data/AF3_PREDICTIONS/OPERONS"
    
    af_json_files <- list.files(af_dir, full.names = TRUE, recursive = TRUE, pattern = "job_request") %>% str_subset("(?i)dld") %>% str_subset("\\.json$")
    
    for (af_json_file in af_json_files) {
      
      # af_json_file = af_json_files[1] # for debugging
      json <- fromJSON(af_json_file)
      
      af3_name     <- json[[1]][[1]]
      af3_seed     <- as.numeric(json[[2]][[1]])
      af3_seq      <- json[[3]][[1]]$proteinChain$sequence
      af3_n_chains <- json[[3]][[1]]$proteinChain$count
      
      summary_confidences_files <- list.files(af_dir, pattern = "summary_confidences", full.names = TRUE, recursive = TRUE) %>%
        str_subset(paste0("(?i)", gsub("-", "_", af3_name))) %>%
        str_subset("\\.json$")
      
      # Initialize a list to store the extracted information from each summary_confidences file
      sum_conf_extract_list <- list()
      
      for (summary_confidences_file in summary_confidences_files) {
        
        # summary_confidences_file <- summary_confidences_files[1] # for debugging
        af3_model               <- basename(summary_confidences_file)
        
        json <- fromJSON(summary_confidences_file)
        
        af3_chain_iptm          <- json[["chain_iptm"]]
        af3_chain_pair_iptm     <- json[["chain_pair_iptm"]]
        af3_chain_pair_pae_min  <- json[["chain_pair_pae_min"]]
        af3_chain_ptm           <- json[["chain_ptm"]]
        af3_fraction_disordered <- json[["fraction_disordered"]]
        af3_has_clash           <- json[["has_clash"]]
        ##########################################################################
        af3_iptm                <- json[["iptm"]]
        af3_n_recycles          <- json[["num_recycles"]]
        af3_ptm                 <- json[["ptm"]]
        af3_ranking_score       <- json[["ranking_score"]]
        
        sum_conf_extract <- cbind(af3_name,
                                  af3_seed,
                                  af3_seq,
                                  af3_n_chains,
                                  af3_model,
                                  af3_chain_iptm,
                                  af3_chain_pair_pae_min,
                                  af3_chain_ptm,
                                  af3_fraction_disordered,
                                  af3_has_clash,
                                  af3_iptm,
                                  af3_n_recycles,
                                  af3_ptm,
                                  af3_ranking_score
        ) %>%  as.data.frame()
        
        # Add the data frame to the list
        sum_conf_extract_list <- append(sum_conf_extract_list, list(sum_conf_extract))
        
        # Combine the af_json_extract data frame with all the summary_confidences data frames
        large_af3_df <- bind_rows(sum_conf_extract_list)
        
        # Add the large data frame to the list
        large_af3_dfs <- append(large_af3_dfs, list(large_af3_df))
      }
      
      # Combine all large data frames into one
      mega_af3_df <- bind_rows(large_af3_dfs) %>%
        mutate(ID = str_extract(af3_name, "(?i)^(?:[^_]*_)?([A-Z]+[12]?DLD\\d?)_\\d{2}"),
               ID = str_replace(ID, "(DLD|BDLD[12])_(\\d{2})", "\\1_\\2"))
      
      mega_af3_df$model = as.numeric(unlist(lapply(strsplit(gsub(".json", "", x = mega_af3_df$af3_model), "_", fixed=TRUE), function(x) return(x[9]))))
      
      # Subset the data to save a mini version with unique rows based on specific columns
      mini_af3_df <- mega_af3_df %>%
        select(
          ID,
          af3_name,
          af3_seed,
          af3_seq,
          af3_n_chains,
          af3_model,
          model,
          af3_iptm,
          af3_n_recycles,
          af3_ptm,
          af3_ranking_score
        ) %>% distinct(af3_model, .keep_all = T)
      
      
      mini_af3_df$af3_ranking_score <- as.numeric(mini_af3_df$af3_ranking_score)
      mini_af3_df$af3_iptm          <- as.numeric(mini_af3_df$af3_iptm)
      
      # filter highest ranked model (0)
      # AF3_SUMMARY <- mini_af3_df %>% group_by(af3_name) %>% filter(grepl("summary_confidences_0", af3_model))
      
      plot(operon_plot <- mini_af3_df %>%
             ggplot(aes(af3_ranking_score, af3_iptm, label = model)) +
             geom_point(show.legend = F, size = 3, col = "#009E73", alpha = 0.6) +
             geom_label_repel() +
             theme(axis.text = element_text(size = 15),
                   axis.title = element_text(size = 15))) +
        lims(x = c(0, 1), y = c(0, 1))
    }
  }
  
}

if (GENERATE_LCI_PLOTS) {
  
  NAME_KEY        <- fread("ELISA_CL_KEY.csv", header = T, na.strings = c("", "NA", NA)) %>% mutate(ID = str_replace(ID, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2"), CL_NAME_ON_PLOT = str_replace(CL_NAME_ON_PLOT, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2")) %>% arrange(ID)
  
  cell_table      <- fread(file.path(temp_dir, "lci_cell_table.csv.gz"))
  mean_cell       <- fread(file.path(temp_dir, "lci_mean_cell.csv.gz"))       %>% left_join(NAME_KEY[, c("COHORT", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOT_ID", "PLOTTING_COLOR")], by = "COHORT", relationship = "many-to-many") %>% unique()
  mean_replicates <- fread(file.path(temp_dir, "lci_mean_replicates.csv.gz")) %>% left_join(NAME_KEY[, c("COHORT", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOT_ID", "PLOTTING_COLOR")], by = "COHORT", relationship = "many-to-many") %>% unique()
  mean_total      <- fread(file.path(temp_dir, "lci_mean_total.csv.gz"))      %>% left_join(NAME_KEY[, c("COHORT", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOT_ID", "PLOTTING_COLOR")], by = "COHORT", relationship = "many-to-many") %>% unique()
  
  mean_cell$CL_NAME_ON_PLOT       <- factor(reorder(mean_cell$CL_NAME_ON_PLOT ,       -mean_cell$ORDER_NO))
  mean_replicates$CL_NAME_ON_PLOT <- factor(reorder(mean_replicates$CL_NAME_ON_PLOT , -mean_replicates$ORDER_NO))
  mean_total$CL_NAME_ON_PLOT      <- factor(reorder(mean_total$CL_NAME_ON_PLOT ,      -mean_total$ORDER_NO))

  
  # EXTRA FILTERS
  mean_cell       <- mean_cell %>% filter(CL_NAME_ON_PLOT != "CHARMS-DHF91 (L)" & CL_NAME_ON_PLOT != "CHARMS-TIR (L)")
  mean_replicates <- mean_replicates %>% filter(CL_NAME_ON_PLOT != "CHARMS-DHF91 (L)" & CL_NAME_ON_PLOT != "CHARMS-TIR (L)")
  mean_total      <- mean_total %>% filter(CL_NAME_ON_PLOT != "CHARMS-DHF91 (L)" & CL_NAME_ON_PLOT != "CHARMS-TIR (L)")
    
  if (GENERATE_TRAF6_VIOLIN_PLOTS) {
    
    LOW_CAT     <- "< 4 s"
    MEDIUM_CAT  <- "4-40 s"
    HIGH_CAT    <- "≥ 40 s"
    
    # Violin plot for long-lived TRAF6 HIGH_CAT ----
    plot(traf6_pct_over_40s <-
           ggplot(mean_cell, aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT, fill = PLOTTING_COLOR)) +
           geom_violin(data   = mean_cell       %>% filter(CATEGORY_DWELL_TIME == HIGH_CAT), scale = "width", width = 0.5, alpha = 0.3) +
           geom_jitter(data   = mean_replicates %>% filter(CATEGORY_DWELL_TIME == HIGH_CAT), size = 5, shape = 21, position = position_jitter(width = 0.2, height = 0, seed = 1000)) +
           geom_crossbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == HIGH_CAT), aes(ymin = PCT_RECRUITMENT, ymax = PCT_RECRUITMENT), width = 0.2, fatten = 1) +
           geom_errorbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == HIGH_CAT), aes(ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT, ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT), width = 0.15) +
           coord_flip() +
           scale_color_identity() +
           scale_fill_identity() +
           labs(x = "", y = "Percentage of long-lived TRAF6 on CHARMS", subtitle = "Category ≥40 s") +
           theme_cowplot() +
           theme(legend.position ="none"))
    
    # Violin plot for medium-lived TRAF6 MEDIUM_CAT ----
    plot(traf6_pct_4_to_40s <-
           ggplot(mean_cell, aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT, fill = PLOTTING_COLOR)) +
           geom_violin(data   = mean_cell       %>% filter(CATEGORY_DWELL_TIME == MEDIUM_CAT), scale = "width", width = 0.5, alpha = 0.3) +
           geom_jitter(data   = mean_replicates %>% filter(CATEGORY_DWELL_TIME == MEDIUM_CAT), size = 5, shape = 21, position = position_jitter(width = 0.2, height = 0, seed = 1000)) +
           geom_crossbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == MEDIUM_CAT), aes(ymin = PCT_RECRUITMENT, ymax = PCT_RECRUITMENT), width = 0.2, fatten = 1) +
           geom_errorbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == MEDIUM_CAT), aes(ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT, ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT), width = 0.15) +
           coord_flip() +
           scale_color_identity() +
           scale_fill_identity() +
           labs(x = "", y = "Percentage of medium-lived TRAF6 on CHARMS", subtitle = "Category 4-40 s") +
           theme_cowplot() +
           theme(legend.position ="none"))
    
    # Violin plot for short-lived TRAF6 LOW_CAT ----
    plot(traf6_pct_under_4s <-
           ggplot(mean_cell, aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT, fill = PLOTTING_COLOR)) +
           geom_violin(data   = mean_cell       %>% filter(CATEGORY_DWELL_TIME == LOW_CAT), scale = "width", width = 0.5, alpha = 0.3) +
           geom_jitter(data   = mean_replicates %>% filter(CATEGORY_DWELL_TIME == LOW_CAT), size = 5, shape = 21, position = position_jitter(width = 0.2, height = 0, seed = 1000)) +
           geom_crossbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == LOW_CAT), aes(ymin = PCT_RECRUITMENT, ymax = PCT_RECRUITMENT), width = 0.2, fatten = 1) +
           geom_errorbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == LOW_CAT), aes(ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT, ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT), width = 0.15) +
           coord_flip() +
           scale_color_identity() +
           scale_fill_identity() +
           labs(x = "", y = "Percentage of short-lived TRAF6 on CHARMS", subtitle = "Category 0-4 s") +
           theme_cowplot() +
           theme(legend.position ="none"))
    
    
    if (OVERWRITE) {
      ggsave(file.path(temp_dir, "figures/lci_traf6_pct_over_40s.svg"), plot = traf6_pct_over_40s, device = "svg", width = 12, height = 7, dpi = 800)
      ggsave(file.path(temp_dir, "figures/lci_traf6_pct_4_to_40s.svg"), plot = traf6_pct_4_to_40s, device = "svg", width = 12, height = 7, dpi = 800)
      ggsave(file.path(temp_dir, "figures/lci_traf6_pct_under_4s.svg"), plot = traf6_pct_under_4s, device = "svg", width = 12, height = 7, dpi = 800)
    }
  }
  
  if (GENERATE_TRAF6_LT_PLOTS) {
    
    cell_table      <- fread(file.path(temp_dir, "lci_cell_table.csv.gz")) %>% filter(COHORT != "MyD88-DHF91-TRAF6-BD TRAF6" & COHORT != "MyD88-TIR-TRAF6-BD TRAF6")
    
    NAME_KEY <- fread("ELISA_CL_KEY.csv", header = T) %>% select(CL_NAME_ON_PLOT, COHORT, PLOTTING_COLOR)
    ORDER_NO <- fread("ELISA_CL_KEY.csv", header = T) %>% select(CL_NAME_ON_PLOT, COHORT, ORDER_NO)
    
    cell_table$CATEGORY_DWELL_TIME <- factor(cell_table$CATEGORY_DWELL_TIME, levels = c(LOW_CAT, MEDIUM_CAT, HIGH_CAT))
    
    Mean_LT <- cell_table %>%
      filter(!is.na(COHORT)) %>%
      group_by(COHORT) %>% 
      summarise(LT_TRAF6 = mean(DWELL_TIME), 
                SEM_LT_TRAF6 = sem(DWELL_TIME)) %>%
      left_join(NAME_KEY)
    
    Mean_Total <- cell_table %>%
      group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
      count(CATEGORY_DWELL_TIME, COHORT, name = "N_CATEGORY_DWELL_TIME", .drop = FALSE) %>% 
      group_by(COHORT) %>% 
      mutate(PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME)) %>% 
      group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
      summarise(PCT_RECRUITMENT = mean(PCT_RECRUITMENT)) %>% 
      as.data.table() %>%
      left_join(NAME_KEY, relationship = "many-to-many") %>%
      unique()
    
    
    # Reorder levels of COHORT based on decreasing PCT_RECRUITMENT of HIGH_CAT
    ordered_cohorts <- Mean_Total %>%
      filter(CATEGORY_DWELL_TIME == HIGH_CAT) %>%
      arrange(desc(PCT_RECRUITMENT)) %>%
      select(CL_NAME_ON_PLOT) %>%
      pull() %>%
      factor(levels = unique(.)) %>%
      unique()
    
    # Adjust the factor levels in your Mean_Total dataframe:
    Mean_Total$CL_NAME_ON_PLOT   <- factor(Mean_Total$CL_NAME_ON_PLOT, levels = ordered_cohorts)
    
    Mean_Total_pivot <- Mean_Total %>% 
      pivot_wider(id_cols = c(CL_NAME_ON_PLOT), 
                  names_from = CATEGORY_DWELL_TIME, 
                  values_from = PCT_RECRUITMENT) %>%
      arrange(desc(HIGH_CAT))
    
    plotting_data <- Mean_Total
    
    plotting_data <- left_join(plotting_data, ORDER_NO, relationship = "many-to-many") %>% unique()
    plotting_data$CL_NAME_ON_PLOT <- reorder(plotting_data$CL_NAME_ON_PLOT, -plotting_data$ORDER_NO)
    plotting_data$BINNING_COLOR <- case_match(plotting_data$CATEGORY_DWELL_TIME, LOW_CAT    ~ "gray90", MEDIUM_CAT ~ "#FDEAA6", HIGH_CAT   ~ "#0072B2")
    
    # Create the ggplot
    plot(figure_TRAF6_LT <- ggplot(data = plotting_data, aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT * 100, 
                                                            # fill = CATEGORY_DWELL_TIME, 
                                                            fill = BINNING_COLOR,
                                                            group = CATEGORY_DWELL_TIME)) +
      geom_col(width = 0.7, size = 0.75, 
               color = "black",
               linewidth = 0.75, alpha = 1) +
      scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_data$BINNING_COLOR,  breaks = plotting_data$BINNING_COLOR, labels = case_match(plotting_data$BINNING_COLOR, "#fff3e5" ~ paste0(LOW_CAT), "#FDEAA6" ~ paste0(MEDIUM_CAT), "#b41f24" ~ paste0(HIGH_CAT))) +
      # scale_fill_scico_d(palette = 'nuuk', direction = -1) +
      labs(y = "% of total recruitments", x = "") +
      theme_classic(base_size = 20) +
      theme(legend.key.width = unit(65, "mm"),
            legend.key.height = unit(20, "mm"),
            legend.position = "top",
            legend.direction = "horizontal",
            text = element_text(family = "Helvetica"),
            plot.margin = margin(l = 0, t = 140, r = 0, b = 0, unit = "pt")) +
      guides(color = "none", fill = guide_legend(title = "TRAF6 Lifetime", reverse = TRUE)) +
      coord_flip())
    

    if (OVERWRITE) {
      ggsave(file.path(temp_dir, "figures/lci_TRAF6_LT.svg"), plot = figure_TRAF6_LT, device = "svg", width = 7, height = 14, dpi = 800)
    }
  }
}


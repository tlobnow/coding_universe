### SOTA ANALYSIS
SETTINGS = TRUE

if (SETTINGS) {
  
  LOAD_LIBS = TRUE
  LOAD_DIRS = TRUE
  LOAD_COLS = TRUE
  
  COMPARE_MSA_ALIGNERS = TRUE
  
  GENERATE_CORR_PLOTS  = TRUE
  GENERATE_AF3_PLOTS   = TRUE
  GENERATE_LCI_PLOTS   = TRUE
  
  if (LOAD_LIBS) {library(pacman)
    pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, 
                   readxl, scales, knitr, tidyr, visdat, purrr, fs, jsonlite, ggrepel)}
  
  if (LOAD_DIRS) {
    setwd("/Users/u_lobnow/Documents/Github/coding_universe/")
    main       <- "SOTA"
    input_dir  <- file.path(main, "01_raw_data")
    temp_dir   <- file.path(main, "02_temp_data")
    output_dir <- file.path(main, "03_output_data")
    
    source("scripts/SOTA_FX.R")
    source("scripts/SOTA_COLNAMES.R") ; ls(pattern = "_cols")
  }
}

# sota <- fread(file.path(output_dir, "SOTA.csv"), na.strings = c("", "NA", NA)) 
# ncbi <- fread(file.path(output_dir, "genome_info.csv"), na.strings = c("", "NA", NA)) 
# 
# af3 <- sota %>% 
#   select(any_of(c(
#     basics_cols,
#     taxonomy_cols,
#     # accession_cols,
#     # basics_extended_cols,
#     environment_cols,
#     sequence_cols,
#     annotation_cols,
#     alphaFold_cols#,
#     # gBlock_cols,
#     # cl_gen_cols,
#     # elisa_cols,
#     # lci_cols,
#     # miscellaneous,
#     # ncbi_genome_cols
#     )))
# 
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
# if (GENERATE_CORR_PLOTS) {
#   # start putting together final outcome df that provides a good evaluation per assay
#   
#   # vis_cor(dplyr::select_if(sota, is.numeric))
#   # vis_cor(dplyr::select_if(af3, is.numeric))
#   # vis_cor(dplyr::select_if(elisa, is.numeric))
# }
# 
# if (GENERATE_AF3_PLOTS) {
#   af3_plt_subset <- af3 %>% 
#     filter(af3_n_chains == 10) %>% 
#     distinct(ID, .keep_all = T)
#   
#   plt_2 <- af3 %>% 
#     select(ID, af3_iptm, af3_ptm, af3_ranking_score, af_confidence)
#     # select(any_of(alphaFold_cols))
#   plot(plt_2)
# 
#   af3 %>%
#     ggplot(aes(af3_ranking_score, af3_ptm, 
#                col = af_confidence,
#                label = ID)) +
#     geom_point(show.legend = F) +
#     geom_label_repel(max.overlaps = 50)
#   
# }
# 
# if (COMPARE_MSA_ALIGNERS) {
#   # 
#   # # Load required libraries
#   # library(msa)        # install.packages("msa")
#   # library(ips)        # install.packages("ips")
#   # library(Biostrings) # install.packages("Biostrings")
#   # library(seqinr)     # install.packages("seqinr")
#   # library(bioseq)
#   # 
#   # # Load sequences from a FASTA file
#   # seqs_aa <- readAAStringSet(file.path(temp_dir, "DLD_sequences.fasta")) ; head(seqs_aa)
#   # 
#   
#   # alignment_mafft   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_mafft.fasta"))
#   # 
#   # # Align sequences using MUSCLE
#   # alignment_muscle <- msa(seqs_aa, method = "Muscle")
#   # 
#   # # Align sequences using Clustal Omega
#   # alignment_clustalo <- msa(seqs_aa, method = "ClustalOmega")
#   # 
#   # 
#   # # Convert alignments to character vectors
#   # sequences_mafft <- as.character(alignment_mafft)
#   # sequences_muscle <- as.character(alignment_muscle)
#   # sequences_clustalo <- as.character(alignment_clustalo)
#   # 
#   # # Calculate conserved sites for each alignment
#   # conserved_mafft    <- calculate_conserved_sites(sequences_mafft)
#   # conserved_muscle   <- calculate_conserved_sites(sequences_muscle)
#   # conserved_clustalo <- calculate_conserved_sites(sequences_clustalo)
#   # 
#   # # Identify specific conserved sites for each alignment
#   # specific_mafft    <- identify_specific_conserved_sites(sequences_mafft)
#   # specific_muscle   <- identify_specific_conserved_sites(sequences_muscle)
#   # specific_clustalo <- identify_specific_conserved_sites(sequences_clustalo)
# 
# }
# 
# 
# # aln     <- AlignSeqs(seqs_aa, useStructures = TRUE) # BrowseSeqs(aln, highlight=0)
# # 
# # writeXStringSet(aln, file = file.path(temp_dir, "DLD_sequences_aln.fasta"))
# 
# ConsensusSequence(aln, threshold = 0.5)
# 
# 
# ### Generate Alignments with different aligment tools
# # Load required libraries
# library(DECIPHER) # install.packages("DECIPHER")
# 
# seqs_aa <- readAAStringSet(file.path(temp_dir, "DLD_sequences.fasta")) ; print(seqs_aa) # BrowseSeqs(seqs_aa, highlight=0)
# 
# # Align sequences using MAFFT
# if (file.exists(file.path(temp_dir, "DLD_sequences_aligned_mafft.fasta"))) { 
#   alignment_mafft   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_mafft.fasta"))
# } else {
#   if (OVERWRITE) {system(paste0("mafft --auto --inputorder --thread 4 ", file.path(temp_dir, "DLD_sequences.fasta"), " > ", file.path(temp_dir, "DLD_sequences_aligned_mafft.fasta")))}
#   alignment_mafft   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_mafft.fasta"))
# }
# 
# # Align sequences using MUSCLE
# if (file.exists(file.path(temp_dir, "DLD_sequences_aligned_muscle.fasta"))) { 
#   alignment_muscle   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_muscle.fasta"))
# } else {
#   alignment_muscle <- msa(seqs_aa, method = "Muscle")
#   if (OVERWRITE) {
#     writeAAMultipleAlignment(alignment_muscle, file.path(temp_dir, "DLD_sequences_aligned_muscle.fasta"))
#   }
#   alignment_muscle   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_muscle.fasta"))
# }
# 
# # Align sequences using Clustal Omega
# if (file.exists(file.path(temp_dir, "DLD_sequences_aligned_clustalo.fasta"))) { 
#   alignment_clustalo   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_clustalo.fasta"))
# } else {
#   if (OVERWRITE) {system(paste0("clustalo -i ", file.path(temp_dir, "DLD_sequences.fasta"), " -o ", file.path(temp_dir, "DLD_sequences_aligned_clustalo.fasta")))}
#   alignment_clustalo   <- readAAStringSet(file.path(temp_dir, "DLD_sequences_aligned_clustalo.fasta"))
# }
# 
# # Align sequences using Clustal W
# if (file.exists(file.path(temp_dir, "DLD_sequences_aligned_clustalw.fasta"))) { 
#   alignment_clustalw   <- readAAMultipleAlignment(file.path(temp_dir, "DLD_sequences_aligned_clustalw.fasta"))
# } else {
#   if (OVERWRITE) {system(paste0("clustalw -infile=", file.path(temp_dir, "DLD_sequences.fasta"), " -outfile=", file.path(temp_dir, "DLD_sequences_aligned_clustalw.fasta")))}
#   alignment_clustalw   <- readAAStringSet(file.path(temp_dir, "DLD_sequences_aligned_clustalw.fasta"))
# }
# 
# # Get consensus sequence from AAMultipleAlignment object
# print(consensus_mafft    <- ConsensusSequence(alignment_mafft, threshold = 0.5))
# print(consensus_muscle   <- ConsensusSequence(alignment_muscle, threshold = 0.5))
# print(consensus_clustalo <- ConsensusSequence(alignment_clustalo, threshold = 0.5))
# 
# # Function to calculate conserved sites
# calculate_conserved_sites <- function(sequences) {
#   # Create a matrix to store the number of conserved sites
#   conserved_sites <- matrix(0, nrow = ncol(sequences), ncol = nrow(sequences))
#   
#   # Loop through each pair of sequences
#   for (i in 1:(ncol(sequences) - 1)) {
#     for (j in (i + 1):ncol(sequences)) {
#       # Calculate the number of conserved sites between the two sequences
#       conserved_sites[i, j] <- sum(sequences[, i] == sequences[, j])
#     }
#   }
#   
#   # Return the matrix of conserved sites
#   return(conserved_sites)
# }


if (GENERATE_LCI_PLOTS) {
  
  NAME_KEY        <- fread("ELISA_CL_KEY.csv", header = T, na.strings = c("", "NA", NA)) %>% mutate(ID = str_replace(ID, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2"), CL_NAME_ON_PLOT = str_replace(CL_NAME_ON_PLOT, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2")) %>% arrange(ID)
  
  cell_table      <- fread(file.path(temp_dir, "lci_cell_table.csv.gz"))
  mean_cell       <- fread(file.path(temp_dir, "lci_mean_cell.csv.gz"))       %>% left_join(NAME_KEY[, c("COHORT", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOT_ID")], by = "COHORT", relationship = "many-to-many") %>% unique()
  mean_replicates <- fread(file.path(temp_dir, "lci_mean_replicates.csv.gz")) %>% left_join(NAME_KEY[, c("COHORT", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOT_ID")], by = "COHORT", relationship = "many-to-many") %>% unique()
  mean_total      <- fread(file.path(temp_dir, "lci_mean_total.csv.gz"))      %>% left_join(NAME_KEY[, c("COHORT", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOT_ID")], by = "COHORT", relationship = "many-to-many") %>% unique()
  
  mean_cell$CL_NAME_ON_PLOT       <- factor(reorder(mean_cell$CL_NAME_ON_PLOT , -mean_cell$ORDER_NO))
  mean_replicates$CL_NAME_ON_PLOT <- factor(reorder(mean_replicates$CL_NAME_ON_PLOT , -mean_replicates$ORDER_NO))
  mean_total$CL_NAME_ON_PLOT      <- factor(reorder(mean_total$CL_NAME_ON_PLOT , -mean_total$ORDER_NO))
  
  # Violin plot for long-lived TRAF6 ">=40 s" ----
  ggplot(mean_cell, aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT)) +
    geom_violin(data   = mean_cell       %>% filter(CATEGORY_DWELL_TIME == ">=40 s"), fill = "gray", scale = "width", width = 0.5, alpha = 0.3) +
    geom_jitter(data   = mean_replicates %>% filter(CATEGORY_DWELL_TIME == ">=40 s"), fill = "gray", size = 5, shape = 21, position = position_jitter(width = 0.2, height = 0, seed = 1000)) +
    geom_crossbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == ">=40 s"), aes(ymin = PCT_RECRUITMENT, ymax = PCT_RECRUITMENT), width = 0.2, fatten = 1) +
    geom_errorbar(data = mean_total      %>% filter(CATEGORY_DWELL_TIME == ">=40 s"), aes(ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT, ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT), width = 0.15) +
    coord_flip() +
    labs(x = "", y = "Percentage of long-lived TRAF6 on CHARMS") +
    theme_cowplot() +
    theme(legend.position ="none") 
  
  ggsave(
    "Pct of TRAF6 on CHARMS_greater than equal to 40s.pdf",
    #scale = 1,
    units = "cm",
    family = "Helvetica",
    height = 9,
    width = 7.2
  )
  
  # Violin plot for TRAF6 "<4 s" ----
  ggplot() +
    geom_violin(
      data = mean_cell %>% filter(CATEGORY_DWELL_TIME == "0 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT
      ),
      color = "black",
      fill = "#fff3e5",
      scale = "width",
      width = 0.5,
      alpha = 0.3
    ) +
    geom_jitter(
      data = Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == "0 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT
      ),
      color = "black",
      fill = "#fff3e5",
      size = 5,
      shape = 21,
      position = position_jitter(width = 0.2, height = 0, seed = 1150)
    ) +
    geom_crossbar(
      data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "0 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT,
        ymin = PCT_RECRUITMENT,
        ymax = PCT_RECRUITMENT
      ),
      width = 0.2,
      color = "black",
      fatten = 1
    ) +
    geom_errorbar(
      data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "0 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT,
        ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT,
        ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT
      ),
      width = 0.15,
      color = "black"
    ) +
    coord_flip() +
    scale_x_discrete(
      labels = c("1", "2", "3","4","5")
    ) +
    labs(
      x = "",
      y = "% of TRAF6 on CHARMS"
    ) +
    theme_classic() +
    theme(legend.position ="none") 
  
  ggsave(
    "Pct of TRAF6 on CHARMS_less than 4s.pdf",
    #scale = 1,
    units = "cm",
    family = "Helvetica",
    height = 12,
    width = 6.5
  )
  
  # Violin plot for TRAF6 "<=4 x < 40 s" ----
  ggplot() +
    geom_violin(
      data = mean_cell %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT
      ),
      color = "black",
      fill = "#FDEAA6",
      scale = "width",
      width = 0.5,
      alpha = 0.3
    ) +
    geom_jitter(
      data = Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT
      ),
      color = "black",
      fill = "#FDEAA6",
      size = 5,
      shape = 21,
      position = position_jitter(width = 0.2, height = 0, seed = 1228)
    ) +
    geom_crossbar(
      data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT,
        ymin = PCT_RECRUITMENT,
        ymax = PCT_RECRUITMENT
      ),
      width = 0.2,
      color = "black",
      fatten = 1
    ) +
    geom_errorbar(
      data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
      aes(
        x = COHORT,
        y = PCT_RECRUITMENT,
        ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT,
        ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT
      ),
      width = 0.15,
      color = "black"
    ) +
    coord_flip() +
    scale_x_discrete(
      labels = c("1", "2", "3","4","5")
    ) +
    labs(
      x = "",
      y = "% of TRAF6 on CHARMS"
    ) +
    theme_classic() +
    theme(legend.position ="none") 
  
  ggsave(
    "Pct of TRAF6 on CHARMS_from 4 to less than 40.pdf",
    #scale = 1,
    units = "cm",
    family = "Helvetica",
    height = 12,
    width = 6.5
  )
}


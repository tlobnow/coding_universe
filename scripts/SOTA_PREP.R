### Generating a State Of The Art (SOTA) Table on BDLD summary
LOAD_SETTINGS_AND_PATHS = TRUE
FEEDBACK                = TRUE # Set to TRUE to get feedback on the progress of the script

### TODO: Add feedback functionality to the script (instead of automatic print / head calls)

if (LOAD_SETTINGS_AND_PATHS) {
  setwd("/Users/u_lobnow/Documents/Github/coding_universe")
  LOAD_PATHS   <- T
  GET_OVERVIEW <- T
  OPEN_DATA    <- T
  OVERWRITE    <- T
  CLEAN        <- T
  JOIN         <- T
  RM_CLEANED_COLS         <- T
  GENERATE_GENOME_TABLE   <- T
  
  # GENERATE_ELISA_TABLE    <- T
  GENERATE_ELISA_TABLE    <- F
  
  GENERATE_AF3_JSON_TABLE <- T
  GENERATE_MSA_TABLE      <- T
  GENERATE_OPERON_TABLE   <- T
  DOWNLOAD_AND_PROCESS_GENOME_FOLDERS <- F
  EXTRACT_OPERONS              <- T
  GENERATE_SOTA_TABLE     <- T
  GENERATE_LCI_TABLE      <- F
  
  library(pacman)
  pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, 
                 readxl, scales, knitr, tidyr, visdat, purrr, fs, jsonlite,
                 rentrez, taxize, pwalign, Biostrings, plotrix, ape, rtracklayer,
                 read.gb) ; unloadNamespace("msa")
  
  # BiocManager::install("pwalign")
  # BiocManager::install("rtracklayer")
  
  reduce <- purrr::reduce
  select <- dplyr::select
  
  source("scripts/SOTA_FX.R")
  source("scripts/SOTA_PATHS.R")
  source("scripts/SOTA_COLNAMES.R") ; ls(pattern = "_cols")
}

##############################################################
##############################################################
##############################################################

if (OPEN_DATA) {
  # Load DD Library Data
  sheets         <- lapply(excel_files, read_excel_allsheets)         # Read all sheets from all Excel files
  sheets         <- purrr::flatten(sheets)                            # Flatten the list of lists into a single list of data frames
  sheets_with_id <- keep(sheets, ~ "ID" %in% colnames(.x))            # Filter out sheets that do not contain the 'ID' column
  dd_lib_tmp     <- purrr::reduce(sheets_with_id, full_join) %>% arrange(ID) # Join all data frames by the 'ID' column
  dd_lib_tmp     <- dd_lib_tmp %>% arrange(ID) %>% group_by(ID) %>% 
    fill(c(everything()), .direction = "downup") %>% 
    ungroup() %>% 
    filter(!is.na(taxend)) %>%
    distinct(ID, .keep_all = T)
  
  # Load DD Order Data
  dd_order_list <- lapply(order_files, function(file) {fread(order_files[grepl(file, order_files, ignore.case = TRUE)], na.strings = c("", "NA", NA))})
  dd_order_tmp  <- rbindlist(dd_order_list, fill = TRUE) %>% 
                    mutate(ID = str_replace(ID, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2")) %>% # Clean up the ID column
                    select(-c("n", "Image", "ID_LONG", "FASTA_bp_optimized"))
    
  # Load Experimental Data
  dates         <- c("20230313", "20230313_2", "20230912", "20240516")
  dd_lab_list   <- lapply(dates, function(date) {fread(csv_files[grepl(paste0(date, "_LAB.csv"), csv_files, ignore.case = TRUE)], na.strings = c("", "NA", "-", NA))})
  dd_lab_tmp    <- rbindlist(dd_lab_list, fill = TRUE) %>% 
                    mutate(ID = str_replace(ID, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2")) %>% # Clean up the ID column
                    select(-c("NAME_2", "Sequence Name")) %>% filter(grepl("DLD", ID)) %>% arrange(ID) %>% group_by(ID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(ID, col, .keep_all = T)
  
  # Load Cell Line Key
  NAME_KEY      <- fread("ELISA_CL_KEY.csv", header = T, na.strings = c("", "NA", NA)) %>%
                   mutate(ID = str_replace(ID, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2"),
                          CL_NAME_ON_PLOT = str_replace(CL_NAME_ON_PLOT, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2")) %>% arrange(ID)
  
  RH_DamFRET    <- fread(csv_files[grepl(paste0("_RH_DamFRET.csv"), csv_files, ignore.case = TRUE)], na.strings = c("", "NA", "-", NA))
 
  # declutter environment and remove all objects that are no longer needed
  rm(sheets, sheets_with_id, dd_order_list, dates, dd_lab_list)
}

if (CLEAN) {
  
  # Join  DD Library and Experimental Data
  SOTA_to_clean <- full_join(dd_lab_tmp, dd_lib_tmp, relationship = "many-to-many") %>%
    full_join(dd_order_tmp) %>% filter(grepl("DLD", ID)) %>% arrange(ID) %>% 
    arrange(ID) %>% group_by(ID) %>% fill(c(everything()), .direction = "downup") %>% 
    distinct(ID, col, .keep_all = T) %>%
    ungroup() %>%
    select(-c(taxid, FASTA_bp_optimized, BP_LEN)) %>%
    
    dplyr::rename(
      gene              = Abbr,
      genbank_acc       = NCBI_Abbr,
      genbank_acc_2     = NCBI_Abbr_2,
      aa_start          = AA_START,
      aa_end            = AA_END,
      architecture      = ARCHITECTURE,
      genbank_gca       = GenBank_GCA,
      genbank_gca_2     = GenBank_GCA_new,
      refseq_gcf        = RefSeq_GCF,
      refseq_gcf_2      = RefSeq_GCF_new,
      abs_260_280       = `260/280`,
      af_confidence     = Confidence,
      annotated_domains = `ANNOTATED_DOMAINS NCBI [UniProt]`,
      bp_len            = `Sequence Length`,
      calc_mol_wt       = `Calculated Molecular Weight`,
      quantity          = Quantity,
      tirf_coalesce     = COALESCE,
      colony            = col,
      conc_ng_ul        = `Conc. [ng/µL]`,
      conc_ug_ul        = `Conc. [µg/µL]`,
      elisa_stim_4      = `ELISA #4 DATE`,
      elisa_stim_5      = `ELISA #5 DATE`,
      elisa_stim_6      = `ELISA #6 DATE`,
      elisa_assay       = `ELISA ASSAY DATE`,
      elisa_assay_2     = `ELISA ASSAY DATE 2`,
      fm_ng             = `fm/ng`,
      gblock_seq_nt     = `Sequence (nt)`,
      gblock_sales_no   = `Sales Order Number`,
      gblock_man_id     = `Manufacturing ID`,
      seq_DLD_nt        = FASTA_bp_stripped,
      env_temp          = `Likely Growth Properties`,
      operon            = OPERON,
      reference         = Reference,
      taxid             = TAXID,
      ref_number        = `Reference Number`,
      terminus          = Terminus,
      potential_fct     = POTENTIAL_FCT,
      tirf_lci_1        = `TIRF LCI DATE #1`,
      tirf_lci_2        = `TIRF LCI DATE #2`,
      tirf_lci_3        = `TIRF LCI DATE #3`,
      frozen_date       = `FROZEN DATE`,
      frozen_date_us    = `UNSORTED FROZEN DATE`,
      frap_date         = `FRAP DATE`,
      facs_date         = `FACS Date`,
      facs_yield        = `FACS Yield`,
      tirf_rela_date    = `TIRF RelA DATE`,
      tirf_traf6_puncta = `TRAF6 puncta`,
      tirf_gfp_puncta   = `GFP puncta`,
      plasmid_volume    = `Volume needed (need 6.3 total µg plasmid DNA + 2.1 µg Plasmid)`,
      td_cl204          = `EL-4 cl204 Transduction`,
      td_date           = `Transduction Date`,
      gene_suppr        = Suppressed_Gene,
      genome_suppr      = Suppressed_Genome
      
      ) %>% group_by(ID) %>% arrange(ID) %>%
    
    mutate(uniprot_acc    = coalesce(UNIPROT_ACCESSION, UNIPROT_ID, Accession, ACCESSION),
           info           = coalesce(Info, INFO),
           pmjtl          = coalesce(pMJTL, pMJTL_no),
           genbank_acc    = case_when(is.na(genbank_acc) ~ GenBank, T ~ genbank_acc), 
           uniprot_acc    = case_when(uniprot_acc == genbank_acc ~ NA, T ~ uniprot_acc),
           species        = coalesce(Species, SPECIES),
           seq_DLD_aa     = coalesce(BDLD_SEQUENCE, Sequence, `Full Sequence (aa, not prepped as gblock)`, FASTA_aa_stripped),
           seq_FL_aa      = coalesce(Sequence_FL, FASTA_aa_FL),
           elisa_stim_1   = coalesce(`ELISA #1 DATE`, `ELISA STIM #1`),
           elisa_stim_2   = coalesce(`ELISA #2 DATE`, `ELISA STIM #2`),
           elisa_stim_3   = coalesce(`ELISA #3 DATE`, `ELISA STIM #3`),
           note           = coalesce(Note, NOTE),
           env            = coalesce(Environment, `Environment isolated from`),
           tf_date        = coalesce(`Transfection Date`, `HEK Transfection`),
           mag            = as.logical(MAG),
           gblock_norm    = "250ng in 25µl of RNase-Free Water"
           ) %>%
    select(-c("UNIPROT_ACCESSION", "UNIPROT_ID", "Accession", "ACCESSION", 
              "AA_LEN", "Abbr2", "...26", "BP_LEN_2", "Info", "INFO", "FASTA", "ELISA_STATUS", "STATUS",
              "GenBank", "pMJTL", "pMJTL_no", "DALI", "FILE", "NAME", "NAME_PASTED", "Sequence Name",
              "ORDER_seq_pasted", "ORDER_seq_pasted_raw", "Organism", "PDB_ID", "Product", "Randal Halfmann",
              "GENOME_AVAILABLE", "GENOME_LEVEL", "GENOME_NOTE", "MAG", "Normalization", "PROTEIN_ANNOTATION",
              "BDLD_SEQUENCE", "Sequence", "Full Sequence (aa, not prepped as gblock)", "FASTA_aa_stripped",
              "ELISA #1 DATE", "ELISA STIM #1", "ELISA #2 DATE", "ELISA STIM #2", "ELISA #3 DATE", "ELISA STIM #3",
              "Construct", "taxend2", "Species", "SPECIES", "Sequence_FL", "FASTA_aa_FL", "Note", "NOTE",
              "NCBI_ACCESSION", "NCBI_ACC", "FASTA_aa", "FASTA_bp", "BDLD_ID", "FoldSeek_Downloaded", "Expression",
              "Environment", "Environment isolated from", "Transfection Date", "HEK Transfection")) %>%
    distinct(ID, .keep_all = T) %>% ungroup()
  
  cl_name_df    <- NAME_KEY      %>% select(ID, EXPRESSION_LVL, CL_NUMBER) %>% filter(!is.na(ID)) %>% unique() %>% pivot_wider(id_cols = ID, names_from = EXPRESSION_LVL, values_from = CL_NUMBER) %>% select(ID, ENDO, HIGH) ; setnames(cl_name_df, c("ENDO", "HIGH"), c("cl_name_endo", "cl_name_high")); SOTA_to_clean <- SOTA_to_clean %>% select(-c(colnames(SOTA_to_clean[grepl("cell", colnames(SOTA_to_clean))]), `CELL LINE`))
  taxid_df_init <- SOTA_to_clean %>% select(ID, taxid) ; if (file.exists(file.path(temp_dir, "taxend.csv"))) { taxend_df <- fread(file.path(temp_dir, "taxend.csv"), na.strings = c("", "NA", NA)) } else { taxend_df <- taxid_df_init %>% rowwise() %>% mutate(tax_info = list(safe_get_taxonomy_df(taxid)$result)) %>% unnest_wider(c(tax_info)) %>% ungroup() %>% unique() }
  taxend_df     <- taxend_df %>% unnest(cols = colnames(taxend_df)) ; fwrite(taxend_df, file.path(temp_dir, "taxend.csv"))
  taxid_df      <- taxid_df_init %>% left_join(taxend_df); rm(taxid_df_init, taxend_df)
  genbank_ncbi_df     <- SOTA_to_clean %>% select(ID, genbank_acc) %>% left_join(if (file.exists(file.path(temp_dir, "genbank.csv"))) {gb_name_df <- fread(file.path(temp_dir, "genbank.csv"), na.strings = c("", "NA", NA))} else {gb_name_df <- genbank_ncbi_df %>% mutate(gene_name = sapply(genbank_acc, get_name_from_ncbi)) %>% select(genbank_acc, gene_name) %>% unique() ; fwrite(gb_name_df, file.path(temp_dir, "genbank.csv"))})
  damFRET_df          <- SOTA_to_clean %>% select(ID, genbank_acc) %>% left_join(RH_DamFRET, relationship = "many-to-many") ; fwrite(damFRET_df, file.path(temp_dir, "taylor_x_halfmann.csv")) ; damFRET_df <- damFRET_df %>% mutate(damfret_rh = case_when(!is.na(expression_clone) ~ TRUE, TRUE ~ FALSE)) %>% select(ID, damfret_rh) %>% distinct(ID, .keep_all = T)
  
  # initiate clean SOTA DF
  dfs_to_join <- list(cl_name_df, taxid_df, genbank_ncbi_df, damFRET_df)
  
  # Create the initial data frame
  SOTA_final <- SOTA_to_clean %>% 
    select(any_of(c(
      basics_cols,
      basics_extended_cols,
      accession_cols,
      environment_cols,
      sequence_cols,
      annotation_cols,
      alphaFold_cols,
      gBlock_cols,
      cl_gen_cols,
      elisa_cols,
      lci_cols,
      miscellaneous,
      ncbi_genome_cols
    ))) %>% 
    mutate(aa_len = aa_end - aa_start + 1,
           bp_len = case_when(!is.na(bp_len) ~ bp_len, TRUE ~ aa_len*3)) %>%
    unique()
  
  SOTA_final <- reduce(dfs_to_join, left_join, .init = SOTA_final) %>% 
    select(any_of(unique(c(basics_cols, basics_extended_cols, accession_cols, environment_cols, sequence_cols, annotation_cols, 
                           alphaFold_cols, gBlock_cols, cl_gen_cols, elisa_cols, lci_cols, miscellaneous))), everything())
  
  # remove added data frames to de-clutter environment
  for (var in c("cl_name_df, taxid_df")) {if (exists(var)) {rm(list = var)}}
  
  if (GET_OVERVIEW) {
    
    vis_miss(SOTA_to_clean[, setdiff(colnames(SOTA_to_clean), colnames(SOTA_final)), drop = FALSE])
    
    # nms_SOTA_to_clean <- as.data.frame(sort(names(SOTA_to_clean))) %>%
    #   unique() %>%
    #   rename(Column_Name = 1) %>%
    #   filter(!Column_Name %in% names(SOTA_final))
    # 
    # nms_SOTA_cleaned <- as.data.frame(sort(names(SOTA_final))) %>%
    #   unique() %>%
    #   rename(Column_Name = 1)
  }
  
  # Save files
  if (OVERWRITE) {
    
    # save SOTA_to_clean as CSV
    fwrite(SOTA_to_clean, file = file.path(temp_dir, "SOTA_to_clean.csv"), append = F, row.names = F, col.names = T)
    print(paste0("Overwritten: ", file.path(temp_dir, "SOTA_to_clean.csv")))
    SOTA_to_clean <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA_to_clean.csv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))

    # save SOTA_final as CSV
    fwrite(SOTA_final, file = file.path(temp_dir, "SOTA.csv"), append = F, row.names = F, col.names = T)
    print(paste0("Overwritten: ", file.path(temp_dir, "SOTA.csv")))
    SOTA_final <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA.csv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))

  } else {
    # or simply load pre-existing data if you don't want to overwrite
    print("Not overwriting existing CSV file. Set OVERWRITE = TRUE to generate new SOTA CSV summary files!")
    SOTA_to_clean <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA_to_clean.csv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
    SOTA_final    <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA.csv$",          ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
  }
  
  # declutter environment and remove all objects that are no longer needed
  # rm(dd_lab_tmp, dd_lib_tmp, dd_order_tmp, cl_name_df, taxid_df)
  # charms_vars <- grep("^CHARMS", ls(), value = TRUE) ; rm(list = charms_vars)
}

if (GENERATE_GENOME_TABLE) {
  # Load NCBI Genome Info Data
  ncbi_genome_data <- fread(file.path(list.files(input_dir, pattern = "\\NCBI_GENOMES.tsv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
  setnames(ncbi_genome_data, "Assembly Accession", "assembly")
  
  # Join Library Data with ncbi genome information
  genome_info_tbl <- SOTA_to_clean[,c("genbank_gca", "genbank_gca_2", "refseq_gcf", "refseq_gcf_2", "ID")] %>% 
    pivot_longer(cols = 1:4, names_to = "ORIGINAL_COL", values_to = "assembly") %>%
    filter(!is.na(assembly), assembly != "NA") %>%
    select(-ORIGINAL_COL) %>%
    left_join(ncbi_genome_data) %>%
    select_if(~ !all(is.na(.))) %>%
    unique()
  
  if (OVERWRITE) {
    # save NCBI Genome Info Data as CSV
    fwrite(genome_info_tbl, file = file.path(output_dir, "genome_info.csv"), append = F, row.names = F, col.names = T)
  } else {print("Not overwriting existing CSV file. Set OVERWRITE = TRUE to generate new SOTA CSV summary files!")}
  
  rm(ncbi_genome_data)
}

if (GENERATE_ELISA_TABLE) {
  run_settings_and_prep     <- TRUE
  run_processing_and_subset <- TRUE
  prep_and_process          <- TRUE
  plot_elisa                <- TRUE
  join_elisa_df             <- TRUE

  if (run_settings_and_prep) {
    
    # GENERAL SETTINGS
    figure <- "SOTA"
    SAVE   <- TRUE
    
    RELATIVE_SECRETION      <- TRUE
    FOLD_CHANGE             <- TRUE
    
    # PLOT SETTINGS
    FONT   <- "Helvetica"
    SIZE   <- 25
    POINTS <- 6
    TEXT   <- 8
    
    
    # GATHER DATA & FUNCTIONS
    Input_Directory <- elisa_dir
    
    # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
    plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory = temp_dir)
    plate_data     <- left_join(plate_data_raw, NAME_KEY, relationship = "many-to-many")
    plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
    
    # ensure correct column type assignment
    plate_data$CONDITION <- as.factor(plate_data$CONDITION)
    
    # streamline ID (BDLD_6 --> BDLD_06)
    plate_data <- plate_data %>% mutate(CL_NAME_ON_PLOT = str_replace(CL_NAME_ON_PLOT, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2"))
  }
  
  if (prep_and_process) {
    
    
    # Applying the function to each set of parameters
    pre_CHARMS        <- plate_data %>% filter(Date == "2022-06-23" |  Date == "2022-07-01") %>% mutate(ID = "CHARMS")
    CHARMS_norm_3E10  <- process_ELISA_data(DF = pre_CHARMS, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "3E10") %>%
                          filter(CL_NAME_ON_PLOT %in% c("CHARMS", "3E10"), Date == "2022-06-23") %>%
                          mutate(DF = "pre_CHARMS", NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "3E10")
    
    CHARMS_norm_cl069 <- process_ELISA_data(DF = pre_CHARMS, NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069") %>%
      filter(CL_NAME_ON_PLOT %in% c("CHARMS", "cl069"), Date == "2022-06-23") %>%
      mutate(DF = "pre_CHARMS", NEGATIVE_CTRL = "tKO_EL4", POSITIVE_CTRL = "cl069")
    
    rm(pre_CHARMS)
    
    
    CHARMS_BDD_02_H_norm   <- process_plate_data(plate_data, "2023-06-21", TRUE,                 "BDLD_02", "pre_BDLD_02_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_06_H_norm   <- process_plate_data(plate_data, "2023-07-04", Plate == 3,           "BDLD_06", "pre_BDLD_06_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_07_H_norm   <- process_plate_data(plate_data, "2023-07-04", Plate != 3,           "BDLD_07", "pre_BDLD_07_H",   "tKO_EL4", "3E10")
    CHARMS_RANDOM_norm     <- process_plate_data(plate_data, "2023-06-21", TRUE,                 "RANDOM",  "pre_RANDOM",      "tKO_EL4", "3E10")
    CHARMS_BDD_11_H_norm   <- process_plate_data(plate_data, "2024-03-04", Plate != 3,           "BDLD_11", "pre_BDLD_11_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_13_H_norm   <- process_plate_data(plate_data, "2023-05-29", TRUE,                 "BDLD_13", "pre_BDLD_13_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_14_E_norm   <- process_plate_data(plate_data, "2023-05-04", TRUE,                 "BDLD_14", "pre_BDLD_14_E",   "tKO_EL4", "3E10")
    CHARMS_BDD_14_H_norm   <- process_plate_data(plate_data, "2023-05-12", TRUE,                 "BDLD_14", "pre_BDLD_14_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_16_H_norm   <- process_plate_data(plate_data, "2024-03-04", Plate != 3,           "BDLD_16", "pre_BDLD_16_H",   "tKO_EL4", "3E10")
    # CHARMS_BDD_27_E_norm   <- process_plate_data(plate_data, "2023-05-04", TRUE,                 "BDLD_27", "pre_BDLD_27_E",   "tKO_EL4", "3E10")
    CHARMS_BDD_27_H_norm   <- process_plate_data(plate_data, "2023-05-12", TRUE,                 "BDLD_27", "pre_BDLD_27_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_38_H_norm   <- process_plate_data(plate_data, "2024-03-04", Plate != 3,           "BDLD_38", "pre_BDLD_38_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_46_H_norm   <- process_plate_data(plate_data, "2023-07-04", Plate != 3,           "BDLD_46", "pre_BDLD_46_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_50_E_norm   <- process_plate_data(plate_data, "2023-05-04", TRUE,                 "BDLD_50", "pre_BDLD_50_E",   "tKO_EL4", "3E10")
    CHARMS_BDD_50_H_norm   <- process_plate_data(plate_data, "2023-05-12", TRUE,                 "BDLD_50", "pre_BDLD_50_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_57_H_norm   <- process_plate_data(plate_data, "2023-06-21", TRUE,                 "BDLD_57", "pre_BDLD_57_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_57_H_2_norm <- process_plate_data(plate_data, "2024-03-04", Plate == 3,           "BDLD_57", "pre_BDLD_57_H_2", "tKO_EL4", "cl069")
    CHARMS_BDD_58_H_norm   <- process_plate_data(plate_data, "2024-03-04", Plate != 3,           "BDLD_58", "pre_BDLD_58_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_62_H_norm   <- process_plate_data(plate_data, "2023-06-21", TRUE,                 "BDLD_62", "pre_BDLD_62_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_67_H_norm   <- process_plate_data(plate_data, "2023-07-04", Plate == 3,           "BDLD_67", "pre_BDLD_67_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_69_H_norm   <- process_plate_data(plate_data, "2023-07-04", Plate != 3,           "BDLD_69", "pre_BDLD_69_H",   "tKO_EL4", "3E10")
    CHARMS_BDD_73_H_norm   <- process_plate_data(plate_data, "2023-07-04", Plate != 3,           "BDLD_73", "pre_BDLD_73_H",   "tKO_EL4", "3E10")
    CHARMS_EDD_01_E_norm   <- process_plate_data(plate_data, "2024-01-31", Plate == 3,           "EDLD_01", "pre_EDLD_01_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_01_E_2_norm <- process_plate_data(plate_data, "2024-02-08", Plate == 2,           "EDLD_01", "pre_EDLD_01_E_2", "tKO_EL4", "cl069")
    CHARMS_EDD_03_E_norm   <- process_plate_data(plate_data, "2024-02-08", Plate == 2,           "EDLD_03", "pre_EDLD_03_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_06_E_norm   <- process_plate_data(plate_data, "2024-01-31", Plate == 3,           "EDLD_06", "pre_EDLD_06_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_09_E_norm   <- process_plate_data(plate_data, "2024-01-31", Plate == 3,           "EDLD_09", "pre_EDLD_09_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_25_E_norm   <- process_plate_data(plate_data, "2024-05-17", CELL_LINE != "cl456", "EDLD_25", "pre_EDLD_25_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_25_H_norm   <- process_plate_data(plate_data, "2024-05-17", CELL_LINE != "cl455", "EDLD_25", "pre_EDLD_25_H",   "tKO_EL4", "cl069")
    CHARMS_EDD_26_E_norm   <- process_plate_data(plate_data, "2024-05-17", CELL_LINE != "cl458", "EDLD_26", "pre_EDLD_26_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_26_H_norm   <- process_plate_data(plate_data, "2024-05-17", CELL_LINE != "cl457", "EDLD_26", "pre_EDLD_26_H",   "tKO_EL4", "cl069")
    CHARMS_EDD_27_E_norm   <- process_plate_data(plate_data, "2024-05-17", CELL_LINE != "cl460", "EDLD_27", "pre_EDLD_27_E",   "tKO_EL4", "cl069")
    CHARMS_EDD_27_H_norm   <- process_plate_data(plate_data, "2024-05-17", CELL_LINE != "cl459", "EDLD_27", "pre_EDLD_27_H",   "tKO_EL4", "cl069")
    
    # Combine all normalized datasets into one
    plotting_data <- bind_rows(
      CHARMS_norm_3E10, CHARMS_norm_cl069,
      CHARMS_BDD_02_H_norm, CHARMS_BDD_06_H_norm, CHARMS_BDD_07_H_norm, CHARMS_RANDOM_norm,
      CHARMS_BDD_11_H_norm, CHARMS_BDD_13_H_norm, CHARMS_BDD_14_E_norm, CHARMS_BDD_14_H_norm,
      CHARMS_BDD_16_H_norm, #CHARMS_BDD_27_E_norm, 
      CHARMS_BDD_27_H_norm, CHARMS_BDD_38_H_norm,
      CHARMS_BDD_46_H_norm, CHARMS_BDD_50_E_norm, CHARMS_BDD_50_H_norm, CHARMS_BDD_57_H_norm,
      CHARMS_BDD_57_H_2_norm, CHARMS_BDD_58_H_norm, CHARMS_BDD_62_H_norm, CHARMS_BDD_67_H_norm,
      CHARMS_BDD_69_H_norm, CHARMS_BDD_73_H_norm, CHARMS_EDD_01_E_norm, CHARMS_EDD_01_E_2_norm,
      CHARMS_EDD_03_E_norm, CHARMS_EDD_06_E_norm, CHARMS_EDD_09_E_norm, CHARMS_EDD_25_E_norm,
      CHARMS_EDD_25_H_norm, CHARMS_EDD_26_E_norm, CHARMS_EDD_26_H_norm, CHARMS_EDD_27_E_norm,
      CHARMS_EDD_27_H_norm
    ) %>% left_join(NAME_KEY)
    
    charms_vars <- grep("^CHARMS", ls(), value = TRUE) ; rm(list = charms_vars)
  }
  
  controls <- c("cl069", "3E10")
  for (control in controls) {
    POS_CTRL = control
    plotting_data_main <- process_data_for_plot(plotting_data) %>% filter(POSITIVE_CTRL == POS_CTRL)
    plotting_means     <- prepare_plotting_means(data = plotting_data_main, group_var = c("PLOT_ID", "CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "Date")) %>% left_join(NAME_KEY) %>% 
      fill_columns_by_id(id = "PLOT_ID", target_cols = c("CL_NUMBER", "ORIGIN", "ID", "EXPRESSION_LVL", "PURPOSE")) %>%
      mutate(POSITIVE_CTRL = POS_CTRL)
    
    stat_significance_dt_rltv <- process_statistical_analysis(data = plotting_means, group_var = "PLOT_ID", value_var = "Relative_Intensity_mean")
    setnames(stat_significance_dt_rltv, "CL_NAME_ON_PLOT", "PLOT_ID")
    stat_significance_dt_rltv <- stat_significance_dt_rltv %>%
      left_join(NAME_KEY) %>%
      left_join(plotting_means[colnames(plotting_means) %in% c("PLOT_ID", "Date", "POSITIVE_CTRL")], 
                relationship = "many-to-many") %>%
      unique()
    
    stat_significance_dt_real <- process_statistical_analysis(data = plotting_means, group_var = "PLOT_ID", value_var = "IL2_concentration_Dilution_Factor_mean")
    setnames(stat_significance_dt_real, "CL_NAME_ON_PLOT", "PLOT_ID")
    stat_significance_dt_real <- stat_significance_dt_real %>%
      left_join(NAME_KEY) %>%
      left_join(plotting_means[colnames(plotting_means) %in% c("PLOT_ID", "Date", "POSITIVE_CTRL")], 
                relationship = "many-to-many") %>%
      unique()
    
    
    plotting_stats_rltv <- prepare_plotting_stats(plotting_data_main, stat_significance_dt_rltv, group_var = c("PLOT_ID", "CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "Date", "POSITIVE_CTRL")) %>%
      select(-c("ID", "CL_NUMBER", "PURPOSE", "INFO", "ORIGIN", "EXPRESSION_LVL")) %>%
      left_join(NAME_KEY) %>%
      fill_columns_by_id(id = "PLOT_ID", target_cols = c("ID", "CL_NUMBER", "PURPOSE", "INFO", "ORIGIN", "EXPRESSION_LVL", "p_value", "significance"))
    
    plotting_stats_real <- prepare_plotting_stats(plotting_data_main, stat_significance_dt_real, group_var = c("PLOT_ID", "CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "Date", "POSITIVE_CTRL")) %>%
      select(-c("ID", "CL_NUMBER", "PURPOSE", "INFO", "ORIGIN", "EXPRESSION_LVL")) %>%
      left_join(NAME_KEY) %>%
      fill_columns_by_id(id = "PLOT_ID", target_cols = c("ID", "CL_NUMBER", "PURPOSE", "INFO", "ORIGIN", "EXPRESSION_LVL", "p_value", "significance"))
    
    
    # reordering
    plotting_means$PLOT_ID <- factor(reorder(plotting_means$PLOT_ID , -plotting_means$ORDER_NO))
    plotting_stats_rltv$PLOT_ID <- factor(reorder(plotting_stats_rltv$PLOT_ID , -plotting_stats_rltv$ORDER_NO))
    plotting_stats_real$PLOT_ID <- factor(reorder(plotting_stats_real$PLOT_ID , -plotting_stats_real$ORDER_NO))
    
    ################################################################################
    
    # plotting
    if (plot_elisa) {
      figure_ELISA_rltv <- ggplot(data = plotting_stats_rltv, aes(x = Relative_Intensity_mean, y = PLOT_ID, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
        geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
        geom_errorbar(aes(y = PLOT_ID,
                          xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                          xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                      linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
        geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = PLOT_ID, fill = PLOTTING_COLOR), col = "black",
                   shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
        geom_text(data = plotting_stats_rltv, aes(x = max(Relative_Intensity_mean) + 0.5, y = PLOT_ID, label = significance_rltv),
                  hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
        scale_x_continuous(breaks = seq(from = 0, to = 1.2, by = 0.5), position = "top") +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_manual(name  = "PLOT_ID", values = plotting_stats_rltv$PLOTTING_COLOR, breaks = plotting_stats_rltv$PLOTTING_COLOR, labels = ifelse(plotting_stats_rltv$CONDITION == "UNSTIM", paste0("- ", plotting_stats_rltv$STIMULANT), paste0("+ ", plotting_stats_rltv$STIMULANT))) +
        scale_color_manual(name = "PLOT_ID", values = plotting_stats_rltv$PLOTTING_COLOR, breaks = plotting_stats_rltv$PLOTTING_COLOR, labels = ifelse(plotting_stats_rltv$CONDITION == "UNSTIM", paste0("- ", plotting_stats_rltv$STIMULANT), paste0("+ ", plotting_stats_rltv$STIMULANT))) +
        labs(x = "Relative IL-2 secretion", y = "", subtitle = paste0("Normalized to: ", unique(plotting_stats_rltv$POSITIVE_CTRL))) +
        guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
        theme_cowplot(font_size = SIZE, font_family = FONT) +
        theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
              axis.title.y      = element_blank(),
              legend.position   = "bottom",
              legend.title      = element_blank(),
              legend.text       = element_text(size = SIZE),
              legend.key.size   = unit(9, "mm"))
      
      figure_ELISA_rltv
      
      figure_ELISA_real <- ggplot(data = plotting_stats_real, aes(x = IL2_concentration_Dilution_Factor_mean, y = PLOT_ID, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
        geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
        geom_errorbar(aes(y = PLOT_ID,
                          xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                          xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                      linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
        geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = PLOT_ID, fill = PLOTTING_COLOR), col = "black",
                   shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
        geom_text(data = plotting_stats_real, aes(x = max(IL2_concentration_Dilution_Factor_mean) + 0.5, y = PLOT_ID, label = significance_real),
                  hjust = .5, vjust = 1, size = TEXT, angle = 90) + 
        scale_x_continuous(breaks = seq(from = 0, to = 1.2, by = 0.5), position = "top") +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_manual(name  = "PLOT_ID", values = plotting_stats_real$PLOTTING_COLOR, breaks = plotting_stats_real$PLOTTING_COLOR, labels = ifelse(plotting_stats_real$CONDITION == "UNSTIM", paste0("- ", plotting_stats_real$STIMULANT), paste0("+ ", plotting_stats_real$STIMULANT))) +
        scale_color_manual(name = "PLOT_ID", values = plotting_stats_real$PLOTTING_COLOR, breaks = plotting_stats_real$PLOTTING_COLOR, labels = ifelse(plotting_stats_real$CONDITION == "UNSTIM", paste0("- ", plotting_stats_real$STIMULANT), paste0("+ ", plotting_stats_real$STIMULANT))) +
        labs(x = "Real IL-2 secretion", y = "", subtitle = paste0("Normalized to: ", unique(plotting_stats_real$POSITIVE_CTRL))) +
        guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
        theme_cowplot(font_size = SIZE, font_family = FONT) +
        theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
              axis.title.y      = element_blank(),
              legend.position   = "bottom",
              legend.title      = element_blank(),
              legend.text       = element_text(size = SIZE),
              legend.key.size   = unit(9, "mm"))
      
      figure_ELISA_real
      
      }
    
    ################################################################################
    
    # rename 
    setnames(plotting_stats_rltv, c("p_value", "significance"), c("p_value_rltv", "significance_rltv"))
    setnames(plotting_stats_real, c("p_value", "significance"), c("p_value_real", "significance_real"))
    
    # generate a condensed (pivoted) table of plotting stats per DLD:
    piv_stats <- plotting_stats_rltv %>% left_join(plotting_stats_real) %>% 
      select(-c(PLOT_ID, CELL_LINE, CL_NAME_ON_PLOT, PLOTTING_COLOR, ORDER_NO, PURPOSE, INFO)) %>%
      # filter(!is.na(ID)) %>%
      pivot_wider(names_from  = c(CONDITION), 
                  values_from = c(IL2_concentration_Dilution_Factor_mean, IL2_concentration_Dilution_Factor_sem, Relative_Intensity_mean, Relative_Intensity_sem)) %>%
      select(ID, CL_NUMBER, ORIGIN, EXPRESSION_LVL, POSITIVE_CTRL, Date, 
             p_value_rltv, p_value_real, significance_rltv, significance_real,
             IL2_concentration_Dilution_Factor_mean_UNSTIM, IL2_concentration_Dilution_Factor_mean_STIM,
             IL2_concentration_Dilution_Factor_sem_UNSTIM, IL2_concentration_Dilution_Factor_sem_STIM,
             Relative_Intensity_mean_UNSTIM, Relative_Intensity_mean_STIM,
             Relative_Intensity_sem_UNSTIM, Relative_Intensity_sem_STIM,
             PATHWAY, STIMULANT, STIM_CONCENTRATION)
    
    if (OVERWRITE) {
      
      if (dir.exists(temp_dir))   {print(paste0("Files will be saved to ", temp_dir))}   else {dir.create(temp_dir,   recursive = T); print(paste0("Files will be saved to ", temp_dir))}
      if (dir.exists(output_dir)) {print(paste0("Files will be saved to ", output_dir))} else {dir.create(output_dir, recursive = T); print(paste0("Files will be saved to ", output_dir))}
      
      # save figure
      if (plot_elisa) {
        ggsave(paste0(temp_dir, "/FIGURES/", "ELISA_rltv_secretion_", control,"_normalized.svg"), plot = figure_ELISA_rltv, device = "svg", width = 12, height = 24)
        ggsave(paste0(temp_dir, "/FIGURES/", "ELISA_real_secretion_", control,"_normalized.svg"), plot = figure_ELISA_real, device = "svg", width = 12, height = 24)
      }
      
      # save tables
      fwrite(plate_data,                paste0(temp_dir, "/", "plate_data_",                control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "plate_data_",                control,"_normalized.csv"))
      fwrite(plotting_data,             paste0(temp_dir, "/", "plotting_data_",             control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "plotting_data_",             control,"_normalized.csv"))
      fwrite(plotting_means,            paste0(temp_dir, "/", "plotting_means_",            control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "plotting_means_",            control,"_normalized.csv"))
      fwrite(plotting_stats_rltv,       paste0(temp_dir, "/", "plotting_stats_rltv_",       control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "plotting_stats_rltv_",       control,"_normalized.csv"))
      fwrite(plotting_stats_real,       paste0(temp_dir, "/", "plotting_stats_real_",       control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "plotting_stats_real_",       control,"_normalized.csv"))
      fwrite(stat_significance_dt_rltv, paste0(temp_dir, "/", "stat_significance_dt_rltv_", control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "stat_significance_dt_rltv_", control,"_normalized.csv"))
      fwrite(stat_significance_dt_real, paste0(temp_dir, "/", "stat_significance_dt_real_", control,"_normalized.csv"), append = F) ; print(paste0("Overwritten: ", temp_dir, "/", "stat_significance_dt_real_", control,"_normalized.csv"))
      
      fwrite(piv_stats,                 paste0(temp_dir, "/", "ELISA_SUMMARY_", control,"_normalized.csv"), append = F)
    }
  }
  
  if (FOLD_CHANGE) {
    # Select Data
    COHORT_DATA <- plotting_data
    
    # Set negative values to 1
    COHORT_DATA$MEASUREMENT <- ifelse(COHORT_DATA$MEASUREMENT < 0, yes = 0, COHORT_DATA$MEASUREMENT)
    
    # Normalization for fold change from UNSTIM to STIM
    NORMALIZED_TO_CONTROL <- COHORT_DATA %>%
      group_by(Date, STIM_DAY, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PLOTTING_COLOR, PLOT_ID) %>%
      summarise(mean_per_day = mean(MEASUREMENT),
                Concentration = mean(Concentration))
    
    # First, calculate the mean for the UNSTIM condition separately.
    unstim_means <- COHORT_DATA %>%
      filter(CONDITION == "UNSTIM") %>%
      group_by(CELL_LINE) %>%
      summarise(mean_unstim = ifelse(mean(MEASUREMENT) == 0, 1, mean(MEASUREMENT)))
    
    # Now, join this back to the main dataset.
    COHORT_SUBSET <- COHORT_DATA %>%
      left_join(unstim_means, by = "CELL_LINE")
    
    # Compute the fold change. 
    # The `case_when` logic ensures that the fold change is calculated only for the STIM condition.
    NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
      mutate(fold_change = case_when(CONDITION == "STIM" ~ MEASUREMENT / mean_unstim, TRUE ~ NA_real_ )) %>% # set NA for non-STIM conditions
      ungroup() %>%
      group_by(CELL_LINE, CONDITION, PLOTTING_COLOR) %>%
      mutate(trip_mean = mean(fold_change),
             fold_change_sd = sd(fold_change)) %>%
      unique()
    
    # Helper function to perform a t-test and return p-value and annotation
    perform_ttest <- function(data) {
      ttest_result <- t.test(data$MEASUREMENT[data$CONDITION == "STIM"],
                             data$MEASUREMENT[data$CONDITION == "UNSTIM"],
                             paired = F)
      p_value <- ttest_result$p.value
      annotation <- ifelse(p_value < 0.0001, '****', 
                           ifelse(p_value < 0.001, '***', 
                                  ifelse(p_value < 0.01, '**',
                                         ifelse(p_value < 0.05, '*', 'ns'))))
      return(list(p_value_fc = p_value, fc_annotation = annotation))
    }
    
    results <- COHORT_SUBSET %>% group_by(PLOT_ID, Date) %>%
      do({
        data <- .
        fold_change_data <- data %>%
          summarise(fold_change = round(mean(MEASUREMENT[CONDITION == "STIM"]) / mean(MEASUREMENT[CONDITION == "UNSTIM"]), digits = 2))
        ttest_results <- perform_ttest(data)
        cbind(fold_change_data, ttest_results)
      }) %>%
      ungroup() %>%
      left_join(NAME_KEY[, c("PLOT_ID", "CL_NUMBER", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")], relationship = "many-to-many") %>%
      distinct(PLOT_ID, Date, .keep_all = T)
    
    # Extract annotations for plotting
    annotations        <- results$fc_annotation
    names(annotations) <- results$PLOT_ID
    
    # Reorder
    results$PLOT_ID = reorder(results$PLOT_ID, -results$ORDER_NO)
    
    ## WHY IS THIS AN INFINITE VALUE??? --> INVESTIGATE FURTHER..
    ## THE AMOUNT OF IL-2 RELEASE IS SUPER LOW.. FOLD CHANGE IS NOT AN APPROPRIATE MEASURE HERE.
    results$fold_change[results$PLOT_ID == "CHARMS-BDLD_14_HIGH"] <- NaN
    
    ################################################################################
    
    # Plotting the fold changes
    fc_plot_elisa <- ggplot(results, aes(x = PLOT_ID, y = fold_change)) +
      geom_col(aes(col = PLOTTING_COLOR, fill = PLOTTING_COLOR), position = position_dodge(width = 0.8), alpha = 0.3) +
      geom_text(aes(label = annotations), size = TEXT, hjust = -0.3) +
      geom_hline(yintercept = 2, linetype = "longdash") +
      scale_color_identity() +
      scale_fill_identity() +
      scale_y_continuous(limits = c(0, max(results$fold_change, na.rm = T))) +
      labs(x = "", y = "fold change") +
      coord_flip() +
      theme_cowplot(font_size = SIZE, font_family = FONT) +
      theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
            axis.title.y      = element_blank(),
            legend.position.inside = c(x = 0.4, y = 0.73),
            legend.title      = element_blank(),
            legend.text       = element_text(size = SIZE),
            legend.key.size   = unit(9, "mm")) 
    
    fc_plot_elisa
    
    if (OVERWRITE) {
      
      if (dir.exists(temp_dir))   {print(paste0("Files will be saved to ", temp_dir))}   else {dir.create(temp_dir,   recursive = T); print(paste0("Files will be saved to ", temp_dir))}
      if (dir.exists(output_dir)) {print(paste0("Files will be saved to ", output_dir))} else {dir.create(output_dir, recursive = T); print(paste0("Files will be saved to ", output_dir))}
      
      # save figure
      if (plot_elisa) {
        ggsave(paste0(temp_dir, "/FIGURES/", "ELISA_fold_change.svg"), plot = fc_plot_elisa, device = "svg", width = 12, height = 24); print(paste0("Overwritten: ", temp_dir, "/FIGURES/", "ELISA_fold_change.svg"))
        }
      
      # save tables
      fwrite(results, paste0(temp_dir, "/", "ELISA_SUMMARY_fc.csv"), append = F); print(paste0("Overwritten: ", temp_dir, "/", "ELISA_SUMMARY_fc.csv"))
    }
    
    # rm(COHORT_DATA, NORMALIZED_TO_CONTROL, unstim_means, COHORT_SUBSET, results, annotations, fc_plot_elisa)
    
  }
  
  if (join_elisa_df) {
    
    elisa_fc_df     <- fread(paste0(temp_dir, "/", "ELISA_SUMMARY_fc.csv"))
    elisa_3E10_df   <- fread(paste0(temp_dir, "/", "ELISA_SUMMARY_3E10_normalized.csv")) %>% mutate(POSITIVE_CTRL = "3E10")
    elisa_cl069_df  <- fread(paste0(temp_dir, "/", "ELISA_SUMMARY_cl069_normalized.csv"))
    elisa_rltv_real <- rbind(elisa_3E10_df, elisa_cl069_df)
    elisa_main      <- elisa_rltv_real %>% left_join(elisa_fc_df) %>% left_join(NAME_KEY, relationship = "many-to-many") %>% filter(PURPOSE == "SAMPLE")
    elisa_main_cln  <- elisa_main %>% arrange(ID, POSITIVE_CTRL) %>% group_by(ID, POSITIVE_CTRL) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>%  distinct(ID, POSITIVE_CTRL, .keep_all = T)
    rm(elisa_fc_df, elisa_3E10_df, elisa_cl069_df, elisa_rltv_real, elisa_main)
    
    if (OVERWRITE) {
      fwrite(elisa_main_cln, file.path(output_dir, "ELISA_SUMMARY.csv"), append = F); print(paste0("Overwritten: ", output_dir, "/", "ELISA_SUMMARY.csv"))
    }
  }
}

### TODO: finish AF3 extraction script
if (GENERATE_AF3_JSON_TABLE) {
  
  mega_af3_df <- extract_af3_json_files(af_json_files = af_json_files)[[1]]
  mini_af3_df <- extract_af3_json_files(af_json_files = af_json_files)[[2]]
  AF3_SUMMARY <- extract_af3_json_files(af_json_files = af_json_files)[[3]]
  
  if (OVERWRITE) {
    
    fwrite(mega_af3_df, file.path(temp_dir,   "AF3_SUMMARY_ALL_MODELS.csv"),  append = F, row.names = F, col.names = T)
    fwrite(mini_af3_df, file.path(temp_dir,   "AF3_SUMMARY_BEST_MODELS.csv"), append = F, row.names = F, col.names = T)
    fwrite(AF3_SUMMARY, file.path(output_dir, "AF3_SUMMARY.csv"),             append = F, row.names = F, col.names = T)
    
    print("Overwriting AF3 files: ")
    print(file.path(temp_dir,     "AF3_SUMMARY_ALL_MODELS.csv"))
    print(file.path(temp_dir,     "AF3_SUMMARY_BEST_MODELS.csv"))
    print(file.path(output_dir,   "AF3_SUMMARY.csv"))
  
    AF3_SUMMARY <- fread(file.path(output_dir, "AF3_SUMMARY.csv"))
    
  } else {
    print("Not overwriting AF3 files.")
    AF3_SUMMARY <- fread(file.path(output_dir, "AF3_SUMMARY.csv"))
  }
  
  # declutter environment and remove all objects that are no longer needed
  rm(large_af3_dfs, large_af3_df, json, mega_af3_df, mini_af3_df, sum_conf_extract_list,
     af3_name, af3_seed, af3_seq, af3_n_chains, summary_confidences_files, af3_model,
     af3_chain_iptm, af3_chain_pair_iptm, af3_chain_pair_pae_min, af3_chain_ptm,
     af3_fraction_disordered, af3_has_clash, af3_iptm, af3_n_recycles, af3_ptm,
     af3_ranking_score, sum_conf_extract)
}

### TODO: Add the new tables to SOTA
if (GENERATE_SOTA_TABLE) {
  SOTA_SUMMARY  <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA.csv$",   ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA)) %>% unique()
  # SOTA_final[duplicated(SOTA_final$seq_DLD_aa)]$ID
  ELISA_SUMMARY <- fread(file.path(output_dir, "ELISA_SUMMARY.csv")) ; elisa_cols     <- colnames(ELISA_SUMMARY)
  AF3_SUMMARY   <- fread(file.path(output_dir, "AF3_SUMMARY.csv"))   ; alphaFold_cols <- colnames(AF3_SUMMARY)
  
  sota <- 
    left_join(SOTA_SUMMARY, ELISA_SUMMARY, relationship = "many-to-many") %>% 
    left_join(AF3_SUMMARY, relationship = "many-to-many") %>% 
    unique()
  
  if (OVERWRITE) {fwrite(sota, file = file.path(output_dir, "SOTA.csv"), append = F, row.names = F, col.names = T)}
  sota <- fread(file.path(output_dir, "SOTA.csv"))
  # rm(SOTA_SUMMARY, ELISA_SUMMARY, AF3_SUMMARY)
  
}

### TODO: perform correlation analysis
if (GENERATE_MSA_TABLE) {
  
  library(msa)

  AF3_SUMMARY <- fread(file.path(output_dir, "AF3_SUMMARY.csv"))
  sota <- fread(file.path(output_dir, "SOTA.csv")) %>% select(ID, species) %>% unique()
  
  # write a FASTA file of all ordered sequences
  df_for_msa <- AF3_SUMMARY %>% distinct(ID, af3_seq, .keep_all = F) %>% filter(!is.na(af3_seq))
  write_fasta(df_for_msa, "ID", "af3_seq", file.path(temp_dir, "SEQS_AND_ALNS/DLD_sequences.fasta")); rm(df_for_msa)
  msa_sota <- readAAStringSet(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences.fasta")) ; head(msa_sota)
  
  methods <- c("ClustalW", "ClustalOmega", "Muscle", "MAFFT", "kalign")
  for (method in methods) {
    output_file <- paste0(temp_dir, "/SEQS_AND_ALNS/DLD_sequences_aligned_", tolower(method), ".fasta")
    if (method %in% c("ClustalW", "ClustalOmega", "Muscle")) {
      aln_msa_sota <- msa(msa_sota, method = method)
      print(aln_msa_sota, show = "complete")
      if (OVERWRITE) {
        writeXStringSet(as(aln_msa_sota, "AAStringSet"), filepath = output_file)
      }
    } else if (method == "MAFFT") {
      system(paste0("mafft --auto --inputorder --thread 4 --reorder ", file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences.fasta"), " > ", output_file))
    } else if (method == "kalign") {
      # run in command line - kalign is a bit meh running directly from R
      # kalign -i SOTA/02_temp_data/SEQS_AND_ALNS/DLD_sequences.fasta -o SOTA/02_temp_data/SEQS_AND_ALNS/DLD_sequences_aligned_kalign.fasta
    }
  }
  
  ## create a function to process all alignments in a streamlined way
  run_iterative_msa_masking = FALSE
  if (run_iterative_msa_masking) {
    for (file_path in list.files(file.path(temp_dir, "SEQS_AND_ALNS"), pattern = "^DLD_sequences_aligned_.*\\.fa(sta)?$", full.names = TRUE)) {
      iterative_msa_masking(file_path, temp_dir, sota,
                            min_fraction_start_end = seq(0.1, 0.6, by = 0.1),
                            min_block_width_start_end = c(1, 2, 3))
    }
  }
  
  # Initialize an empty data frame to store all consensus sequences
  consensus_sequences <- data.frame()
  
  # Loop through each consensus sequence file
  for (cons_seq_file in list.files(file.path(temp_dir, "SEQS_AND_ALNS/iterative_msa_masking/"), pattern = paste0("(", paste(tolower(methods), collapse = "|"), ")*_cons_string\\.csv$"), full.names = TRUE)) {
    print(cons_seq_file)
    consensus_sequences <- rbind(consensus_sequences, fread(cons_seq_file))
  }
  
  # Initialize counts vectors
  num_sequences  <- nrow(consensus_sequences)
  count_hash     <- integer(num_sequences)
  count_dash     <- integer(num_sequences)
  count_question <- integer(num_sequences)
  count_letters  <- integer(num_sequences)
  consensus_seqs <- character(num_sequences)
  
  # Iterate through each consensus sequence
  for (i in seq_len(num_sequences)) {
    seq <- as.character(consensus_sequences[i, ])  # Ensure seq is treated as character vector
    
    # Count characters
    count_hash[i]     <- sum(str_count(seq, "#"))
    count_dash[i]     <- sum(str_count(seq, "-"))
    count_question[i] <- sum(str_count(seq, "\\?"))
    count_letters[i]  <- sum(str_count(seq, "[A-Za-z]"))
    consensus_seqs[i] <- seq
  }
  
  # Create a data frame to store the counts
  consensus_counts <- data.frame(
    Count_Hash     = count_hash,
    Count_Dash     = count_dash,
    Count_Question = count_question,
    Count_Letters  = count_letters,
    Consensus_seq  = consensus_seqs
  )
  
  # Optionally, set row names to the sequences
  rownames(consensus_counts) <- basename(list.files(file.path(temp_dir, "SEQS_AND_ALNS/iterative_msa_masking/"), pattern = "cons_string.csv$", full.names = TRUE))
  rownames(consensus_counts) <- str_replace_all(string = rownames(consensus_counts), "DLD_sequences_aligned_", "")
  rownames(consensus_counts) <- str_replace_all(string = rownames(consensus_counts), "_cons_string.csv", "")
  consensus_counts$min_block_width <- as.numeric(str_extract(rownames(consensus_counts), "(?<=_min_block_width_)\\d+"))
  consensus_counts$min_fraction    <- as.numeric(str_extract(rownames(consensus_counts), "(?<=_min_fraction_)\\d+\\.\\d+"))
  consensus_counts$method          <- str_extract(rownames(consensus_counts), "^[^_]+")
   
  # Print the resulting data frame
  print(consensus_counts)
  
  # Function to extract amino acids only from consensus sequences
  extract_amino_acids <- function(sequence) {gsub("[^A-Z]", "", sequence)}

  # Apply function to extract amino acids
  amino_acids_only <- sapply(consensus_counts$Consensus_seq, extract_amino_acids)

  # Create a matrix representation of sequences
  max_length <- max(nchar(amino_acids_only))
  alignment_matrix <- matrix("", nrow = length(amino_acids_only), ncol = max_length)
  for (i in seq_along(amino_acids_only)) {
    alignment_matrix[i, 1:nchar(amino_acids_only[i])] <- strsplit(amino_acids_only[i], "")[[1]]
  }
  
  # Print matrix representation
  print(alignment_matrix)
  
  alignment <- msa(consensus_counts$Consensus_seq, type = "protein")
  
  # Access the aligned sequences and their names
  aligned_sequences <- as.character(alignment)
  aligned_names <- rownames(consensus_counts)
  
  # Print or use the aligned sequences and their names as needed
  print(aligned_sequences)
  print(aligned_names)
  
  # Save the alignment with names
  if (OVERWRITE) {
    # Prepare aligned sequences with names
    aligned_seqs_with_names <- as(alignment, "AAStringSet")
    names(aligned_seqs_with_names) <- aligned_names
    
    # Write aligned sequences to a FASTA file
    writeXStringSet(aligned_seqs_with_names, file = file.path(temp_dir, "SEQS_AND_ALNS", "DLD_aligned_consensus_across_methods.fasta"))
  }
  
  #### Generate the consensus of the conensus
  origMAlign <- Biostrings::readAAMultipleAlignment(file.path(temp_dir, "SEQS_AND_ALNS", "DLD_aligned_consensus_across_methods.fasta"))
  
  ## list the names of the sequences in the alignment
  rownames(origMAlign)
  
  # ## automatically masks (removes or marks as unreliable) columns in the alignment that have a high proportion of gaps.
  # ## Improves the reliability of subsequent analyses by excluding poorly aligned regions.
  autoMasked <- maskGaps(origMAlign, min.fraction=0.5, min.block.width=2)
  colmask(autoMasked)
  
  ## Calculates the frequency of each amino acid in the masked alignment.
  ## Provides insights into the composition of the sequences
  alph_freq <- Biostrings::alphabetFrequency(autoMasked)
  
  ## Provides insights into the consensus sequence across the alignment.
  cons_mat  <- Biostrings::consensusMatrix(autoMasked, baseOnly=TRUE)[, 1:nrow(autoMasked)]
  
  ## Generates a consensus sequence from the alignment
  ## "#" represents a masked column
  ## "-" represents a gap
  ## "?" represents an ambiguous position
  cons_string <- consensusString(autoMasked)
  cons_view   <- consensusViews(autoMasked)
  sdist       <- stringDist(as(autoMasked,"AAStringSet"), method="hamming")
  clust       <- hclust(sdist, method = "single")
  
  if (OVERWRITE) {
    fwrite(as.data.frame(cons_string), file.path(temp_dir, "SEQS_AND_ALNS", "DLD_aligned_consensus_across_methods_cons_string.csv"))
    print(paste0("Overwritten: ", temp_dir, "/", "DLD_aligned_consensus_across_methods_cons_string.csv"))
  }
  
  ## Add the MSA to SOTA table
  sota <- fread(file.path(output_dir, "SOTA.csv"))
  msa_cons_string <- fread(file.path(temp_dir, "SEQS_AND_ALNS", "DLD_aligned_consensus_across_methods_cons_string.csv"))
  msa_sota <- readAAStringSet(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences.fasta")) ; head(msa_sota)
  
  kalign   <- readAAMultipleAlignment(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences_aligned_kalign.fasta")) ; kalign
  # mafft    <- readAAMultipleAlignment(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences_aligned_mafft.fasta"))  ; mafft
  # muscle   <- readAAMultipleAlignment(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences_aligned_muscle.fasta")) ; muscle
  # clustalw <- readAAMultipleAlignment(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences_aligned_clustalw.fasta")) ; clustalw
  # clustalo <- readAAMultipleAlignment(file.path(temp_dir, "/SEQS_AND_ALNS/DLD_sequences_aligned_clustalomega.fasta")) ; clustalo
  
  # Split the MSA into a data frame with one column per position
  msa_split_df <- data.frame(ID = rownames(kalign), seq = kalign@unmasked) %>% 
    separate(seq, into = paste0("pos_", sprintf("%02d", 1:max(nchar(kalign@unmasked)))), sep = "")
  
  # write msa_split_df to a CSV file
  if (OVERWRITE) {
    fwrite(msa_split_df, file.path(output_dir, "MSA_SUMMARY.csv"), col.names = T)
    print(paste0("Overwritten: ", output_dir, "/", "MSA_SUMMARY.csv"))
  }

  # One-liner to apply the consensus function and convert the result to a single string
  # curly brackets allow you to treat the result of the pipe as an argument in the paste function.
  # The unlist function is used to convert the result to a vector before pasting it together.
  consensus_strings <- data.frame() ; max_pct_stored    <- data.frame()

  for (max_pct in seq(0.5, 1.0, by = 0.1)) {
    print(consensus_string <- msa_split_df %>% 
      summarise(across(starts_with("pos_"), ~consensus_char(.x, min_pct = 0.2, max_pct = max_pct))) %>% 
      {paste(unlist(.), collapse = "")})
    print(max_pct_current <- max_pct)
    
    consensus_strings <- rbind(consensus_strings, consensus_string)
    max_pct_stored <- rbind(max_pct_stored, max_pct_current)
  }
  
  # Save the consensus strings to a CSV file
  if (OVERWRITE) {
    fwrite(cbind(max_pct_stored, consensus_strings), file.path(temp_dir, "DLD_consensus_seq_iterative.csv"), col.names = F)
    print(paste0("Overwritten: ", temp_dir, "/", "DLD_consensus_seq_iterative.csv"))
  }
  
  # declutter environment and remove all objects that are no longer needed
  rm(consensus_strings, max_pct_stored, msa_split_df, msa_cons_string, sota, msa_sota, kalign)
}

if (GENERATE_LCI_TABLE) {
  
  GENERATE_CELL_SUMMARY = FALSE
  
  if (GENERATE_CELL_SUMMARY) {
    
    source("scripts/SOTA_LCI.R")
    
    lci_table <- lci_data %>% 
      select(PROTEIN,COHORT,IMAGE,CELL,UNIVERSAL_TRACK_ID,TIME,TIME_ADJUSTED,FRAME,FRAMES_ADJUSTED,
             FRAMES_SINCE_LANDING,LIFETIME,NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,
             STARTING_NORMALIZED_INTENSITY,COMPLEMENTARY_PROTEIN_1,COMPLEMENTARY_NORMALIZED_INTENSITY_1) %>%
      arrange(UNIVERSAL_TRACK_ID) %>%
      as.data.table()
    
    lci_table$COHORT <- as.factor(lci_table$COHORT) ; rm(lci_data)
    
    LOW_CAT     <- "< 4 s"
    MEDIUM_CAT  <- "4-40 s"
    HIGH_CAT    <- "≥ 40 s"
    
    # Clean up table, and Measure dwell frames/time of CHARMS associated TRAF6 ----
    cell_table <- lci_table %>% 
      # Step 1: Filter out the noise
      filter(PROTEIN == "MyD88", MAX_NORMALIZED_INTENSITY >= 1) %>% 
      # Step 2: Group the data by 'UNIVERSAL_TRACK_ID'
      group_by(UNIVERSAL_TRACK_ID) %>% 
      # Step 3: Filter tracks with at least two time points based on 'FRAMES_ADJUSTED'
      filter(max(FRAMES_ADJUSTED) >= 2) %>%
      # Step 4: Calculate colocalization and create a streak identifier ('TRAF6_Streak')
      mutate(COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1) %>%
      mutate(TRAF6_Streak = cumsum(!COLOCALIZATION)) %>%
      # Step 5: Filter rows with colocalization events
      filter(COLOCALIZATION == 1) %>%
      # Step 6: Group by multiple variables for further summarization
      group_by(COHORT, IMAGE, CELL, UNIVERSAL_TRACK_ID, TRAF6_Streak) %>%
      # Step 7: Summarize the colocalization events
      summarise(DWELL_FRAMES = sum(COLOCALIZATION), # Number of frames that complementary protein is above threshold in a continuous stretch
                DWELL_TIME = (sum(COLOCALIZATION) - 1) * 4) %>%
      # Step 8: Categorize dwell times
      mutate(CATEGORY_DWELL_TIME = fcase(DWELL_TIME ==  0, LOW_CAT,
                                         DWELL_TIME <  40 & DWELL_TIME != 0, MEDIUM_CAT,
                                         DWELL_TIME >= 40, HIGH_CAT)) %>%
      # Step 9: Convert the result to a data.table
      as.data.table()
    
    if (OVERWRITE) {
      fwrite(cell_table, file.path(temp_dir, "/lci_cell_table.csv"), row.names = F, quote = F, sep = ",")
      system(paste0("gzip -f ",    temp_dir, "/lci_cell_table.csv"))
      
      cell_table <- fread(file.path(temp_dir, "/lci_cell_table.csv.gz"))
    }
    
    # Measure Pct of TRAF6 CATEGORY_DWELL_TIME ----
    mean_cell <-
      cell_table %>% 
      group_by(COHORT, IMAGE, CELL, CATEGORY_DWELL_TIME) %>% 
      # Values to count
      count(CATEGORY_DWELL_TIME, COHORT, IMAGE, CELL, 
            # Name of the Column with the counts
            name = "N_CATEGORY_DWELL_TIME", 
            # .drop =F ensures that output of 0 will be kept
            .drop = FALSE ) %>% 
      group_by(COHORT, IMAGE, CELL) %>% 
      # more than 3 recruitment events per cell
      filter(sum(N_CATEGORY_DWELL_TIME) >= 3) %>%
      mutate(PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME)*100,
             #DATE = strsplit(IMAGE, " ")[[1]][1]
      ) %>%
      arrange(CATEGORY_DWELL_TIME) %>%
      as.data.table()
    
    mean_replicates <- mean_cell %>%
      group_by(COHORT, IMAGE, CATEGORY_DWELL_TIME) %>%
      summarise(PCT_RECRUITMENT = mean(PCT_RECRUITMENT)) %>%
      arrange(CATEGORY_DWELL_TIME) %>%
      as.data.table()
    
    mean_total <- mean_replicates %>% 
      group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
      summarise(SE_PCT_RECRUITMENT = std.error(PCT_RECRUITMENT),
                PCT_RECRUITMENT = mean(PCT_RECRUITMENT)) %>% 
      arrange(CATEGORY_DWELL_TIME) %>%
      as.data.table()
    
    if (OVERWRITE) {
      write.csv(mean_cell,       file.path(temp_dir, "/lci_mean_cell.csv"),       row.names = F)
      write.csv(mean_replicates, file.path(temp_dir, "/lci_mean_replicates.csv"), row.names = F)
      write.csv(mean_total,      file.path(temp_dir, "/lci_mean_total.csv"),      row.names = F)
      
      if (file.exists(file.path(temp_dir, "/lci_mean_cell.csv.gz"))) {file.remove(file.path(temp_dir, "/lci_mean_cell.csv.gz"))}
      if (file.exists(file.path(temp_dir, "/lci_mean_replicates.csv.gz"))) {file.remove(file.path(temp_dir, "/lci_mean_replicates.csv.gz"))}
      if (file.exists(file.path(temp_dir, "/lci_mean_total.csv.gz"))) {file.remove(file.path(temp_dir, "/lci_mean_total.csv.gz"))}
      
      system(paste0("gzip -f ", temp_dir, "/lci_mean_cell.csv"))
      system(paste0("gzip -f ", temp_dir, "/lci_mean_replicates.csv"))
      system(paste0("gzip -f ", temp_dir, "/lci_mean_total.csv"))
      
    }
    
    mean_repl_long   <- mean_replicates %>% filter(CATEGORY_DWELL_TIME == HIGH_CAT)   %>% group_by(COHORT) %>% mutate(num_replicates = n())
    mean_repl_medium <- mean_replicates %>% filter(CATEGORY_DWELL_TIME == MEDIUM_CAT) %>% group_by(COHORT) %>% mutate(num_replicates = n())
    mean_repl_short  <- mean_replicates %>% filter(CATEGORY_DWELL_TIME == LOW_CAT)    %>% group_by(COHORT) %>% mutate(num_replicates = n())
    
    # Perform statistics using t.test (only possible for cell lines with 2 or more replicates)
    stat.test_pct_long_lived_TRAF6   <- mean_repl_long   %>% filter(num_replicates >= 2) %>% ggpubr::compare_means(data = . , PCT_RECRUITMENT ~ COHORT, method = "t.test")
    stat.test_pct_medium_lived_TRAF6 <- mean_repl_medium %>% filter(num_replicates >= 2) %>% ggpubr::compare_means(data = . , PCT_RECRUITMENT ~ COHORT, method = "t.test")
    stat.test_pct_short_lived_TRAF6  <- mean_repl_short  %>% filter(num_replicates >= 2) %>% ggpubr::compare_means(data = . , PCT_RECRUITMENT ~ COHORT, method = "t.test")
    
    if (OVERWRITE) {
      write.csv(stat.test_pct_long_lived_TRAF6,   file.path(temp_dir, "lci_stat.test_pct_long_lived_TRAF6.csv"),   row.names = F)
      write.csv(stat.test_pct_medium_lived_TRAF6, file.path(temp_dir, "lci_stat.test_pct_medium_lived_TRAF6.csv"), row.names = F)
      write.csv(stat.test_pct_short_lived_TRAF6,  file.path(temp_dir, "lci_stat.test_pct_short_lived_TRAF6.csv"),  row.names = F) 
    }
  }
  
  cell_table <- fread(file.path(temp_dir, "lci_cell_table.csv.gz"))
  mean_cell <- fread(file.path(temp_dir, "/lci_mean_cell.csv.gz"))
  mean_replicates <- fread(file.path(temp_dir, "/lci_mean_replicates.csv.gz"))
  mean_total <- fread(file.path(temp_dir, "/lci_mean_total.csv.gz"))
  
}

if (GENERATE_OPERON_TABLE) {

  sota      <- fread(file.path(output_dir, "SOTA.csv"), na.strings = c("", NA, "NA"), strip.white = TRUE)
  operon_df <- fread(file.path(output_dir, "SOTA.csv"), na.strings = c("", NA, "NA"), strip.white = TRUE) %>%
    # create a subset df that only contains the columns we need
    select(ID, genbank_acc, uniprot_acc, operon, architecture, genbank_gca, genbank_gca_2, refseq_gcf, refseq_gcf_2) %>%
    # remove rows with missing operon information
    filter(!is.na(operon)) %>%
    # remove duplicates
    distinct(ID, .keep_all = TRUE)
  
  # Apply the function to each operon in operon_df and bind rows
  processed_operons   <- mapply(interpret_operon_direction, ID = operon_df$ID, operon = operon_df$operon, SIMPLIFY = FALSE) %>% bind_rows()
  
  # Merge the processed operons with the original data frame
  operon_df_processed <- operon_df %>% left_join(processed_operons) %>% filter(!is.na(part), part != c("")) %>% 
    # trim white space from the 'part' column
    mutate(part = str_trim(part)) %>%
    group_by(ID) %>%
    mutate(n_parts = n(),
           # grab the rows per ID that have the part that contains the pattern 'DLD' and count them
           n_DLD = sum(str_detect(part, "DLD")),
           # Identify DLD positions using case_when
           DLD_position = case_when(
             # Beginning (bDLD+X+Y+..)
             str_detect(part, "bDLD") & !str_detect(part, "\\+.*bDLD") ~ "NTD",
             # End (X+Y+Z+bDLD)
             str_detect(part, "bDLD") & !str_detect(part, "bDLD.*\\+") ~ "CTD",  
             # Middle (X+bDLD+Z)
             str_detect(part, "bDLD") ~ "MID",
             # Single (bDLD)
             str_detect(part, "bDLD") ~ "SINGLE"),
           part_pos_count     = row_number(),
           DLD_position_count = ifelse(str_detect(part, "bDLD"), row_number(), NA_integer_),
           DLD_pos_cnt_ratio  = ifelse(!is.na(DLD_position_count), paste0(DLD_position_count, "/", n_parts), NA_integer_),
           n_upstream         = DLD_position_count - 1,
           n_downstream       = n_parts - DLD_position_count) %>% 
    ungroup()
  
  genomes_to_download <- sota %>% select(genbank_gca, genbank_gca_2) %>% 
    pivot_longer(cols = c(genbank_gca, genbank_gca_2), names_to = "assembly_type", values_to = "assembly") %>% 
    filter(!is.na(assembly)) %>% 
    distinct(assembly)
  
  
  if (OVERWRITE) {
    fwrite(operon_df,           file.path(temp_dir,  "operon_df.csv"))
    fwrite(operon_df_processed, file.path(temp_dir,  "operon_df_processed.csv"))
    fwrite(genomes_to_download, file.path(input_dir, "genome_download_list.csv"), append = F, row.names = F, col.names = F)
    
    # declutter environment and remove all objects that are no longer needed
    rm(operon_df, operon_df_processed, processed_operons, genomes_to_download)
  }
}
  
if (DOWNLOAD_AND_PROCESS_GENOME_FOLDERS) {
    ### Download genomes (script #1)
    # run this directly in the command line - connection cannot be established via R CLI
    # go to coding_universe > scripts > run SOTA_GENOME_DOWNLOAD.sh
    # ************************************************************************ #
    
    ### Flatten the folder structure (script #2) and add the genome prefix to all folder files (script #3)
    system("bash ./scripts/SOTA_GENOME_FLATTEN.sh") ; system("bash ./scripts/SOTA_GENOME_RENAME.sh") ; system("sleep 3")
    
    # check if the genomes have been downloaded and processed by comparing files in file.path(input_dir, "GENOMES") to genomes_to_download
    # genome_files <- list.files(file.path(input_dir, "GENOMES"), pattern = "_extracted$", full.names = FALSE) %>% stringr::str_remove("_extracted$") %>% na.omit()
    # genomes_to_download <- genomes_to_download %>% pull(assembly)
    # setdiff(genomes_to_download, genome_files)
  }
  
  # # remotes::install_github("thackl/gggenomes")
  # # library(gggenomes)
  # 
  # operon_df_processed <- fread(file.path(temp_dir, "operon_df_processed.csv"))
  # sota                <- fread(file.path(output_dir, "SOTA.csv"), na.strings = c("", NA, "NA"), strip.white = TRUE)
  # genome_files        <- list.files(file.path(input_dir, "GENOMES"), pattern = "_extracted$", full.names = FALSE) %>% stringr::str_remove("_extracted$") %>% na.omit()
  # 
  # # genome_file <- "GCA_002154725.1" # NOSTOC SP. 106C
  # # genome_file <- "GCF_963875125.1" # CALOTHRIX NIES-2098 NEW
  # # genome_file <- genome_files[3]
  # 
  # genome_list <- list.files(file.path(input_dir, "GENOMES"), pattern = genome_file, full.names = TRUE, recursive = TRUE)
  # genome      <- readDNAStringSet(genome_list[grepl(".fna", genome_list)])
  # genome_gff  <- read.gff(genome_list[grepl(".gff", genome_list)])
  # # genome_gb   <- read_gbk(genome_list[grepl(".gb", genome_list)])
  # 
  # # genome_gb_CDS <- genome_gb %>% filter(type == "CDS")
  # 
  # gb_names <- search_target_by_pattern_col(sota, genome_file, "genbank_gca", "refseq_gcf", "genbank_acc") %>% unique()
  
  
  if (EXTRACT_OPERONS) {
    
    operon_df_processed <- fread(file.path(temp_dir, "operon_df_processed.csv"))
    
    # list genomes (fna files)
    genome_files <- list.files(file.path(input_dir, "GENOMES"), pattern = "_extracted$", full.names = FALSE) %>% stringr::str_remove("_extracted$") %>% na.omit()
    
    # genome_file <- genome_files[3] # for debugging
    # genome_file <- "GCA_000522425.1"
    # genome_file <- "GCA_002154725.1"
    # genome_file <- "GCA_001610855.1"
    for (genome_file in genome_files) {
      if (genome_file == "") next  # Skip if the genome_file is empty
      tryCatch({
        
        ### (1) Load raw files (genome fasta, gff, maybe gb) ###################
        genome_list <- list.files(file.path(input_dir, "GENOMES"), pattern = genome_file, full.names = TRUE, recursive = TRUE)
        genome      <- readDNAStringSet(genome_list[grepl(".fna", genome_list)])
        genome_gff  <- read.gff(genome_list[grepl(".gff", genome_list)])
        
        ####### Fetch relevant gene names and genbank accessions
        gene_names <- search_target_by_pattern_col(sota, genome_file, "assembly", "assembly_2", "gene")
        gb_names   <- search_target_by_pattern_col(sota, genome_file, "assembly", "assembly_2", "genbank_acc")
        
        ####### Initialize a list to collect data
        
        data_list <- list(
          genbank_acc = character(),
          n_DLD = numeric(),
          n_parts = numeric(),
          gene_base_start = numeric(),
          gene_base_end = numeric(),
          operon_base_start = numeric(),
          operon_base_end = numeric(),
          operon_components = character()
        )
        
        ### (2) Process each genbank accession #################################
        
        for (i in seq_along(gb_names)) {
          cat("************************************************************* \n")
          gb_name <- gb_names[i]
          cat("*** Starting ", gb_name, " *** \n")
          
          # Find the gene or gb accession in the GFF file
          if (!is.na(gb_name) && nchar(gb_name) > 0) {
            n_parts <- operon_df_processed$n_parts[operon_df_processed$genbank_acc == gb_name] %>% unique()
            bDD_parts <- operon_df_processed$DLD_position_count[operon_df_processed$genbank_acc == gb_name]
            bDD_parts <- bDD_parts[!is.na(bDD_parts)]
            
            for (j in seq_along(bDD_parts)) {
              bDD_part <- bDD_parts[j]
              cat("*** bDD position: ", bDD_part, "/", n_parts, " *** \n")
              
              get_n_upstream   <- as.numeric(bDD_part - 1)
              get_n_downstream <- as.numeric(n_parts   - bDD_part)
              
              # Get gene start and end positions
              gene_start <- genome_gff$start[grepl(gb_name, genome_gff$attributes)] %>% unique()
              gene_end <- genome_gff$end[grepl(gb_name, genome_gff$attributes)] %>% unique()
              
              if (length(gene_start) == 0 || length(gene_end) == 0) {
                cat("Gene start or end not found for ", gb_name, "\n")
                next
              }
              
              if (get_n_downstream == 0 & get_n_upstream == 0) {
                upstream_genes    <- NA
                downstream_genes  <- NA
                operon_start      <- gene_start
                operon_end        <- gene_end
                operon_components <- gb_name
                
              } else if (get_n_downstream == 0 & get_n_upstream > 0) {
                upstream_genes    <- get_upstream_genes(gb_name, get_n_upstream)
                downstream_genes  <- NA
                operon_start      <- genome_gff$start[grepl(upstream_genes, genome_gff$attributes)] %>% unique()
                operon_end        <- gene_end
                operon_components <- c(upstream_genes, gb_name)
                
              } else if (get_n_downstream > 0 & get_n_upstream == 0) {
                upstream_genes    <- NA
                downstream_genes  <- get_downstream_genes(gb_name, get_n_downstream)
                operon_start      <- gene_start
                operon_end        <- genome_gff$end[grepl(downstream_genes[length(downstream_genes)], genome_gff$attributes)] %>% unique()
                operon_components <- c(gb_name, downstream_genes)
                
              } else {
                upstream_genes    <- get_upstream_genes(gb_name, get_n_upstream)
                downstream_genes  <- get_downstream_genes(gb_name, get_n_downstream)
                operon_start      <- genome_gff$start[grepl(upstream_genes[1], genome_gff$attributes)] %>% unique()
                operon_end        <- genome_gff$end[grepl(downstream_genes[length(downstream_genes)], genome_gff$attributes)] %>% unique()
                operon_components <- c(upstream_genes, gb_name, downstream_genes)
              }
              
              # Ensure consistency in length
              data_list$genbank_acc       <- c(data_list$gb_name, gb_name)
              data_list$n_DLD             <- c(data_list$bDD_part, bDD_part)
              data_list$n_parts           <- c(data_list$n_parts, n_parts)
              data_list$gene_base_start   <- c(data_list$gene_base_start, ifelse(length(gene_start) > 0, gene_start, NA))
              data_list$gene_base_end     <- c(data_list$gene_base_end, ifelse(length(gene_end) > 0, gene_end, NA))
              data_list$operon_base_start <- c(data_list$operon_base_start, ifelse(length(operon_start) > 0, operon_start, NA))
              data_list$operon_base_end   <- c(data_list$operon_base_end, ifelse(length(operon_end) > 0, operon_end, NA))
              data_list$operon_components <- c(data_list$operon_components, paste(operon_components, collapse=", "))
            }
          } else {
            warning(paste("No genbank accession found for", genome_file, "in SOTA or matched in the genome_gff table. Skipping."))
            next
          }
        }
        
        # Combine lists into a data frame
        operon_df <- data.frame(
          genbank_acc       = data_list$genbank_acc,
          n_DLD             = data_list$n_DLD,
          n_parts           = data_list$n_parts,
          gene_base_start   = data_list$gene_base_start,
          gene_base_end     = data_list$gene_base_end,
          operon_base_start = data_list$operon_base_start,
          operon_base_end   = data_list$operon_base_end,
          operon_components = data_list$operon_components,
          stringsAsFactors  = FALSE
        )
        
        operon_df_processed_final <- left_join(operon_df_processed, operon_df)
        
        opf <- operon_df_processed_final %>%
          group_by(ID) %>%
          select(ID, part_pos_count, n_parts, operon_components) %>%
          separate_rows(operon_components, sep = ",\\s*") %>%
          group_by(ID) %>%
          select(ID, operon_components) %>%
          distinct() %>%
          mutate(part_pos_count = row_number()) %>%
          filter(!is.na(operon_components)) %>%
          ungroup()
        
        # Step 3: Merge the processed operons with the original data frame
        operon_df_processed_final <- operon_df_processed_final %>% select(-operon_components) %>% left_join(opf)
        
      }, error = function(e) {
        message("Error processing genome file: ", genome_file, "\n", e)
      })
    }
    
    if (OVERWRITE) {
      fwrite(operon_df_processed, file.path(temp_dir, "operon_df_processed.csv"))
    }
  }
    
    
    
    
   
    
    
    # extract bDD-containing gene names from sota
    # bDD_genes <- sota$genbank_acc %>% str_extract(".*(?=\\.)") %>% unique()
    
    
  # TODO: extract operon information from the genomes
  # TODO: Load genome data and annotation data
  # TODO: Extract operons based on up/downstream genes
  # TODO: Extract operons based on gene clustering?
  # TODO: Extract operons based on distance (extract +/- 10kb from the bDD-containing gene)
  # TODO: Set up a line plot for comparing the different approaches
  
  




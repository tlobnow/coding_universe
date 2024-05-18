### Generating a State Of The Art (SOTA) Table on BDLD summary

LOAD_FX    <- T
LOAD_COLS  <- T
LOAD_PATHS <- T

OPEN_DATA  <- T

RM_CLEANED_COLS <- T

GENERATE_GENOME_TABLE <- T
OVERWRITE <- T
CLEAN     <- T

##############################################################
##############################################################
##############################################################

# PACKAGES
library(pacman) ; pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, 
                                 dplyr, cowplot, readxl, scales, knitr, tidyr, visdat, purrr)

# FX
if (LOAD_FX) {
  source("~/Documents/Github/Analysis_Synthethic_Myddosome/FigureScripts/functions.R")
  
  read_excel_allsheets <- function(filename, tibble = FALSE) {
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
  }
  
  fill_column_by_id <- function(data, id, target_col) {
    # Fill the target column based on id groups
    data <- data %>%
      group_by(!!sym(id)) %>%
      fill(!!sym(target_col), .direction = "downup") %>%
      ungroup()
    
    return(data)
  }
  
  fill_subset <- function(data, id, cols_to_join, new_col_name) {
    temp_data <- data %>%
      select(!!sym(id), all_of(cols_to_join)) %>%
      pivot_longer(names_to = "Temp", values_to = new_col_name, cols = all_of(cols_to_join)) %>%
      arrange(!!sym(id)) %>%
      group_by(!!sym(id)) %>%
      fill(everything(), .direction = "downup") %>%
      ungroup() %>%
      select(!!sym(id), !!sym(new_col_name)) %>%
      distinct(!!sym(id), .keep_all = TRUE)
    
    return(temp_data)
  }
  
  join_subset_to_main <- function(data, temp_data, id, cols_to_join) {
    data <- data %>%
      full_join(temp_data, by = id) %>%
      select(-all_of(cols_to_join)) %>%
      arrange(!!sym(id)) %>%
      group_by(!!sym(id)) %>%
      fill(everything(), .direction = "downup") %>%
      ungroup() %>%
      distinct(!!sym(id), .keep_all = TRUE)
    
    return(data)
  }
  
  join_columns <- function(data, id, cols_to_join, new_col_name) {
    temp_data <- fill_subset(data, id, cols_to_join, new_col_name)
    data <- join_subset_to_main(data, temp_data, id, cols_to_join)
    
    return(data)
  }
  
  
}

##############################################################
##############################################################
##############################################################

# SEE COLUMN DOCUMENTATION --> ALL THE WAY AT THE END OF THE FILE!

#### Select Columns ############################################################
if (LOAD_COLS) {
  basics_cols           <- c("ID", "SPECIES", "ASSEMBLY", 
                             "PMJTL", "GENE", "NCBI_ACC", "UNIPROT_ACC", "SEQ_AA")
  
  basics_extended_cols  <- c("ID", "SPECIES", "ASSEMBLY", "ASSEMBLY_new", 
                             "PMJTL", "GENE", "NCBI_ACC", "UNIPROT_ACC", "SEQ_AA", 
                             "INFO")
  
  annotation_cols       <- c("ID","NCBI_ACC", "UNIPROT_ACC", "GENE", "ANNOTATED_DOMAINS")
  
  gBlock_cols           <- c("ID","")
  
  cl_gen_cols           <- c("ID","CELL_LINE", "COLONY", "ELISA ASSAY DATE")
  
  lci_cols              <- c("")
  
  elisa_cols            <- c("ID", 
                             "ELISA ASSAY DATE", "ELISA_ #1 DATE", "ELISA_ #2 DATE", "ELISA_ #3 DATE", "ELISA #4 DATE")
  
  ncbi_genome_cols      <- c("ID", "ASSEMBLY", "Assembly Name", "Organism Name", "Organism Infraspecific Names Strain", 
                             "Organism Infraspecific Names Isolate", "Annotation Name", "Assembly Stats Total Sequence Length", 
                             "Assembly Level", "Assembly Release Date", "WGS project accession", "Assembly Stats Contig N50", 
                             "Assembly Stats Scaffold N50", "Assembly Sequencing Tech", "Assembly Submitter", "Assembly BioProject Accession", 
                             "Assembly BioSample Accession", "Annotation Count Gene Total", "Annotation Count Gene Protein-coding", 
                             "Annotation Count Gene Pseudogene", "Type Material Display Text", 
                             "CheckM marker set", "CheckM completeness", "CheckM contamination")
}


if (LOAD_PATHS) {
  input_dir  <- "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/01_raw_data"
  temp_dir   <- "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/02_temp_data"
  output_dir <- "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/03_output_data"
  
  input_files <- list.files(input_dir,  full.names = TRUE)
  temp_files  <- list.files(temp_dir,   full.names = TRUE)
  out_files   <- list.files(output_dir, full.names = TRUE)
  
  # group file lists
  excel_files <- input_files[grepl("\\.xlsx$", input_files,  ignore.case = TRUE)]
  csv_files   <- input_files[grepl("\\.csv$",  input_files,  ignore.case = TRUE)]
  tsv_files   <- input_files[grepl("\\.tsv$",  input_files,  ignore.case = TRUE)]
}

if (OPEN_DATA) {
  # Load DD Library Data
  sheets <- lapply(excel_files, read_excel_allsheets)      # Read all sheets from all Excel files
  sheets <- flatten(sheets)                                # Flatten the list of lists into a single list of data frames
  sheets_with_id <- keep(sheets, ~ "ID" %in% colnames(.x)) # Filter out sheets that do not contain the 'ID' column
  dd_lib_tmp     <- reduce(sheets_with_id, full_join, by = "ID") %>% arrange(ID) # Join all data frames by the 'ID' column
  
  # Load Experimental Data
  dates         <- c("20230313", "20230313_2", "20230612", "20230912", "20240516")
  dd_lab_list <- lapply(dates, function(date) {fread(csv_files[grepl(paste0(date, "_LAB.csv"), csv_files, ignore.case = TRUE)], na.strings = c("", "NA", "-", "-------------", NA))})
  dd_lab_tmp      <- rbindlist(dd_lab_list, fill = TRUE) %>% 
    mutate(ID = str_replace(ID, "(DLD|BDLD[12])_(\\d)(?!\\d)", "\\1_0\\2")) %>% # Clean up the ID column
    select(-c("NAME_2", "Sequence Name"))
}



if (CLEAN) {
  
  # Join  DD Library and Experimental Data
  sota_tmp <- full_join(dd_lib_tmp, dd_lab_tmp) %>% filter(grepl("DLD", ID)) %>% arrange(ID)
  
  # fix smaller column issues
  species_df      <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("Species", "SPECIES.x", "Species.x", "Species.x.x", "SPECIES.y", "Species.y", "Species.y.y"), new_col_name = "SPECIES")
  aa_start_df     <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("AA_START.x", "AA_START.y"), new_col_name = "AA_START")
  aa_end_df       <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("AA_END.x", "AA_END.y"), new_col_name = "AA_END")
  gene_df         <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("Abbr", "Abbr.x", "Abbr.x.x", "Abbr.y", "Abbr.y.y"), new_col_name = "GENE")
  info_df         <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("Info", "INFO"), new_col_name = "INFO")
  pmjtl_df        <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("pMJTL", "pMJTL_no", "pMJTL.x", "pMJTL.y"), new_col_name = "PMJTL") %>% filter(!is.na(ID))
  construct_df    <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("construct"), new_col_name = "CONSTRUCT")
  anno_domains_df <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("ANNOTATED_DOMAINS NCBI [UniProt].x", "ANNOTATED_DOMAINS NCBI [UniProt].y"), new_col_name = "ANNOTATED_DOMAINS")
  
  
  # Accession Chaos
  up_df   <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("UNIPROT_ACCESSION.x", "UNIPROT_ACCESSION.y", "UNIPROT_ID", "Accession", "ACCESSION", "ACCESSION.x", "ACCESSION.y", "ACCESSION.x.x", "ACCESSION.y.y"), new_col_name = "UNIPROT_ACC")
  ncbi_df <- sota_tmp %>% fill_subset(id = "ID", cols_to_join = c("NCBI_Abbr", "NCBI_ACCESSION.x", "NCBI_ACCESSION.y"), new_col_name = "NCBI_ACC")
  acc_df  <- sota_tmp %>% group_by(ID) %>% arrange(ID) %>% select(ID, GenBank) %>%
    left_join(ncbi_df) %>%
    left_join(up_df) %>%
    mutate(NCBI_ACC    = case_when(is.na(NCBI_ACC) ~ GenBank, T ~ NCBI_ACC),
           UNIPROT_ACC = case_when(UNIPROT_ACC == NCBI_ACC ~ NA, T ~ UNIPROT_ACC)) %>% 
    select(-GenBank) %>%
    distinct(ID, .keep_all = T)
  
  
  # initiate clean SOTA DF
  sota_final <- data.frame(ID = sota_tmp$ID) %>% 
    unique() %>% 
    left_join(species_df) %>%
    left_join(aa_start_df) %>% left_join(aa_end_df) %>% 
    left_join(gene_df) %>%
    left_join(acc_df) %>%
    left_join(anno_domains_df) %>%
    
    arrange(ID)
  
  if (RM_CLEANED_COLS) {
    
    sota_tmp <- sota_tmp %>% 
      select(-c("Species", "SPECIES.x", "Species.x", "Species.x.x", "SPECIES.y", "Species.y", "Species.y.y",
                "AA_START.x", "AA_START.y", "AA_END.x", "AA_END.y", "AA_LEN.x", "AA_LEN.y",
                "Abbr", "Abbr.x", "Abbr.x.x", "Abbr.y", "Abbr.y.y",
                "UNIPROT_ACCESSION.x", "UNIPROT_ACCESSION.y", "UNIPROT_ID", 
                "Accession", "ACCESSION.x", "ACCESSION", "ACCESSION.y", "ACCESSION.x.x", "ACCESSION.y.y",
                "NCBI_Abbr", "NCBI_ACCESSION.x", "NCBI_ACCESSION.y",
                "ANNOTATED_DOMAINS NCBI [UniProt].x", "ANNOTATED_DOMAINS NCBI [UniProt].y",
                "pMJTL", "pMJTL_no", "pMJTL.x", "pMJTL.y",
                "construct",
                
      ))
    
    nms_sota_tmp <- as.data.frame(sort(names(sota_tmp))) %>% unique()
    nms_sota_fnl <- as.data.frame(sort(names(sota_final))) %>% unique()
  }
  
  if (GENERATE_GENOME_TABLE) {
    # Load NCBI Genome Info Data
    ncbi_genome_data <- fread(file.path(list.files(input_dir, pattern = "\\NCBI_GENOMES.tsv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
    setnames(ncbi_genome_data, "Assembly Accession", "ASSEMBLY")
    
    # Join Library Data with ncbi genome information
    genome_info_tbl <- sota_tmp[,c("ASSEMBLY", "ASSEMBLY_new", "ID")] %>% 
      pivot_longer(cols = 1:2, names_to = "ORIGINAL_COL", values_to = "ASSEMBLY") %>%
      filter(!is.na(ASSEMBLY)) %>%
      select(-ORIGINAL_COL) %>%
      left_join(ncbi_genome_data) %>%
      # remove all cols with NA only
      select_if(~ !all(is.na(.)))
  }
  
  if (OVERWRITE) {
    
    # save sota_tmp as CSV
    fwrite(sota_tmp, file = file.path(temp_dir, "SOTA_tmp.csv"), append = F, row.names = F, col.names = T)
    print(paste("Overwriting existing CSV file:", file.path(temp_dir, "SOTA_tmp.csv")))
    sota_tmp <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA_tmp.csv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
    
    # save sota_final as CSV
    fwrite(sota_final, file = file.path(output_dir, "SOTA.csv"), append = F, row.names = F, col.names = T)
    print(paste("Overwriting existing CSV file:", file.path(output_dir, "SOTA.csv")))
    sota_final <- fread(file.path(list.files(output_dir, pattern = "\\SOTA.csv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
    
    # save NCBI Genome Info Data as CSV
    fwrite(genome_info_tbl, file = file.path(output_dir, "genome_info.csv"), append = F, row.names = F, col.names = T)
    
  } else {
    # or simply load pre-existing data if you don't want to overwrite
    print("Not overwriting existing CSV file. Set OVERWRITE = TRUE to generate new SOTA CSV summary files!")
    sota_tmp   <- fread(file.path(list.files(temp_dir, pattern = "\\SOTA_tmp.csv$", ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
    sota_final <- fread(file.path(list.files(output_dir, pattern = "\\SOTA.csv$",   ignore.case = TRUE, full.names = TRUE)), na.strings = c("", "NA", NA))
  }
}

vis_miss(sota_final)




################################################################################
##### COLUMN DOCUMENTATION #####################################################
################################################################################

# ID         .. Name
# PMJTL      .. Plasmid ID              --> pMJTL, pMTJL_no, transduced
# CONSTRUCT  .. Plasmid Name            --> construct
# UNIPROT    .. UniProt Accession       --> Accession, ACCESSION
# LOCUS      .. NCBI Locus Tag          --> Abbr
# NCBI       .. NCBI Accession          --> NCBI_Abbr
# GENBANK    .. GenBank Accession       --> Genbank
# COMMON     .. Common Name             --> Info, INFO
# SPECIES    .. Species Name            --> Species
# TAXID      .. Taxonomy ID             --> taxid
# SALES_NO   .. gBlocks Sales Number    --> `Sales Order Number`
# REFERENCE  .. gBlocks Ref No          --> `Reference Number`
# MAN_ID     .. gBlocks manufact. ID    --> `Manufacturing ID`
# MOL_WT     .. gBlocks calc mol weight --> `Calculated Molecular Weight`
# NORM       .. gBlocks normalization   --> Normalization
# QUANTITY   .. gBlocks mass [ng]       --> Quantity
# SEQ_LEN_NT .. gBlocks Sequence Length --> `Sequence Length`
# FL_AA      .. FL AA sequence          --> `Full Sequence (aa, not prepped as gblock)`
# GBLOCK_NT  .. gBlocks Sequence        --> `Sequence (nt)`
# AF2C_CONF  .. Confidence in AF2C Pred --> `Confidence`
# ENV        .. Natural Environment     --> `Environment isolated from`
# GROWTH_TMP .. Env. Growth Temp.       --> `Likely Growth Properties`
# MULTICELL  .. Multicellularity known  --> Multicellularity
# RH_DamFRET .. RH DamFRET Assay performed --> `Randal Halfmann`


# GENOME_AVAILABLE .. LOGICAL TRUE/FALSE   --> GENOME_AVAILABLE
# GENOME_LEVEL     .. FACTOR contig/scaffold/chromosome/complete

################################################################################
### COLUMN TABLE GENOME INFO TABLE #############################################
################################################################################

# "ID"
# "ASSEMBLY"
# "Assembly Name"                       
# "Organism Name"                        
# "Organism Infraspecific Names Strain"  
# "Organism Infraspecific Names Isolate"
# "Annotation Name"                      
# "Assembly Stats Total Sequence Length" "Assembly Level"                      
# "Assembly Release Date"                
# "WGS project accession"                
# "Assembly Stats Contig N50"           
# "Assembly Stats Scaffold N50"          
# "Assembly Sequencing Tech"             
# "Assembly Submitter"                  
# "Assembly BioProject Accession"        
# "Assembly BioSample Accession"
# "Annotation Count Gene Total"         
# "Annotation Count Gene Protein-coding"
# "Annotation Count Gene Pseudogene"
# "Type Material Display Text"          
# "CheckM marker set"
# "CheckM completeness"
# "CheckM contamination"



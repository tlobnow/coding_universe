### SPLITTING 20241029_PF20690_architecture_bacterial-Death-like-domain-3.tsv ###

pacman::p_load(tidyverse, data.table, stringr, Biostrings, seqinr, bioseq, glue)
source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")

RUN_FROM_SCRATCH = T

if (RUN_FROM_SCRATCH) {
  
  ARCHITECTURAL_DIGEST = F
  FRESH_FASTA          = F
  FILL_INFO_GAPS       = F
  SPLIT_DOMAINS        = F
  PASTE_DOMAINS        = F
  EXTEND_DOMAINS       = F
  PULL_IT_TOGETHER     = F
  CREATE_FASTA_FILES   = F
  JOIN_WITH_SOTA       = F
  
  if (ARCHITECTURAL_DIGEST) {
    representative_architectures <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/01_input/20241029_PF20690_architecture_bacterial-Death-like-domain-3.tsv")
    architecture_groups          <- list.files("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/01_input", pattern = "architecture_group", full.names = TRUE)
    architecture_groups_df       <- data.frame()
    
    for (file in architecture_groups) {
      temp <- fread(file)
      temp$filename <- file
      architecture_groups_df <- rbind(architecture_groups_df, temp, fill = TRUE)
      architecture_groups_df <- architecture_groups_df[architecture_groups_df$Accession != TRUE, ]
      architecture_groups_df$group_representative_acc <- gsub(".*architecture_group_\\d{1,2}_([A-Za-z0-9]+)\\.tsv$", "\\1", architecture_groups_df$filename)
    }
    
    architecture_groups_extended_df <- architecture_groups_df %>% 
      full_join(representative_architectures, by = c("group_representative_acc" = "Representative Accession")) %>%
      mutate(Accession = case_when(is.na(Accession) ~ group_representative_acc, TRUE ~ Accession),
             Length = case_when(is.na(Length) ~ `Representative Length`, TRUE ~ Length)) %>%
      select(Accession, Name, `Tax ID`, `Tax Name`, Length, group_representative_acc, `Unique Proteins`, `Representative Domains`)
    
    # replace spaces from colnames with underscores
    colnames(architecture_groups_extended_df) <- gsub(" ", "_", colnames(architecture_groups_extended_df))
    colnames(architecture_groups_extended_df) <- tolower(colnames(architecture_groups_extended_df))
    
    # write the data frame to a file
    fwrite(architecture_groups_extended_df, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/01_architecture_groups.csv")
  }
  
  # load outcome that should have been created by code chunk above
  architecture_groups_extended_df <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/01_architecture_groups.csv", na.strings = c("NA", "NA_character_", ""))
  
  if (FRESH_FASTA) {
    fasta_accessions <- architecture_groups_extended_df %>% select(accession) %>% distinct()
    
    for (i in 1:nrow(fasta_accessions)) {
      fasta_accessions$fasta_FL[i] <- ""
      fasta_accessions$Length[i] <- 0
      
      url <- paste0("https://www.uniprot.org/uniprot/", fasta_accessions$accession[i], ".fasta")
      temp <- readLines(url)
      
      # save first row as description
      fasta_accessions$description[i] <- temp[1]
      
      # remove the header
      temp <- temp[-1]
      
      # concatenate the rest of the lines
      fasta_accessions$fasta_FL[i] <- paste(temp, collapse = "")
      
      # count the number of amino acids
      fasta_accessions$length[i] <- nchar(fasta_accessions$fasta_FL[i]) ; fasta_accessions$length[i]
    }
    
    fasta_accessions$description <- gsub("*>tr\\|", "", fasta_accessions$description)
    
    fasta_accessions <- fasta_accessions %>% select(-Length)
    
    fwrite(fasta_accessions, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_df.csv")
  }
  
  # load outcome that should have been created by code chunk above
  fasta_accessions <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_df.csv", na.strings = c("NA", "NA_character_", ""))
  
  if (FILL_INFO_GAPS) {
    
    main <- left_join(architecture_groups_extended_df, fasta_accessions)
    
    main_filled <- main %>%
      mutate(name     = case_when(is.na(name)     ~ sub(".*\\|[^|]+\\|([^OX=]+)|\\([^OS=]+)\\s.*", "\\1", description),  TRUE ~ name),
             tax_name = case_when(is.na(tax_name) ~ sub(".*OS=([^OX=]+)\\sOX=.*",   "\\1",                description),  TRUE ~ tax_name),
             tax_id   = case_when(is.na(tax_id)   ~ as.numeric(sub(".*OX=(\\d+).*", "\\1",                description)), TRUE ~ tax_id)
             )
    
    fwrite(main_filled, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/03_main_df.csv")
    
  }
  
  # load outcome that should have been created by code chunk above
  main_filled <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/03_main_df.csv")

  if (SPLIT_DOMAINS) {
    
    domain_data <- data.frame(
      Accession     = character(),
      pfam_id       = character(),
      pfam_name     = character(),
      interpro_id   = character(),
      interpro_name = character(),
      start         = integer(),
      end           = integer(),
      length        = integer(),
      stringsAsFactors = FALSE
    )
    
    # Split file into individualrepresentative domains
    rep_doms <- strsplit(main_filled$representative_domains, ",")
    
    for (i in 1:length(rep_doms)) {
      
      acc <- main_filled$accession[[i]] ; acc
      
      # rep_dom <- rep_doms[[1]]
      for (rep_dom in rep_doms[[i]]) {
        
        # Extract PFAM ID, InterPro ID and respective names
        parts   <- unlist(strsplit(rep_dom, "\\{|\\}:|\\}\\[|\\]|\\:")) ; parts
        pid     <- parts[1] ; pid
        pname   <- parts[2] ; pname
        ip_id   <- parts[3] ; ip_id
        ip_name <- parts[4] ; ip_name
        domain_boundary <- as.numeric(unlist(strsplit(parts[5], "\\-"))) ; domain_boundary
        domain_start    <- domain_boundary[1] ; domain_start
        domain_end      <- domain_boundary[2] ; domain_end

        # Add to data frame
        domain_data <- rbind(domain_data, data.frame(
          accession     = acc,
          pfam_id       = pid,
          pfam_name     = pname,
          interpro_id   = ip_id,
          interpro_name = ip_name,
          start         = domain_start,
          end           = domain_end,
          domain_length = domain_end - domain_start + 1,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    domain_data <- domain_data %>% arrange(accession, start)
    
    # fix minor poopie
    domain_data_fixed <- domain_data
    domain_data_fixed$start[domain_data_fixed$accession         == "A0A836SDX0" & domain_data_fixed$pfam_name == "Trypsin_2"] <- 29
    domain_data_fixed$end[domain_data_fixed$accession           == "A0A836SDX0" & domain_data_fixed$pfam_name == "Trypsin_2"] <- 167
    domain_data_fixed$interpro_name[domain_data_fixed$accession == "A0A836SDX0" & domain_data_fixed$pfam_name == "Trypsin_2"] <- NA
    domain_data_fixed$interpro_id[domain_data_fixed$accession   == "A0A836SDX0" & domain_data_fixed$pfam_name == "Trypsin_2"] <- NA
    
    fwrite(domain_data_fixed, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/04_domain_data.csv")
    
  }
  
  # load outcome that should have been created by code chunk above
  domain_data_fixed <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/04_domain_data.csv", na.strings = c("NA", "NA_character_", ""))

  if (PASTE_DOMAINS) {
    pasted_domains <- domain_data_fixed %>%
      group_by(accession) %>%
      mutate(architecture = paste0(pfam_name, collapse = "___")) %>% 
      rowwise() %>%
      # add information about domain location (C-terminal, N-terminal, internal) depending on pfam_name and Architecture
      mutate(domain_position = case_when(
        grepl(paste0("^",   pfam_name, "___"), architecture) ~ "N",
        grepl(paste0("___", pfam_name, "___"), architecture) ~ "i",
        grepl(paste0("___", pfam_name, "$"),   architecture) ~ "C",
        TRUE ~ NA)) %>% 
      ungroup()
    
    # write the data frame to a file
    fwrite(pasted_domains, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/05_domain_data.csv")
  }
  
  # load outcome that should have been created by code chunk above
  pasted_domains <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/05_domain_data.csv", na.strings = c("NA", "NA_character_", ""))

  if (EXTEND_DOMAINS) {
    
    limit_df <- pasted_domains %>% group_by(accession) %>% 
      filter(!is.na(pfam_name)) %>%
      left_join(fasta_accessions)
    
    # Use reported domain boundaries to extract the corresponding amino acid sequence of the domain
    limit_df$fasta_domain <- substr(limit_df$fasta_FL, limit_df$start, limit_df$end)
    
    # Extend domain boundaries based on surrounding domains (to include the surrounding unstructured bits)
    extended_df <- extend_boundaries(limit_df)
    
    extended_df_reordered <- extended_df %>%
      mutate(length_fasta_domain = end - start + 1,
             length_fasta_extended = extended_end - extended_start + 1) %>%
      select(accession, pfam_id, pfam_name, architecture, length, 
             domain_position, start, end, length_fasta_domain,
             extended_start, extended_end, length_fasta_extended,
             fasta_domain, fasta_extended, fasta_FL, description)
    
    # write the data frame to a file
    fwrite(extended_df_reordered, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/06_extended_domains.csv")
  }
  
  extended_df <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/06_extended_domains.csv", na.strings = c("NA", "NA_character_", ""))
  
  if (PULL_IT_TOGETHER) {
    
    final_df <- extended_df %>%
      left_join(main_filled %>% select(-c(representative_domains))) %>%
      select(accession, tax_id, tax_name,
             pfam_id, pfam_name, architecture, length, name,
             domain_position, start, end, length_fasta_domain,
             extended_start, extended_end, length_fasta_extended,
             fasta_domain, fasta_extended, fasta_FL, description,
             group_representative_acc, unique_proteins)
    
    # Apply backtranslation to each amino acid sequence in the dataframe
    final_df_bt <- final_df %>%
      mutate(nt_sequence_domain   = sapply(fasta_domain, backtranslate, codon_table = ecoli_codon_table),
             nt_sequence_extended = sapply(fasta_extended, backtranslate, codon_table = ecoli_codon_table),
             nt_sequence_FL       = sapply(fasta_FL, backtranslate, codon_table = ecoli_codon_table))
    
    fwrite(final_df_bt, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/07_final_df.csv")
  }
  
  final_df <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/07_final_df.csv", na.strings = c("NA", "NA_character_", ""))
  
  if (CREATE_FASTA_FILES) {
    
    domain = "bDLD3"
    date <- Sys.Date() %>% gsub("-", "", .)
    
    fasta_data <- prepare_fasta_data(final_df, domain_name = domain)  # using the previously created extended_df
    
    fasta_data[["strict"]] <- fasta_data[["strict"]] %>% filter(!is.na(start))
    fasta_data[["extended"]] <- fasta_data[["extended"]] %>% filter(!is.na(extended_start))
    
    write_domain_fasta(fasta_data,
                       strict_filename   = paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/", date, "_", domain, "_strict_sequences.fasta"),
                       extended_filename = paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/", date, "_", domain, "_extended_sequences.fasta")
    )
  }
  
}



 
# domain <- "bDLD3"
# 
# # Prepare data and write to files
# fasta_data <- prepare_fasta_data(extended_df, domain_name = domain)  # using the previously created extended_df
# 
# write_domain_fasta(fasta_data,
#                    strict_filename = paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/",   domain, "_strict_sequences.fasta"),
#                    extended_filename = paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/", domain, "_extended_sequences.fasta")
# )
# 
# 
# tsv_ag_merged_out %>% 
#   filter(pfam_name %in% domain) %>%
#   mutate(fasta_domain = case_when(is.na(fasta_domain) ~ fasta_FL, TRUE ~ fasta_domain)) %>%
#   rowwise() %>% mutate(fasta_entry = paste0(">", description, "\n", fasta_domain)) %>% pull(fasta_entry) %>%
#   writeLines("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/20241030_PF20690.fasta")
# 
# # read the fasta
# fasta <- readAAStringSet("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/20241030_PF20690.fasta")
# 
# # perform a MSA with the sequences
# aln_msa <- msa(fasta)
# print(aln_msa, show = "complete")
# 
# # write the MSA to a file
# writeXStringSet(as(aln_msa, "AAStringSet"), filepath = "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/20241030_PF20690.aln")

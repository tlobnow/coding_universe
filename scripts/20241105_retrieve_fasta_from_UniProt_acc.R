### DOWNLOAD FASTA FILES

pacman::p_load(tidyverse, data.table, stringr, Biostrings, seqinr, bioseq, glue)
source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")

architecture_groups_extended_df <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/01_architecture_groups.csv", na.strings = c("NA", "NA_character_", ""))

fasta_accessions <- architecture_groups_extended_df %>% select(accession) %>% distinct()
taxids <- architecture_groups_extended_df %>% select(tax_id) %>% distinct()

fwrite(fasta_accessions, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_files/fasta_accessions.csv", col.names = F, row.names = F)
fwrite(taxids, "/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_files/taxids.txt", col.names = F, row.names = F)

for (i in 1:nrow(fasta_accessions)) {
  fasta_accessions$fasta_FL[i] <- ""
  fasta_accessions$Length[i] <- 0
  
  url <- paste0("https://www.uniprot.org/uniprot/", fasta_accessions$accession[i], ".fasta")
  
  # download fasta file
  fasta_file <- readLines(url) ; fasta_file
  
  # create fasta directory in 02_temp
  dir.create("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_files", showWarnings = FALSE)
  
  # write fasta file
  writeLines(fasta_file, paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_files/", fasta_accessions$accession[i], ".fasta"))
  
}

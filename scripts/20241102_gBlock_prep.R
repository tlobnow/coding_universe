### gBlock Prep ###

pacman::p_load(tidyverse, data.table, stringr, Biostrings, seqinr, bioseq, glue)

source("/Users/u_lobnow/Documents/Github/coding_universe/scripts/SOTA_FX.R")

PREP_GBLOCKS <- T
  
if (PREP_GBLOCKS) {
  
  main_filled      <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/03_main_df.csv")
  pasted_domains   <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/05_domain_data.csv", na.strings = c("NA", "NA_character_", ""))
  fasta_accessions <- fread("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/02_temp/02_fasta_df.csv", na.strings = c("NA", "NA_character_", ""))
  sota             <- fread(file.path("/Users/u_lobnow/Documents/Github/coding_universe/SOTA/03_output_data/SOTA.csv"), na.strings = c("", "NA", NA)) %>% filter(!is.na(gblock_man_id))
  
  PREP_REPRESENTATIVES <- T 
  
  
  if (PREP_REPRESENTATIVES) {
    
    df              <- left_join(main_filled, pasted_domains)
    df$fasta_domain <- substr(df$fasta_FL, df$start, df$end)
    
    extended_df <- extend_boundaries(df)
    extended_df <- extended_df %>%
      filter(pfam_name == "bDLD3",
             accession %in% unique(extended_df$group_representative_acc)) %>%
      mutate(length_fasta_domain = end - start + 1,
             length_fasta_medium = medium_end - medium_start + 1,
             length_fasta_extended = extended_end - extended_start + 1) %>%
      select(accession, pfam_id, pfam_name, architecture, length, 
             domain_position, start, end, length_fasta_domain,
             medium_start, medium_end, length_fasta_medium,
             extended_start, extended_end, length_fasta_extended,
             fasta_domain, fasta_medium, fasta_extended, fasta_FL, description)
    
    fwrite(extended_df, paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/", gsub(x = Sys.Date(), pattern = "-", replacement = ""), "_representative_architecture_badds.csv"))
    
    rm(df, main_filled, pasted_domains, fasta_accessions)
  }
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  
  # Find the latest file that was created
  extended_df <- fread(tail(sort(list.files("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/",
                                            pattern = "_representative_architecture_badds.csv$", 
                                            full.names = TRUE)), 1))

  # Prepare the gBlocks for BaDDS and CHARMS
  badds_sequences_prep <- extended_df %>% 
    filter(pfam_name == "bDLD3") %>% 
    select(accession, architecture, fasta_domain, fasta_medium, fasta_extended)
  
  # the bDLD3-only group is more diverese.. 
  # The representative domain architecture of I9KD35 (NCBI Accession EIV91905.1 or FraQA3DRAFT_1390) 
  # is probably more accurately described as "TCAD9___bDLD3"
  # since this architecture is not reflected in the main representation set, we can include it as a test case.
  
  # BDLD_57 and BDLD_06 also have a domain at the start that is not captured by the PFAM search:
  #   EAD9 (Kaur et al. 2021) aka
  #   COG4913 (CDD annotation: Uncharacterized conserved protein, contains a C-terminal ATPase domain [Function unknown])
  
  badds_sequences <- list()
  
  for (acc in badds_sequences_prep$accession) {
    
    data      <- badds_sequences_prep %>% filter(accession == acc) #print(data)
    
    badds_sequences[[acc]] <- list(
      fasta_domain   = data$fasta_domain,
      fasta_medium   = data$fasta_medium,
      fasta_extended = data$fasta_extended
    )
  }
  
  badds_sequences[["BDLD_57"]] <- list(
    fasta_domain   = "QTKIFICRRLTQDWQDLADYFEIQPHERAGFKPGREPHSIWEWLEQRNRLGELESALIEIGREDLAQELKKN",
    fasta_medium   = "PGQTKIFICRRLTQDWQDLADYFEIQPHERAGFKPGREPHSIWEWLEQRNRLGELESALIEIGREDLAQELKKN",
    fasta_extended = "NPNIIDKTPNLPVVELPANKQIPGQTKIFICRRLTQDWQDLADYFEIQPHERAGFKPGREPHSIWEWLEQRNRLGELESALIEIGREDLAQELKKN"
  )
  
  badds_sequences[["BDLD_27"]] <- list(
    fasta_domain   = "KIKISICNKLVRDWEDLADYFDIKLHERETFDKGKEARRVWEWLEQRTKLHELEDAFIAIGREDLVEELNK",
    fasta_medium   = "SGKIKISICNKLVRDWEDLADYFDIKLHERETFDKGKEARRVWEWLEQRTKLHELEDAFIAIGREDLVEELNK",
    fasta_extended = "MNQEALPQSNLKPQNLVNENPVKATSTVTTKNYSGKIKISICNKLVRDWEDLADYFDIKLHERETFDKGKEARRVWEWLEQRTKLHELEDAFIAIGREDLVEELNK"
  )
  
  badds_sequences[["MYD88_MOUSE_DD"]] = list(
    #MSAGDPRVGSGSLDSFMFSIPLVALNVGVRRRLSLFLNPRTPVAADWTLLAEEMGFEYLEIRELETRPDPTRSLLDAWQGRSGASVGRLLELLALLDREDILKELKSRIEEDCQKYLGKQQNQESEKPLQVARVESSVPQTKELGGITTLDDPLGQTPELFDAFICYCPNDIEFVQEMIRQLEQTDYRLKLCVSDRDVLPGTCVWSIASELIEKRCRRMVVVVSDDYLQSKECDFQTKFALSLSPGVQQKRLIPIKYKAMKKDFPSILRFITICDYTNPCTKSWFWTRLAKALSLP
    fasta_domain   = "MGFEYLEIRELETRPDPTRSLLDAWQGRSGASVGRLLELLALLDREDILKELKSRI",
    # same as boundaries of 3MOP human DD from Hao Wu paper
    fasta_medium   = "IPLVALNVGVRRRLSLFLNPRTPVAADWTLLAEEMGFEYLEIRELETRPDPTRSLLDAWQGRSGASVGRLLELLALLDREDILKELKSRIEEDCQKYLG",
    fasta_extended = "MSAGDPRVGSGSLDSFMFSIPLVALNVGVRRRLSLFLNPRTPVAADWTLLAEEMGFEYLEIRELETRPDPTRSLLDAWQGRSGASVGRLLELLALLDREDILKELKSRIEEDCQKYLGKQQNQE"
  )
  
  badds_sequences[["IRAK4_MOUSE_DD"]] <- list(
    fasta_domain   = "RKLSDFIDPQEGWKKLAVAIKKPSGDDRYNQFHIRRFEALLQTGKSPTCELLFDWGTTNCTVGDLVDLLVQIELFAPATLLLPDA",
    fasta_medium   = "MNKPLTPSTYIRNLNVGILRKLSDFIDPQEGWKKLAVAIKKPSGDDRYNQFHIRRFEALLQTGKSPTCELLFDWGTTNCTVGDLVDLLVQIELFAPATLLLPDA",
    fasta_extended = "MNKPLTPSTYIRNLNVGILRKLSDFIDPQEGWKKLAVAIKKPSGDDRYNQFHIRRFEALLQTGKSPTCELLFDWGTTNCTVGDLVDLLVQIELFAPATLLLPDAVPQTVKSLPPREAATVAQTHGPCQEKDRTSVMPMPKLEHSCEPPDSSSPDNRSVESSDTRFHSFS" #SAGGNRMGEGGFGVVYKGCVNNTIVAVKKLGAMVEISTEELKQQFDQEIKVMATCQHENLVELLGFSSDSDNLCLVYAYMPNGSLLDRLSCLDGTPPLSWHTRCKVAQGTANGIRFLHENHHIHRDIKSANILLDKDFTAKISDFGLARASARLAQTVMTSRIVGTTAYMAPEALRGEITPKSDIYSFGVVLLELITGLAAVDENREPQLLLDIKEEIEDEEKTIEDYTDEKMSDADPASVEAMYSAASQCLHEKKNRRPDIAKVQQLLQEMSA
  )
  
  badds_sequences[["IRAK1_MOUSE_DD"]] <- list(
    fasta_domain   = "MCRFYKVMDALEPADWCQFAALIVRDQTELRLCERSEQRTASVLWPWINRNARVADLVHILTHLQLLRARDIITAWHPPA",
    fasta_medium   = "MAGGPGPGEPVVPGAQHFLYEVPPWVMCRFYKVMDALEPADWCQFAALIVRDQTELRLCERSEQRTASVLWPWINRNARVADLVHILTHLQLLRARDIITAWHPPA",
    fasta_extended = "MAGGPGPGEPVVPGAQHFLYEVPPWVMCRFYKVMDALEPADWCQFAALIVRDQTELRLCERSEQRTASVLWPWINRNARVADLVHILTHLQLLRARDIITAWHPPAPVVPPSTAAPRPSSISAGSEAGDWSPRKLQSSASTFLSPAFPGSQTHSESELLQVPLPVSLGPPLPSSAPSSTKSSPESPVSGLQRAHPSPFCWPFCEISQGTCN" # FSEELRIGEGGFGCVYRAVMRNTTYAVKRLKEEADLEWTMVKQSFLTEVEQLSRFRHPNIVDFAGYCAESGLYCLVYGFLPNGSLEDQLHLQTQACSPLSWPQRLDILLGTARAIQFLHQDSPSLIHGDIKSSNVLLDERLMPKLGDFGLARFSRFAGAKASQSSTVARTSTVRGTLAYLPEEYIKTGRLAVDTDTFSFGVVILETLAGQRAVRTQGAKTKYLKDLIEDEAEEAGVTLKSTQPTLWVGVATDAWAAPIAAQIYKKHLDSRPGPCPPQLGLALAQLACCCMHRRAKKRPPMTQVYKRLEGLQAGPPWELEVAGHGSPSPQENSYMSTTGSAQSGDEPWQPLVVTTRAPAQAAQQLQRSPNQPVESDESVPGLSATLHSWHLTPGSHPSPASFREASCTQGGTTRESSVRSSPGFQPTTMEGSPTGSSSLLSSEPPQIIINPARQKMVQKLALYEEGVLDSLQLLSSGFFPGLDLEPEKSQGPEESDEFQS
  )
  
  badds_sequences[["bCARD_1"]] <- list(
    fasta_domain   = "MVQADYSRAEALAAWTRLSDEFIGNCYVSVRPRHAPAWEVVVASAAGSLRLEAFKRAHDHDFLDRLAVAIGNWEQKAQRPDHEIAQMLDQVG",
    fasta_medium   = NA,
    fasta_extended = NA
  )
  
  charms_sequences <- list(
    "MYD88_MOUSE_DD" = badds_sequences[["MYD88_MOUSE_DD"]],
    "IRAK4_MOUSE_DD" = badds_sequences[["IRAK4_MOUSE_DD"]],
    "IRAK1_MOUSE_DD" = badds_sequences[["IRAK1_MOUSE_DD"]],
    "BDLD_57"        = badds_sequences[["BDLD_57"]]
  )
  
  # gBlock generation function
  generate_gblock <- function(sequence, codon_table, cut_site, kozak, linker, sticky_end) {
    if (is.na(sequence)) return(NA)  # Skip gBlock generation if sequence is NA
    nt_sequence <- backtranslate(sequence, codon_table = codon_table)
    start_codon <- ifelse(substr(nt_sequence, 1, 3) == "ATG", "", "atg")
    gblock <- paste0(cut_site, kozak, start_codon, nt_sequence, linker, sticky_end)
    return(gblock)
  }
  
  # Define constants for gBlock components
  gblock_components <- list(
    badds  = list(cut_site = "cgatcggagcgcggttctaga", kozak = "gccgccacc", linker = "GGAGGTAGTGGCGGGTCC", sticky_end = "gattatggattaatgcaaggtatgacaaaccctgataaag"),
    charms = list(cut_site = "ttggaggcctaggctacgcgt", kozak = "gccgccacc", linker = "GGAGGTAGTGGCGGGTCC", sticky_end = "gagaagcctttacaggtggccaga")
  )
  
  # Generate gBlocks for BaDDS
  badds_gblocks <- tibble::tibble(
    ID = names(badds_sequences),
    fasta_domain = map_chr(names(badds_sequences),    ~ badds_sequences[[.x]]$fasta_domain),
    fasta_medium = map_chr(names(badds_sequences),    ~ badds_sequences[[.x]]$fasta_medium),
    fasta_extended = map_chr(names(badds_sequences),  ~ badds_sequences[[.x]]$fasta_extended),
    
    strict_gblock   = map_chr(names(badds_sequences), ~ generate_gblock(badds_sequences[[.x]]$fasta_domain,   ecoli_codon_table, gblock_components$badds$cut_site, gblock_components$badds$kozak, gblock_components$badds$linker, gblock_components$badds$sticky_end)),
    medium_gblock   = map_chr(names(badds_sequences), ~ generate_gblock(badds_sequences[[.x]]$fasta_medium,   ecoli_codon_table, gblock_components$badds$cut_site, gblock_components$badds$kozak, gblock_components$badds$linker, gblock_components$badds$sticky_end)),
    extended_gblock = map_chr(names(badds_sequences), ~ generate_gblock(badds_sequences[[.x]]$fasta_extended, ecoli_codon_table, gblock_components$badds$cut_site, gblock_components$badds$kozak, gblock_components$badds$linker, gblock_components$badds$sticky_end)),
    ORIGIN = "BaDDS"
  )
  
  # Generate gBlocks for CHARMS
  charms_gblocks <- tibble::tibble(
    ID = names(charms_sequences),
    fasta_domain = map_chr(names(charms_sequences),    ~ charms_sequences[[.x]]$fasta_domain),
    fasta_medium = map_chr(names(charms_sequences),    ~ charms_sequences[[.x]]$fasta_medium),
    fasta_extended = map_chr(names(charms_sequences),  ~ charms_sequences[[.x]]$fasta_extended),
    
    strict_gblock = map_chr(names(charms_sequences),   ~ generate_gblock(charms_sequences[[.x]]$fasta_domain, mouse_codon_table, gblock_components$charms$cut_site, gblock_components$charms$kozak, gblock_components$charms$linker, gblock_components$charms$sticky_end)),
    medium_gblock = map_chr(names(charms_sequences),   ~ generate_gblock(charms_sequences[[.x]]$fasta_medium, mouse_codon_table, gblock_components$charms$cut_site, gblock_components$charms$kozak, gblock_components$charms$linker, gblock_components$charms$sticky_end)),
    extended_gblock = map_chr(names(charms_sequences), ~ generate_gblock(charms_sequences[[.x]]$fasta_extended, mouse_codon_table, gblock_components$charms$cut_site, gblock_components$charms$kozak, gblock_components$charms$linker, gblock_components$charms$sticky_end)),
    ORIGIN = "CHARMS"
  )
  
  # paste both data frames and assign either BaDDS or charms to the ORIGIN column
  gblocks <- rbind(badds_gblocks, charms_gblocks)
  setnames(gblocks, "ID", "accession")
  gblocks <- left_join(gblocks, extended_df)
  
  # transpose fasta_domain	fasta_medium	fasta_extended to column fasta
  gblocks_final <- gblocks %>%
    pivot_longer(cols = c(strict_gblock, medium_gblock, extended_gblock), names_to = "gBlock", values_to = "sequence") %>%
    mutate(gBlock = case_when(
      gBlock == "strict_gblock"   ~ "S",
      gBlock == "medium_gblock"   ~ "M",
      gBlock == "extended_gblock" ~ "L"),
      # create the ORDER_ID by pasting origin, number of the origin and the accession
      ORDER_ID = glue("{ORIGIN}_{gBlock}_{accession}_")
      ) %>%
    select(ORDER_ID, gBlock, sequence, architecture) %>%
    filter(!is.na(sequence))
  
  fwrite(gblocks, paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/", gsub(x = Sys.Date(), pattern = "-", replacement = ""), "_gblocks.csv"))
  fwrite(gblocks_final, paste0("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/", gsub(x = Sys.Date(), pattern = "-", replacement = ""), "_gblocks_final.csv"))
  
}

gblocks <- fread(tail(sort(list.files("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/",
                                    pattern = "_gblocks.csv$", 
                                    full.names = TRUE)), 1))

gblocks_final <- fread(tail(sort(list.files("/Users/u_lobnow/Documents/etc/PhD/DD_Manuscript/20241029_PFAM_InterPro/03_output/",
                                          pattern = "_gblocks_final.csv$", 
                                          full.names = TRUE)), 1))

View(gblocks_final)



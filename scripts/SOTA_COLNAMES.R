
################################################################################
### COLNAMES FOR SOTA ##########################################################
################################################################################


basics_cols           <- c("ID", "species", "cl_name_endo", "cl_name_high")

taxonomy_cols         <- c("ID", "species", "genus", "family", "order", "class", "phylum", "clade", "superkingdom")

accession_cols        <- c("assembly", "gene", "genbank_acc", "uniprot_acc")

basics_extended_cols  <- c("ID", "assembly", "assembly_2", "pmjtl", "gene", "gene_name", "genbank_acc", "uniprot_acc", "info", "taxid", "taxend")

environment_cols      <- c("aquatic", "terrestrial", "saltwater", "freshwater", "soil", "hotspring", "env", "env_temp")

sequence_cols         <- c("seq_FL_aa", "seq_FL_nt", "seq_DLD_aa", "aa_start", "aa_end", "aa_len", "bp_len")

annotation_cols       <- c("genbank_acc", "uniprot_acc", "gene", "annotated_domains", "architecture", "operon", "gene_suppr", "genome_suppr")

alphaFold_cols        <- c("af_confidence", "af3_name", "af3_seq", "af3_n_chains", "af3_model", "af3_iptm", "af3_n_recycles", "af3_ptm", "af3_ranking_score", "af3_seed")

gBlock_cols           <- c("pmjtl", "seq_DLD_aa", "calc_mol_wt", "quantity", "ref_number", 
                           "gblock_norm", "gblock_man_id", "gblock_sales_no", "gblock_seq_nt")

cl_gen_cols           <- c("pmjtl", "cell_line_endo", "cl_line_high", "cell_line", "colony", "abs_260_280", 
                           "conc_ng_ul", "conc_ug_ul", "fm_ng", "use", "plasmid_volume",
                           "tf_date", "td_date", "td_cl204", "facs_date", "facs_yield", "frozen_date", "frozen_date_us")

elisa_cols            <- c("elisa_stim_1", "elisa_stim_2", "elisa_stim_3", "elisa_stim_4", "elisa_stim_5", "elisa_stim_6", "elisa_assay", "elisa_assay_2",
                           "ID", "ORIGIN", "EXPRESSION_LVL",  "POSITIVE_CTRL", "Date", "p_value_rltv",
                           "p_value_real", "significance_rltv", "significance_real",
                           "IL2_concentration_Dilution_Factor_mean_UNSTIM", "IL2_concentration_Dilution_Factor_mean_STIM",
                           "IL2_concentration_Dilution_Factor_sem_UNSTIM", "IL2_concentration_Dilution_Factor_sem_STIM",
                           "Relative_Intensity_mean_UNSTIM", "Relative_Intensity_mean_STIM", 
                           "Relative_Intensity_sem_UNSTIM", "Relative_Intensity_sem_STIM",
                           "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOT_ID", "fold_change", "p_value_fc",
                           "fc_annotation", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR", "CELL_LINE",
                           "PURPOSE", "INFO")


lci_cols              <- c("tirf_coalesce", "tirf_lci_1", "tirf_lci_2","tirf_traf6_puncta", "tirf_lci_3", "tirf_rela_date", "frap_date")

miscellaneous         <- c("reference", "nr", "terminus", "note", "mag")

ncbi_genome_cols      <- c("ID", "ASSEMBLY", "Assembly Name", "Organism Name", "Organism Infraspecific Names Strain", 
                           "Organism Infraspecific Names Isolate", "Annotation Name", "Assembly Stats Total Sequence Length", 
                           "Assembly Level", "Assembly Release Date", "WGS project accession", "Assembly Stats Contig N50", 
                           "Assembly Stats Scaffold N50", "Assembly Sequencing Tech", "Assembly Submitter", "Assembly BioProject Accession", 
                           "Assembly BioSample Accession", "Annotation Count Gene Total", "Annotation Count Gene Protein-coding", 
                           "Annotation Count Gene Pseudogene", "Type Material Display Text", 
                           "CheckM marker set", "CheckM completeness", "CheckM contamination")


# ls(pattern = "_cols")

#     basics_cols,
#     accession_cols,
#     basics_extended_cols,
#     taxonomy_cols,
#     environment_cols,
#     sequence_cols,
#     annotation_cols,
#     alphaFold_cols,
#     gBlock_cols,
#     cl_gen_cols,
#     elisa_cols,
#     lci_cols,
#     miscellaneous,
#     ncbi_genome_cols


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

library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, ggh4x)

REGENERATE_DATA <- FALSE

if (REGENERATE_DATA) {
  # Analysis setup ----
  TablePaths <- c(
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/1_069-1R_MyD88 TRAF6/20230413 2.5nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/1_069-1R_MyD88 TRAF6/20231102 2nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/1_069-1R_MyD88 TRAF6/20231114 2nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/1_069-1R_MyD88 TRAF6/20240307 3nM_cl069_MyD88_TRAF6 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/2_232_MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/2_232_MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/2_232_MyD88-TRAF6-BD TRAF6/20231114 2nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/2_232_MyD88-TRAF6-BD TRAF6/20231114 2nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/3_240_MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/3_240_MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/3_240_MyD88-TIR-TRAF6-BD TRAF6/20230405 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/3_240_MyD88-TIR-TRAF6-BD TRAF6/20231114 2nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/3_240_MyD88-TIR-TRAF6-BD TRAF6/20231115 2nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/3_240_MyD88-TIR-TRAF6-BD TRAF6/20231121 2nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/5_236_MyD88-DHF91-TRAF6-BD TRAF6/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/5_236_MyD88-DHF91-TRAF6-BD TRAF6/20221207 4nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/5_236_MyD88-DHF91-TRAF6-BD TRAF6/20230413 2.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/5_236_MyD88-DHF91-TRAF6-BD TRAF6/20240307 3nM nM_cl236_MyD88_TRAF6_DHF91 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/5_236_MyD88-DHF91-TRAF6-BD TRAF6/20240307 3nM nM_cl236_MyD88_TRAF6_DHF91 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/6_321_BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl321-BDLD57H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/6_321_BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230623 3nM_cl321-BDLD57H_TRAF6_MyD88 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/1_Image analysis/6_321_BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231102 2nM_cl321-BDLD57H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230616_batch/Output/BDLD_2H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl322-BDLD2H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_6H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl328-BDLD6H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230623_batch/Output/BDLD_6H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230623 3nM_cl328-BDLD6H_TRAF6_MyD88 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_7H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl325-BDLD7H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline_2/Analysis/new_pipeline/pending_processing/20231201_batch_1/Output/BDLD_11H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231201 3nM_cl332-BDLD11H_TRAF6_MyD88 003/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230602_batch/Output/BDLD_13H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230602 6nM_cl319-BDLD13H_TRAF6_MyD88_001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230519_batch/Output/BDLD_14H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230519 3nM_cl315-BDLD14H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline_2/Analysis/new_pipeline/pending_processing/20231201_batch_2/Output/BDLD_16H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231201 3nM_cl335-BDLD16H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_22H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl329-BDLD22H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_23H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl330-BDLD23H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_24H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl331-BDLD24H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230519_batch/Output/BDLD_27H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230519 3nM_cl311-BDLD27H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline_2/Analysis/new_pipeline/pending_processing/20231201_batch_2/Output/BDLD_38H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231201 3nM_cl338-BDLD38H_TRAF6_MyD88 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230616_batch/Output/BDLD_46H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl324-BDLD46H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline_2/Analysis/new_pipeline/pending_processing/20231201_batch_1/Output/BDLD_58H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231201 3nM_cl333-BDLD58H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230602_batch/Output/BDLD_62H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230602 6nM_cl318-BDLD62H_TRAF6_MyD88_001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230616_batch/Output/BDLD_62H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl318-BDLD62H_TRAF6_MyD88 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230616_batch/Output/BDLD_67H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl327-BDLD67H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_67H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl327-BDLD67H_TRAF6_MyD88 003/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230623_batch/Output/BDLD_67H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230623 3nM_cl327-BDLD67H_TRAF6_MyD88 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_69H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl326-BDLD69H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230616_batch/Output/BDLD_73H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl323-BDLD73H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline/pending_processing/20230619_batch/Output/BDLD_73H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230619 3nM_cl323-BDLD73H_TRAF6_MyD88 002/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline_2/Analysis/new_pipeline/pending_processing/20231201_batch_2/Output/BDLD_74H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231201 3nM_cl337-BDLD74H_TRAF6_MyD88 001/Essential.csv.gz",
    "/Volumes/TAYLOR-LAB/Finn_v2/01_PIPELINE/new_pipeline_2/Analysis/new_pipeline/pending_processing/20231201_batch_1/Output/BDLD_82H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20231201 3nM_cl334-BDLD82H_TRAF6_MyD88 001/Essential.csv.gz"
  )
  
  file.exists(TablePaths)
  
  # Load the data
  TableList        <- lapply(TablePaths, function(path) {fread(file = path)})
  essential_Table  <- rbindlist(TableList, fill=TRUE) 
  
  # fill = TRUE:  if not all images are two color
  # fill = FALSE: if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA
  
  # rename cohorts if necessary
  # essential_Table$COHORT <- gsub("*-MyD88-TIR-TRAF6-BD-GFP TRAF6", replacement = "-TIR-T6BD TRAF6", x = essential_Table$COHORT)
  # essential_Table$COHORT <- gsub("MyD88-DHF91", replacement = "DHF91-TIR", x = essential_Table$COHORT)
  # essential_Table$COHORT <- gsub("*TRAF6-BD", replacement = "T6BD", x = essential_Table$COHORT)
  
  # get all the different cohorts in Table
  unique(essential_Table$COHORT)
  
  unique(essential_Table$LIGAND_DENSITY_CAT)
  # [1] 32 10
  
  fwrite(x = essential_Table,
         file = "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/02_temp_data/lci_essential_table.csv",
         quote = F,
         sep = ",")
  system("gzip -f /Users/u_lobnow/Documents/Github/coding_universe/SOTA/02_temp_data/lci_essential_table.csv")
}

# Load the data
lci_data <- fread(file = "/Users/u_lobnow/Documents/Github/coding_universe/SOTA/02_temp_data/lci_essential_table.csv.gz")


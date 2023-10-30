#         .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^

load_dtConnOrLocal <- function(PATH_DT, PATH_LOCAL, FILE) {
  if ( DATA_TAY_CONNECTION == F ) { # source from local folder
    read_xlsx(path = paste0(PATH_LOCAL, FILE), sheet = 1, na = c("", " ", "NA", "NaN"))
  } else { # source from Data-Tay
    read_xlsx(path = paste0(PATH_DT, FILE), sheet = 1, na = c("", " ", "NA", "NaN"))
  }}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

retrieveAccessionIDs <- function(DF,OUT="retrievedAccessionIDs.txt") {
  # Select UniProt Accession Numbers and retain unique values
  DF_uniq <- DF %>% select(Protein.IDs) %>% unique()
  
  # filter out contaminants
  DF_fil <- DF_uniq %>% filter(!str_detect(Protein.IDs, paste("CON__")))
  
  # separate the joined protein IDs (sometimes multiple per row, separated by ";")
  DF_sep <- unlist(lapply(strsplit(DF_fil$Protein.IDs, ";", fixed=TRUE), function(x) return(x[1:50]))) %>%
    unique()
  
  # filter out NAs
  DF_sep <- as.data.frame(DF_sep) %>% drop_na()
  
  # filter out invalid Accession IDs ("REV__")
  DF_sep <- DF_sep %>% filter(!str_detect(DF_sep, paste("REV__")))
  
  # write text file
  write.table(x = DF_sep, 
              file = OUT, 
              quote = F, sep = "\t", row.names = F, col.names = F)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


# jsonExtract <- function(JSON, OUT, FILE) {
#   # REMOVE ALL "INFINITY" STRINGS IN YOUR JSON FILES USING:
#   # grep -rl Infinity . | xargs sed -i 's/Infinity/9999/g'
#   json <- fromJSON(JSON)
#   # EXTRACT FILE
#   #FILE_A  <- unlist(lapply(strsplit(json[["chains"]][["A"]], "_", fixed=TRUE), function(x) return(x[1])))
#   #FILE_B  <- unlist(lapply(strsplit(json[["chains"]][["B"]], "_", fixed=TRUE), function(x) return(x[1])))
#   #FILE    <- paste(FILE_A, FILE_B, "x1", sep = "_")
#   FILE    <- FILE
#   MODEL   <- unlist(lapply(strsplit(as.data.frame(json$order)[1,], "_", fixed=TRUE), function(x) return(x[2])))
#   TOL     <- as.data.frame(json$tol_values) %>% pivot_longer(names_to = "RECYCLE", values_to = "TOL", cols = 1:ncol(.))
#   pLDDT   <- as.data.frame(json$plddts) %>% pivot_longer(names_to = "RECYCLE", values_to = "pLDDT", cols = 1:ncol(.))
#   pTM     <- as.data.frame(json$ptms) %>% pivot_longer(names_to = "RECYCLE", values_to = "pTM", cols = 1:ncol(.))
#   piTM    <- as.data.frame(json$pitms) %>% pivot_longer(names_to = "RECYCLE", values_to = "piTM", cols = 1:ncol(.))
#   iScore  <- as.data.frame(json$`interface score`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iScore", cols = 1:ncol(.))
#   iRes    <- as.data.frame(json$`interfacial residue number`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iRes", cols = 1:ncol(.))
#   iCnt    <- as.data.frame(json$`interficial contact number`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iCnt", cols = 1:ncol(.))
#   FILE_MODEL    <- paste(FILE, MODEL, sep = "_")
#   NUM_CLUSTERS  <- json[["clusters"]][[iScore$RECYCLE[1]]][["num_clusters"]]
#   N_MONOMERS    <- length(json[["chains"]])
#   # JOIN, REMOVE DUPLICATES, WRITE TO CSV
#   EXTRACT   <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
#   EXTRACT   <- EXTRACT[, !duplicated(colnames(EXTRACT))]
#   write.table(EXTRACT, file = paste0(OUT,"_fromJSON.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
# }

# jsonExtract <- function(JSON, OUT, FILE) {
#   json       <- fromJSON(JSON)
#   MODEL      <- strsplit(json$order[[1]], "_", fixed = TRUE)[[1]][2]
#   TOL        <- json$tol_values %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "TOL") 
#   pLDDT      <- json$plddts %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pLDDT") 
#   pTM        <- json$ptms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pTM") 
#   piTM       <- json$pitms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "piTM") 
#   iScore     <- json$`interface score` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iScore") 
#   iRes       <- json$`interfacial residue number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iRes") 
#   iCnt       <- json$`interficial contact number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iCnt") 
#   FILE_MODEL <- paste(FILE, MODEL, sep = "_")
#   NUM_CLUSTERS <- json$clusters[[iScore$RECYCLE[1]]]$num_clusters
#   N_MONOMERS <- length(json$chains)
#   DATE       <- sub(".*_(\\d+)_.*", "\\1", names(json$clusters))
#   
#   EXTRACT <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
#   EXTRACT$DATE = lubridate::as_date(DATE)
#   EXTRACT <- EXTRACT[, !duplicated(colnames(EXTRACT))]
#   
#   write.table(EXTRACT, file = paste0(OUT,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
#   
#   EXTRACT_noRecycle <- EXTRACT %>% filter(!str_detect(RECYCLE, "_recycled_"))
#   write.table(EXTRACT_noRecycle, file = paste0(OUT,"_noRecycle.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
# }

jsonExtract <- function(JSON, OUT, FILE) {
  json       <- fromJSON(JSON)
  MODEL      <- strsplit(json$order[[1]], "_", fixed = TRUE)[[1]][2]
  TOL        <- json$tol_values %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "TOL") 
  pLDDT      <- json$plddts %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pLDDT") 
  pTM        <- json$ptms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pTM") 
  piTM       <- json$pitms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "piTM") 
  iScore     <- json$`interface score` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iScore") 
  iRes       <- json$`interfacial residue number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iRes") 
  iCnt       <- json$`interficial contact number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iCnt") 
  FILE_MODEL <- paste(FILE, "MODEL", MODEL, sep = "_")
  NUM_CLUSTERS <- json$clusters[[iScore$RECYCLE[1]]]$num_clusters
  N_MONOMERS <- length(json$chains)
  
  EXTRACT <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT <- EXTRACT[, !duplicated(colnames(EXTRACT))]
  
  write.table(EXTRACT, file = paste0(OUT,"_withRecycles.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
  
  EXTRACT_noRecycle <- EXTRACT %>% filter(!str_detect(RECYCLE, "_recycled_"))
  write.table(EXTRACT_noRecycle, file = paste0(OUT,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

slurmExtract <- function(SLURM, OUT) {
  SLURM        <- fread(SLURM, sep = "\t", header = F)
  SLURM        <- SLURM %>% mutate(RECYCLED = F, TOP = F, INFO = F, GRAB = F, N_MON = F, MODEL = F)
  SLURM$RECYCLED[grep("*_recycled_*",      SLURM$V1)] <- T
  SLURM$INFO[grep("*Info:*",               SLURM$V1)] <- T
  SLURM$GRAB[grep("Info: num_clusters*",   SLURM$V1)] <- T
  SLURM$N_MON[grep("* to model *",         SLURM$V1)] <- T
  SLURM$MODEL[grep("*pLDDT*",              SLURM$V1)] <- T
  
  MODEL_INFO        <- SLURM %>% filter(MODEL == T & RECYCLED == F) %>% select(V1)
  MODEL_INFO        <- separate(data = MODEL_INFO, col = V1,   sep = ",",   into = c("NAME", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes"), convert = T)
  MODEL_INFO        <- separate(data = MODEL_INFO, col = NAME, sep = " ",   into = c("INFO", "FILE", "MODEL", "PERFORMED", "X", "CYCLE"), convert = T)
  MODEL_INFO$TOL    <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$TOL,   "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$pLDDT  <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$pLDDT, "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$pTM    <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$pTM,   "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$piTM   <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$piTM,  "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$iScore <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$iScore,"= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO        <- separate(data = MODEL_INFO, col = iRes, sep =  "iCnt = ", into = c("iRes", "iCnt"), convert = T)
  MODEL_INFO$iRes   <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$iRes,  "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$iCnt   <- as.numeric(MODEL_INFO$iCnt)
  MODEL_INFO$Clash_Indicator <- MODEL_INFO$iRes / MODEL_INFO$iCnt
  
  CLUSTER_INFO               <- SLURM %>% filter(GRAB == T) %>% select(V1)
  CLUSTER_INFO               <- separate(data = CLUSTER_INFO, col = V1,   sep = " = ", into = c("V1", "NUM_CLUSTERS", "CLUSTER_SIZES", "CLUSTERS"), convert = T, remove = T)
  CLUSTER_INFO$NUM_CLUSTERS  <- as.numeric(unlist(lapply(strsplit(CLUSTER_INFO$NUM_CLUSTERS, ", ", fixed=TRUE), function(x) return(x[1]))))
  CLUSTER_INFO$CLUSTER_SIZES <- unlist(lapply(strsplit(CLUSTER_INFO$CLUSTER_SIZES, ",  clusters", fixed=TRUE), function(x) return(x[1])))
  CLUSTER_INFO$CLUSTER_SIZES <- str_replace(CLUSTER_INFO$CLUSTER_SIZES, ", ", "_")
  CLUSTER_INFO$CLUSTERS      <- str_replace(CLUSTER_INFO$CLUSTERS, ", ", "/")
  
  N_MONOMERS_INFO <- SLURM %>% filter(N_MON == T) %>% select(V1)
  N_MONOMERS_INFO <- separate(data = N_MONOMERS_INFO, col = V1,   sep = c("chain"), into = c("N_MONOMERS", "TRASH"), convert = T, remove = T) %>% unique()
  N_MONOMERS_INFO$N_MONOMERS <- unlist(lapply(strsplit(N_MONOMERS_INFO$N_MONOMERS, ": ", fixed=TRUE), function(x) return(x[2]))) %>% as.numeric()
  N_MONOMERS_INFO <- N_MONOMERS_INFO %>% select(N_MONOMERS)
  
  EXTRACT       <- bind_cols(MODEL_INFO, CLUSTER_INFO)
  EXTRACT       <- EXTRACT %>% 
    mutate(FILE_MODEL = paste(FILE, MODEL, sep = "_"), N_MONOMERS = N_MONOMERS_INFO$N_MONOMERS) %>% 
    select(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT$MODEL <- unlist(lapply(strsplit(EXTRACT$MODEL, "_multimer", fixed=TRUE), function(x) return(x[1])))
  EXTRACT       <- unique(EXTRACT)
  
  write.table(EXTRACT, file = paste0(OUT ,"_fromSLURM.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}

################################################################################
################################################################################
################################################################################

join_timepoints <- function(DF1, DF2, DF3, OUT) {
  # Join the initial dataframes derived from protein X pull down
  DF <- full_join(DF1, DF2)
  DF <- full_join(DF, DF3)
}

################################################################################
################################################################################
################################################################################

ELISA_Fx <- function(Input_Directory, Output_Directory) {
  # Initialize as an empty data frame
  All_plates_data = data.frame()
  
  # Get the list of subdirectories matching the pattern "Plate_"
  subdirs <- list.files(Input_Directory, recursive = FALSE, full.names = TRUE, pattern = "^Plate_\\d+_\\d{8}$")
  
  if (length(subdirs) > 0) {
    print("Plates exist!")
    for (input_plate_dir in subdirs) {
      
      Input_plate <- input_plate_dir
      
      #Reading Plate Treatment 
      MEASUREMENTS <- fread(paste0(Input_plate, "/MEASUREMENTS.csv"), header = F)
      CELL_LINES   <- fread(paste0(Input_plate, "/CELL_LINES.csv"), header = F)
      CONDITIONS   <- fread(paste0(Input_plate, "/CONDITIONS.csv"), header = F)
      STIM_DAYS    <- fread(paste0(Input_plate, "/STIM_DAYS.csv"), header = F)
      
      #Converting tables into vector for to make a single table
      MEASUREMENTS <- as.vector(as.matrix(MEASUREMENTS))
      CELL_LINES   <- as.vector(as.matrix(CELL_LINES))
      CONDITIONS   <- as.vector(as.matrix(CONDITIONS))
      STIM_DAYS    <- as.vector(as.matrix(STIM_DAYS))
      
      #Creating Table containing all plate Information
      Plate <- NULL
      Plate$MEASUREMENTS <- MEASUREMENTS
      Plate$CELL_LINES <- CELL_LINES
      Plate$CONDITIONS <- CONDITIONS
      Plate$STIM_DAYS <- STIM_DAYS
      
      rm(MEASUREMENTS, CELL_LINES, CONDITIONS, STIM_DAYS)
      
      Plate <- Plate %>% as.data.table()
      
      #Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINES != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate %>% 
        filter(CONDITIONS == "CALIBRATION") %>% 
        group_by(CELL_LINES) %>% 
        summarise(MEASUREMENTS_mean = mean(MEASUREMENTS)) %>%  #, 
        #MEASUREMENTS_median = median(MEASUREMENTS)) %>%  
        mutate(CELL_LINES = as.numeric(CELL_LINES),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}"))) %>% 
        arrange(CELL_LINES)
      
      # Plate_Date = Plate_Standards$Date %>% unique()
      # 
      # if (Plate_Date > "2023-01-01" & Plate_Date <= "2023-05-09") {
      #   Fit <- lm(CELL_LINES ~ MEASUREMENTS_mean*5 - 1, data = Plate_Standards) #linear model of the Standard curve. -1 omits the intercept
      # } else if (Plate_Date < "2023-01-01" | Plate_Date > "2023-05-09") {
      #   Fit <- lm(CELL_LINES ~ MEASUREMENTS_mean - 1, data = Plate_Standards) #linear model of the Standard curve. -1 omits the intercept
      # }
      
      Fit <- lm(CELL_LINES ~ MEASUREMENTS_mean - 1, data = Plate_Standards) #linear model of the Standard curve. -1 omits the intercept
      
      R       <- summary(Fit)$r.squared
      Rsquare <- signif(R, digits = 4)
      
      print(paste0("IL2-Amount = slope*Intensity"))
      print(paste0("IL2-Amount = ", Fit$coefficients[1],"*Intensity"))
      
      Plate_Standards <- Plate_Standards %>% 
        mutate(Fit_Test = (Fit$coefficients[1]*MEASUREMENTS_mean))
      
      # Plotting Standard Curve
      p <- ggplot(data = Plate_Standards) +
        geom_point(aes(x = MEASUREMENTS_mean, y = CELL_LINES), size = 5) +
        geom_line(aes(x = MEASUREMENTS_mean, y = Fit_Test), linetype = "dashed") +
        annotate('text', x = 0.15, y = 700, label = paste0("R^2 = ", Rsquare), size = 10) +
        annotate('text', 
                 x = max(Plate_Standards$MEASUREMENTS_mean) - (0.25 * max(Plate_Standards$MEASUREMENTS_mean)),
                 y = 150, label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        labs(x = "Measured Values",
             y = "IL-Concentration (pg/mL)") +
        ggtitle(label = paste0(basename(Input_plate)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 20))
      
      # Saving the plot
      Save_Name <- paste0(basename(Input_plate), "_Standard_Curve.pdf")
      Save_Name <- file.path(Output_Directory, Save_Name)
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      
      # Fitting Data To Standarad Curve ----------------------------------------
      Plate <- Plate %>% 
        filter(CONDITIONS != "CALIBRATION") %>% 
        mutate(Plate = as.numeric(gsub("^Plate_(\\d+)_\\d{8}$", "\\1", basename(input_plate_dir))),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}")),
               Dilution_Factor = case_when(#Date <= "2023-01-01" ~ DILUTION_FACTOR_2,
                 Date <= "2023-05-09" ~ DILUTION_FACTOR_5,
                 Date >  "2023-05-09" ~ DILUTION_FACTOR_10),
               MEASUREMENTS = as.numeric(MEASUREMENTS),
               IL2_concentration = (Fit$coefficients[1]*MEASUREMENTS),
               IL2_concentration_DILUTION_FACTOR = IL2_concentration*Dilution_Factor)
      
      All_plates_data <- rbind(All_plates_data, Plate)
    }
  } else {
    print("No plates found!")
  }
  return(All_plates_data)
}

################################################################################
################################################################################
################################################################################

correct_colnames <- function(DF) {
  colnames(DF) <- gsub(" ", ".", colnames(DF))
}

################################################################################
################################################################################
################################################################################

add2Summary <- function(NEW_DF, EXISTING_DF) {
  # read dataframes
  existing_df <- fread(EXISTING_DF)
  new_df      <- fread(NEW_DF)
  
  # Combine the two data frames
  combined_df <- rbind(existing_df, new_df, fill = T)
  
  unq <- combined_df %>% distinct(RECYCLE, .keep_all = TRUE)
  
  # Write to the existing CSV
  fwrite(combined_df, EXISTING_DF, row.names = F, sep = ",")
}


################################################################################
################################################################################
################################################################################

#' Extracts and Processes AlphaFold Results
#'
#' This function extracts data from JSON files, processes it, and appends it to summary CSV files.
#' It's designed to work specifically with AlphaFold results.
#'
#' @param LOC A string representing the location of the target data folder.
#' @param MAIN A string representing the main directory where the summaries are stored. Defaults to "~/Desktop/SUMMARIES/LOC/".
#' @param SUMMARY_FOLDER A string representing the directory where summary files are saved. Defaults to "~/Desktop/SUMMARIES/".
#' @param ADD_2_EXISTING_DF A logical value indicating whether to append the extracted data to an existing dataframe. Defaults to FALSE.
#' @param EXISTING_DF A string representing the path of the existing dataframe to which the data should be appended, if ADD_2_EXISTING_DF is TRUE.
#'
#' @return None. This function primarily performs side effects, such as modifying files in place.
#' 
#' @examples
#' \dontrun{
#' run_extraction(LOC = "my_location")
#' run_extraction(LOC = "my_location", ADD_2_EXISTING_DF = TRUE, EXISTING_DF = "~/path_to_existing_df.csv")
#' }
#'

run_extraction <- function(LOC, MAIN = NULL, SUMMARY_FOLDER = NULL, ADD_2_EXISTING_DF = F, EXISTING_DF = NULL) {
  pacman::p_load(dplyr, tidyr, stringr, fs, jsonlite, purrr, utils, data.table)
  if (is.null(MAIN)) {
    MAIN <- paste0("~/Desktop/SUMMARIES/", LOC, "/")
  }
  
  if (is.null(SUMMARY_FOLDER)) {
    SUMMARY_FOLDER <- "~/Desktop/SUMMARIES/"
  }
  
  # Check for the summary directory and the RECYCLES folder inside
  if (!file.exists(file.path(SUMMARY_FOLDER, "/RECYCLES"))) {
    # If it doesn't exist, create it
    dir.create(file.path(SUMMARY_FOLDER, "/RECYCLES"), recursive = TRUE)
  }
  
  if (isTRUE(ADD_2_EXISTING_DF) & is.null(EXISTING_DF)) {
    print("If you wish to add the extracted data to an existing DF, please provide the path for EXISTING_DF.")
  }
  
  summaryWithRecycles = paste0(SUMMARY_FOLDER, "/RECYCLES/", LOC, "_summaryWithRecycles.csv")
  summary             = paste0(SUMMARY_FOLDER, LOC, ".csv")
  if (file.exists(summaryWithRecycles))  {file.remove(summaryWithRecycles)}
  if (file.exists(summary))              {file.remove(summary)}
  
  COUNTER <- 0
  FILES   <- list.files(MAIN, pattern = "_")
  LEN     <- as.numeric(length(FILES))
  
  # if you wish to extract a single file, unhash below and provide file name
  #FILES = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"
  #FILE = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"
  
  # Loop through each FILE in FILES
  for (FILE in FILES) {
    
    tryCatch({
      
      FOLDER = paste0(MAIN, FILE, "/")
      dir_create(FOLDER, "CSV")
      csv_file <- paste0(FOLDER, "CSV/", FILE, "_withRecycles.csv")
      csv_file2 <- paste0(FOLDER, "CSV/", FILE, ".csv")
      
      if (file.exists(csv_file))  {file.remove(csv_file)}
      if (file.exists(csv_file2)) {file.remove(csv_file2)}
      
      json_folder <- ifelse(dir.exists(file.path(FOLDER, "JSON")), yes = file.path(FOLDER, "JSON"), file.path(FOLDER))
      maxJSON <- list.files(json_folder, pattern = ".json")
      
      if (is_empty(maxJSON)) {
        next
      }
      
      json_files <- dir_ls(json_folder, regexp = "\\.json$", recurse = TRUE)
      
      for (file in json_files) {
        json <- readLines(file)
        json <- str_replace_all(json, "Infinity", "9999")
        writeLines(json, file)
      }
      
      if (length(maxJSON) > 0) {
        for (i in 1:length(maxJSON)) {
          JSON <- file.path(json_folder, maxJSON[i])
          OUT <- paste(FOLDER, "CSV", FILE, sep = "/")
          jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
        }
      }
      
      CSV_FILES <- c(csv_file, csv_file2)
      
      for (CSV_FILE in CSV_FILES) {
        JSON_EXTRACT <- data.table::fread(CSV_FILE, header = FALSE) %>%
          dplyr::mutate(ORIGIN = FILE) %>%
          dplyr::rename(FILE = V1, MODEL = V2, RECYCLE = V3, TOL = V4, pLDDT = V5, pTM = V6, piTM = V7, iScore = V8, iRes = V9, iCnt = V10, FILE_MODEL = V11, NUM_CLUSTERS = V12, N_MONOMERS = V13)
        
        JE <- JSON_EXTRACT %>%
          dplyr::mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_RECYCLE_", RECYCLE), RANK = NA) %>%
          dplyr::distinct(FILE_RECYCLE, .keep_all = TRUE) %>%
          dplyr::group_by(FILE) %>%
          dplyr::mutate(RANK = frank(desc(iScore), ties.method = "min"))
        
        data.table::fwrite(JE, CSV_FILE, row.names = FALSE)
        data.table::fwrite(JE, summaryWithRecycles, row.names = FALSE, append = T)
      }
      
      # Increment counter after each successful round
      COUNTER <- COUNTER + 1
      
      # Print current count out of total
      cat("\r", paste(COUNTER, "/", LEN, " done"), fill = F)
      
    }, error = function(e) {
      # Handle error
      cat("An error occurred for FILE:", FILE, "\n")
      cat("Error message:", e$message, "\n")
    })
    
  }
  
  bigboy      <- fread(summaryWithRecycles, fill = T) %>% distinct(RECYCLE, .keep_all = T)
  bigboy$DATE <- ifelse(test = str_detect(bigboy$RECYCLE, pattern = "multimer"), 
                        yes = sapply(strsplit(bigboy$RECYCLE, "_"), function(x) x[6]),
                        no = sapply(strsplit(bigboy$RECYCLE, "_"), function(x) x[5]))
  bigboy$DATE <- lubridate::as_date(bigboy$DATE)
  data.table::fwrite(bigboy, summaryWithRecycles, row.names = FALSE, append = F)
  
  if (isTRUE(ADD_2_EXISTING_DF) &  !is.null(EXISTING_DF)) {
    add2Summary(NEW_DF = summaryWithRecycles, EXISTING_DF = EXISTING_DF)
  }
  
  smolboy <- bigboy %>% filter(!str_detect(RECYCLE, "_recycled_")) %>% distinct(RECYCLE, .keep_all = T)
  data.table::fwrite(smolboy, summary, row.names = FALSE, append = F)
  
  if (isTRUE(ADD_2_EXISTING_DF) &  !is.null(EXISTING_DF)) {
    add2Summary(NEW_DF = summary, EXISTING_DF = EXISTING_DF)
  }
}

################################################################################
################################################################################
################################################################################

#' Plot AlphaFold Results
#'
#' This function takes AlphaFold results from a specified location and plots them.
#' It categorizes the `iScore` into different confidence levels and allows filtering based on a given pattern.
#'
#' @param LOC A character string specifying the location of the data. This will be appended to the SUMMARY_FOLDER path to generate the full path.
#' @param SUMMARY_FOLDER The folder where summary results are stored. By default, it is set to "~Desktop/SUMMARIES/".
#' @param xlab A character string for the x-axis label. Default is "iScore".
#' @param ylab A character string for the y-axis label. Default is "piTM".
#' @param plot_title A character string for the plot title. If NULL (the default), a default title is generated based on the LOC value.
#' @param pattern A character string. If provided, the function will filter the dataset to include only rows where the FILE column matches this pattern.
#' @param best_only A logical value. If TRUE, only the best results (as determined by iScore) will be plotted.
#'
#' @return A ggplot object showing AlphaFold results.
#' 
#' @examples
#' plot_alphafold_results(LOC = "ExampleLocation")
#' plot_alphafold_results(LOC = "ExampleLocation", pattern = "specific_pattern", best_only = TRUE)
#'
#' @importFrom ggplot2 ggplot annotate geom_abline geom_point labs expand_limits
#' @importFrom dplyr mutate group_by ungroup filter select arrange
#' @importFrom data.table fread
#' @importFrom stringr str_detect
#' @importFrom ggalt geom_encircle
#' @importFrom ggrepel geom_label_repel
#' 

plot_alphafold_results <- function(LOC, SUMMARY_FOLDER = NULL, xlab = "iScore", ylab = "piTM", 
                                   plot_title = NULL, pattern = NULL, best_only = FALSE) {
  pacman::p_load(ggplot2, dplyr, data.table, ggrepel, ggalt, stringr)
  
  if (is.null(SUMMARY_FOLDER)) {
    SUMMARY_FOLDER <- "~/Desktop/SUMMARIES/"
  }
  
  file_path <- paste0(SUMMARY_FOLDER, LOC, ".csv")
  
  if (!file.exists(file_path)) {
    stop("The specified file does not exist.")
  }
  
  DF <- data.table::fread(file_path)
  
  if (!all(c("iScore", "piTM", "FILE", "RECYCLE") %in% names(DF))) {
    stop("The data does not have the expected columns.")
  }
  
  DF <- DF %>%
    mutate(Confidence = case_when(iScore < 0.4 ~ "Low",
                                  iScore >= 0.4 & iScore < 0.5 ~ "Medium",
                                  iScore >= 0.5 & iScore < 0.7 ~ "High",
                                  TRUE ~ "Very High"),
           Clash_Indicator = iRes/iCnt)
  
  # Adjust the levels of the Confidence factor variable
  DF$Confidence <- factor(DF$Confidence, levels = c("Low", "Medium", "High", "Very High"))
  
  if (!is.null(pattern)) {  
    DF <- DF %>% filter(stringr::str_detect(FILE, fixed(pattern)))
    
    # ensure that the DF isn't empty
    if (nrow(DF) == 0) {
      stop("The file you want to filter for was not found. Is it spelled correctly?")
    }
  }
  
  
  max_iScore <- DF %>%
    group_by(as.factor(FILE)) %>%
    slice_max(order_by = iScore, n = 1) %>%
    ungroup()
  
  if (isTRUE(best_only)) {  
    DF <- max_iScore
  }
  
  mean_labeling <- DF %>%
    group_by(FILE) %>%
    mutate(mean_iScore = mean(iScore),
           mean_piTM   = mean(piTM)) %>%
    distinct(FILE, .keep_all = T) %>%
    ungroup()
  
  
  print(max_iScore %>%
          dplyr::select(FILE, RECYCLE, iScore, piTM) %>% 
          dplyr::arrange(desc(iScore)))
  
  if (is.null(plot_title)) {
    plot_title = paste0("AlphaFold Results for ", LOC)
  }
  
  # Base plot
  plot <- ggplot(DF) +
    annotate("rect", xmin = 0, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.3) +
    annotate("rect", xmin = 0.4, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "gray40", alpha = 0.3) +
    annotate("rect", xmin = 0.5, xmax = 0.7, ymin = -Inf, ymax = Inf, fill = "cornflowerblue", alpha = 0.3) +
    annotate("rect", xmin = 0.7, xmax = 1, ymin = -Inf, ymax = Inf, fill = "lightgreen", alpha = 0.3) +
    annotate("text", x = 0.41, y = 0.99, col = "black", label = "medium confidence", angle = 90, hjust = 1) +
    annotate("text", x = 0.51, y = 0.99, col = "blue", label = "high confidence", angle = 90, hjust = 1) +
    annotate("text", x = 0.71, y = 0.99, col = "darkgreen", label = "very high confidence", angle = 90, hjust = 1) +
    geom_abline(col = "gray")
  
  # Add encircle only if best_only is FALSE
  if (!best_only & length(unique(DF$FILE)) <= 10) {
    plot <- plot + geom_encircle(aes(iScore, piTM, fill = FILE), alpha = 0.1)
  }
  
  # Check number of unique FILE names and continue adding layers based on condition
  if (length(unique(DF$FILE)) > 10) {
    plot <- plot +
      geom_point(aes(iScore, piTM, col = Confidence), size = 3) +
      geom_point(data = max_iScore, aes(iScore, piTM, col = Confidence), size = 5) +
      scale_color_manual(name = "Confidence",
                         values = c("Low" = "gray80",
                                    "Medium" = "gray40",
                                    "High" = "cornflowerblue",
                                    "Very High" = "lightgreen"))
  } else {
    plot <- plot +
      geom_point(aes(iScore, piTM, color = FILE)) +
      geom_point(data = max_iScore, aes(iScore, piTM, color = FILE), size = 4) +
      geom_label_repel(data = mean_labeling, aes(mean_iScore, mean_piTM, label = FILE, color = FILE))
  }
  
  # Add remaining layers and return
  plot <- plot + 
    expand_limits(x=c(0,1), y=c(0,1)) +
    labs(x = xlab, y = ylab, title = plot_title)
  
  return(plot)
}

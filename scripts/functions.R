# FUNCTIONS NEW VERSION 20231125

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
  # FILE_MODEL <- paste(FILE, "MODEL", MODEL, sep = "_")
  NUM_CLUSTERS <- json$clusters[[iScore$RECYCLE[1]]]$num_clusters
  N_MONOMERS <- length(json$chains)
  
  EXTRACT <- cbind(FILE, #MODEL, 
                   TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, #FILE_MODEL, 
                   NUM_CLUSTERS, N_MONOMERS)
  EXTRACT <- EXTRACT[, !duplicated(colnames(EXTRACT))]
  
  write.table(EXTRACT, file = paste0(OUT,"_withRecycles.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
  
  EXTRACT_noRecycle <- EXTRACT %>% filter(!str_detect(RECYCLE, "_recycled_"))
  write.table(EXTRACT_noRecycle, file = paste0(OUT,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
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
  # Input_Directory = "~/Desktop/ALL"
  # Output_Directory = "~/Desktop/ALL"
  
  # Get the list of subdirectories matching the pattern "Plate_"
  subdirs <- list.files(Input_Directory, recursive = FALSE, full.names = TRUE, pattern = "Plate_\\d+_\\d{8}$")
  # input_plate_dir = "~/Desktop/ALL/Plate_1_20220609"
  
  if (length(subdirs) > 0) {
    print("Plates exist!")
    for (input_plate_dir in subdirs) {
      
      Input_plate <- input_plate_dir
      
      #Reading Plate Treatment 
      MEASUREMENTS <- fread(paste0(Input_plate, "/MEASUREMENTS.csv"), header = F)
      CELL_LINES   <- fread(paste0(Input_plate, "/CELL_LINES.csv"),   header = F)
      CONDITIONS   <- fread(paste0(Input_plate, "/CONDITIONS.csv"),   header = F)
      STIM_DAYS    <- fread(paste0(Input_plate, "/STIM_DAYS.csv"),    header = F)
      DILUTIONS    <- fread(paste0(Input_plate, "/DILUTIONS.csv"),    header = F)
      
      #Converting tables into vector for to make a single table
      MEASUREMENT <- as.vector(as.matrix(MEASUREMENTS))
      CELL_LINE   <- as.vector(as.matrix(CELL_LINES))
      CONDITION   <- as.vector(as.matrix(CONDITIONS))
      STIM_DAY    <- as.vector(as.matrix(STIM_DAYS))
      DILUTION    <- as.vector(as.matrix(DILUTIONS))
      
      #Creating Table containing all plate Information
      Plate <- NULL
      Plate$MEASUREMENT <- MEASUREMENT
      Plate$CELL_LINE   <- CELL_LINE
      Plate$CONDITION   <- CONDITION
      Plate$STIM_DAY    <- STIM_DAY
      
      if (exists("DILUTION")) {
        Plate$DILUTION <- DILUTION
      } else {
        Plate$DILUTION <- ifelse(Plate$CONDITION != "CALIBRATION", 5, NA_real_)
      }
      Plate$DILUTION <- as.numeric(Plate$DILUTION)
      
      rm(MEASUREMENT, CELL_LINE, CONDITION, STIM_DAY, DILUTION)
      
      Plate <- Plate %>% as.data.table()
      
      #Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINE != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate[Plate$CONDITION == "CALIBRATION"]
      Plate_Standards$CELL_LINE <- as.numeric(Plate_Standards$CELL_LINE)
      
      Plate_Standards <- Plate_Standards %>%
        group_by(CELL_LINE) %>% 
        summarise(MEASUREMENT_mean = mean(MEASUREMENT)) %>%
        mutate(CELL_LINE = as.numeric(CELL_LINE),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}"))) %>% 
        arrange(CELL_LINE)
      
      # linear model of the Standard curve. -1 omits the intercept
      # We will only use standard curve values of 1 and below (machine is optimized to measure absorptions between 0 and 1.1)
      Fit <- lm(CELL_LINE ~ MEASUREMENT_mean - 1, data = Plate_Standards[Plate_Standards$MEASUREMENT_mean <= 1.1, ])
      
      R       <- summary(Fit)$r.squared
      Rsquare <- signif(R, digits = 4)
      
      print(paste0("Secretion = slope*Intensity"))
      print(paste0("Secretion = ", Fit$coefficients[1],"*Intensity"))
      Plate_Standards <- Plate_Standards %>% mutate(Fit_Test = (Fit$coefficients[1]*MEASUREMENT_mean))
      
      # Plotting Standard Curve
      p <- ggplot(data = Plate_Standards) +
        geom_point(aes(x = MEASUREMENT_mean, y = CELL_LINE, col = MEASUREMENT_mean <= 1.5), size = 5) +
        geom_line(aes(x = MEASUREMENT_mean, y = Fit_Test), linetype = "dashed") +
        annotate('text', x = 0.15, y = 700, label = paste0("R^2 = ", Rsquare), size = 10) +
        annotate('text', 
                 x = max(Plate_Standards$MEASUREMENT_mean) - (0.25 * max(Plate_Standards$MEASUREMENT_mean)),
                 y = 150, label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        labs(x = "Measured Values",
             y = "IL-Concentration (pg/mL)") +
        ggtitle(label = paste0(basename(Input_plate)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 20)) +
        theme(legend.position = "none")
      
      # Saving the plot
      Save_Name <- paste0(basename(Input_plate), "_Standard_Curve.pdf")
      Save_Name <- file.path(Output_Directory, Save_Name)
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      
      # Fitting Data To Standard Curve ----------------------------------------
      Plate <- Plate %>% 
        filter(CONDITION != "CALIBRATION") %>% 
        mutate(Plate = as.numeric(gsub("(DR_)?Plate_(\\d+)_\\d{8}$", "\\2", basename(input_plate_dir))),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}")),
               MEASUREMENT = as.numeric(MEASUREMENT),
               
               # METHOD I: Correct negative values to zero, multiply measurements by dilution factor and THEN extrapolate values based on SC
               # ########  20231206 @Fakun: We should adjust for the dilution factor BEFORE extrapolating values based on the Standard Curve..
               # ########
               # MEASUREMENT_DIL_ADJ = (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT)*DILUTION),
               # Concentration = (Fit$coefficients[1]*MEASUREMENT_DIL_ADJ),
               
               # METHOD II: Extrapolate values based on SC and THEN multiply measurements by dilution factor
               # ########  20231206 @Finn: set machine measurement errors to zero (raw values below zero should be set to zero before extrapolating etc)
               # ########  20231206 @Finn:  We remove the column name *Concentration_DILUTION_FACTOR* from downstream analysis and simply stick to *Concentration*!
               # ########
               # Concentration = (Fit$coefficients[1]*MEASUREMENT),
               # Concentration_DILUTION_FACTOR = Concentration*DILUTION,
               Concentration = (Fit$coefficients[1] * (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT))),
               Concentration = Concentration * DILUTION,
               
               Is_Dose_Response = ifelse(str_detect(basename(input_plate_dir), "^DR_"), TRUE, FALSE)
        )
      
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
    MAIN <- paste0("~/Desktop/", LOC, "/")
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
          dplyr::rename(FILE = V1, 
                        # MODEL = V2, 
                        RECYCLE = V2, 
                        TOL = V3, 
                        pLDDT = V4, 
                        pTM = V5, 
                        piTM = V6, 
                        iScore = V7, 
                        iRes = V8, 
                        iCnt = V9, 
                        # FILE_MODEL = V11,
                        NUM_CLUSTERS = V10, 
                        N_MONOMERS = V11)
        
        JE <- JSON_EXTRACT %>%
          # dplyr::mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_RECYCLE_", RECYCLE), RANK = NA) %>%
          dplyr::mutate(FILE_RECYCLE = paste0(FILE, "_RECYCLE_", RECYCLE), 
                        # RANK = NA
          ) %>%
          dplyr::distinct(FILE_RECYCLE, .keep_all = TRUE)# %>%
        # dplyr::group_by(FILE) %>%
        # dplyr::mutate(RANK = frank(desc(iScore), ties.method = "min"))
        
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
                                   plot_title = NULL, pattern = NULL, best_only = FALSE,
                                   plot_interactive = FALSE) {
  pacman::p_load(ggplot2, dplyr, data.table, ggrepel, ggalt, stringr, plotly)
  
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
  
  
  if (plot_interactive) {
    plotly_plot <- ggplot(DF, aes(label = FILE)) +
      geom_abline(col = "gray")
    # Add encircle only if best_only is FALSE
    if (!best_only & length(unique(DF$FILE)) <= 10) {
      plotly_plot <- plotly_plot + geom_encircle(aes(iScore, piTM, fill = FILE), alpha = 0.1)
    }
    # Check number of unique FILE names and continue adding layers based on condition
    if (length(unique(DF$FILE)) > 10) {
      plotly_plot <- plotly_plot +
        geom_point(aes(iScore, piTM, col = Confidence), size = 3) +
        geom_point(data = max_iScore, aes(iScore, piTM, col = Confidence), size = 5) +
        scale_color_manual(name = "Confidence",
                           values = c("Low" = "gray80",
                                      "Medium" = "gray40",
                                      "High" = "cornflowerblue",
                                      "Very High" = "lightgreen"))
    } else {
      plotly_plot <- plotly_plot +
        geom_point(aes(iScore, piTM, color = FILE)) +
        geom_point(data = max_iScore, aes(iScore, piTM, color = FILE), size = 4) +
        geom_label_repel(data = mean_labeling, aes(mean_iScore, mean_piTM, label = FILE, color = FILE))
    }
    # Add remaining layers and return
    plotly_plot <- plotly_plot + 
      expand_limits(x=c(0,1), y=c(0,1)) +
      labs(x = xlab, y = ylab, title = plot_title)
    return(ggplotly(plotly_plot))
    
  } else return(plot)
  
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Define the function to save plots with specified dimensions
save_plots <- function(folder_name, plots, path = save_to, width = 12, height = 8) {
  # Create the full path for the folder
  folder_path <- file.path(save_to, folder_name)
  
  # Check if the folder exists, and if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # Loop through the plots and save them
  for (i in seq_along(plots)) {
    file_base <- file.path(folder_path, paste0(folder_name, "_plot_", i))
    
    # Save as SVG
    ggsave(paste0(file_base, ".svg"), plot = plots[[i]], device = "svg", width = width, height = height)
    
    # Save as PNG
    ggsave(paste0(file_base, ".png"), plot = plots[[i]], device = "png", width = width, height = height)
  }
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Vectorized function for pattern matching
matches_any_pattern_vec <- Vectorize(function(string, patterns) {
  any(sapply(patterns, function(pattern) grepl(pattern, string)))
}, vectorize.args = "string")

# FILTER_VALUES = c("MyD88-GFP-synTRAF6-BD-1x", "MyD88-GFP-synTRAF6-BD-3x", "MyD88-GFP-synTRAF6-BD-5x")
# FILTER_TYPE = "COHORT"        # Should be "COHORT" or "DAY"
# COLOR = "salmon"
# SEED = 600
# plot_results = TRUE
# plot_pval = TRUE
# run_anova = FALSE
# plot_faceted_by_date = FALSE
# save_to = "~/Desktop"
# POSITIVE_CTRL = c("WT")
# NEGATIVE_CTRL = 0
# DATA = All_plates_data

filter_data <- function(DATA, FILTER_VALUES, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL) {
  filter_pattern_func <- if (FILTER_TYPE %in% "COHORT") {
    function(value, df) matches_any_pattern_vec(df$CELL_LINE, value) | matches_any_pattern_vec(df$CL_NAME_ON_PLOT, value)
  } else if (FILTER_TYPE %in% "DAY") {
    function(value, df) matches_any_pattern_vec(df$Date, value)
  } else {
    stop("Invalid FILTER_TYPE. Must be either 'COHORT' or 'DAY'.")
  }
  
  # Generate lists for plates, dates, and stim_days based on the filter type
  plates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Plate) %>%
      pull(Plate)
  })
  
  dates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Date) %>%
      pull(Date)
  })
  
  stim_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(STIM_DAY) %>%
      pull(STIM_DAY)
  })
  
  condition_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(CONDITION) %>%
      pull(CONDITION)
  })
  
  # Adding names to the list elements
  names(plates_list)    <- FILTER_VALUES
  names(dates_list)     <- FILTER_VALUES
  names(stim_list)      <- FILTER_VALUES
  names(condition_list) <- FILTER_VALUES
  
  # Create subset based on the lists
  subset_list <- lapply(names(plates_list), function(value) {
    plates <- plates_list[[value]]
    dates  <- dates_list[[value]]
    stim   <- stim_list[[value]]
    cond   <- condition_list[[value]]
    
    DATA %>%
      filter(Date %in% dates,
             Plate %in% plates,
             STIM_DAY %in% stim,
             CONDITION %in% cond,
             filter_pattern_func(c(POSITIVE_CTRL, NEGATIVE_CTRL, value), DATA))
  })
  
  # Combine the subsets and return
  FILTERED_SUBSET <- bind_rows(subset_list)
  return(FILTERED_SUBSET)
}

# DATA = FILTERED_SUBSET
calculate_baseline_and_control <- function(DATA, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL) {
  # Determine the group_by and filter parameters based on FILTER_TYPE
  if (FILTER_TYPE == "COHORT") {
    group_vars <- c("Date", "STIM_DAY")
    filter_var <- "CELL_LINE"
  } else if (FILTER_TYPE == "DAY") {
    group_vars <- c("Date", "STIM_DAY", "Plate")
    filter_var <- "CELL_LINE"  # Assuming the correct variable for DAY filter
  } else {
    stop("Invalid FILTER_TYPE. Must be either 'cohort' or 'day'.")
  }
  
  # Debugging information
  # print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  # print(paste("Filtering using variable:", filter_var))
  
  # Calculate the baseline control value
  if (is.character(NEGATIVE_CTRL)) {
    baseline <- DATA %>%
      group_by(!!!syms(group_vars)) %>%
      filter(matches_any_pattern_vec(!!sym(filter_var), NEGATIVE_CTRL), CONDITION == "UNSTIM") %>% 
      summarise(baseline_control_value = mean(Concentration))
    
    # Join the calculated values with the dataset
    data <- left_join(DATA, baseline, by = group_vars) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value,
                                                               TRUE ~ Concentration))
  } else {
    # baseline <- DATA %>%
    #   group_by(!!!syms(group_vars)) %>%
    #   summarise(baseline_control_value = min(Concentration))
    baseline <- data %>%
      group_by(Date) %>%
      filter(Date == FILTER_VALUE) %>%
      summarise(baseline_control_value = min(Concentration))
    
    # Join the calculated values with the dataset
    data <- left_join(DATA, baseline) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value,
                                                               TRUE ~ Concentration))
  }
  
  
  if (is.character(POSITIVE_CTRL)) {
    control_mean_per_day <- data %>%
      filter(matches_any_pattern_vec(!!sym(filter_var), POSITIVE_CTRL), CONDITION == "STIM") %>% 
      group_by(!!!syms(group_vars)) %>%
      summarise(control_mean_MEASUREMENT = case_when(mean(Concentration_REDUCED) > 0 ~ mean(Concentration_REDUCED), TRUE ~ -Inf))
    
    # Join the calculated control means
    data <- left_join(data, control_mean_per_day, by = group_vars)
    
    # Perform normalization
    DATA_NORMALIZED <- data %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_mean_MEASUREMENT < 0 ~ 0,
                                                                  TRUE ~ Concentration_REDUCED / control_mean_MEASUREMENT),
             triplicate_mean_per_day = mean(Concentration_NORMALIZED)) %>%
      ungroup()
    
  } else {
    # Calculate control mean using an alternative method, e.g., mean value per group
    # control_mean_per_day <- data %>%
    #   group_by(!!!syms(group_vars)) %>%
    #   summarise(control_mean_MEASUREMENT = mean(Concentration_REDUCED))
    control_mean_per_day <- data %>%
      group_by(STIM_DAY) %>%
      reframe(Concentration_REDUCED = Concentration_REDUCED) %>%
      ungroup() %>%
      group_by(STIM_DAY) %>%
      reframe(control_MEASUREMENT = max(Concentration_REDUCED),
              STIM_DAY = STIM_DAY) %>%
      unique()
    
    # Join the calculated control means
    data <- left_join(data, control_mean_per_day)
    
    # Perform normalization
    DATA_NORMALIZED <- data %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_MEASUREMENT < 0 ~ 0,
                                                                  TRUE ~ Concentration_REDUCED / control_MEASUREMENT),
             triplicate_mean_per_day = mean(Concentration_NORMALIZED)) %>%
      ungroup()
  }
  
  return(DATA_NORMALIZED)
}

# DATA = DATA_NORMALIZED
perform_statistical_analysis <- function(DATA, GROUP_BY_COLUMN) {
  # Internal function for pairwise t-test
  pairwise_ttest <- function(DATA, return_annotation = FALSE) {
    p_values <- pairwise.t.test(DATA$triplicate_mean_per_day, DATA$CONDITION, p.adjust.method = "none")$p.value
    
    if (return_annotation) {
      p_annotation <- ifelse(p_values < 0.001, '***', 
                             ifelse(p_values < 0.01, '**', 
                                    ifelse(p_values < 0.05, '*', 
                                           'ns')))
      return(p_annotation)
    } else {
      return(formatC(p_values, format = "e", digits = 3))
    }
  }
  
  # Calculate statistical significance using a t-test for each group
  annotations <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), pairwise_ttest)
  p_values    <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), function(DATA) pairwise_ttest(DATA, return_annotation = TRUE))
  
  return(list(annotations = annotations, p_values = p_values))
}
# Example usage
# results <- perform_statistical_analysis(MEANS, "CL_NAME_ON_PLOT")

create_plot <- function(FILTER_VALUES, FILTER_TYPE, MEANS, MOM_SUBSET, STATISTICAL_RESULTS, COLOR = "salmon", 
                        x_label = xlabel, y_label = y_label, plot_title = plot_title,
                        subtitle = paste0("Assay performed for ", FILTER_VALUES, ", based on filter type ", FILTER_TYPE),
                        SEED = 600, plot_pval = T, plot_faceted_by_date = F) {
  plot <- ggplot(MEANS, aes(x = CL_NAME_ON_PLOT)) +
    geom_col(     data = MOM_SUBSET,  aes(y = triplicate_mean_per_day, fill = CONDITION), position = position_dodge(width = 1), alpha = 0.5) +
    geom_point(   data = MEANS,       aes(y = triplicate_mean_per_day, group = CONDITION, shape = STIM_DAY),               position = position_jitterdodge(jitter.height = 0, jitter.width = 1.2, seed = SEED), col = "white", size = 4) +
    geom_point(   data = MEANS,       aes(y = triplicate_mean_per_day, group = CONDITION, shape = STIM_DAY), position = position_jitterdodge(jitter.height = 0, jitter.width = 1.2, seed = SEED), 
                  col = "black",
                  size = 3) +
    geom_errorbar(data = MOM_SUBSET, aes(ymin = triplicate_mean_per_day - triplicate_sd_per_day, 
                                         ymax = triplicate_mean_per_day + triplicate_sd_per_day, group = CONDITION), width = 0.25, position = position_dodge(width = 1)) +
    labs(x = x_label,
         y = y_label) +
    scale_fill_manual(values = c("UNSTIM" = "gray50", "STIM" = COLOR)) +
    ggtitle(plot_title) +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (FILTER_TYPE == "DAY") {
    plot_faceted_by_date = T
    plot_pval = F
  }
  
  if (plot_pval == T & plot_faceted_by_date == F) {
    plot <- plot +
      geom_text(data = data.frame(CL_NAME_ON_PLOT = names(STATISTICAL_RESULTS[["annotations"]]), y = max(MEANS$triplicate_mean_per_day) + 0.1, label = STATISTICAL_RESULTS[["annotations"]]), aes(y = y, label = label), col = "gray40", size = 4) +
      geom_text(data = data.frame(CL_NAME_ON_PLOT = names(STATISTICAL_RESULTS[["p_values"]]), y = max(MEANS$triplicate_mean_per_day) + 0.2, label = STATISTICAL_RESULTS[["p_values"]]), aes(y = y, label = label), size = 7)
  } else if (plot_pval == F & plot_faceted_by_date == T) {
    plot <- plot  + 
      facet_wrap(~Date, scales = "free_x")
  } else {
    plot_pval = F
    plot <- plot  + 
      facet_wrap(~Date, scales = "free_x")
    print("Overriding p-value plotting. Option not available for date-faceted plots.")
  }
  
  return(plot)
}

# I am providing an example cohort of the All_plates_data data frame to debug the functions:
# FILTER_VALUES = c("BDLD_57", "BDLD_6-")  # Example cohort names
# FILTER_TYPE = "COHORT"        # Should be "COHORT" or "DAY"
# save_to = "~/Desktop"
# POSITIVE_CTRL = c("3E10_GFP", "MyD88_GFP")
# NEGATIVE_CTRL = c("204_TRIPLE_KO", "tKO_EL4")

plot_ELISA <- function(FILTER_VALUES, FILTER_TYPE, POSITIVE_CTRL = 0, NEGATIVE_CTRL = 0,
                       COLOR = "salmon", SEED = 600, 
                       plot_pval = TRUE,  plot_results = TRUE, 
                       run_anova = FALSE, plot_faceted_by_date = FALSE, 
                       x_label = "", y_label = "relative IL-2 conc.", plot_title = "IL-2 ELISA",
                       subtitle = paste0("Assay performed for ", FILTER_VALUES, ", based on filter type ", FILTER_TYPE)) {
  
  ALL_COHORT_DATA <- lapply(FILTER_VALUES, function(FILTER_VALUE) {
    # Step 1: Filter & Subset the main data into lists
    COHORT_DATA <- filter_data(All_plates_data, FILTER_VALUE, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL)
    # Step 2 and Step 3
    NORMALIZED_DATA <- calculate_baseline_and_control(COHORT_DATA, FILTER_TYPE, POSITIVE_CTRL, NEGATIVE_CTRL)
    return(NORMALIZED_DATA)
  })
  
  # Step 4: Combine the subset lists into one data frame
  COMBINED_DATA <- bind_rows(ALL_COHORT_DATA)
  # print(COMBINED_DATA)  # Debugging line
  
  # Step 5: Calculate the means per cohort
  MEANS <- COMBINED_DATA %>%
    group_by(CELL_LINE, CONDITION, Date) %>%
    distinct(triplicate_mean_per_day, STIM_DAY, .keep_all = TRUE) %>%
    ungroup()
  # print(MEANS)  # Debugging line
  
  # Step 6: Calculate the mean of means (MOM) per cohort
  MOM_SUBSET <- MEANS %>%
    group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION, Date) %>%
    summarise(triplicate_sd_per_day = sd(triplicate_mean_per_day),
              triplicate_mean_per_day = mean(triplicate_mean_per_day)) %>%
    distinct(.keep_all = TRUE)
  # print(MOM_SUBSET)  # Debugging line
  
  # Step 7: Perform statistics per cohort
  STATISTICAL_RESULTS <- perform_statistical_analysis(MEANS, "CL_NAME_ON_PLOT")
  # print(STATISTICAL_RESULTS)  # Debugging line
  
  # Step 8: Plot the results
  if (plot_results) {
    if (all(POSITIVE_CTRL == 0) & all(NEGATIVE_CTRL == 0)) {
      PLOT <- ggplot(MEANS, aes(x = CL_NAME_ON_PLOT)) +
        geom_col(     data = MOM_SUBSET,  aes(y = triplicate_mean_per_day, fill = CONDITION, group = CONDITION), position = position_dodge(width = 1), alpha = 0.5) +
        labs(x = x_label,
             y = y_label) +
        scale_fill_manual(values = c("UNSTIM" = "gray50", "STIM" = COLOR)) +
        labs(x = x_label, y = y_label) +
        ggtitle(label = plot_title, subtitle = subtitle) +
        theme_cowplot() +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~STIM_DAY)
      
      PLOT
      
    } else {
      PLOT <- create_plot(FILTER_VALUES, FILTER_TYPE, MEANS, MOM_SUBSET, STATISTICAL_RESULTS, COLOR = "salmon", SEED = 600,
                          x_label = x_label, y_label = y_label, plot_title = plot_title, subtitle = subtitle)
    }
    list(PLOT, COMBINED_DATA, MEANS, MOM_SUBSET, MEANS)
  }
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

plot_dose_response_ELISA <- function(DATA, FILTER_VALUES, 
                                     COLOR_PALETTE_KEY = color_palette_key,
                                     SEED = 600, 
                                     plot_pval = TRUE, plot_results = TRUE, 
                                     run_anova = FALSE, plot_faceted_by_date = FALSE, 
                                     x_label = "Stimulant Conc.", 
                                     y_label = "Relative Response", 
                                     plot_title = "Dose Response IL-2 ELISA",
                                     subtitle = paste0("Assay performed on ", FILTER_VALUES),
                                     without_scientific_notation = FALSE) {
  
  group_vars <- c("CELL_LINE", "CL_NAME_ON_PLOT", "CONDITION")
  extra_var <- "STIM_DAY" 
  cat("group_vars used:\n")
  cat(group_vars, "\n")
  cat(paste("extra_var used:", extra_var, "\n"))
  
  ALL_COHORT_DATA <- lapply(FILTER_VALUES, function(FILTER_VALUE) {
    
    cat("Step 1: Filter & Subset the main data into lists\n")
    COHORT_DATA <- All_plates_data %>%
      group_by(CELL_LINE, STIM_DAY, CONDITION, Date, Plate) %>%
      filter(Date == FILTER_VALUE)
    
    cat("Step 2: Calculate baseline values (lowest value per plate)\n")
    baseline <- COHORT_DATA %>%
      group_by(Plate) %>%
      filter(Date == FILTER_VALUE) %>%
      summarise(baseline_control_value = min(Concentration))
    
    cat("Step 3: Subtract the lowest value from all measurements to normalize per plate (Machine Normalization).\n")
    COHORT_DATA <- left_join(COHORT_DATA, baseline) %>%
      mutate(Concentration_REDUCED    = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value,
                                                                  TRUE ~ Concentration))
    cat("Step 4: Baseline-corrected per cell line, per condition, per stimulation and retrieve the max values to normalize for day-to-day differences.\n")
    control_mean_per_day <- COHORT_DATA %>%
      group_by(CELL_LINE, CONDITION, Date, STIM_DAY, Plate) %>%
      reframe(Concentration_REDUCED = Concentration_REDUCED) %>%
      ungroup() %>%
      group_by(STIM_DAY) %>%
      reframe(control_MEASUREMENT = max(Concentration_REDUCED))
    
    cat("Step 5: Normalize the baseline-corrected values by dividing each value by the max value per stimulation day.\n")
    COHORT_DATA <- left_join(COHORT_DATA, control_mean_per_day) %>%
      group_by(CELL_LINE, CONDITION, Date, STIM_DAY) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_MEASUREMENT < 0 ~ 0, TRUE ~ Concentration_REDUCED / control_MEASUREMENT),
             triplicate_mean_per_day = mean(Concentration_NORMALIZED)) %>%
      ungroup()
    
    # Make sure that the conditions are now read as numbers for plotting!
    COHORT_DATA$CONDITION <- as.numeric(as.character(COHORT_DATA$CONDITION))
    
    # Additional processing steps can be added here
    return(COHORT_DATA)
  })
  
  COHORT_DATA <- ALL_COHORT_DATA[[1]]
  
  group_and_summarize <- function(COHORT_DATA, group_vars = c("CELL_LINE", "CL_NAME_ON_PLOT", "CONDITION"), extra_var = "STIM_DAY") {
    
    cat("Step 6: Calculate the normalized means per grouped variables (group_vars + extra_vars)\n")
    MEANS <- COHORT_DATA %>%
      group_by_at(c(group_vars, extra_var)) %>%
      distinct(triplicate_mean_per_day, STIM_DAY, .keep_all = TRUE) %>%
      ungroup()
    # print(MEANS)  # Debugging line
    
    cat("Step 7: Calculate the mean of means (MOM) per cell line and per condition.\n")
    MOM_SUBSET <- COHORT_DATA %>%
      group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION) %>%
      summarise(triplicate_sd_per_day = sd(triplicate_mean_per_day),
                triplicate_mean_per_day = mean(triplicate_mean_per_day)) %>%
      distinct(.keep_all = TRUE)
    # print(MOM_SUBSET)  # Debugging line
    
    return(list(MEANS = MEANS, MOM_SUBSET = MOM_SUBSET))
  }
  
  # Group, Summarize, and Extract Results from group_and_summarize
  result_list <- group_and_summarize(COHORT_DATA, group_vars, extra_var)
  MEANS       <- result_list$MEANS
  MOM_SUBSET  <- result_list$MOM_SUBSET
  
  cat("Step 8: Plot the results.\n")
  
  if (plot_results) {
    create_dose_response_plot <- function(MEANS = MEANS, MOM_SUBSET = MOM_SUBSET, plot_title = plot_title, subtitle = subtitle, x_label = x_label, y_label = y_label, COLOR_PALETTE_KEY = COLOR_PALETTE_KEY) {
      
      # Deriving breaks from unique CONDITION values
      condition_breaks <- sort(unique(MOM_SUBSET$CONDITION))
      
      PLOT <- ggplot(data = MOM_SUBSET, aes(x = CONDITION, y = triplicate_mean_per_day, group = CELL_LINE, color = CELL_LINE)) +
        geom_path(size = 1.5) +
        color_palette(palette = COLOR_PALETTE_KEY) +
        fill_palette(palette  = COLOR_PALETTE_KEY) +
        geom_point(data = MEANS, aes(x = CONDITION, y = triplicate_mean_per_day, fill = CELL_LINE), size = 3, color = "black", shape = 21) +
        geom_errorbar(data = MOM_SUBSET, col = "black", linewidth = .4, width = 0.25,
                      aes(x = CONDITION, y = triplicate_mean_per_day, 
                          ymin = triplicate_mean_per_day - triplicate_sd_per_day, 
                          ymax = triplicate_mean_per_day + triplicate_sd_per_day)) +
        labs(x = x_label, y = y_label) +
        ggtitle(label = plot_title, subtitle = subtitle) +
        theme_cowplot() +
        theme(legend.position = "bottom", axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
      
      if (without_scientific_notation) {
        PLOT <- PLOT +
          scale_x_continuous(trans = "log10", breaks = condition_breaks,
                             # removes the scientific notation
                             labels = scales::label_comma())
      } else {
        PLOT <- PLOT +
          scale_x_continuous(trans = "log10", breaks = condition_breaks)
      }
      
      return(PLOT)
    }
    
    DR_PLOT <- create_dose_response_plot(MEANS = MEANS, MOM_SUBSET = MOM_SUBSET, plot_title = plot_title, subtitle = subtitle, x_label = x_label, y_label = y_label, COLOR_PALETTE_KEY = COLOR_PALETTE_KEY)
    list(DR_PLOT, MOM_SUBSET, MEANS, COHORT_DATA)
  }
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# Plot for a specific day
# FILTER_VALUES = c("2023-11-21")
# FILTER_VALUES = c("2023-11-16")
# FILTER_TYPE = "DAY"
# POSITIVE_CTRL = 0
# NEGATIVE_CTRL = 0
# DATA = All_plates_data

plot_dose_response_ELISA_2 <- function(FILTER_VALUES, FILTER_TYPE = "DAY", DATA = All_plates_data, POSITIVE_CTRL = 0, NEGATIVE_CTRL = 0,
                                       x_label = "X-Axis", y_label = "Y-Axis", 
                                       plot_title = "Dose Response Plot", 
                                       subtitle = paste0("Assay performed on ", FILTER_VALUES),
                                       COLOR_PALETTE_KEY = color_palette_key) {
  cat("This function is optimized for FILTER_TYPE = 'DAY'. Using this as default for subsetting data.\n")
  
  filter_pattern_func <- if (FILTER_TYPE %in% "COHORT") {
    # function(value, df) matches_any_pattern_vec(df$CELL_LINE, value) | matches_any_pattern_vec(df$CL_NAME_ON_PLOT, value)
    print("plot_dose_response_ELISA_2() is currently optizimed for 'DAY' inputs only. Please adjust your input accordingly to use this function.")
  } else if (FILTER_TYPE %in% "DAY") {
    function(value, df) matches_any_pattern_vec(df$Date, value)
    # cat(paste0("Filter Type is 'DAY'. Subsetting for the Date supplied in FILTER_VALUES: ", FILTER_VALUES, "\n"))
  } else {
    stop("Invalid FILTER_TYPE. Must be either 'COHORT' or 'DAY'.")
  }
  
  # Generate lists for plates, dates, and stim_days based on the filter type
  plates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Plate) %>%
      pull(Plate)
  })
  
  dates_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(Date) %>%
      pull(Date)
  })
  
  stim_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(STIM_DAY) %>%
      pull(STIM_DAY)
  })
  
  condition_list <- lapply(FILTER_VALUES, function(value) {
    DATA %>%
      filter(filter_pattern_func(value, DATA)) %>%
      distinct(CONDITION) %>%
      pull(CONDITION)
  })
  
  # Adding names to the list elements
  names(plates_list)    <- FILTER_VALUES
  names(dates_list)     <- FILTER_VALUES
  names(stim_list)      <- FILTER_VALUES
  names(condition_list) <- FILTER_VALUES
  
  # Create subset based on the lists
  subset_list <- lapply(names(plates_list), function(value) {
    plates <- plates_list[[value]]
    dates  <- dates_list[[value]]
    stim   <- stim_list[[value]]
    cond   <- condition_list[[value]]
    
    DATA %>%
      filter(Date %in% dates,
             Plate %in% plates,
             STIM_DAY %in% stim,
             CONDITION %in% cond,
             filter_pattern_func(c(POSITIVE_CTRL, NEGATIVE_CTRL, value), DATA))
  })
  
  # Combine the subsets and return
  FILTERED_SUBSET <- bind_rows(subset_list)
  # return(FILTERED_SUBSET)
  cat("Filtered Data.\n")
  
  # DATA = FILTERED_SUBSET
  
  group_vars <- c("Date", "STIM_DAY", "Plate")
  filter_var <- "CELL_LINE"  # Assuming the correct variable for DAY filter
  
  # Debugging information
  # print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  # print(paste("Filtering using variable:", filter_var))
  
  # DATA_NORMALIZED <- DATA %>% # Debugging line
  DATA_NORMALIZED <- FILTERED_SUBSET %>%
    group_by(Date, Plate) %>%
    filter(Date == FILTER_VALUES) %>%
    mutate(baseline_control_value =  min(Concentration),
           Concentration_REDUCED = Concentration - baseline_control_value,
           Concentration_REDUCED = case_when(Concentration_REDUCED < 0 ~ 0, TRUE ~ Concentration_REDUCED)) %>%
    ungroup() %>%
    group_by(Plate) %>%
    mutate(max_control_value = max(Concentration),
           Concentration_REDUCED_NORMALIZED = Concentration_REDUCED/max_control_value)
  cat("Extracted min value per plate to perform baseline-normalization.\n")
  cat("Extracted max value per plate to perform upper limit normalization for relative secretion.\n")
  
  DATA_NORMALIZED <- DATA_NORMALIZED %>%
    group_by(CELL_LINE, CONDITION, STIM_DAY) %>%
    mutate(triplicate_mean_per_day = mean(Concentration_REDUCED_NORMALIZED),
           triplicate_sd_per_day = sd(Concentration_REDUCED_NORMALIZED))
  # return(DATA_NORMALIZED)
  cat("Calculated means per cell line, condition, and stimulation.\n")
  
  # Set the dodge width
  dodge_width  <- 0.9
  jitter_width <- 0.3
  
  cat("Plotting results:\n")
  DR2_PLOT_1 <- ggplot(DATA_NORMALIZED, aes(x = CONDITION, y = triplicate_mean_per_day, fill = CELL_LINE, group = CELL_LINE)) +
    geom_col(position = position_dodge(width = dodge_width)) +
    geom_point(aes(y = Concentration_REDUCED_NORMALIZED, shape = as.factor(DILUTION)), 
               position = position_jitterdodge(dodge.width = dodge_width, jitter.width = jitter_width), 
               size = 3,
               show.legend = T) +
    geom_errorbar(aes(ymin = triplicate_mean_per_day - triplicate_sd_per_day,
                      ymax = triplicate_mean_per_day + triplicate_sd_per_day),
                  position = position_dodge(width = dodge_width)) +
    facet_wrap(~paste0(STIM_DAY, "h"), nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  DR2_PLOT_2 <- ggplot(DATA_NORMALIZED, aes(x = CONDITION, y = Concentration_REDUCED_NORMALIZED, fill = CELL_LINE, group = CELL_LINE)) +
    geom_col(position = position_dodge(width = dodge_width)) +
    facet_wrap(~paste0(STIM_DAY, "h")+DILUTION, nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  DR2_PLOT_3 <- ggplot(DATA_NORMALIZED, aes(x = CONDITION, y = triplicate_mean_per_day, fill = STIM_DAY, group = STIM_DAY)) +
    geom_col(position = position_dodge(width = 1)) +
    facet_wrap(~CELL_LINE, nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  DR2_PLOT_4 <- ggplot(DATA_NORMALIZED, aes(x = STIM_DAY, y = triplicate_mean_per_day, fill = CONDITION, group = CONDITION)) +
    geom_col(position = position_dodge(width = 1)) +
    facet_wrap(~CELL_LINE+CONDITION, nrow = 1) +
    labs(x = x_label, y = y_label) +
    ggtitle(label = plot_title, subtitle = subtitle) +
    theme_cowplot() +
    theme(legend.position = "bottom")
  
  
  return(list(FILTERED_SUBSET = FILTERED_SUBSET, 
              DATA_NORMALIZED = DATA_NORMALIZED, 
              DR2_PLOT_1 = DR2_PLOT_1, 
              DR2_PLOT_2 = DR2_PLOT_2, 
              DR2_PLOT_3 = DR2_PLOT_3, 
              DR2_PLOT_4 = DR2_PLOT_4))
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
# DALI Z-Score Plot

#Load libraries
library(pacman)
pacman::p_load(data.table)

extract_dali_scores <- function(data_path) {
  dir_path  <- dirname(data_path)
  file_name <- sub("\\..*$", "", basename(data_path)) 
  
  # read the input file
  lines     <- readLines(data_path, warn = FALSE)
  
  # identify the lines of interest (all lines that contain the "MOLECULE:" pattern)
  molecule_lines <- grep("MOLECULE:", lines, value = TRUE)
  
  # 1) Split Lines into Parts: Split each line into parts. For lines with a long description, the "%id" column will be part of the description.
  # 2) Extract and Concatenate Description: Identify and concatenate the description parts correctly. 
  #    This involves splitting the line at "MOLECULE:" and handling the parts before and after this keyword separately.
  # 3) Reconstruct the Parts: Reconstruct the line with the correct number of columns, ensuring that the "%id" and "PDB" columns are correctly identified.
  data_list <- lapply(molecule_lines, function(line) {
    # Split the line by spaces, but not when spaces are part of the description
    parts <- strsplit(line, "\\s+(?![^MOLECULE:]*\")", perl = TRUE)[[1]]
    # Adjust based on the expected number of columns before the description
    if (length(parts) > 8) {
      # Concatenate the description parts back together
      description <- paste(parts[9:length(parts)], collapse = " ")
      parts       <- parts[1:9]
      parts[9]    <- description
    }
    return(parts)
  })
  
  # Convert to a DataFrame
  data_df <- do.call(rbind, data_list)
  data_df <- as.data.frame(data_df, stringsAsFactors = FALSE)
  
  # Set column names
  colnames(data_df) <- c("empty", "No", "PDB_Chain", "Z", "RMSD", "Lali", "Nres", "%ID", "Description")
  
  # Remove 'MOLECULE:' pattern
  data_df$Description <- gsub("MOLECULE:", "", data_df$Description)
  
  # Replace underscores with spaces
  data_df$Description <- gsub("_", " ", data_df$Description)
  
  # Remove the semicolons
  data_df$Description <- gsub(";", "", data_df$Description)
  
  # Split the 'Chain' Column 
  split_chain <- strsplit(data_df$PDB_Chain, "-", fixed = TRUE)
  
  # Create New Columns
  data_df$PDB    <- sapply(split_chain, function(x) ifelse(length(x) > 1, x[1], NA))
  data_df$Chain  <- sapply(split_chain, function(x) ifelse(length(x) > 1, x[2], x[1]))
  
  # Reset the row names to normal numbers
  rownames(data_df) <- seq(nrow(data_df))
  
  data_df_final         <- data_df[, c("PDB", "Chain", "Z", "RMSD", "Lali", "Nres", "%ID", "Description", "PDB_Chain")]
  data_df_final$Origin  <- file_name
  data_df_final$ID      <- unlist(lapply(strsplit(data_df_final$Origin,   "_", fixed=TRUE), function(x) return(x[2])))
  data_df_final$Origin  <- unlist(lapply(strsplit(data_df_final$Origin,   "_", fixed=TRUE), function(x) return(x[3])))
  
  # Write the final data frame to a csv file
  fwrite(x = data_df_final, file = paste0(dir_path, "/", file_name, "_extracted.csv"))
}

# Function to process all .txt files in a directory and combine CSVs
process_directory <- function(directory_path) {
  # List all .txt files in the directory
  txt_files <- list.files(directory_path, pattern = "\\.txt$", full.names = TRUE)
  
  # Extract data from each .txt file and create a CSV summary
  for (txt_file in txt_files) {
    extract_dali_scores(data_path = txt_file)
  }
  
  # List all created CSV files in the directory
  csv_files <- list.files(directory_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Read each CSV file and combine them into one data frame
  combined_df <- rbindlist(lapply(csv_files, fread), fill = TRUE)
  
  # Write the combined DataFrame to a CSV file
  fwrite(combined_df, file = file.path(directory_path, "combined_summary.csv"))
}

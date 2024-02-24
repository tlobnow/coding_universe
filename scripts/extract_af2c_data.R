# Load required libraries (if you don't have pacman installed yet, please unhash the installation line and run it)
#install.packages("pacman")
library(pacman); p_load(data.table, knitr, dplyr, fs)

source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/coding_universe/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/coding_universe/main/scripts/functions.R",
                     no  =  "~/Documents/Github/coding_universe/scripts/functions.R"))

OVERWRITE = T
COPY_BEST_MODELS_TO_EXTRA_FOLDER = T

# Provide the name of the folder where the run folders are located
# LOC = "my_run_files"
# LOC = "CONTROLS"
# LOC = "ADLD"
# LOC = "MAJER"
# LOC = "TEST"
# LOC = "EDLD"
# LOC = "MYD88"
# LOC = "IRAK4"
# LOC = "IRAK1"
# LOC = "TUBULIN"
# LOC = "project_2xTIR"
LOC = "20240222"


# Provide the path where this folder lies
# PATH = "~/Documents/etc/PhD/2xTIR/"
# PATH = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/"
# PATH = "~/Desktop/output/"
# PATH = "/Volumes/TAYLOR-LAB/Finn/NiranjanThanksYou/"
PATH = "~/Desktop/"

# Provide the path to where summaries should be saved
# SUMMARY_FOLDER = "~/Documents/etc/PhD/2xTIR/SUMMARIES/"
# SUMMARY_FOLDER = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/SUMMARIES/"
SUMMARY_FOLDER = file.path(PATH, "SUMMARIES/")
if (dir.exists(SUMMARY_FOLDER)) {print(paste0(basename(SUMMARY_FOLDER), " directory exists."))} else {dir.create(SUMMARY_FOLDER)}

# Where to copy the best models to?
base_path = file.path(PATH, "BOTB")
if (dir.exists(base_path)) {print(paste0(basename(base_path), " directory exists."))} else {dir.create(base_path)}

# Automatically generates the main folder variable
MAIN      = file.path(PATH, LOC, "/")


# Extract all JSON files from the provided location automatically and store the results in two files (with/without recycles)
if (file.exists(paste0(SUMMARY_FOLDER, LOC, ".csv")) & OVERWRITE == F) {
  dataframe = fread(paste0(SUMMARY_FOLDER, LOC, ".csv"))
  # dataframe_with_recycles = fread(paste0(SUMMARY_FOLDER, LOC, "_summaryWithRecycles.csv"))
} else {
  # run_extraction(LOC = LOC, MAIN = MAIN, SUMMARY_FOLDER = SUMMARY_FOLDER)
  ### dataframe_with_recycles = fread(paste0(SUMMARY_FOLDER, LOC, "_summaryWithRecycles.csv"))
  run_extraction_wrapper(LOC = LOC, MAIN = MAIN, SUMMARY_FOLDER = SUMMARY_FOLDER)
  dataframe = fread(paste0(SUMMARY_FOLDER, LOC, ".csv"))  %>% filter(!is.na(pTM))
}

dataframe <- dataframe %>%
  mutate(FILE = str_replace_all(FILE, "IL1R_MOUSE-TIR_x1_|_IL1RAP_MOUSE-TIR_1|IL1RAP_MOUSE-TIR_x1_", ""))


if (COPY_BEST_MODELS_TO_EXTRA_FOLDER) {
  copy_max_score_files(
    dataframe = dataframe,
    main_dir = MAIN,
    base_path = base_path,
    file_col = "FILE",       # column name for file identifier
    origin_col = "ORIGIN",   # column name for origin identifier
    recycle_col = "RECYCLE"  # column name for recycle identifier
  )
}

max_iScores_df <- dataframe %>%
  group_by(as.factor(FILE)) %>%
  slice_max(order_by = iScore, n = 1) %>%
  select(FILE, RECYCLE, iScore, piTM) %>% 
  arrange(desc(iScore)) %>%
  ungroup()

top30 <- head(max_iScores_df[,2:5], 30)
kable(top30)  

plot_alphafold_results(LOC = LOC, SUMMARY_FOLDER = SUMMARY_FOLDER, plot_interactive = T)
plot_alphafold_results(LOC = LOC, SUMMARY_FOLDER = SUMMARY_FOLDER, plot_interactive = T, best_only = T)
plot_alphafold_results(LOC = LOC, SUMMARY_FOLDER = SUMMARY_FOLDER, 
                       best_only = TRUE,
                       # plot_interactive = T,
                       # pattern = "MYD88_MOUSE-DD"
                       # pattern = "DC7",
                       #plot_labels = TRUE
                       )
plot_alphafold_results(LOC = LOC, SUMMARY_FOLDER = SUMMARY_FOLDER) 
                       


################################################################################
################ ADDITIONAL DOCUMENTATION OF THE FUNCTIONS USED ################
################################################################################

#' Extracts and Processes AlphaFold Results
#'
#' This function extracts data from JSON files, processes it, and appends it to summary CSV files.
#' It's designed to work specifically with results from AF2C.
#' @param LOC A string representing the location of the target data folder.
#' @param MAIN A string representing the main directory where the summaries are stored. Defaults to "~/Desktop/SUMMARIES/LOC/".
#' @param SUMMARY_FOLDER A string representing the directory where summary files are saved. Defaults to "~/Desktop/SUMMARIES/".
#' @param ADD_2_EXISTING_DF A logical value indicating whether to append the extracted data to an existing dataframe. Defaults to FALSE.
#' @param EXISTING_DF A string representing the path of the existing dataframe to which the data should be appended, if ADD_2_EXISTING_DF is TRUE.
#' \dontrun{
#' run_extraction(LOC = "my_location")
#' run_extraction(LOC = "my_location", ADD_2_EXISTING_DF = TRUE, EXISTING_DF = "~/path_to_existing_df.csv")
#' }
#' 

#' Plot AlphaFold Results
#'
#' This function takes AlphaFold results from a specified location and plots them.
#' It categorizes the `iScore` into different confidence levels and allows filtering based on a given pattern.
#' @param LOC A character string specifying the location of the data. This will be appended to the SUMMARY_FOLDER path to generate the full path.
#' @param SUMMARY_FOLDER The folder where summary results are stored. By default, it is set to "~/Desktop/SUMMARIES/".
#' @param xlab A character string for the x-axis label. Default is "iScore".
#' @param ylab A character string for the y-axis label. Default is "piTM".
#' @param plot_title A character string for the plot title. If NULL (the default), a default title is generated based on the LOC value.
#' @param pattern A character string. If provided, the function will filter the dataset to include only rows where the FILE column matches this pattern.
#' @param best_only A logical value. If TRUE, only the best results (as determined by iScore) will be plotted.
#' @examples
#' plot_alphafold_results(LOC = LOC, SUMMARY_FOLDER = SUMMARY_FOLDER)
#' plot_alphafold_results(LOC = LOC, SUMMARY_FOLDER = SUMMARY_FOLDER, pattern = "specific_pattern", best_only = TRUE)

################################################################################
### SOTA PATHS #################################################################
################################################################################

input_dir  <- "SOTA/01_raw_data"
temp_dir   <- "SOTA/02_temp_data"
output_dir <- "SOTA/03_output_data"
elisa_dir  <- file.path(input_dir, "ELISA_PLATES")
af_dir     <- file.path(input_dir, "AF3_PREDICTIONS")

input_files   <- list.files(input_dir,  full.names = TRUE)
temp_files    <- list.files(temp_dir,   full.names = TRUE)
out_files     <- list.files(output_dir, full.names = TRUE)
elisa_plates  <- list.files(elisa_dir,  full.names = TRUE)
af_json_files <- list.files(af_dir,     full.names = TRUE, recursive = TRUE, pattern = "job_request") %>% str_subset("(?i)dld") %>% str_subset("\\.json$")

# group file lists
excel_files <- input_files[grepl("\\.xlsx$",  input_files,  ignore.case = TRUE)]
csv_files   <- input_files[grepl("\\.csv$",   input_files,  ignore.case = TRUE)]
tsv_files   <- input_files[grepl("\\.tsv$",   input_files,  ignore.case = TRUE)]
elisa_files <- elisa_plates[grepl("\\.xlsx$", elisa_plates, ignore.case = TRUE)]
order_files <- csv_files[grepl("ORDER",       csv_files,    ignore.case = TRUE)]
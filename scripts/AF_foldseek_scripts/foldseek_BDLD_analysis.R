### FoldSeek Results of BDLDs ###

library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr)

main <- "/Users/u_lobnow/Documents/Github/phd/foldseek_results/"

head(folder_names <- list.dirs(main)[-1]) # the [-1] omits the first row: [1] "/Users/u_lobnow/Documents/Github/phd/foldseek_results/"

### FUNCTIONS ###

# Function to read files and add ID and db columns
read_and_add_ID_and_db <- function(folder_path, file_names, db_names) {
  ID <- unlist(lapply(strsplit(basename(folder_path), "_mmseqs", fixed=TRUE), function(x) return(x[1])))
  file_paths <- file.path(folder_path, file_names)
  data_list <- lapply(seq_along(file_paths), function(i) {
    data <- fread(file_paths[i], sep = "\t")
    data$ID <- ID
    data$db <- db_names[i]
    return(data)
  })
  return(data_list)
}

# Define database abbreviations
db_names <- c("proteome", "swissprot", "afdb50", "gmgcl_id", "mgnify_esm30", "pdb100")


### READ IN DATA ###


# Initialize an empty list to store data frames
df_list <- list()

# Iterate over folder_names and read files
for (folder in folder_names) {
  data_list <- read_and_add_ID_and_db(folder, c("alis_afdb-proteome.m8", 
                                                "alis_afdb-swissprot.m8", 
                                                "alis_afdb50.m8", 
                                                "alis_gmgcl_id.m8", 
                                                "alis_mgnify_esm30.m8", 
                                                "alis_pdb100.m8"),
                                      db_names)
  df_list <- c(df_list, data_list)
}

# Combine all data frames into a single data frame
df <- rbindlist(df_list, fill = TRUE)

colnames(df) <- c("query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
                  "tstart", "tend", "prob", "evalue", "V13", "V14", "V15",
                  "qaln", "taln", "bits", "tseq", "taxid", "taxname", "ID", "db")

### FILTERING ###

df_filtered <- df %>%
  filter(!is.na(query),   # remove empty rows
         !is.na(taxid),   # remove rows with unknown taxonomy match
         evalue <= 0.05   # filter for significant matches
         )

# BDLD_table <- table(as.factor(df_filtered$taxname), as.factor(df_filtered$ID))


### PLOTTING ###

heat_prep <- df_filtered %>% select(ID, taxname)
heat_prep$ID      <- as.factor(heat_prep$ID)
heat_prep$taxname <- as.factor(heat_prep$taxname)

heat_prep2 <- heat_prep %>%
  group_by(ID, taxname) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  complete(ID, taxname, fill = list(count = 0))


p <- ggplot(data = heat_prep2, aes(x = ID, y = taxname, fill = count)) +
  geom_tile(color = "white",
            lwd = 0.3,
            linetype = 1) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  coord_fixed()

p

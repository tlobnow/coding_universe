# Bio3D

pacman::p_load(bio3d,tidyverse,msa,stringr,ggrepel)

PATH    = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/PDB_BEST_iSCORE/BOTB"
SAVE_TO = "~/Documents/Github/coding_universe/bio3d/raw_pdbs"
ids     = paste0(basename.pdb(list.files("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/PDB_BEST_iSCORE/BOTB", pattern = ".pdb")), "_A")

TESTED <-  c("BDLD_67", "BDLD_14", "BDLD_7", "BDLD_69",                       # 0 PASSES ~ LOW
             "BDLD_10", #<--- negative control! 
             "BDLD_50", "BDLD_2",                                             # 1/2 PASSES ~ MEDIUM HIGH
             "BDLD_46", "BDLD_57", "BDLD_62", "BDLD_23",                      # 1/2 PASSES ~ MEDIUM
             "BDLD_27", "BDLD_73", "BDLD_13", "BDLD_6", "BDLD_22", "BDLD_24") # 3 PASSES ~ HIGH

# pdb <- read.pdb("3MOP")
# ids <- "3MOP_A"
# 
# raw_files <- get.pdb(ids = ids, path = SAVE_TO)
# files     <- pdbsplit(raw_files, ids, path = paste0(SAVE_TO, "/split_chain"))
# pdb       <- read.pdb(files)
# print(pdb)
# head(pdb$atom)
# print(pdb$xyz)
# 
# # Perform a Normal Mode Analysis (NMA for short) to predict the flexibility of this particular structure
# modes <- nma(pdb)
# plot(modes)
# plot(modes, sse=pdb)
# plot.bio3d(modes$fluctuations, sse=pdb, sheet.col="orange", helix.col="purple", typ="l", lwd=3, ylab="Fluctuations from NMA (custom label)")

# # List all PDB files in the specified directory
# raw_files <- list.files(PATH, pattern = ".pdb", full.names = TRUE)

# raw_files <- get.pdb(c("6I3N", "3MOP"))

# # Create ids for each PDB file
# ids <- paste0(basename.pdb(raw_files), "_A")
# 
# # Split PDB files based on chain
# files <- pdbsplit(raw_files, ids, path = paste0(SAVE_TO, "/split_chain"))

# pdb = read.pdb("/Users/u_lobnow/Documents/Github/coding_universe/P22366.pdb")
# modes <- nma(pdb)
# plot(modes)
# plot.bio3d(modes$fluctuations, sse=pdb, sheet.col="orange", helix.col="purple", typ="l", lwd=3, ylab="Fluctuations from NMA (custom label)")


# Function to extract data for each file
extract_data <- function(file) {
  pdb <- read.pdb(paste0(SAVE_TO, "/split_chain/", file))
  modes <- nma(pdb)
  data.frame(
    File = basename.pdb(file),
    X = 1:length(modes$fluctuations),
    Y = modes$fluctuations
  )
}

# MYD88_MOUSE <- extract_data("MYD88_MOUSE.pdb")

# Trim the PDB to include only residues 32 to 108
modes_MYD88_MOUSE <- nma(trim.pdb(read.pdb("/Users/u_lobnow/Documents/Github/coding_universe/bio3d/raw_pdbs/P22366.pdb"), resno=c(32:108)))
modes_MYD88_HUMAN <- nma(trim.pdb(read.pdb("/Users/u_lobnow/Documents/Github/coding_universe/bio3d/raw_pdbs/Q99836.pdb"), resno=c(32:108)))

MYD88_MOUSE <- data.frame(
  File = "MYD88_MOUSE",
  X = 1:length(modes_MYD88_MOUSE$fluctuations),
  Y = modes_MYD88_MOUSE$fluctuations)

MYD88_HUMAN <- data.frame(
  File = "MYD88_HUMAN",
  X = 1:length(modes_MYD88_HUMAN$fluctuations),
  Y = modes_MYD88_HUMAN$fluctuations)

# Extract data for all files
files    <- list.files(paste0(SAVE_TO, "/split_chain"))
all_data_raw <- lapply(files, extract_data) %>% bind_rows() %>% as.data.frame()
all_data_raw <- rbind(all_data_raw, MYD88_MOUSE, MYD88_HUMAN)

all_data <- all_data_raw %>%
  mutate(ID = case_when(!is.na(str_extract(File, "BDLD_\\d{1,2}")) ~ str_extract(File, "BDLD_\\d{1,2}"), TRUE ~ File),
         X = case_when(File == "6I3N_A" ~ X - 29,
                       File == "3MOP_A" ~ X - 29,
                       File == "MYD88_MOUSE" ~ X - 2,
                       File == "MYD88_HUMAN" ~ X - 2,
                       TRUE ~ X)) #%>% filter(X < 101)

rescue_data <- all_data %>%
  filter(ID %in% c("3MOP_A", "6I3N_A", "MYD88_MOUSE", "MYD88_HUMAN", TESTED)) %>%
  mutate(ORIGIN = ifelse(ID %in% TESTED, "BDLD", "METAZOAN"))

max_data          <- all_data %>% group_by(File) %>% filter(Y < 101) %>% filter(Y == max(Y))
filtered_max_data <- max_data %>% filter(ID %in% unique(rescue_data$ID))
SELECTION         <- c("BDLD_57", "BDLD_6", "3MOP_A", "6I3N_A", "MYD88_MOUSE", "MYD88_HUMAN", 
                       "BDLD_27", "BDLD_14")

label_data <- rescue_data %>%
  group_by(ID) %>%
  filter(X <= 80, X >= 0, Y <= 1) %>%
  filter(Y == max(Y)) %>%
  ungroup()

label_data <- rescue_data %>%
  filter(ID %in% SELECTION, 
         X <= 80, X >= 0, Y <= 1) %>%
  group_by(ID) %>%
  filter(Y == max(Y)) %>%
  ungroup()

rescue_data %>%
  filter(X <= 100, Y <= 1) %>%
  ggplot(aes(x = X, y = Y, group = ID)) +
  geom_line(aes(color = ID), alpha = 0.1) +
  geom_line(data = filter(rescue_data, ID %in% SELECTION), aes(color = ID), alpha = 1) +
  geom_label_repel(data = label_data, aes(label = ID, color = ID)) +
  ylim(0, 1) +
  xlim(0, 80) +
  theme_minimal() +
  theme(legend.position = 0) +
  facet_wrap(~ORIGIN, ncol = 1)

# all_data %>%
#   ggplot(aes(x = X, y = Y, group = ID)) +
#   geom_line(aes(color = ID)) +
#   geom_label(data = max_data, aes(label = ID, color = ID), nudge_y = 0.05) +
#   ylim(0, 100) +
#   theme_minimal() +
#   theme(legend.position = 0)

# for (file in files) {
#   pdb <- read.pdb(paste0(SAVE_TO, "/split_chain/", file))
#   modes <- nma(pdb)
#   png(paste0(SAVE_TO, "/figures/", basename.pdb(file),".png"), width=800, height=600)
#   plot <- plot.bio3d(modes$fluctuations, sse=pdb, sheet.col="orange", helix.col="purple", typ="l", lwd=3, ylab="Fluctuations from NMA (custom label)")
#   dev.off()
# }
# 
# # Define the PNG path and initialize the graphics device
# png(paste0(SAVE_TO, "/figures/combined_plot.png"), width=800, height=600)
# 
# # Read the first pdb file and set up an initial plot
# pdb <- read.pdb(paste0(SAVE_TO, "/split_chain/", files[1]))
# modes <- nma(pdb)
# plot(modes$fluctuations, type="l", lwd=3, ylab="Fluctuations from NMA (custom label)", ylim = c(0,10))
# 
# # Now loop through the rest of the files and add to the plot
# for (i in 2:length(files)) {
#   pdb <- read.pdb(paste0(SAVE_TO, "/split_chain/", files[i]))
#   modes <- nma(pdb)
#   lines(modes$fluctuations, col=i)  # Change color for each line
# }
# 
# # Close the graphics device
# dev.off()

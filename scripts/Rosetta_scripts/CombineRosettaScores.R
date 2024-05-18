library(pacman)
p_load(tidyr, data.table, ggplot2, biomaRt, dplyr, scales, cowplot, ggrepel, ggalt, stringr)
select <- dplyr::select
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# # get list of all available info
# filters <- listFilters(ensembl)
# attributes <- listAttributes(ensembl)
# Symbol_fx <- function (x) {
#   getBM(attributes=c('uniprot_gn_symbol','description'), 
#         filters = 'pdb', 
#         values = x, 
#         mart = ensembl)
# }

# Step 1: List and read all files with the pattern "*.sc"
filenames  <- list.files("~/Documents/Github/master_thesis/rosetta_scores/scores/", pattern = "\\.sc$", full.names = TRUE)
filenames2 <- list.files("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Rosetta_results/merged_scores_chain_A", pattern = "\\.sc$", full.names = TRUE)

data_list  <- lapply(filenames,  read.csv, header = TRUE, sep = "")
data_list2 <- lapply(filenames2, read.csv, header = TRUE, sep = "")

# Step 2: Apply extend_fx to add columns for all possibilities
extend_fx <- function(x) {
  for (col in paste("dG_A", LETTERS[2:26], sep = "")) {
    NAMES <- colnames(x)[-c(1:4)]
    if (!(col %in% NAMES)) {
      x[[col]] <- NA
    }
  }
  return(x)
}

data_list_extended <- lapply(data_list, extend_fx)
data_list_extended2 <- lapply(data_list2, extend_fx)

data_list_extended2 <- data_list_extended2[-20]
data_list_extended2 <- data_list_extended2[-14]
data_list_extended2 <- data_list_extended2[-8]

# Step 3: Combine the prepared lists with rbind() and save as "score_all"
score_all1 <- do.call(rbind, data_list_extended)
score_all2 <- do.call(rbind, data_list_extended2)
score_all <- rbind(score_all1, score_all2)


# Step 4: Mutate the "name" column by splitting the "description" column
score_all <- score_all %>%
  mutate(name = strsplit(as.character(description), split = "_startconstraints") %>% sapply(`[[`, 1))

# Step 5: Calculate the number of interfaces per row
score_all$n_rows <- rowSums(!is.na(score_all[, 6:length(score_all)]))

# Step 6: Remove columns with only NA values
# score_all_clean <- score_all %>% discard(~all(is.na(.) | . == ""))
# Remove columns with only NA values
score_all_clean <- score_all %>%
  select(-where(~all(is.na(.))))

# Step 7: Calculate the interface_min using the min() function for available columns with the dG_A* pattern
score_all_clean$interface_min <- apply(score_all_clean[, grep("^dG_A", colnames(score_all_clean))], 1, min, na.rm = TRUE)

# Step 8: Calculate the mean for available columns with the dG_A* pattern
score_all_clean$mean_dG <- rowMeans(score_all_clean[, grep("^dG_A", colnames(score_all_clean))], na.rm = TRUE)
score_all_clean$median_dG <- apply(score_all_clean[, grep("^dG_A", colnames(score_all_clean))], 1, median, na.rm = TRUE)

# Step 9: Calculate the sum of all values in the columns that match the pattern "^dG_A" and subtract the value of the "interface_min" column from the calculated sum
score_all_clean$sum_sec_interface <- sum(score_all_clean[, grep("^dG_A", colnames(score_all_clean))], na.rm = TRUE) - score_all_clean$interface_min

# Step 10: add cols
score_all_prepped <- score_all_clean %>%
  mutate(BDLD_ID = gsub(".*(BDLD_[0-9]{1,2}).*", "\\1", name),
         NUM_MONOMERS = gsub(".*_x([0-9]+).*", "\\1", name),
         MODEL_NUMBER = gsub(".*model_([0-9]+).*", "\\1", name),
         type = fcase(substr(name, start = 1, stop = 4) == "BDLD", "BDLD",
                      substr(name, start = 1, stop = 5) == "CHIMY", "CHIMY_BDLD",
                      substr(name, start = 1, stop = 3) == "6E9", "DHF",
                      substr(name, start = 1, stop = 3) == "HA1", "HA13",
                      substr(name, start = 1, stop = 3) != "6E9" & substr(name, start = 1, stop = 3) != "HA1", "DD")) %>%
  filter(!is.na(dG_AB)) %>%
  group_by(BDLD_ID, NUM_MONOMERS, type) %>%
  summarize(mean_dG_all       = mean(mean_dG, na.rm = TRUE),
            median_dG_all     = median(mean_dG, na.rm = TRUE),
            interface_min_all = min(interface_min),
            interface_min     = interface_min,
            median_dG         = median_dG,
            mean_dG           = mean_dG,
            BDLD_ID           = BDLD_ID,
            type              = type,
            sum_sec_interface = sum_sec_interface,
            name              = name)


# Output the final data frame
#score_all_clean



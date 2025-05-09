suppressPackageStartupMessages({
library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(zoo)
library(tidyr)
library(data.table)
library(R.utils)
library(argparse)
})

# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for running EWAS")
parser$add_argument("--infile",
                    required=TRUE,
                    help="Path to the metabolomics results excel file.")
parser$add_argument("--outfile",
                    required=TRUE,
                    help="Path for the output file")


# parse arguments
args <- parser$parse_args()
raw_results <- args$infile
outfile_path <- args$outfile

# Read and transpose the header rows
header_rows <- read_excel(raw_results, range = "A1:CX2", col_names = FALSE)
header_t <- t(as.matrix(header_rows))
header_df <- as.data.frame(header_t, stringsAsFactors = FALSE)
colnames(header_df) <- c("group", "label")

# Forward-fill group before dropping empty columns
header_df <- header_df %>%
  mutate(group = zoo::na.locf(group, na.rm = FALSE),
         label = replace_na(label, ""),
         combined = str_trim(paste(group, label)))

# Drop the first two columns (empty columns)
header_df <- header_df[-c(1, 2), ]

# Read in data without the headers
data <- read_excel(raw_results, skip = 2, col_names = FALSE)

# Align and apply headers
header_df <- header_df[seq_len(ncol(data)), ]
colnames(data) <- header_df$combined

# Filter and subset data
filtered_data <- data %>%
  select(
    `Sample Name`,
    `Sample Type`,
    contains("Area")
  ) %>%
  filter(`Sample Type` %in% c("Cal", "Sample")) %>% # select only the standards and samples
  select(-`Sample Type`) %>%
  mutate(across(`Sample Name`, str_replace, "(^std[^_]*)(.+)", "\\1")) %>% # Remove anything after the first '_' in the standards names
  column_to_rownames(var="Sample Name")

# Compute the peak area ratio for each metabolite ===
area_cols <- colnames(filtered_data) # Identify all area columns
# Extract the base metabolite names and abbreviations
metabolite_matches <- str_match(area_cols, "(^[^(]*)\\s\\(([^()]*).+")

metabolite_info <- tibble(
  full_name = area_cols,
  name      = metabolite_matches[, 2],
  abbrev    = metabolite_matches[, 3]
) %>%
  filter(!is.na(abbrev)) %>%
  mutate(
    base_name = str_replace(name, "-d[3|4]$", ""),
    is_heavy = str_detect(name, "-d[3|4]") #create a column that makes it easy to filter sample vs heavy standard
  )

# Compute peak area ratios
ratio_df <- purrr::map_dfc(unique(metabolite_info$base_name), function(base) {
  rows <- metabolite_info %>% filter(base_name == base)

  heavy_row     <- rows %>% filter(is_heavy)
  non_heavy_row <- rows %>% filter(!is_heavy)

  if (nrow(heavy_row) == 1 && nrow(non_heavy_row) == 1) {
    num <- filtered_data[[non_heavy_row$full_name]]
    den <- filtered_data[[heavy_row$full_name]]
    tibble(!!base := num / den) # Use the base names as the columns for the new dataframe with the peak area ratio values
  } else {
    NULL
  }
})

# Add sample names back as rownames
ratio_df <- as.data.frame(ratio_df) # convert to dataframe because adding rownames to a tibble is depreciated
rownames(ratio_df) <- rownames(filtered_data)

# Transpose dataframe to make the wide format expected by SCALiR
ratio_df <- as.data.frame(t(ratio_df))

# Add rownames back to dataframe as a column with the label expected by SCALiR
ratio_df <- ratio_df %>%
  rownames_to_column(var="peak_label") %>%
  mutate(peak_label = str_to_lower(peak_label))

# Make sure the order of the metabolites is arranged alphabetically
ratio_df <- ratio_df %>% arrange(peak_label)

# Export file
fwrite(ratio_df, file = outfile_path)
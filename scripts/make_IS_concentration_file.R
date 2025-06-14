suppressPackageStartupMessages({
library(readxl)
library(dplyr)
library(data.table)
library(stringr)
library(tibble)
library(tidyr)
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
standards_path <- args$infile
outfile_path <- args$outfile

# Read in file
suppressMessages(standards_df <- read_excel(standards_path))

# Add a column for base_name
standards_df <- standards_df %>%
  mutate(
    Compound = str_to_lower(Compound),
    base_name = str_replace(Compound, "^d[3|4]-", ""))

# Make a dataframe of the internal standard stock concentration as a ratio of the working standard/heavy standard
std_ratio_df <- purrr::map_dfc(unique(standards_df$base_name), function(base) {
  rows <- standards_df %>% filter(base_name == base)
  heavy_row     <- rows %>% filter(Type=="Heavy (IS)") %>% select(Compound, `Final Concentration (ng/uL)`)
  non_heavy_row <- rows %>% filter(Type=="Light (WS)") %>% select(Compound, `Final Concentration (ng/uL)`)
  if (nrow(heavy_row) == 1 && nrow(non_heavy_row) == 1) {
    num <- non_heavy_row$`Final Concentration (ng/uL)`
    den <- heavy_row$`Final Concentration (ng/uL)`
    tibble(!!base := num / den)
  } else {
    NULL
  }
})

# Transpose dataframe
std_ratio_df <- as.data.frame(std_ratio_df)
rownames(std_ratio_df) <- "std-1"
std_ratio_df <- as.data.frame(t(std_ratio_df))

# Extrapolate other standard concentration levels
std_ratio_df <- std_ratio_df %>%
  mutate(
    "std-0" = 0,
    "std-0.5" = `std-1`*0.5,
    "std-2.5" = `std-1`*2.5,
    "std-5" = `std-1`*5,
    "std-10" = `std-1`*10,
    "std-25" = `std-1`*25,
    "std-50" = `std-1`*50,
    "std-100" = `std-1`*100,
    "unit" = "ng/uL",
    "internal_standard_injection_concentration" = 0.25
  ) %>%
  # removed std-2.5 due to consistent poor accuracy in measurement
  select("std-0", "std-0.5", "std-1", "std-5",
         "std-10", "std-25", "std-50", "std-100", unit,
         internal_standard_injection_concentration) %>%
  rownames_to_column(var="peak_label")

std_ratio_df <- std_ratio_df %>% arrange(peak_label)

fwrite(std_ratio_df, file = outfile_path)
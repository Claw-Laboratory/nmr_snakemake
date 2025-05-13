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
parser$add_argument("--qc-file",
                    required=TRUE,
                    help="Path to the metabolomics results excel file.")
parser$add_argument("--batch-file",
                    required=TRUE,
                    help="Path to the metabolomics results excel file.")
parser$add_argument("--outfile",
                    required=TRUE,
                    help="Path for the output file")


# parse arguments
args <- parser$parse_args()
qc_file <- args$qc_file
batch_file <- args$batch_file
outfile_path <- args$outfile

qc_file <- "data/QCSolutions_002_20220630.xlsx"
batch_file <- "results/NIC_MET_QUANT_RESULTS_02.14.23/NIC_MET_QUANT_RESULTS_02.14.23_concentrations.csv"

# Read in QC stock solution file (skip the first 20 rows with information we don't need)
suppressMessages(qc_df <- read_excel(qc_file, skip = 20, col_names = TRUE))
qc_df <- qc_df %>% mutate(Compound = str_to_lower(Compound))
glimpse(qc_df)

# Read in SCALiR predicted concentrations file
batch_df <- fread(batch_file)
glimpse(batch_df)

# Subset batch file to only QC Samples
batch_test <- batch_df %>%
  filter(grepl('QC', ms_file)) %>%
  mutate(qc_level = str_replace(ms_file, "^QC_(Low|Med|High)_.+", "\\1")) %>%
  left_join(qc_df, by=c("peak_label"="Compound")) %>%
  mutate(dif_from_nominal = )
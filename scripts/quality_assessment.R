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
library(ggplot2)
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
# Wrangle data into a format that can be merged with the predicted sample data
qc_df <- qc_df %>%
  mutate(Compound = str_to_lower(Compound)) %>%
  select(Compound, `QC Low Mass (ng)`, `QC Med Mass (ng)`, `QC High Mass (ng)`) %>%
  rename("Low" = `QC Low Mass (ng)`,
         "Med" = `QC Med Mass (ng)`,
         "High" = `QC High Mass (ng)`) %>%
  gather(Low, Med, High, key="qc_level", value = "nominal_conc")

# Read in SCALiR predicted concentrations file
batch_df <- fread(batch_file)
glimpse(batch_df)

# Combine QC sample nominal and predicted concentration data
batch_df <- batch_df %>%
  filter(grepl('QC', ms_file)) %>%
  mutate(qc_level = str_replace(ms_file, "^QC_(Low|Med|High)_.+", "\\1")) %>%
  left_join(qc_df, by=c("peak_label"="Compound", "qc_level")) %>%
  mutate(dif_from_nominal = abs(pred_conc - nominal_conc),
         percent_dif = (dif_from_nominal/nominal_conc)*100)

# Precision & Accuracy estimates
batch_stats <- batch_df %>%
  select(ms_file, peak_label, pred_conc, qc_level, dif_from_nominal, percent_dif) %>%
  group_by(peak_label, qc_level) %>%
  summarize(pred_conc_sd = sd(pred_conc),
            pred_conc_avg = mean(pred_conc),
            variance = pred_conc_sd/pred_conc_avg,
            percent_variance = round(precision*100, digits = 2),
            avg_dif = mean(dif_from_nominal),
            avg_percent_dif = round(mean(percent_dif), digits = 2))

# Plot Precision/Variance estimate
batch_df %>%
  mutate(qc_level = factor(qc_level, levels = c("Low", "Med", "High"))) %>%
  ggplot(aes(x=peak_label, y=((pred_conc)*(50/1000)), col=qc_level)) +
    geom_boxplot(varwidth = T) +
    facet_wrap(~qc_level, nrow=3, scales="free_y") +
    theme_minimal(base_size = 14) +
    labs(title = "Predicted Concentration Variance",
         x = "Metabolite",
         y = "Predicted Concentration (ng/mL)") +
    guides(color=guide_legend(title="QC Level"))

# Plot Accuracy estimate
batch_df %>%
  mutate(qc_level = factor(qc_level, levels = c("Low", "Med", "High"))) %>%
  ggplot(aes(x=peak_label, y=percent_dif, col=qc_level)) +
    geom_boxplot() +
    facet_wrap(~qc_level, nrow=3, scales="free_y") +
    theme_minimal(base_size = 14) +
    labs(title = "Predicted Concentration Accuracy",
         x = "Metabolite",
         y = "Percent Different from Nominal") +
    guides(color=guide_legend(title="QC Level"))
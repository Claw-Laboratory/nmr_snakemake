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
                    help="Path to the QC sample stock solutions file.")
parser$add_argument("--batch-file",
                    required=TRUE,
                    help="Path to the SCALiR predicted concentrations file.")
parser$add_argument("--batch-name",
                    required=TRUE,
                    help="Name of batch")
parser$add_argument("--out-dir",
                    required=TRUE,
                    help="Path for the output files")


# parse arguments
args <- parser$parse_args()
qc_file <- args$qc_file
batch_file <- args$batch_file
batch_name <- args$batch_name
out_dir <- args$out_dir

# Define variables from arguments
outfile_prefix <- paste0(out_dir, "/", batch_name)

# Read in QC stock solution file (skip the first 20 rows with information we don't need)
suppressMessages(qc_df <- read_excel(qc_file, skip = 20, col_names = TRUE))

# Wrangle data into a format that can be merged with the sample data
qc_df <- qc_df %>%
  mutate(Compound = str_to_lower(Compound)) %>%
  select(Compound, `QC Low Mass (ng)`, `QC Med Mass (ng)`, `QC High Mass (ng)`) %>%
  rename("Low" = `QC Low Mass (ng)`,
         "Med" = `QC Med Mass (ng)`,
         "High" = `QC High Mass (ng)`) %>%
  gather(Low, Med, High, key="qc_level", value = "nominal_conc")

# Read in SCALiR predicted concentrations file
batch_df <- fread(batch_file)

# Combine QC sample nominal and predicted concentration data
batch_df <- batch_df %>%
  filter(grepl('QC', ms_file)) %>%
  mutate(qc_level = str_replace(ms_file, "^QC(L|M|H).+|^QC_(Low|Med|High)_.+", "\\1\\2")) %>%
  mutate(qc_level = case_when(
           qc_level == "L" ~ "Low",
           qc_level == "M" ~ "Med",
           qc_level == "H" ~ "High",
           .default = qc_level
         )) %>%
  left_join(qc_df, by=c("peak_label"="Compound", "qc_level")) %>%
  mutate(dif_from_nominal = abs(pred_conc - nominal_conc),
         percent_dif = (dif_from_nominal/nominal_conc)*100,
         smape = dif_from_nominal/((abs(pred_conc)+abs(nominal_conc))/2)*100)

# Precision & Accuracy estimates
batch_stats <- batch_df %>%
  select(ms_file, peak_label, pred_conc, qc_level, in_range, dif_from_nominal,
         percent_dif, nominal_conc, smape) %>%
  group_by(peak_label, qc_level) %>%
  summarize(in_range = mean(in_range),
            pred_conc_sd = sd(pred_conc),
            pred_conc_avg = mean(pred_conc),
            lower_limit = pred_conc_avg - 2*pred_conc_sd,
            upper_limit = pred_conc_avg + 2*pred_conc_sd,
            variance = pred_conc_sd/pred_conc_avg,
            percent_variance = round(variance*100, digits = 2),
            avg_dif = mean(dif_from_nominal),
            avg_percent_dif = round(mean(percent_dif), digits = 2),
            avg_smape = round(mean(smape), digits = 2),
            accuracy_check = (lower_limit <= nominal_conc) & (nominal_conc <= upper_limit))

# Make tables pretty for output
pretty_batch_df <- batch_df %>%
  select(ms_file, peak_label, qc_level, value, in_range,
         pred_conc, nominal_conc, dif_from_nominal, percent_dif, smape) %>%
  rename("Sample Name" = "ms_file",
         "Metabolite" = "peak_label",
         "Peak Ratio" = "value",
         "Predicted Concentration" = "pred_conc",
         "Valid" = "in_range",
         "QC Level" = "qc_level",
         "Nominal Concentration" = "nominal_conc",
         "Error" = "dif_from_nominal",
         "Percent Error" = "percent_dif",
         "SMAPE" = "smape")
fwrite(pretty_batch_df, file = paste0(outfile_prefix, "_qc_samples.csv"))

pretty_batch_stats <- batch_stats %>%
  select(peak_label, qc_level, in_range, pred_conc_avg, pred_conc_sd,
         variance, percent_variance, avg_dif, avg_percent_dif, avg_smape, accuracy_check) %>%
  rename("Metabolite" = "peak_label",
         "QC Level" = "qc_level",
         "Predicted Concentration Average" = "pred_conc_avg",
         "Predicted Concentration Standard Deviation" = "pred_conc_sd",
         "Variance" = "variance",
         "Percent Variance" = "percent_variance",
         "Average Error" = "avg_dif",
         "Average Percent Error" = "avg_percent_dif",
         "Average SMAPE" = "avg_smape",
         "Accuracy Check" = "accuracy_check")

fwrite(pretty_batch_stats, file = paste0(outfile_prefix, "_qc_statistics.csv"))

# Plot Precision/Variance estimate
var_plot <- batch_df %>%
  mutate(qc_level = factor(qc_level, levels = c("Low", "Med", "High"))) %>%
  ggplot(aes(x=((pred_conc)/(50/1000)), y=peak_label, fill=peak_label)) +
    geom_boxplot(varwidth = T) +
    #geom_point(aes(x=(nominal_conc*(50/1000)), y=peak_label)) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~qc_level, nrow=3, scales="free", strip.position = "top") +
    theme_bw(base_size = 14) +
    labs(title = "Predicted Concentration Variance",
         y = "Metabolite",
         x = "Predicted Concentration (ng/mL)") +
    theme(axis.title.x = element_text(margin = margin(b=15), vjust=-3),
          axis.title.y = element_text(margin = margin(l=15), vjust=6),
          legend.position = "none")

ggsave(paste0(outfile_prefix, "_precision_plot.jpg"),
       plot = var_plot,
       width=21.5,
       height=20,
       units="cm")

# Plot Accuracy estimate
acc_plot <- batch_df %>%
  mutate(qc_level = factor(qc_level, levels = c("Low", "Med", "High"))) %>%
  ggplot(aes(y=peak_label, x=smape, fill=peak_label)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~qc_level, nrow=3, scales = "free") +
    theme_bw(base_size = 14) +
    labs(title = "Predicted Concentration Accuracy",
         y = "Metabolite",
         x = "Symmetric Percent-Error") +
    theme(axis.title.x = element_text(margin = margin(b=15), vjust=-3),
          axis.title.y = element_text(margin = margin(l=15), vjust=6),
          legend.position = "none")

ggsave(paste0(outfile_prefix, "_accuracy_plot.jpg"),
       plot = acc_plot,
       width=16,
       height=20,
       units="cm")

# Plot of predicted concentrations with nominal concentrations
conc1 <- batch_df %>%
  mutate(qc_level = factor(qc_level, levels = c("Low", "Med", "High"))) %>%
  ggplot(aes(x=((pred_conc)/(50/1000)), y=peak_label, fill=peak_label)) +
    geom_boxplot(varwidth = T) +
    geom_point(aes(x=(nominal_conc/(50/1000)), y=peak_label)) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~qc_level, ncol=1, scales="free", strip.position = "top") +
    theme_bw(base_size = 14) +
    labs(title = "Predicted & Nominal Concentrations",
         y = "Metabolite",
         x = "Predicted Concentration (ng/mL)") +
    theme(axis.title.x = element_text(margin = margin(b=15), vjust=-3),
          axis.title.y = element_text(margin = margin(l=15), vjust=6),
          plot.margin = margin(10, 12, 10, 10),
          legend.position = "none")

ggsave(paste0(outfile_prefix, "_qc_plot_1.jpg"),
       plot = conc1,
       width=18,
       height=20,
       units="cm")

conc2 <- batch_df %>%
  mutate(qc_level = factor(qc_level, levels = c("Low", "Med", "High"))) %>%
  ggplot(aes(x=((pred_conc)/(50/1000)), y=qc_level, fill=peak_label)) +
    geom_boxplot(varwidth = T) +
    geom_point(aes(x=(nominal_conc/(50/1000)), y=qc_level)) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(~peak_label, nrow=3, scales="free",strip.position = "right") +
    theme_bw(base_size = 14) +
    labs(title = "Predicted & Nominal Concentrations",
         y = "QC Level",
         x = "Predicted Concentration (ng/mL)") +
    theme(axis.title.x = element_text(margin = margin(b=15), vjust=-3),
          axis.title.y = element_text(margin = margin(l=15), vjust=6),
          plot.margin = margin(10, 10, 10, 10),
          legend.position = "none")

ggsave(paste0(outfile_prefix, "_qc_plot_2.jpg"),
       plot = conc2,
       width=21.5,
       height=20,
       units="cm")

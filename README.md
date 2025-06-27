# SCALiR LCMS Nicotine Metabolomics Workflow

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Configuration](#configuration)
  - [Snakemake Configuration (](#snakemake-configuration-configconfigyaml)[`config/config.yaml`](#snakemake-configuration-configconfigyaml)[)](#snakemake-configuration-configconfigyaml)
  - [Default Profile (](#default-profile-profilesdefaultconfigv8yaml)[`profiles/default/config.v8+.yaml`](#default-profile-profilesdefaultconfigv8yaml)[)](#default-profile-profilesdefaultconfigv8yaml)
- [Usage](#usage)
  - [1. Clone the Repository](#1-clone-the-repository)
  - [2. Prepare Input Data & Configuration](#2-prepare-input-data--configuration)
  - [3. Install Miniforge (Conda)](#3-install-miniforge-conda)
  - [4. Create & Activate the Snakemake Environment](#4-create--activate-the-snakemake-environment)
  - [5. Running the Workflow](#5-running-the-workflow)
  - [Example Run](#example-run)
- [Workflow Rule Graph](#workflow-rule-graph)
- [Input Data Format](#input-data-format)
  - [Batch Results Files](#batch-results-files)
  - [QC Solutions File](#qc-solutions-file)
  - [Stock Solutions File](#stock-solutions-file)
- [Outputs Summary](#outputs-summary)
- [Citation](#citation)
- [Troubleshooting](#troubleshooting)
- [Issues & Collaborations](#issues--collaborations)

## Overview

**SCALiR** (Standard Curve Application for determining Linear Ranges) is an open-source, vendor-agnostic tool from the Lewis Research Group’s [ms-conc](https://github.com/LewisResearchGroup/ms-conc) package (see DOI: [10.1186/1471-2105-15-102](https://pmc.ncbi.nlm.nih.gov/articles/PMC10461962/)). Designed to eliminate the manual workload of metabolite quantification, SCALiR automatically fits standard curves from known concentration mixtures, computes the equation of best fit, and transforms LC–MS signal intensities into absolute concentrations. Key features include:

- **Rapid throughput**: processes over 75 compounds in under 10 seconds.
- **High accuracy**: achieves R² = 0.99 against expected concentrations in a 77-metabolite mix.
- **Reproducibility**: average coefficient of variation (CV) of 0.091 across ten analytical batches.
- **Calibration flexibility**: supports both external and internal standard correction for matrix effects.
- **Deterministic limits**: evaluates quantification ranges automatically.

In this workflow, the original web application has been adapted into a CLI (`scalir_cli.py`) for seamless integration into Snakemake.

The pipeline automates the processing of raw LCMS data from the CU Anschutz Medical Campus Core Facilities to quantify nicotine metabolites. It:

- Converts raw Excel outputs into SCALiR-compatible input files.
- Builds calibration curves and predicts concentrations using the SCALiR CLI tool.
- Performs quality assessment based on QC samples by batch.
- Aggregates and filters results across batches.
- Calculates the Nicotine Metabolite Ratio (NMR) for valid samples.
- **Combined Results**: If more than one `model_fit` is specified, the workflow outputs additional “combined” results. These identify the best-fit model for each batch × metabolite × concentration level (low, medium, high) based on SMAPE, and report the optimal-fit concentration for every sample.

## Repository Structure

```text
├── config
│   └── config.yaml
├── data
│   ├── raw_results_batch1.xlsx
│   ├── raw_results_batch2.xlsx
│   ├── raw_results_batch3.xlsx
│   ├── qc_solutions.xlsx
│   └── stock_solutions.xlsx
├── envs
│   ├── python.yaml
│   ├── r_env.yaml
│   └── scalir_env.yaml
├── profiles
│   └── default
│       └── config.v8+.yaml
├── README.md
├── results
│   ├── combined
│   ├── fixed
│   ├── interval
│   ├── run_config.yaml
│   └── wide
├── scripts
│   ├── calculate_nmr.py
│   ├── combine_results.py
│   ├── make_IS_concentration_file.R
│   ├── make_peakslist_file.R
│   ├── ms_conc/
│   ├── quality_assessment.R
│   ├── report_template.rmd
│   ├── results_comparison.rmd
│   ├── save_yaml.py
│   ├── scalir_cli.py
│   └── splice_optimal_results.py
├── Snakefile
└── utils.py
```

## Requirements

- **Conda** (e.g., Miniconda or Anaconda)
- **Snakemake** (≥9.3.3)

> All other software dependencies are installed automatically via the Conda environment YAMLs in the `envs/` directory (see `envs/*.yaml`).

## Configuration

### Snakemake Configuration (`config/config.yaml`)

This Snakemake config file is a YAML document where you define workflow parameters—such as input/output paths, model options, and QC criteria—separately from the Snakefile. Snakemake reads these settings at runtime to generate rules, manage dependencies, and produce outputs without editing pipeline code.

Key parameters:

- `standards_stock_file`: Path to the stock solutions file (Excel) containing known standard concentrations for calibration. fileciteturn3file1
- `lcms_results`: List of raw LCMS result files (Excel) to be processed. fileciteturn3file1
- `qc_file`: Path to the QC sample concentration file (Excel) used in the QC assessment step. fileciteturn3file1
- `model_fit`: Regression mode for calibration curves. Options:
  - `fixed`: simple linear fit
  - `interval`: confidence interval bands
  - `wide`: prediction interval bands fileciteturn3file1
- `percent_error_threshold`: SMAPE threshold (Symmetric Mean Absolute Percentage Error) for retaining valid samples. SMAPE calculates the average, symmetric, relative error between observed and predicted concentrations. Commonly used thresholds are **20%** (stricter, higher accuracy but fewer samples) or **25%** (more lenient, more samples retained). fileciteturn3file1
- `out_dir`: Base output directory where all SCALiR results and QC summaries will be written. fileciteturn3file1

### Default Profile (`profiles/default/config.v8+.yaml`)

This profile is configured for local execution and automatically applies key Snakemake options without needing to specify them on the command line.

```yaml
software-deployment-method: conda
configfile: config/config.yaml
rerun-incomplete: True
use-conda: True
conda-prefix: ../.snakemake
cores: 3
```

- **software-deployment-method: conda**\
  Use Conda to create and manage environments for each rule.
- **configfile**\
  Path to the Snakemake config file (`config/config.yaml`).
- **rerun-incomplete: True**\
  Automatically rerun any jobs that were left incomplete in previous runs.
- **use-conda: True**\
  Enable Conda environment creation by default (equivalent to `--use-conda`).
- **conda-prefix**\
  Directory where Conda environments will be stored (relative to the repo root).
- **cores**\
  Number of parallel jobs to run locally. **Adjust this value** to match the number of CPU cores you have available.

## Usage

### 1. Clone the Repository

```bash
git clone https://github.com/krferrier/nmr_snakemake.git
cd nmr_snakemake
```

### 2. Prepare Input Data & Configuration

Copy all necessary raw input files (Excel) into the `data/` directory. Then open `config/config.yaml` and modify the following fields to match your run:

```yaml
lcms_results:
  - data/your_batch1.xlsx
  - data/your_batch2.xlsx
qc_file: data/your_qc_solutions.xlsx
standards_stock_file: data/your_stock_solutions.xlsx
model_fit:
  - fixed
  # or other fits: interval, wide
percent_error_threshold: 25
out_dir: results
```

### 3. Install Miniforge (Conda)

If you don’t already have Conda, install **Miniforge** for **Linux x86\_64**. For other platforms, see the official guide: [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
```

### 4. Create & Activate the Snakemake Environment

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake==9.3.3 python=3.8
conda activate snakemake
```

### 5. Running the Workflow

Simply run the pipeline using the default profile:

```bash
snakemake
```

To preview the execution plan without running jobs (dry-run):

```bash
snakemake -n
```

### Example Run

```bash
$ snakemake
Using workflow specific profile profiles/default for setting default command line arguments.
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 3
Rules claiming more threads will be scaled down.
Job stats:
job                                    count
-----------------------------------  -------
all                                        1
calculate_nmr                              3
calculate_nmr_optimal_results              1
combine_optimal_results                    1
combine_results                            3
create_scalir_IS_concentration_file       33
create_scalir_peaklist_file               33
run_quality_assessment                    33
run_scalir                                33
save_run_config_and_paths_df               1
total                                    142

Select jobs to execute...
Execute 3 jobs...

[Date Time]
localrule create_scalir_IS_concentration_file:
    input: data/stock_solutions.xlsx
    output: results/fixed/raw_results_batch2/scalir_input/internal_standards_concentration_file.csv
    jobid: 17
    reason: Missing output files: results/fixed/raw_results_batch2/scalir_input/internal_standards_concentration_file.csv
    wildcards: fit=fixed, batch_name=raw_results_batch2
    resources: tmpdir=/tmp

[Date Time]
localrule create_scalir_peaklist_file:
    input: data/raw_results_batch3.xlsx
    output: results/wide/raw_results_batch3/scalir_input/raw_results_batch3_peaklist.csv
    jobid: 72
    reason: Missing output files: results/wide/raw_results_batch3/scalir_input/raw_results_batch3_peaklist.csv
    wildcards: fit=wide, batch_name=raw_results_batch3
    resources: tmpdir=/tmp

# ... rest of jobs execute similarly ...
[Date Time]
All jobs done.
```

### Workflow Rule Graph



## Input Data Format

### Batch Results Files

- **File naming**: any Excel file in `data/` containing the word `results` and a date (e.g., `NIC_MET_QUANT_RESULTS_YYYYMMDD.xlsx`).
- **Sheet**: first sheet of the workbook.
- **Headers**: two header rows:
  1. **Analyte columns**: for each compound, a pair of “Method” and “Results” columns, plus corresponding “Qualifier” columns for quantifier and qualifier transitions.
  2. **Metadata columns**: `Sample`, `Name`, `Data File`, `Type`, `Level`, `Acq. Date-Time`, `Pos.`, followed by repeating groups of `Precursor Ion`, `Product Ion`, `Area`, `RT`, `Ratio` for each analyte.

**Example Batch Results (first sheet)**



### QC Solutions File

- **Filename**: one file matching `QCSolutions_*.xlsx` in `data/`.
- **Sheet**: first sheet.
- **QC table**: begins at line 20 and includes columns:
  - `Compound`
  - `Abbreviation`
  - `Plasma Amount (uL)`
  - `QC Low Conc. (ng/mL)`
  - `QC Med Conc. (ng/mL)`
  - `QC High Conc. (ng/mL)`
  - `QC Low Mass (ng)`
  - `QC Med Mass (ng)`
  - `QC High Mass (ng)`

**Anonymized example** of the QC table (starting at row 20):

| Compound          | Abbreviation | Plasma Amount (uL) | QC Low Conc. (ng/mL) | QC Med Conc. (ng/mL) | QC High Conc. (ng/mL) | QC Low Mass (ng) | QC Med Mass (ng) | QC High Mass (ng) |
| ----------------- | ------------ | ------------------ | -------------------- | -------------------- | --------------------- | ---------------- | ---------------- | ----------------- |
| nicotine          | NIC          | 50                 | 4.99                 | 45.82                | 252.00                | 0.25             | 2.29             | 12.60             |
| nornicotine       | NORN         | 50                 | 5.16                 | 47.35                | 260.40                | 0.26             | 2.37             | 13.02             |
| 3-hydroxycotinine | 3HC          | 50                 | 4.60                 | 42.27                | 232.50                | 0.23             | 2.11             | 11.63             |

### Stock Solutions File

- **Filename**: one file matching `STOCKSOLUTION_CONCENTRATIONS_*.xlsx` in `data/`.
- **Sheet**: first sheet.
- **Standards table**: data begins at the top of the sheet with columns including:
  - `Compound`
  - `Type` (e.g., Light/Heavy)
  - `Date Stock Made`
  - `Batch Stock` ID
  - `Volume Made (mL)`
  - `Exact Final Concentration (µg/mL)`
  - `Solvent` (e.g., methanol)
  - `Target Concentration (ng/µL)`

**Anonymized example** of the top of this table:

| Compound  | Type       | Date Stock Made | Batch Stock | Volume Made (mL) | Exact Final Conc. (µg/mL) | Solvent  | Target Conc. (ng/µL) |
| --------- | ---------- | --------------- | ----------- | ---------------- | ------------------------- | -------- | -------------------- |
| compoundA | Light (WS) | 2022-06-03      | 001         | 1                | 1260                      | methanol | 0.5                  |
| compoundB | Light (WS) | 2022-06-03      | 001         | 1                | 2170                      | methanol | 0.5                  |

## Outputs Summary

All files are written to the `out_dir` specified in `config.yaml`:

| File                          | Description                                     |
| ----------------------------- | ----------------------------------------------- |
| `<prefix>_parameters.csv`     | SCALiR regression parameters                    |
| `<prefix>_concentrations.csv` | Raw concentration predictions                   |
| `optimal_conc.tsv`            | Optimal concentrations (pre-grouping)           |
| `all_results.tsv`             | Aggregated results for all (samples + controls) |
| `all_samples.tsv`             | Aggregated target samples                       |
| `valid_samples.tsv`           | In-range samples                                |
| `valid_samples_filtered.tsv`  | In-range, error-filtered samples                |
| `failed_samples.tsv`          | Samples failing QC                              |
| `control_samples.tsv`         | Control sample summaries                        |
| `valid_samples_nmr.tsv`       | NMR and log(NMR) for valid samples              |

## Citation

If you use SCALiR in your research, please cite:

> Bray NL, Pimentel H, Melsted P, Pachter L. **SCALiR: Standard Curve Application for Determining Linear Ranges.** *BMC Bioinformatics*. 2014;15:102. doi:10.1186/1471-2105-15-102

To reference this Snakemake workflow, please cite:

> Ferrier KR. **nmr\_snakemake: A Snakemake Pipeline for SCALiR-based Nicotine Metabolomics.** GitHub repository, 2025. [https://github.com/krferrier/nmr\_snakemake](https://github.com/krferrier/nmr_snakemake)

## Troubleshooting

Below are common issues and suggested resolutions:

- **Conda environment failures**: Ensure `conda` is on your `PATH` and strict channel priority is enabled (`conda config --set channel_priority strict`). Remove any partially created envs under `envs/` and rerun.
- **Missing input files**: Verify that your `config/config.yaml` paths match the filenames in `data/` and that Excel sheets are named correctly.
- **Permission errors**: Check write permissions on `results/` and the `conda-prefix` directory (default: `../.snakemake`).
- **Unexpected rule graph**: Run `snakemake -n` to debug missing dependencies or misnamed wildcards.
- **SCALiR errors**: Examine the logs under each `scalir_input/` folder for detailed error messages from the CLI.

If you encounter other issues, please refer to the **Issues** section below.

## Issues & Collaborations

For bug reports, feature requests, or collaboration inquiries, open an issue on the GitHub repository: [https://github.com/krferrier/nmr\_snakemake/issues](https://github.com/krferrier/nmr_snakemake/issues)

When submitting an issue, please include:

1. **Description** of the problem or request.
2. **Snakemake version** (`snakemake --version`).
3. **Workflow profile** used (contents of `profiles/default/config.v8+.yaml`).
4. **Config file** (`config/config.yaml`) or relevant excerpts.
5. **Exact command** you ran and the **full error message** or log output.
6. **Minimal reproducible example**, if possible (e.g., small dataset and config snippet).

Contributions are welcome! Feel free to fork the repository, submit pull requests, or reach out via GitHub issues to discuss enhancements or collaborations.

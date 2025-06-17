# nmr_snakemake
Snakemake workflow for processing raw LC/MS nicotine metabolite data with SCALiR to get metabolite concentrations.

## Overview

This pipeline processes liquid chromatography-mass spectrometry (LC/MS) data for nicotine metabolites using the SCALiR (Standard Curve and Linearity in R) tool. The workflow processes raw LC/MS data through multiple steps to calculate nicotine metabolite concentrations and the nicotine metabolite ratio (NMR), which is used as a biomarker for nicotine metabolism rate.

## Features

- Automated processing of multiple LC/MS data batches
- Support for different regression model fitting methods (fixed, interval, wide)
- Quality control assessment with accuracy and precision metrics
- Generation of concentration calibration curves and QC plots
- Calculation of nicotine metabolite ratio (NMR) from 3-hydroxycotinine and cotinine concentrations
- Comprehensive filtering of samples based on quality control metrics

## Requirements

### Software Dependencies

- Snakemake (v9.3.3 or higher)

All other dependencies are automatically managed through Conda environments defined in the workflow.

### Input Data

- Raw LC/MS nicotine metabolite data files (.xlsx format)
- Standards concentration file with internal standard information
- Quality control (QC) solutions file

## Workflow Structure

### Main Steps

1. **Data Preparation**:
   - Creation of peaklist files from raw LC/MS data
   - Generation of internal standard concentration files

2. **SCALiR Processing**:
   - Running SCALiR with different fitting methods (fixed, interval, wide)
   - Generation of calibration curves and concentration estimates

3. **Quality Assessment**:
   - Evaluation of accuracy and precision using QC samples
   - Filtering of results based on quality metrics

4. **Results Aggregation**:
   - Combining results from multiple batches
   - Separation of valid and failed samples

5. **NMR Calculation**:
   - Calculating the nicotine metabolite ratio (3-hydroxycotinine/cotinine)
   - Output of final results with NMR values

### Output Files

- Concentration estimates for each batch and fitting method
- Quality control statistics and plots
- Combined results across all batches
- Filtered valid samples with metabolite concentrations
- NMR calculations for valid samples

## Configuration

The workflow is configured using a YAML file (`config/config.yaml`) that specifies:

- Input data files (LC/MS results, standards, QC samples)
- SCALiR parameters including model fitting methods
- Output directory structure
- Quality filtering thresholds

## Usage

1. Configure the input data paths in `config/config.yaml`
2. Run the workflow using Snakemake:

```bash
snakemake --cores [N] --use-conda
```

Where `[N]` is the number of cores to use for parallel processing. The `--use-conda` flag ensures all required dependencies are automatically installed in isolated environments.

3. Results will be organized in the specified output directory with subdirectories for each fitting method.

## Directory Structure

```
├── config/              # Configuration files
├── envs/                # Conda environment files
├── scripts/             # Processing scripts
│   ├── ms_conc/         # SCALiR core modules
│   ├── scalir_cli.py    # SCALiR command-line interface
│   └── ...              # Other processing scripts
├── Snakefile            # Main workflow definition
└── utils.py            # Utility functions for configuration
```

## Citation

[Include citation information for SCALiR and the workflow if applicable]

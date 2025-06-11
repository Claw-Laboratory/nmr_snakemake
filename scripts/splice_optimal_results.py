import pandas as pd
import numpy as np

def main(conc_file_paths, qc_conc_file, qc_stats_file_paths, all_samples_results_path,
         valid_samples_results_path, failed_samples_results_path,
         control_samples_results_path):

    # Create map of breaks from the nominal QC concentrations
    qc_df = pd.read_excel(qc_file, skiprows=20, header=0)
    qc_df['Compound'] = qc_df['Compound'].str.lower()
    qc_df = qc_df[['Compound', 'QC Low Mass (ng)', 'QC Med Mass (ng)', 'QC High Mass (ng)']]
    qc_df = qc_df.rename(columns={
        'QC Low Mass (ng)': 'Low',
        'QC Med Mass (ng)': 'Med',
        'QC High Mass (ng)': 'High'
    })
    qc_breaks = qc_df.copy()
    qc_breaks["b1"] = (qc_breaks["Low"] * qc_breaks["Med"])**0.5
    qc_breaks["b2"] = (qc_breaks["Med"] * qc_breaks["High"])**0.5
    qc_breaks = qc_breaks[["metabolite","b1","b2"]]

    # Build best-model table

    # Merge concentration files

    # Join breaks into sample table

    # Bin each sample to low/med/high by metabolite

    # Join best-model table on metabolite x level x batch

    # Select final concentrations

if __name__ == '__main__':
    main(snakemake.input["conc_files"],
         snakemake.input["qc_conc_file"],
         snakemake.input["qc_stats_files"],
         snakemake.output["all_results_file"],
         snakemake.output["all_samples_file"],
         snakemake.output["valid_samples_file"],
         snakemake.output["failed_samples_file"],
         snakemake.output["control_samples_file"]
         )
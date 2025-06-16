import pandas as pd
import numpy as np
import os

def main(list_of_batch_paths, qc_conc_file, qc_stats_file_paths, threshold, all_results_path, all_samples_results_path,
         valid_samples_results_path,  valid_samples_filtered_results_path, failed_samples_results_path,
         control_samples_results_path):
    # Create map of breaks from the nominal QC concentrations
    qc_df = pd.read_excel(qc_conc_file, skiprows=20, header=0)
    qc_df['Compound'] = qc_df['Compound'].str.lower()
    qc_df = qc_df[['Compound', 'QC Low Mass (ng)', 'QC Med Mass (ng)', 'QC High Mass (ng)']]
    qc_df = qc_df.rename(columns={
        'QC Low Mass (ng)': 'Low',
        'QC Med Mass (ng)': 'Med',
        'QC High Mass (ng)': 'High'
    })
    qc_breaks = qc_df.copy()
    qc_breaks["b1"] = (qc_breaks["Low"] * qc_breaks["Med"]) ** 0.5
    qc_breaks["b2"] = (qc_breaks["Med"] * qc_breaks["High"]) ** 0.5
    qc_breaks = qc_breaks[["Compound", "b1", "b2"]]
    qc_breaks = qc_breaks.rename(columns={"Compound": "peak_label"})
    # Build best-model table
    qc_stats_df = []
    for fp in qc_stats_file_paths:
        _, fit, batch, _, fname = fp.split(os.sep)
        df = pd.read_csv(fp, usecols=["Metabolite", "QC Level", "in_range", "Average SMAPE"])
        df["batch"] = batch
        qc_stats_df.append(df)

    fit_df = pd.concat(qc_stats_df, ignore_index=True)
    fit_df.rename(columns={"QC Level": "level", "Metabolite": "peak_label"}, inplace=True)

    # Merge concentration files
    conc_df = []
    for fp in list_of_batch_paths:
        _, fit, batch, fname = fp.split(os.sep)
        df = pd.read_csv(fp)
        df["batch"] = batch
        conc_df.append(df)

    conc_combined_df = pd.concat(conc_df, ignore_index=True)

    # Join breaks into sample table
    conc_combined_df = conc_combined_df.merge(qc_breaks, how="left", on="peak_label")

    # Bin each sample to low/med/high by metabolite
    def assign_level(row):
        x, b1, b2 = row["pred_conc"], row["b1"], row["b2"]
        if x <= b1:
            return "Low"
        elif x <= b2:
            return "Med"
        else:
            return "High"

    conc_combined_df = conc_combined_df.assign(level=lambda d: d.apply(assign_level, axis=1))

    # Join best-model table on metabolite x level x batch
    fit_df = fit_df.drop(columns=["in_range"])
    conc_combined_df = conc_combined_df.merge(
        fit_df,
        on=["peak_label", "level", "batch"],
        how="left"
    )

    conc_combined_df = conc_combined_df.drop(columns=["b1", "b2"])

    ##### FILTER & FORMAT OUTPUT FILES #####
    filtered_df = conc_combined_df.loc[~conc_combined_df['ms_file'].str.contains('std')]
    filtered_df = filtered_df.loc[~filtered_df['ms_file'].str.contains('QC')]
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'^(\d\d\d\d)BW.+|^(\d\d\d\d)BW' : r'\1\2'}, regex=True)
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'^SHS-(.+)': r'\1'}, regex=True)
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'(.+)_RR$': r'\1'}, regex=True)
    filtered_df['conc_ng_mL'] = filtered_df['pred_conc']/(50/1000) #predicted concentration is in ng/50uL, we want to report ng/mL
    filtered_df['type'] = np.where(filtered_df['ms_file'].str.len() == 4, 'control', 'sample')
    all_grouped = (filtered_df
        .groupby(['ms_file', 'peak_label'], as_index=False)
        .agg(
            sample_id = ('ms_file', 'first'),
            peak_ratio = ('value', 'mean'),
            pred_conc = ('pred_conc', 'mean'),
            conc_ng_mL = ('conc_ng_mL', 'mean'),
            in_range = ('in_range', 'mean'),
            type = ('type', 'first'),
            n = ('ms_file', 'size'),
            level = ('level', 'first'),
            average_smape = ('Average SMAPE', 'mean')
        )
    )
    all_grouped.to_csv(all_results_path, index=False, header=True, sep='\t')

    # Subset to df of samples only
    samples_only = filtered_df.loc[filtered_df['type']=='sample']
    samples_grouped = (samples_only
        .groupby(['ms_file', 'peak_label'], as_index=False)
        .agg(
            sample_id = ('ms_file', 'first'),
            peak_ratio = ('value', 'mean'),
            pred_conc = ('pred_conc', 'mean'),
            conc_ng_mL=('conc_ng_mL', 'mean'),
            in_range = ('in_range', 'mean'),
            type = ('type', 'first'),
            n = ('ms_file', 'size'),
            level=('level', 'first'),
            average_smape=('Average SMAPE', 'mean')
        )
    )
    samples_grouped.to_csv(all_samples_results_path, index=False, header=True, sep='\t')

    # Subset of samples that are out of range
    out_of_range = samples_only.loc[samples_only['in_range']==0]

    # Subset of samples that are in-range
    in_range = samples_only.loc[samples_only['in_range']==1]

    # # Find the intersect of the two dataframes
    # repeated_samples = pd.merge(in_range, out_of_range,
    #                             how ='inner',
    #                             on=['ms_file', 'peak_label'])
    # repeated_samples.to_csv("successfully_repeated_samples.csv", index=False, header=True, sep='\t')

    # Identify samples that are out of range and have not had a successful repeat
    failed_samples = out_of_range[~out_of_range[['ms_file', 'peak_label']].apply(tuple, axis=1).isin(in_range[['ms_file', 'peak_label']].apply(tuple, axis=1))]
    failed_samples_grouped = (failed_samples
        .groupby(['ms_file', 'peak_label'], as_index=False)
        .agg(
            sample_id = ('ms_file', 'first'),
            peak_ratio = ('value', 'mean'),
            pred_conc = ('pred_conc', 'mean'),
            conc_ng_mL=('conc_ng_mL', 'mean'),
            in_range = ('in_range', 'mean'),
            type = ('type', 'first'),
            n = ('ms_file', 'size'),
            level=('level', 'first'),
            average_smape=('Average SMAPE', 'mean')
        )
    )
    failed_samples_grouped.to_csv(failed_samples_results_path, index=False, header=True, sep='\t')

    # If there are multiple in_range samples, report the mean value and pred_conc
    in_range_samples_grouped = (in_range
        .groupby(['ms_file', 'peak_label'], as_index=False)
        .agg(
            sample_id = ('ms_file', 'first'),
            peak_ratio = ('value', 'mean'),
            pred_conc = ('pred_conc', 'mean'),
            conc_ng_mL=('conc_ng_mL', 'mean'),
            in_range = ('in_range', 'mean'),
            type = ('type', 'first'),
            n = ('ms_file', 'size'),
            level=('level', 'first'),
            average_smape=('Average SMAPE', 'mean')
        )
    )
    in_range_samples_grouped.to_csv(valid_samples_results_path, index=False, header=True, sep='\t')

    # Filter in-range samples to those with an average SMAPE below the cutoff threshold
    in_range_filtered_samples_grouped = in_range_samples_grouped.loc[in_range_samples_grouped['average_smape'] <= threshold].copy()
    in_range_filtered_samples_grouped.to_csv(valid_samples_filtered_results_path, index=False, header=True, sep='\t')

    # Subset to df of controls only
    controls_only = filtered_df.loc[filtered_df['type'] == 'control'].copy()
    pos_controls = ["2727", "3819", "3089", "4753", "5487"]
    controls_only['control_type'] = np.where(
        (controls_only['type'] == 'control') & (controls_only['ms_file'].isin(pos_controls)),
        'pos', 'neg')
    controls_grouped = (controls_only
        .groupby(['ms_file', 'peak_label'], as_index=False)
        .agg(
            sample_id = ('ms_file', 'first'),
            peak_ratio = ('value', 'mean'),
            pred_conc = ('pred_conc', 'mean'),
            conc_ng_mL=('conc_ng_mL', 'mean'),
            in_range = ('in_range', 'mean'),
            type = ('type', 'first'),
            control_type = ('control_type', 'first'),
            n = ('ms_file', 'size'),
            level=('level', 'first'),
            average_smape=('Average SMAPE', 'mean')
        )
    )
    controls_grouped.to_csv(control_samples_results_path, index=False, header=True, sep='\t')

if __name__ == '__main__':
    main(snakemake.input["paths"],
         snakemake.input["qc_conc_file"],
         snakemake.input["qc_stats_files"],
         snakemake.params["threshold"],
         snakemake.output["all_results_file"],
         snakemake.output["all_samples_file"],
         snakemake.output["valid_samples_file"],
         snakemake.output["valid_samples_filtered_file"],
         snakemake.output["failed_samples_file"],
         snakemake.output["control_samples_file"]
         )
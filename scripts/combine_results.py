import pandas as pd
import numpy as np
import os

def main(list_of_batch_paths, qc_conc_file, qc_stats_file_paths, threshold, all_results_path, all_samples_results_path,
         valid_samples_results_path,  valid_samples_filtered_results_path, failed_samples_results_path,
         control_samples_results_path):
    """
    Main function for processing QC data and concentration files, performing statistical calculations,
    binning concentration levels, and generating filtered and formatted outputs based on criteria.

    :param list_of_batch_paths: List of file paths containing batch concentration data.
    :type list_of_batch_paths: list[str]
    :param qc_conc_file: Path to the QC concentration file in Excel format.
    :type qc_conc_file: str
    :param qc_stats_file_paths: List of file paths for QC statistics CSV files.
    :type qc_stats_file_paths: list[str]
    :param threshold: SMAPE threshold for filtering in-range samples with acceptable accuracy.
    :type threshold: float
    :param all_results_path: Destination file path for saving all grouped result data.
    :type all_results_path: str
    :param all_samples_results_path: Destination file path for saving all samples grouped result data.
    :type all_samples_results_path: str
    :param valid_samples_results_path: Destination file path for saving valid samples' grouped results data.
    :type valid_samples_results_path: str
    :param valid_samples_filtered_results_path: Destination file path for saving valid samples,
        filtered based on SMAPE threshold.
    :type valid_samples_filtered_results_path: str
    :param failed_samples_results_path: Destination file path for saving grouped results of failed samples.
    :type failed_samples_results_path: str
    :param control_samples_results_path: Destination file path for saving grouped results of control samples.
    :type control_samples_results_path: str
    :return: None
    """
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

    # Build model-fit QC table
    qc_stats_df = []
    for fp in qc_stats_file_paths:
        _, fit, batch, _, fname = fp.split(os.sep)
        df = pd.read_csv(fp, usecols=["Metabolite", "QC Level", "in_range", "Average SMAPE", "Accuracy Check"])
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

    # Join qc_breaks for assigning concentration level categories to the combined dataframe
    conc_combined_df = conc_combined_df.merge(qc_breaks, how="left", on="peak_label")

    # Bin each sample to low/med/high concentration level categories by metabolite
    def assign_level(row):
        x, b1, b2 = row["pred_conc"], row["b1"], row["b2"]
        if x <= b1:
            return "Low"
        elif x <= b2:
            return "Med"
        else:
            return "High"

    # Add a new column called level to the combined dataframe and perform the assign_level function to each row
    conc_combined_df = conc_combined_df.assign(level=lambda d: d.apply(assign_level, axis=1))

    # Join combined dataframe to the dataframe with qc_stats by metabolite x level x batch
    fit_df = fit_df.drop(columns=["in_range"])
    conc_combined_df = conc_combined_df.merge(
        fit_df,
        on=["peak_label", "level", "batch"],
        how="left"
    )
    # Drop the break-point columns
    conc_combined_df = conc_combined_df.drop(columns=["b1", "b2"])

    ##### FILTER & FORMAT OUTPUT FILES #####
    filtered_df = conc_combined_df.loc[~conc_combined_df['ms_file'].str.contains('std')]
    filtered_df = filtered_df.loc[~filtered_df['ms_file'].str.contains('QC')]
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'^(\d\d\d\d)BW.+|^(\d\d\d\d)BW' : r'\1\2'}, regex=True)
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'^SHS-(.+)': r'\1'}, regex=True)
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'(.+)_RR$': r'\1'}, regex=True)
    filtered_df['conc_ng_mL'] = filtered_df['pred_conc']/(50/1000) #predicted concentration is in ng/50uL, we want to report ng/mL
    filtered_df['type'] = np.where(filtered_df['ms_file'].str.len() == 4, 'control', 'sample')

    common_agg = {
        'ms_file': 'first',
        'value': 'mean',
        'pred_conc': 'mean',
        'conc_ng_mL': 'mean',
        'in_range': 'mean',
        'Accuracy Check': 'mean',
        'type': 'first',
        'level': 'first',
        'Average SMAPE': 'mean'
    }

    def group_and_save(df, output_path, include_control_type=False):
        """Helper function to group and save dataframes"""
        # Define aggregation dictionary based on available columns
        agg_dict = {
            'ms_file': 'first',
            'value': 'mean',
            'pred_conc': 'mean',
            'conc_ng_mL': 'mean',
            'in_range': 'mean',
            'Accuracy Check': 'mean',
            'type': 'first',
            'level': 'first',
            'Average SMAPE': 'mean'
        }

        # Add control_type if needed
        if include_control_type and 'control_type' in df.columns:
            agg_dict['control_type'] = 'first'

        grouped = df.groupby(['ms_file', 'peak_label'], as_index=False).agg(agg_dict)

        # Add count column
        grouped['n'] = df.groupby(['ms_file', 'peak_label']).size().values

        # Rename columns to match original output
        grouped = grouped.rename(columns={
            'ms_file': 'sample_id',
            'value': 'peak_ratio',
            'Accuracy Check': 'accuracy_check',
            'Average SMAPE': 'average_smape'
        })

        grouped.to_csv(output_path, index=False, header=True, sep='\t')
        return grouped

    # Process all results including BloodWorks control samples and target samples
    all_grouped = group_and_save(filtered_df, all_results_path)

    # Process target samples
    samples_only = filtered_df.loc[filtered_df['type'] == 'sample']
    samples_grouped = group_and_save(samples_only, all_samples_results_path)

    # Process in-range and out-of-range target samples
    in_range = samples_only.loc[samples_only['in_range'] == 1]
    out_of_range = samples_only.loc[samples_only['in_range'] == 0]

    # Process failed target samples
    failed_samples = out_of_range[~out_of_range[['ms_file', 'peak_label']].apply(tuple, axis=1)
                                .isin(in_range[['ms_file', 'peak_label']].apply(tuple, axis=1))]
    group_and_save(failed_samples, failed_samples_results_path)

    # Process in-range target samples
    in_range_samples_grouped = group_and_save(in_range, valid_samples_results_path)

    # Filter and save target samples below SMAPE threshold
    in_range_filtered = in_range_samples_grouped.loc[in_range_samples_grouped['average_smape'] <= threshold]
    in_range_filtered.to_csv(valid_samples_filtered_results_path, index=False, header=True, sep='\t')

    # Process BloodWorks control samples
    controls_only = filtered_df.loc[filtered_df['type'] == 'control'].copy()
    pos_controls = ["2727", "3819", "3089", "4753", "5487"]
    controls_only['control_type'] = np.where(
        (controls_only['type'] == 'control') & (controls_only['ms_file'].isin(pos_controls)),
        'pos', 'neg')
    group_and_save(controls_only, control_samples_results_path, include_control_type=True)


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
import pandas as pd
import numpy as np
import os

def main(conc_file_paths, qc_conc_file, qc_stats_file_paths, threshold, optimal_conc_file,
         all_results_path, all_samples_results_path, valid_samples_results_path, valid_samples_filtered_results_path,
         failed_samples_results_path, control_samples_results_path
         ):
    """
    Processes and analyzes concentration and statistics data to generate multiple result files, 
    including optimal concentrations, valid and failed sample results, and control sample summaries. 
    This function implements a workflow that uses provided concentration files, QC concentration files, 
    and QC statistics files, applies transformation and filtering, and outputs structured results 
    based on the thresholds specified.

    :param conc_file_paths: List of file paths containing predicted concentration data for analysis.
    :type conc_file_paths: list of str
    :param qc_conc_file: File path of the QC concentrations file in '.xlsx' format for QC levels information.
    :type qc_conc_file: str
    :param qc_stats_file_paths: List of file paths containing quality control statistics data.
    :type qc_stats_file_paths: list of str
    :param threshold: Specified percent error threshold value used for filtering results.
    :type threshold: float
    :param optimal_conc_file: File path to save the computed optimal concentrations as output.
    :type optimal_conc_file: str
    :param all_results_path: File path to save the aggregated results of all samples.
    :type all_results_path: str
    :param all_samples_results_path: File path to save results of all analyzed sample data only.
    :type all_samples_results_path: str
    :param valid_samples_results_path: File path to save results of valid, in-range samples only.
    :type valid_samples_results_path: str
    :param valid_samples_filtered_results_path: File path to save filtered results of valid samples.
    :type valid_samples_filtered_results_path: str
    :param failed_samples_results_path: File path where samples that fail validation will be saved.
    :type failed_samples_results_path: str
    :param control_samples_results_path: File path where control sample summary data will be saved.
    :type control_samples_results_path: str
    :return: Outputs are saved to respective file paths. Does not return a value.
    :rtype: None
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
    qc_breaks["b1"] = (qc_breaks["Low"] * qc_breaks["Med"])**0.5
    qc_breaks["b2"] = (qc_breaks["Med"] * qc_breaks["High"])**0.5
    qc_breaks = qc_breaks[["Compound","b1","b2"]]
    qc_breaks = qc_breaks.rename(columns={"Compound": "peak_label"})
    # Build best-model-fit QC table
    qc_stats_df = []
    for fp in qc_stats_file_paths:
        _, fit, batch, _, fname = fp.split(os.sep)
        df = pd.read_csv(fp, usecols=["Metabolite", "QC Level", "in_range", "Average SMAPE"])
        df["fit"] = fit
        df["batch"] = batch
        qc_stats_df.append(df)

    best_fit_df = pd.concat(qc_stats_df, ignore_index=True)
    # Sort the combined qc dataframe by Average SMAPE (percent error)
    best_fit_df = best_fit_df.sort_values("Average SMAPE")
    # Group the data by batch, metabolite, and qc level, then keep only the first occurence.
    # This will create a dataframe of the lowest percent error/highest accuracy model fit for each batch x metabolite x level
    best_fit_df = (
        best_fit_df
        .drop_duplicates(subset=["batch", "Metabolite", "QC Level"], keep="first")
        .reset_index(drop=True)
    )
    # Rename columns for clarity and ease of merging with the predicted concentration files
    best_fit_df.rename(columns={"fit": "best_fit", "QC Level": "level", "Metabolite": "peak_label"}, inplace=True)

    # Merge concentration files
    conc_df = []
    for fp in conc_file_paths:
        _, fit, batch, fname = fp.split(os.sep)
        df = pd.read_csv(fp)
        df["fit"] = fit
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

    df_levels = (
        conc_combined_df[conc_combined_df["fit"] == "fixed"]
        .copy()
        # Add a new column called level and perform the assign_level function to each row
        .assign(level=lambda d: d.apply(assign_level, axis=1))
        [["ms_file", "peak_label", "batch", "level"]]
    )

    conc_combined_df = conc_combined_df.merge(
        df_levels,
        on=["ms_file", "peak_label", "batch"],
        how="left"
    )

    # Join combined dataframe to the dataframe with best fit model lookup table by metabolite x level x batch
    best_fit_df = best_fit_df.drop(columns=["in_range"])
    conc_filtered_df = conc_combined_df.merge(
        best_fit_df,
        on=["peak_label", "level", "batch"],
        how="left"
    )

    """
    Create a column called score that will be used to filter the optimal concentrations.
    The score will be the sum of two factors: 
        1) if the sample concentration is in-range 
        2) if the model fit of the predicted concentration matches the 'best fit' model
    This ensures that a sample being in-range is prioritized first and then model fit.
    Without this scoring, it was possible for samples to be dropped because the concentration
    from the best fit model was not in-range for that particular sample. 
    """
    conc_filtered_df["score"] = (
        conc_filtered_df["in_range"].astype(int) * 2
      + (conc_filtered_df["fit"] == conc_filtered_df["best_fit"]).astype(int)
    )
    # Group the data by sample, metabolite, and batch, then sort by score and keep only the first occurence.
    # This results in a dataframe of all the samples (including duplicates) with predicted concentrations from
    # the best-fit model for that batch x sample x metabolite
    conc_filtered_df = (
        conc_filtered_df
          .sort_values("score", ascending=False)
          .groupby(["ms_file","peak_label","batch"], as_index=False)
          .first()
    )
    # Drop the break-point, best_fit, and score columns
    conc_filtered_df = conc_filtered_df.drop(columns=["b1", "b2", "best_fit", "score"])
    # While this unfiltered and unformated dataframe is exported by this script, it is
    # currently marked as a temporary file in the Snakefile, so it will not appear in the
    # output results. If you would like to see this file, just remove the temp flag in the
    # snakefile.
    conc_filtered_df.to_csv(optimal_conc_file, index=False, header=True, sep='\t')

    ##### FILTER & FORMAT OUTPUT FILES #####
    filtered_df = conc_filtered_df.loc[~conc_filtered_df['ms_file'].str.contains('std')]
    filtered_df = filtered_df.loc[~filtered_df['ms_file'].str.contains('QC')]
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'^(\d\d\d\d)BW.+|^(\d\d\d\d)BW': r'\1\2'}, regex=True)
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'^SHS-(.+)': r'\1'}, regex=True)
    filtered_df['ms_file'] = filtered_df['ms_file'].replace({r'(.+)_RR$': r'\1'}, regex=True)
    filtered_df['conc_ng_mL'] = filtered_df['pred_conc'] / (
                50 / 1000)  # predicted concentration is in ng/50uL, we want to report ng/mL
    filtered_df['type'] = np.where(filtered_df['ms_file'].str.len() == 4, 'control', 'sample')

    # Define common aggregation dictionary for groupby operations
    agg_dict = {
        'ms_file': 'first',
        'value': 'mean',
        'pred_conc': 'mean',
        'conc_ng_mL': 'mean',
        'in_range': 'mean',
        'type': 'first',
        'level': lambda x: x.unique(),
        'fit': lambda x: x.unique(),
        'Average SMAPE': 'mean'
    }

    def group_and_save(df, output_path, extra_agg=None):
        """Group DataFrame and save to CSV with optional extra aggregations"""
        agg = agg_dict.copy()
        if extra_agg:
            agg.update(extra_agg)
        grouped = df.groupby(['ms_file', 'peak_label'], as_index=False).agg(agg)

        # Add count column
        grouped['n'] = df.groupby(['ms_file', 'peak_label']).size().values

        # Rename columns to match original output
        grouped = grouped.rename(columns={
            'ms_file': 'sample_id',
            'value': 'peak_ratio',
            'Average SMAPE': 'average_smape'
        })

        grouped.to_csv(output_path, index=False, header=True, sep='\t')
        return grouped

    # Process all samples
    all_grouped = group_and_save(filtered_df, all_results_path)

    # Process sample subsets
    samples_only = filtered_df.loc[filtered_df['type'] == 'sample']
    samples_grouped = group_and_save(samples_only, all_samples_results_path)

    # Handle in-range and out-of-range samples
    in_range = samples_only.loc[samples_only['in_range'] == 1]
    out_of_range = samples_only.loc[samples_only['in_range'] == 0]

    # Process failed samples
    failed_samples = out_of_range[~out_of_range[['ms_file', 'peak_label']].apply(tuple, axis=1).isin(
        in_range[['ms_file', 'peak_label']].apply(tuple, axis=1))]
    failed_samples_grouped = group_and_save(failed_samples, failed_samples_results_path)

    # Process in-range samples
    in_range_samples_grouped = group_and_save(in_range, valid_samples_results_path)

    # Filter and save samples below threshold
    in_range_filtered_samples_grouped = in_range_samples_grouped.loc[
        in_range_samples_grouped['average_smape'] <= threshold].copy()
    in_range_filtered_samples_grouped.to_csv(valid_samples_filtered_results_path, index=False, header=True, sep='\t')

    # Process control samples
    controls_only = filtered_df.loc[filtered_df['type'] == 'control'].copy()
    pos_controls = ["2727", "3819", "3089", "4753", "5487"]
    controls_only['control_type'] = np.where(
        (controls_only['type'] == 'control') & (controls_only['ms_file'].isin(pos_controls)),
        'pos', 'neg')
    group_and_save(controls_only, control_samples_results_path,
                   {'control_type': 'first'})


if __name__ == '__main__':
    main(snakemake.input["conc_files"],
         snakemake.input["qc_conc_file"],
         snakemake.input["qc_stats_files"],
         snakemake.params["threshold"],
         snakemake.output["optimal_conc_file"],
         snakemake.output["all_results_file"],
         snakemake.output["all_samples_file"],
         snakemake.output["valid_samples_file"],
         snakemake.output["valid_samples_filtered_file"],
         snakemake.output["failed_samples_file"],
         snakemake.output["control_samples_file"]
         )
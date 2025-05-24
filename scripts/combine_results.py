import pandas as pd
import numpy as np

def main(list_of_batch_paths, all_results_path, all_samples_results_path,
         valid_samples_results_path, failed_samples_results_path,
         control_samples_results_path):
    concat_df = pd.concat(map(pd.read_csv, list_of_batch_paths))
    filtered_df = concat_df.loc[~concat_df['ms_file'].str.contains('std')]
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
            n = ('ms_file', 'size')
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
            n = ('ms_file', 'size')
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
            n = ('ms_file', 'size')
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
            n = ('ms_file', 'size')
        )
    )
    in_range_samples_grouped.to_csv(valid_samples_results_path, index=False, header=True, sep='\t')

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
            n = ('ms_file', 'size')
        )
    )
    controls_grouped.to_csv(control_samples_results_path, index=False, header=True, sep='\t')

if __name__ == '__main__':
    main(snakemake.input["paths"],
         snakemake.output["all_results_file"],
         snakemake.output["all_samples_file"],
         snakemake.output["valid_samples_file"],
         snakemake.output["failed_samples_file"],
         snakemake.output["control_samples_file"]
         )
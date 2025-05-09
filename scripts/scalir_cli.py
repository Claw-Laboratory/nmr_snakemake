import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from ms_conc import calibration_curves as cc
from ms_conc import ConcentrationEstimator as CE

def heav(x):
    return 1 if x > 0.1 else 2

import sys
import datetime

SCALIR_VERSION = "1.0.0"

def main():
    parser = argparse.ArgumentParser(description="SCALiR CLI tool for metabolite concentration estimation")
    parser.add_argument('--standards', required=True, help='Path to standards concentration file (.csv or .xlsx)')
    parser.add_argument('--peaklist', required=True, help='Path to peaklist data file (.csv or .xlsx)')
    parser.add_argument('--program', choices=['Mint', 'Maven'], help='Program used to generate peaklist (required if not internal)')
    parser.add_argument('--mint_table_type', choices=['full results', 'dense peak_max'], help='Mint table type')
    parser.add_argument('--fit_mode', choices=['fixed', 'interval', 'wide'], default='wide', help='Fit mode for regression')
    parser.add_argument('--fit_interval', nargs=2, type=float, metavar=('LOWER', 'UPPER'), default=[0.85, 1.15], help='Slope interval bounds for interval fitting')
    parser.add_argument('--output_prefix', default='scalir_output', help='Prefix for output files')
    parser.add_argument('--plot_dir', help='Directory to save log-log plots for each compound')
    args = parser.parse_args()

    # Prepare log file
    log_file = f"{args.output_prefix}_run.log"
    with open(log_file, 'w') as log:
        log.write(f"SCALiR Version: {SCALIR_VERSION}\n")
        log.write(f"Run Date/Time: {datetime.datetime.now().isoformat()}\n")
        log.write("Command: " + ' '.join(sys.argv) + "\n")
        log.write("\nRun Parameters:\n")
        for arg, val in vars(args).items():
            log.write(f"  {arg}: {val}\n")

    print(f"📝 Log file saved to: {log_file}")

    # Load standards file
    if args.standards.endswith('.csv'):
        std_df = pd.read_csv(args.standards)
    elif args.standards.endswith('.xlsx'):
        std_df = pd.read_excel(args.standards)
    else:
        raise ValueError("Unsupported file type for standards")

    internal_mode = 'internal_standard_injection_concentration' in std_df.columns

    if 'unit' in std_df.columns:
        units = std_df[['peak_label', 'unit']].fillna('')
        std_df = std_df.drop(columns=['unit'])
    else:
        units = std_df[['peak_label']].copy()
        units['unit'] = 'µM'

    if internal_mode:
        internal_stds = std_df[['peak_label', 'internal_standard_injection_concentration']].copy()
        std_df = std_df.drop(columns=['internal_standard_injection_concentration'])

    if args.peaklist.endswith('.csv'):
        raw_df = pd.read_csv(args.peaklist).dropna(how='all')
    elif args.peaklist.endswith('.xlsx'):
        raw_df = pd.read_excel(args.peaklist).dropna(how='all')
    else:
        raise ValueError("Unsupported file type for peaklist")

    if internal_mode:
        by = 'peak_max'
        raw_df = cc.info_from_Mint_dense(raw_df)
    elif args.program == 'Mint':
        if args.mint_table_type == 'full results':
            by = 'peak_max'
            raw_df = cc.info_from_Mint(raw_df, by)
        else:
            raw_df = cc.info_from_Mint_dense(raw_df)
            by = 'peak_max'
    elif args.program == 'Maven':
        by = 'value'
        raw_df = cc.info_from_Maven(raw_df)
    else:
        raise ValueError("Program must be specified for non-internal mode")

    common_peaks = np.intersect1d(raw_df.peak_label.unique(), std_df.peak_label.unique())
    raw_df = raw_df[raw_df.peak_label.isin(common_peaks)]
    std_df = std_df[std_df.peak_label.isin(common_peaks)]

    std_results = cc.setting_from_stdinfo(std_df, raw_df)
    std_results.sort_values(by=['peak_label', 'STD_CONC', by], inplace=True)
    x_train, y_train = cc.training_from_standard_results(std_results, by=by)
    x_train = x_train.rename(columns={by: 'value'})

    ces = CE.ConcentrationEstimator()
    if args.fit_mode == 'interval':
        ces.set_interval(np.array(args.fit_interval))

    ces.fit(x_train, y_train, v_slope=args.fit_mode)

    X = raw_df[['ms_file', 'peak_label', by]].rename(columns={by: 'value'})
    predictions = ces.predict(X)
    X['pred_conc'] = predictions.pred_conc
    X['in_range'] = predictions.in_range

    if internal_mode:
        for metab in X.peak_label.unique():
            factor = internal_stds.loc[internal_stds.peak_label == metab, 'internal_standard_injection_concentration'].values[0]
            X.loc[X.peak_label == metab, 'pred_conc'] *= factor

    param_out = f"{args.output_prefix}_parameters.csv"
    conc_out = f"{args.output_prefix}_results.csv"

    linear_params = ces.params_.sort_values(by='peak_label').drop(['lin_range_min', 'lin_range_max'], axis=1)
    linear_params.rename(columns={'slope': 'log_scale_slope', 'intercept': 'log_scale_intercept'}, inplace=True)
    linear_params['log_scale_slope'] = 1 / linear_params['log_scale_slope']
    linear_params['log_scale_intercept'] = -linear_params['log_scale_intercept'] * linear_params['log_scale_slope']

    linear_params.to_csv(param_out, index=False)
    X.to_csv(conc_out, index=False)

    print(f"✅ Parameters saved to: {param_out}")
    print(f"✅ Concentration predictions saved to: {conc_out}")

    # Optional plotting
    if args.plot_dir:
        os.makedirs(args.plot_dir, exist_ok=True)
        y_train_corrected = cc.train_to_validation(x_train, y_train, ces.params_)
        viz_df = x_train.copy()
        viz_df['pred_conc'] = ces.predict(viz_df).pred_conc
        viz_df['Concentration'] = y_train
        viz_df['Corr_Concentration'] = y_train_corrected
        viz_df['in_range'] = viz_df.Corr_Concentration.apply(heav)
        viz_df = viz_df[viz_df.Concentration > 1e-8]

        for peak in viz_df.peak_label.unique():
            unit = units.loc[units.peak_label == peak, 'unit'].values[0].replace('u', 'μ')
            xlabel = f"{peak} concentration ({unit})"
            ylabel = f"{peak} response ratio" if internal_mode else f"{peak} intensity (AU)"

            subset = viz_df[(viz_df.peak_label == peak) & (viz_df.value > 0)]
            fig = plt.figure(figsize=(4, 4))
            for inr, color, label in zip([2, 1], ['gray', 'black'], ['Out of range', 'In range']):
                xvals = subset.Concentration[subset.in_range == inr]
                yvals = subset.value[subset.in_range == inr]
                plt.plot(xvals, yvals, 'o', color=color, label=label)
            plt.plot(subset.pred_conc[subset.in_range == 1], subset.value[subset.in_range == 1], color='black')
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.xscale('log')
            plt.yscale('log')
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='lightgray')
            plt.legend(fontsize=10)
            plt.tight_layout()
            fig.savefig(os.path.join(args.plot_dir, f"{peak}_standard_curve.png"), dpi=300)
            plt.close(fig)

        print(f"📈 Plots saved to: {args.plot_dir}")

if __name__ == '__main__':
    main()

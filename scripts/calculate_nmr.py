import pandas as pd
import numpy as np

def main(valid_samples, valid_samples_nmr):
    df = pd.read_csv(valid_samples, sep='\t')
    wide_df = df.pivot(index='sample_id', columns='peak_label', values='conc_ng_mL')
    wide_df['NMRcalc'] = wide_df['3-hydroxycotinine']/wide_df['cotinine']
    wide_df['logNMR'] = np.log10(wide_df['NMRcalc'])
    wide_df.to_csv(valid_samples_nmr, index=True, header=True, sep='\t')


if __name__ == '__main__':
    main(snakemake.input["valid_samples"],
         snakemake.output["valid_samples_nmr"],
         )
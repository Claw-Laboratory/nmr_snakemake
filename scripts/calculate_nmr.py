import pandas as pd
import numpy as np

def main(valid_samples, valid_samples_nmr):
    """
    Processes tabular data from a specified input file, reshaping it, performing calculations on
    specific columns, and saving the resulting data to an output file. The script specifically computes
    the Nicotine Metabolite Ratio (NMR) and log(NMR). NMR = 3-hydroxycotinine / and cotinine.

    :param valid_samples: Path to the input csv file containing predicted concentrations for the sample data.
    :type valid_samples: str
    :param valid_samples_nmr: Path to the output csv file.
    :type valid_samples_nmr: str
    :return: None
    """
    df = pd.read_csv(valid_samples, sep='\t')
    wide_df = df.pivot(index='sample_id', columns='peak_label', values='conc_ng_mL')
    wide_df['NMRcalc'] = wide_df['3-hydroxycotinine']/wide_df['cotinine']
    wide_df['logNMR'] = np.log10(wide_df['NMRcalc'])
    wide_df = wide_df.dropna(subset=['logNMR'])
    wide_df.to_csv(valid_samples_nmr, index=True, header=True, sep='\t')


if __name__ == '__main__':
    main(snakemake.input["valid_samples"],
         snakemake.output["valid_samples_nmr"],
         )
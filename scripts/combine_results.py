import pandas as pd

def main(list_of_batch_paths, out_file):
    concat_df = pd.concat(map(pd.read_csv, list_of_batch_paths)).drop_duplicates()
    concat_df.to_csv(out_file, index=False, header=True, sep='\t')


if __name__ == '__main__':
    main(snakemake.input["paths"],
         snakemake.output["combined_results"])
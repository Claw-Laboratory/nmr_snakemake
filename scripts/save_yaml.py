import yaml


def main(run_config, out_yaml_fp):

    with open(out_yaml_fp, 'w', encoding='utf8') as f:
        yaml.dump(run_config, f, allow_unicode=True)



if __name__ == "__main__":
    main(snakemake.params["run_config"],
         snakemake.output["saved_yaml"])
import yaml


def main(run_config, out_yaml_fp):
    """
    Writes a provided configuration dictionary to a specified YAML file.

    This function takes a configuration dictionary and writes it to a YAML file
    at the specified file path. The file is encoded in UTF-8, and the YAML content
    supports Unicode characters.

    :param run_config: The configuration dictionary to be written to the YAML file.
    :param out_yaml_fp: The file path to which the YAML content is to be written.
    :return: None
    """
    with open(out_yaml_fp, 'w', encoding='utf8') as f:
        yaml.dump(run_config, f, allow_unicode=True)



if __name__ == "__main__":
    main(snakemake.params["run_config"],
         snakemake.output["saved_yaml"])
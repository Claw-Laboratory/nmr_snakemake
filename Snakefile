from snakemake.utils import min_version
min_version("8.0")

configfile: "config.yaml"
from utils import ConfigWizard
cf_wizard = ConfigWizard(config)

rule all:
  input:
    expand(cf_wizard.scalir_concentration_file,
      batch_name = cf_wizard.batch_name),
    cf_wizard.run_yaml_file

rule create_scalir_peaklist_file:
  input:
    lcms_results = cf_wizard.batch_file,
    script = "scripts/make_peakslist_file.R"
  output:
    peaklist_file = cf_wizard.scalir_peaklist_file
  conda:
    "envs/r_env.yaml"
  shell:
    """
    Rscript {input.script} \
      --infile {input.lcms_results} \
      --outfile {output.peaklist_file}
    """

rule create_scalir_IS_concentration_file:
  input:
    standards_stock_file = cf_wizard.standards,
    script = "scripts/make_IS_concentration_file.R"
  output:
    standard_conc_file = cf_wizard.scalir_standards_file
  conda:
    "envs/r_env.yaml"
  shell:
    """
    Rscript {input.script} \
      --infile {input.standards_stock_file} \
      --outfile {output.standard_conc_file}
    """

rule run_scalir:
  input:
    scalir_peaklist_file = cf_wizard.scalir_peaklist_file,
    scalir_standard_file = cf_wizard.scalir_standards_file,
    script = "scripts/scalir_cli.py"
  params:
    model_fit = cf_wizard.fit,
    out_prefix = lambda wc: str(cf_wizard.get_batch_results_prefix(wc.batch_name)),
    plot_directory = lambda wc: str(cf_wizard.get_plot_dir(wc.batch_name))
  output:
    sample_conc_file = cf_wizard.scalir_concentration_file
  conda:
    "envs/scalir_env.yaml"
  shell:
    """
    python {input.script} \
      --peaklist {input.scalir_peaklist_file} \
      --standards {input.scalir_standard_file} \
      --fit-mode {params.model_fit} \
      --output-prefix {params.out_prefix} \
      --plot-dir {params.plot_directory}
    """

rule create_batch_concentrations_map:
  input:
    paths=cf_wizard.get_all_concentrations
  output:
    cf_wizard.batch_file.concentrations_map
  conda:
    "envs/python.yaml"
  script:
    "scripts/make_batch_concentrations_map.py"

rule save_run_config_and_paths_df:
    output:
        saved_yaml=cf_wizard.run_yaml_file
    conda:
        "envs/python.yaml"
    params:
        run_config=config
    script:
        "scripts/save_yaml.py"

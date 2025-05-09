from snakemake.utils import min_version
min_version("8.0")

from utils import ConfigWizard

cf_wizard = ConfigWizard(config)

rule all:
  input:
    cf_wizard.final_results

rule create_scalir_peaklist_file:
  input:
    lcms_results = ,
    script = "scripts/make_peakslist_file.R"
  output:
    peakslist_file =
  conda:
    <condayaml>
  shell:
    """
    Rscript {input.script} \
      --infile {input.lcms_results} \
      --outfile {output.peakslist_file}}
    """

rule create_scalir_IS_concentration_file:
  input:
    standards_stock_file = ,
    script = "scripts/make_IS_concentration_file.R"
  output:
    standard_conc_file =
  conda:
    <condayaml>
  shell:
    """
    Rscript {input.script} \
      --infile {input.standards_stock_file} \
      --outfile {output.standard_conc_file} }
    """

rule run_scalir:
  input:
    rules.create_scalir_peaklist_file.output,
    rules.create_scalir_IS_concentration_file.output,
    script = "scripts/scalir_cli.py"
  params:
    model_fit = cf_wizard.model_fit,
    out_directory = ,
    plots_directory = 
  output:
    sample_conc_file =
  conda:
    <condayaml>
  shell:
    """
    python {input.script} \
      --peaklist {input.peaklist_file} \
      --standards {input.standard_conc_file} \
      --fit-mode {params.model_fit} \
      --output-prefix {params.out_directory} \
      --plots-dir {params.plots_directory}
    

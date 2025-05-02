from snakemake.utils import min_version
min_version("8.0")

from utils import ConfigWizard

cf_wizard = ConfigWizard(config)

rule all:
  input:
    cf_wizard.final_results

rule create_scalir_input_files:
  input:
    lcms_results = ,
    core_worksheet = ,
    script = 
  params:
    out_directory = cf_wizard.scalir_infiles_dir
  output:
    peakslist_file = ,
    standard_conc_file_ =
  conda:
    <condayaml>
  shell:
    """
    Rscript {input.script} \
      --lcms {input.lcms_results} \
      --cws {input.core_worksheet} \
      --outdir {params.out_directory}
    """

rule run_scalir:
  input:
    peaklist_file = ,
    standard_conc_file_ = ,
    script = 
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
    

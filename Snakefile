from snakemake.utils import min_version
min_version("9.3.3")

from utils import ConfigWizard
cf_wizard = ConfigWizard(config)

if len(cf_wizard.fit) > 1:
  extra_outputs = [cf_wizard.opt_valid_samples_nmr_file]
else:
  extra_outputs = []

rule all:
  input:
    cf_wizard.run_yaml_file,
    qc_samples= expand(cf_wizard.qc_samples_file,
      batch_name=cf_wizard.batch_name,
      fit=cf_wizard.fit),
    qc_stats= expand(cf_wizard.qc_stats_file,
      batch_name=cf_wizard.batch_name,
      fit=cf_wizard.fit),
    var_plot= expand(cf_wizard.qc_var_plot,
      batch_name=cf_wizard.batch_name,
      fit=cf_wizard.fit),
    acc_plot= expand(cf_wizard.qc_acc_plot,
      batch_name=cf_wizard.batch_name,
      fit=cf_wizard.fit),
    conc_plot_1= expand(cf_wizard.qc_conc_plot_1,
      batch_name=cf_wizard.batch_name,
      fit=cf_wizard.fit),
    conc_plot_2= expand(cf_wizard.qc_conc_plot_2,
      batch_name=cf_wizard.batch_name,
      fit=cf_wizard.fit),
    combined_results = expand(cf_wizard.all_results_file,
      fit=cf_wizard.fit),
    nmr_results = expand(cf_wizard.valid_samples_nmr_file,
      fit=cf_wizard.fit)
    + extra_outputs


rule create_scalir_peaklist_file:
  input:
    lcms_results = cf_wizard.get_batch_file_path
  output:
    peaklist_file = cf_wizard.scalir_peaklist_file
  conda:
    "envs/r_env.yaml"
  shell:
    """
    Rscript scripts/make_peakslist_file.R \
      --infile {input.lcms_results} \
      --outfile {output.peaklist_file}
    """

rule create_scalir_IS_concentration_file:
  input:
    standards_stock_file = cf_wizard.standards
  output:
    standard_conc_file = cf_wizard.scalir_standards_file
  conda:
    "envs/r_env.yaml"
  shell:
    """
    Rscript scripts/make_IS_concentration_file.R \
      --infile {input.standards_stock_file} \
      --outfile {output.standard_conc_file}
    """

rule run_scalir:
  input:
    scalir_peaklist_file = cf_wizard.scalir_peaklist_file,
    scalir_standard_file = cf_wizard.scalir_standards_file
  params:
    model_fit = cf_wizard.get_fit,
    out_prefix = cf_wizard.get_batch_results_prefix,
    plot_directory = cf_wizard.get_plot_dir
  output:
    sample_conc_file = cf_wizard.scalir_concentration_file,
    sample_params_file = cf_wizard.scalir_parameters_file,
    sample_run_file = cf_wizard.scalir_run_file
  conda:
    "envs/scalir_env.yaml"
  shell:
    """
    python scripts/scalir_cli.py \
      --peaklist {input.scalir_peaklist_file} \
      --standards {input.scalir_standard_file} \
      --fit-mode {params.model_fit} \
      --output-prefix {params.out_prefix} \
      --plot-dir {params.plot_directory}
    """

rule run_quality_assessment:
  input:
    qc_file = cf_wizard.qc_file,
    batch_conc_file = cf_wizard.scalir_concentration_file
  params:
    batch_name = cf_wizard.get_batch_name,
    qa_out_dir = cf_wizard.get_qc_out_dir
  output:
    qc_samples = cf_wizard.qc_samples_file,
    qc_stats = cf_wizard.qc_stats_file,
    var_plot = cf_wizard.qc_var_plot,
    acc_plot = cf_wizard.qc_acc_plot,
    conc_plot_1 = cf_wizard.qc_conc_plot_1,
    conc_plot_2 = cf_wizard.qc_conc_plot_2
  conda:
    "envs/r_env.yaml"
  shell:
    """
    Rscript scripts/quality_assessment.R \
      --qc-file {input.qc_file} \
      --batch-file {input.batch_conc_file} \
      --batch-name {params.batch_name} \
      --out-dir {params.qa_out_dir}
    """

rule combine_results:
  input:
    paths = cf_wizard.get_conc_file_paths,
    qc_stats_files= cf_wizard.get_qc_file_paths,
    qc_conc_file= cf_wizard.qc_file
  params:
    threshold = cf_wizard.accuracy_threshold
  output:
    all_results_file = cf_wizard.all_results_file,
    all_samples_file = cf_wizard.all_samples_results_file,
    valid_samples_file = cf_wizard.valid_samples_results_file,
    valid_samples_filtered_file = cf_wizard.valid_samples_filtered_results_file,
    failed_samples_file = cf_wizard.failed_samples_results_file,
    control_samples_file = cf_wizard.control_samples_results_file
  conda:
    "envs/python.yaml"
  script:
    "scripts/combine_results.py"

rule combine_optimal_results:
  input:
    conc_files = cf_wizard.conc_file_paths,
    qc_stats_files = cf_wizard.qc_stats_file_paths,
    qc_conc_file = cf_wizard.qc_file
  params:
    threshold = cf_wizard.accuracy_threshold
  output:
    optimal_conc_file = temp(cf_wizard.optimal_conc_file),
    all_results_file = cf_wizard.opt_all_results_file,
    all_samples_file = cf_wizard.opt_all_samples_results_file,
    valid_samples_file = cf_wizard.opt_valid_samples_results_file,
    valid_samples_filtered_file= cf_wizard.opt_valid_samples_filtered_results_file,
    failed_samples_file = cf_wizard.opt_failed_samples_results_file,
    control_samples_file = cf_wizard.opt_control_samples_results_file
  conda:
    "envs/python.yaml"
  script:
    "scripts/splice_optimal_results.py"

rule calculate_nmr:
  input:
    valid_samples = cf_wizard.valid_samples_filtered_results_file
  output:
    valid_samples_nmr = cf_wizard.valid_samples_nmr_file
  conda:
    "envs/python.yaml"
  script:
    "scripts/calculate_nmr.py"

rule calculate_nmr_optimal_results:
  input:
    valid_samples = cf_wizard.opt_valid_samples_filtered_results_file
  output:
    valid_samples_nmr = cf_wizard.opt_valid_samples_nmr_file
  conda:
    "envs/python.yaml"
  script:
    "scripts/calculate_nmr.py"

rule save_run_config_and_paths_df:
    output:
        saved_yaml = cf_wizard.run_yaml_file
    conda:
        "envs/python.yaml"
    params:
        run_config = config
    script:
        "scripts/save_yaml.py"

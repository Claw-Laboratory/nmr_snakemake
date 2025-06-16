import os
from pathlib import Path
from os.path import expandvars

class ConfigWizard(object):
    __slots__ = ("config", "batch_name_dict", "model_fit")
    def __init__(self, config):
        self.config = config
        self.batch_name_dict =  {os.path.splitext(os.path.basename(f))[0]: f for f in config["lcms_results"]}
        self.model_fit = config["model_fit"]

    @property
    def batch_files(self):
        return self.config["lcms_results"]

    @property
    def batch_name(self):
        return [os.path.splitext(os.path.basename(f))[0] for f in self.batch_files]

    @property
    def standards(self):
        return self.config["standards_stock_file"]

    @property
    def fit(self):
        return self.model_fit

    @property
    def qc_file(self):
        return self.config["qc_file"]

    @property
    def accuracy_threshold(self):
        return self.config["percent_error_threshold"]

    @property
    def out_dir(self):
        return Path(self.config["out_dir"])

    @property
    def batch_results_path(self):
        return self.out_dir.joinpath("{fit}" ,"{batch_name}")

    @property
    def scalir_input_dir(self):
        return self.batch_results_path.joinpath("scalir_input")

    @property
    def scalir_peaklist_file(self):
        return self.scalir_input_dir.joinpath("{batch_name}_peaklist.csv")

    @property
    def scalir_standards_file(self):
        return self.scalir_input_dir.joinpath("internal_standards_concentration_file.csv")

    @property
    def scalir_concentration_file(self):
        return self.batch_results_path.joinpath("{batch_name}_concentrations.csv")

    @property
    def qc_samples_file(self):
        return self.batch_results_path.joinpath("quality_assessment/{batch_name}_qc_samples.csv")

    @property
    def qc_stats_file(self):
        return self.batch_results_path.joinpath("quality_assessment/{batch_name}_qc_statistics.csv")

    @property
    def qc_var_plot(self):
        return self.batch_results_path.joinpath("quality_assessment/{batch_name}_precision_plot.jpg")

    @property
    def qc_acc_plot(self):
        return self.batch_results_path.joinpath("quality_assessment/{batch_name}_accuracy_plot.jpg")

    @property
    def qc_conc_plot_1(self):
        return self.batch_results_path.joinpath("quality_assessment/{batch_name}_qc_plot_1.jpg")

    @property
    def qc_conc_plot_2(self):
        return self.batch_results_path.joinpath("quality_assessment/{batch_name}_qc_plot_2.jpg")

    @property
    def run_yaml_file(self):
        return self.out_dir.joinpath("run_config.yaml")

    @property
    def conc_file_paths(self):
        list_of_files = [str(self.scalir_concentration_file).format(batch_name=batch_name, fit=fit)
                         for batch_name in self.batch_name
                         for fit in self.fit]
        return list_of_files

    @property
    def qc_stats_file_paths(self):
        list_of_files = [str(self.qc_stats_file).format(batch_name=batch_name, fit=fit)
                         for batch_name in self.batch_name
                         for fit in self.fit]
        return list_of_files

    @property
    def all_results_file(self):
        return self.out_dir.joinpath("{fit}", "all_results.csv")

    @property
    def all_samples_results_file(self):
        return self.out_dir.joinpath("{fit}", "samples.csv")

    @property
    def valid_samples_results_file(self):
        return self.out_dir.joinpath("{fit}", "valid_samples.csv")

    @property
    def valid_samples_filtered_results_file(self):
        return self.out_dir.joinpath("{fit}", "valid_filtered_samples.csv")

    @property
    def failed_samples_results_file(self):
        return self.out_dir.joinpath("{fit}", "failed_samples.csv")

    @property
    def control_samples_results_file(self):
        return self.out_dir.joinpath("{fit}", "control_samples.csv")

    @property
    def valid_samples_nmr_file(self):
        return self.out_dir.joinpath("{fit}", "valid_filtered_samples_nmr.csv")

    @property
    def opt_all_results_file(self):
        return self.out_dir.joinpath("combined", "all_results.csv")

    @property
    def opt_all_samples_results_file(self):
        return self.out_dir.joinpath("combined", "samples.csv")

    @property
    def opt_valid_samples_results_file(self):
        return self.out_dir.joinpath("combined", "valid_samples.csv")

    @property
    def opt_valid_samples_filtered_results_file(self):
        return self.out_dir.joinpath("combined", "valid_filtered_samples.csv")

    @property
    def opt_failed_samples_results_file(self):
        return self.out_dir.joinpath("combined", "failed_samples.csv")

    @property
    def opt_control_samples_results_file(self):
        return self.out_dir.joinpath("combined", "control_samples.csv")

    @property
    def opt_valid_samples_nmr_file(self):
        return self.out_dir.joinpath("combined", "valid_filtered_samples_nmr.csv")


    @property
    def optimal_conc_file(self):
        return self.out_dir.joinpath("combined", "optimal_concentrations.csv")

    def get_batch_file_path(self, wildcards):
        return self.batch_name_dict[wildcards.batch_name]

    def get_fit(self, wildcards):
        fit = wildcards.fit
        return fit

    def get_batch_results_prefix(self, wildcards):
        batch_name = wildcards.batch_name
        fit = wildcards.fit
        return self.out_dir.joinpath(fit, batch_name, batch_name)

    def get_plot_dir(self, wildcards):
        batch_name = wildcards.batch_name
        fit = wildcards.fit
        return self.out_dir.joinpath(fit, batch_name, "plots")

    def get_qc_out_dir(self, wildcards):
        batch_name = wildcards.batch_name
        fit = wildcards.fit
        return self.out_dir.joinpath(fit, batch_name, "quality_assessment")

    def get_batch_name(self, wildcards):
        batch_name = wildcards.batch_name
        return batch_name

    def get_conc_file_paths(self, wildcards):
        paths = [str(self.scalir_concentration_file).format(batch_name = b, fit = wildcards.fit) for b in self.batch_name]
        return paths

    def get_qc_file_paths(self, wildcards):
        paths = [str(self.qc_stats_file).format(batch_name = b, fit = wildcards.fit) for b in self.batch_name]
        return paths
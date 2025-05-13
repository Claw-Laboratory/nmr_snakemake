import os
from pathlib import Path
from os.path import expandvars

class ConfigWizard(object):
    __slots__ = ("config")
    def __init__(self, config):
        self.config = config

    @property
    def batch_files(self):
        return self.config.get("lcms_results", 1)

    @property
    def batch_name(self):
        return [os.path.splitext(os.path.basename(f))[0] for f in self.batch_files]


    @property
    def standards(self):
        return self.config.get("standards_stock_file", 1)

    @property
    def fit(self):
        return self.config.get("model_fit", 1)

    @property
    def out_dir(self):
        return Path(self.config["out_dir"])

    @property
    def batch_results_path(self):
        return self.out_dir.joinpath("{batch_name}")

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
    def run_yaml_file(self):
        return self.out_dir.joinpath("run_config.yaml")

    def get_batch_file(self, batch_name):
        for f in self.batch_files:
            return f

    def get_batch_results_prefix(self, batch_name):
        return Path(expandvars(self.config["out_dir"])) / batch_name / batch_name

    def get_plot_dir(self, batch_name):
        return Path(expandvars(self.config["out_dir"])) / batch_name / "plots"
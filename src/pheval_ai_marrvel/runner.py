from dataclasses import dataclass
from pathlib import Path

from pheval.runners.runner import PhEvalRunner

from pheval_ai_marrvel.post_process.post_process import post_process_results
from pheval_ai_marrvel.prepare.prepare import prepare_inputs
from pheval_ai_marrvel.run.run import run_commands
from pheval_ai_marrvel.tool_specific_configuration_options import AIMARRVELConfigurations


@dataclass
class AIMARRVELRunner(PhEvalRunner):
    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        """
        Pre-process phenopackets into tool accepted format.
        """
        print("creating HPO txt files from phenopackets")
        prepare_inputs(testdata_dir=self.testdata_dir)

    def run(self):
        """
        Run AI-MARRVEL to produce the raw output.
        """
        print("running with AI-MARRVEL")
        config = AIMARRVELConfigurations.parse_obj(
            self.input_dir_config.tool_specific_configuration_options
        )
        run_commands(
            tool_input_commands_dir=self.tool_input_commands_dir,
            testdata_dir=self.testdata_dir,
            input_dir=self.input_dir,
            output_dir=self.raw_results_dir,
            environment=config.environment,
        )

    def post_process(self):
        """
        Post-process the raw output into PhEval standardised TSV output.
        """
        print("post processing results to PhEval standardised TSV output.")
        post_process_results(
            raw_results_dir=self.raw_results_dir,
            output_dir=self.output_dir,
            phenopacket_dir=self.testdata_dir.joinpath("phenopackets"),
            variant_analysis=self.input_dir_config.variant_analysis,
            gene_analysis=self.input_dir_config.gene_analysis,
        )

from pathlib import Path

from pheval_ai_marrvel.post_process.post_process_results_format import create_standardised_results


def post_process_results(
    raw_results_dir: Path,
    output_dir: Path,
    phenopacket_dir: Path,
    gene_analysis: bool,
    variant_analysis: bool,
) -> None:
    """
    Post-process AI-MARRVEL raw results and create standardised PhEval TSV results.

    Args:
        raw_results_dir (Path): Path to the raw result directory.
        output_dir (Path): Path to the output directory.
        phenopacket_dir (Path): Path to the directory containing the phenopackets.
        gene_analysis (bool): Whether the generate gene results.
        variant_analysis (bool): Whether the generate variant results.
    """
    create_standardised_results(
        raw_results_dir=raw_results_dir,
        output_dir=output_dir,
        phenopacket_dir=phenopacket_dir,
        gene_analysis=gene_analysis,
        variant_analysis=variant_analysis,
    )

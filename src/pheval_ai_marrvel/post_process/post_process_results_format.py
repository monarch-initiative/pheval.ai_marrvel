import uuid
from pathlib import Path

import polars as pl
from pheval.post_processing.phenopacket_truth_set import calculate_end_pos
from pheval.post_processing.post_processing import (
    SortOrder,
    generate_gene_result,
    generate_variant_result,
)
from pheval.utils.file_utils import all_files


def read_raw_result(raw_result_path: Path) -> pl.DataFrame:
    """
    Read the raw result file.

    Args:
        raw_result_path(Path): Path to the raw result file.

    Returns:
        pl.DataFrame: Contents of the raw result file.
    """
    return (
        pl.read_csv(raw_result_path)
        .rename({"Unnamed: 0": "variant"})
        .select(pl.col(["variant", "predict", "geneSymbol", "geneEnsId"]))
        .unique(["variant", "geneSymbol"], maintain_order=True)
    )


def extract_gene_results(raw_result: pl.DataFrame) -> pl.DataFrame:
    return raw_result.select(
        [
            pl.col("predict").alias("score").cast(pl.Float64),
            pl.col("geneSymbol").alias("gene_symbol").cast(pl.String),
            pl.col("geneEnsId").alias("gene_identifier").cast(pl.String),
            pl.col("variant")
            .map_elements(lambda v: str(uuid.uuid5(uuid.NAMESPACE_DNS, v)), return_dtype=pl.Utf8)
            .alias("grouping_id"),
        ]
    )


def extract_variant_results(raw_result: pl.DataFrame) -> pl.DataFrame:
    return (
        raw_result.unique(subset="variant", maintain_order=True)
        .with_columns(
            [
                pl.col("variant").str.split("-").alias("split_variant"),
            ]
        )
        .select(
            [
                pl.col("split_variant").list.get(0).alias("chrom").cast(pl.String),
                pl.col("split_variant").list.get(1).alias("start").cast(pl.Int64),
                pl.col("split_variant").list.get(2).alias("ref").cast(pl.String),
                pl.col("split_variant").list.get(3).alias("alt").cast(pl.String),
                pl.col("predict").alias("score").cast(pl.Float64),
            ]
        )
        .with_columns(
            pl.struct("start", "ref")
            .map_elements(lambda x: calculate_end_pos(x["start"], x["ref"]))
            .alias("end")
            .cast(pl.String)
        )
    )


def create_standardised_results(
    raw_results_dir: Path,
    output_dir: Path,
    phenopacket_dir: Path,
    gene_analysis: bool,
    variant_analysis: bool,
) -> None:
    """
    Create PhEval gene and variant tsv output from raw results.

    Args:
        raw_results_dir (Path): Path to the raw result directory.
        output_dir (Path): Path to the output directory.
        phenopacket_dir (Path): Path to the phenopacket directory.
        gene_analysis (bool): Whether to generate gene results.
        variant_analysis (bool): Whether to generate variant results.
    """
    raw_results = [file for file in all_files(raw_results_dir) if "_integrated.csv" in file.name]
    for raw_result_path in raw_results:
        raw_result = read_raw_result(raw_result_path)
        if gene_analysis:
            pheval_gene_result = extract_gene_results(raw_result)
            generate_gene_result(
                results=pheval_gene_result,
                sort_order=SortOrder.DESCENDING,
                output_dir=output_dir,
                phenopacket_dir=phenopacket_dir,
                result_path=Path(str(raw_result_path).replace("_integrated", "")),
            )
        if variant_analysis:
            pheval_variant_result = extract_variant_results(raw_result)
            generate_variant_result(
                results=pheval_variant_result,
                sort_order=SortOrder.DESCENDING,
                output_dir=output_dir,
                phenopacket_dir=phenopacket_dir,
                result_path=Path(str(raw_result_path).replace("_integrated", "")),
            )

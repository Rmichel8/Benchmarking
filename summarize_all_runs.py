#!/usr/bin/env python3

import pandas as pd
import argparse


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--input_glob",
        required=True,
        nargs="+",
        help="es: benchmarking/sample*/best/results_gene_gene_metrics.tsv"
    )

    ap.add_argument(
        "--out",
        default="pipeline_summary.tsv"
    )

    args = ap.parse_args()

    files = sorted(args.input_glob)

    if not files:
        raise ValueError("No metric files found!")

    print("Files detected:")
    for f in files:
        print(" ", f)

    # --------------------------------------------------
    # Load all runs
    # --------------------------------------------------

    dfs = []

    for f in files:

        df = pd.read_csv(f, sep="\t")

        required = {
            "sample",
            "pipeline",
            "precision",
            "recall",
            "f1"
        }

        missing = required - set(df.columns)

        if missing:
            raise ValueError(
                f"{f} is missing columns: {sorted(missing)}"
            )

        dfs.append(df)

    all_runs = pd.concat(dfs, ignore_index=True)

    # --------------------------------------------------
    # Check runs
    # --------------------------------------------------

    print("\nRuns detected:")
    print(sorted(all_runs["sample"].unique()))

    print("\nNumber of runs:")
    print(all_runs["sample"].nunique())

    # --------------------------------------------------
    # Expected runs
    # --------------------------------------------------

    expected_runs = [
        "run1",
        "run2",
        "run3",
        "run4",
        "run5",
        "run6"
    ]

    missing_runs = set(expected_runs) - set(all_runs["sample"])

    if missing_runs:
        raise ValueError(
            f"Missing runs: {sorted(missing_runs)}"
        )

    # --------------------------------------------------
    # Expected pipelines
    # --------------------------------------------------

    pipelines = sorted(all_runs["pipeline"].unique())

    print("\nPipelines detected:")
    for p in pipelines:
        print(" ", p)

    # --------------------------------------------------
    # Create complete run × pipeline table
    # --------------------------------------------------

    expected = pd.MultiIndex.from_product(
        [expected_runs, pipelines],
        names=["sample", "pipeline"]
    ).to_frame(index=False)

    actual = all_runs[
        [
            "sample",
            "pipeline",
            "precision",
            "recall",
            "f1"
        ]
    ].copy()

    # Convert metrics to numeric
    for col in ["precision", "recall", "f1"]:
        actual[col] = pd.to_numeric(
            actual[col],
            errors="coerce"
        )

    # --------------------------------------------------
    # Merge
    # --------------------------------------------------

    complete = expected.merge(
        actual,
        on=["sample", "pipeline"],
        how="left"
    )

    # Missing pipeline/run combinations


    complete["recall"] = complete["recall"].fillna(0)

    # Precision and F1 intentionally remain NaN
    # when there are no predictions.

    complete["precision_defined"] = complete["precision"].notna()

    complete["f1_defined"] = complete["f1"].notna()

    # Number of valid runs

    valid_precision = (
        complete
        .groupby("pipeline")["precision"]
        .count()
        .rename("n_precision_runs")
    )

    valid_f1 = (
        complete
        .groupby("pipeline")["f1"]
        .count()
        .rename("n_f1_runs")
    )

    # Verify: every pipeline has 6 recall runs

    print("\nNumber of recall runs per pipeline:")

    recall_counts = (
        complete
        .groupby("pipeline")["recall"]
        .count()
        .sort_values()
    )

    print(recall_counts.to_string())

    print("\nNumber of valid precision runs per pipeline:")

    print(valid_precision.sort_values().to_string())

    print("\nNumber of valid F1 runs per pipeline:")

    print(valid_f1.sort_values().to_string())


    # Print complete matrix

    print("\n=== COMPLETE RUN × PIPELINE MATRIX ===")

    print(
        complete[
            [
                "sample",
                "pipeline",
                "precision",
                "recall",
                "f1"
            ]
        ]
        .sort_values(["pipeline", "sample"])
        .to_string(index=False)
    )

    # Save all runs


    complete.to_csv(
        "ALL_RUNS_gene_metrics.tsv",
        sep="\t",
        index=False
    )


    # Compute summary


    summary = (
        complete
        .groupby("pipeline")
        .agg(
            mean_precision=("precision", "mean"),
            std_precision=("precision", "std"),

            mean_recall=("recall", "mean"),
            std_recall=("recall", "std"),

            mean_f1=("f1", "mean"),
            std_f1=("f1", "std")
        )
        .reset_index()
    )


    # Add number of runs

    summary = summary.merge(
        valid_precision.reset_index(),
        on="pipeline"
    )

    summary = summary.merge(
        valid_f1.reset_index(),
        on="pipeline"
    )

    # Recall is ALWAYS based on 6 runs
    summary["n_recall_runs"] = len(expected_runs)

    # Reorder columns

    summary = summary[
        [
            "pipeline",

            "mean_precision",
            "std_precision",
            "n_precision_runs",

            "mean_recall",
            "std_recall",
            "n_recall_runs",

            "mean_f1",
            "std_f1",
            "n_f1_runs"
        ]
    ]

    # Round

    summary = summary.round(4)


    # Save summary

    summary.to_csv(
        args.out,
        sep="\t",
        index=False
    )

    print("\n=== FINAL PIPELINE SUMMARY ===")

    print(
        summary.to_string(index=False)
    )

    print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()


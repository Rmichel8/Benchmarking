#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os


def plot_metric(summary_df, metric, std_metric, ylabel, output_file,
                pipeline_colors):
    """
    Plot mean metric with standard deviation as error bars.
    """

    df = summary_df.dropna(subset=[metric]).copy()

    pipelines = df["pipeline"].tolist()

    plt.figure(figsize=(12, 7))

    colors = [
        pipeline_colors.get(pipeline, "gray")
        for pipeline in pipelines
    ]

    # Standard deviation is NaN when only one run contributed.
    # In that case, no error bar is shown.
    yerr = None
    if std_metric in df.columns:
        yerr = df[std_metric].fillna(0)

    plt.bar(
        pipelines,
        df[metric],
        yerr=yerr,
        capsize=4,
        color=colors
    )

    plt.ylabel(ylabel)
    plt.xlabel("Pipeline")
    plt.xticks(rotation=45, ha="right")

    plt.ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Written: {output_file}")


def plot_precision_vs_recall(summary_df, output_file,
                             pipeline_colors):
    """
    Scatter plot of mean precision vs mean recall
    with a separate legend for pipeline names.
    """

    df = summary_df.dropna(
        subset=["mean_precision", "mean_recall"]
    ).copy()

    fig, ax = plt.subplots(figsize=(10, 7))

    handles = []

    for _, row in df.iterrows():

        pipeline = row["pipeline"]

        color = pipeline_colors.get(
            pipeline,
            "gray"
        )

        point = ax.scatter(
            row["mean_recall"],
            row["mean_precision"],
            color=color,
            s=90
        )

        handles.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color,
                markersize=8,
                label=pipeline
            )
        )

    ax.set_xlabel("Mean recall")
    ax.set_ylabel("Mean precision")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Legend outside the plotting area
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Written: {output_file}")


def main():

    ap = argparse.ArgumentParser(
        description="Plot gene-level pipeline summary metrics."
    )

    ap.add_argument(
        "--summary",
        default="/nfsd/sysbiobig/Riccardom/benchmarking/pipeline_summarygenedefault.tsv",
        help="Gene-level pipeline summary TSV file"
    )

    ap.add_argument(
        "--out_prefix",
        default="/nfsd/sysbiobig/Riccardom/benchmarking/gene/summary_",
        help="Prefix for output plots"
    )

    args = ap.parse_args()

    # READ SUMMARY

    summary_df = pd.read_csv(
        args.summary,
        sep="\t"
    )

    print("\n=== GENE-LEVEL PIPELINE SUMMARY ===")
    print(summary_df.to_string(index=False))

    # PIPELINE COLORS

    pipeline_colors = {
    "MARTi": "#4C78A8",
    "RUARGpore": "#F58518",
    "Wf-metagenomics_RGI": "#54A24B",
    "Wf-metagenomics_ResFinder": "#E45756",
    "Wf-metagenomics Abricate": "#B279A2",
    "Wf-metagenomics AMRFinderPlus": "#9D755D",
    "NanoBsi_RGI": "#D37295",
    "Bacterial_genomes_RGI": "#79706E",
    "Bacterial_genomes_AMRFinderPlus": "#B5BD61",
    "Bacterial_genomes_Resfinder": "#72B7B2"
    }

    # OUTPUT DIRECTORY

    out_dir = os.path.dirname(args.out_prefix)

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # MEAN PRECISION

    plot_metric(
        summary_df,
        metric="mean_precision",
        std_metric="std_precision",
        ylabel="Mean precision",
        output_file=f"{args.out_prefix}mean_precision.png",
        pipeline_colors=pipeline_colors
    )

    # MEAN RECALL

    plot_metric(
        summary_df,
        metric="mean_recall",
        std_metric="std_recall",
        ylabel="Mean recall",
        output_file=f"{args.out_prefix}mean_recall.png",
        pipeline_colors=pipeline_colors
    )

    # MEAN F1

    plot_metric(
        summary_df,
        metric="mean_f1",
        std_metric="std_f1",
        ylabel="Mean F1",
        output_file=f"{args.out_prefix}mean_f1.png",
        pipeline_colors=pipeline_colors
    )

    # PRECISION VS RECALL

    plot_precision_vs_recall(
        summary_df,
        output_file=f"{args.out_prefix}precision_vs_recall.png",
        pipeline_colors=pipeline_colors
    )

    print("\nAll gene-level summary plots completed.")


if __name__ == "__main__":
    main()
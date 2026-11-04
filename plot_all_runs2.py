#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re
import os
import numpy as np

def plot_metrics_by_run(df: pd.DataFrame, out_prefix: str, colors: dict):
    metrics = ["precision", "recall", "f1"]

    run_mapping = (
        df[["run_num", "run"]]
        .drop_duplicates()
        .sort_values("run_num")
    )
    runs = run_mapping["run_num"].values
    run_labels = run_mapping["run"].values

    pipelines = sorted(df["pipeline"].unique())
    n_pipelines = len(pipelines)
    bar_width = 0.8 / n_pipelines 

    for metric in metrics:
        plt.figure(figsize=(10,6))
        for i, pipeline in enumerate(pipelines):
            if pipeline not in colors:
                continue  # ignora pipeline senza colore
            sub = df[df["pipeline"] == pipeline].sort_values("run_num")
            x_positions = runs + (i - n_pipelines/2) * bar_width + bar_width/2
            plt.bar(
                x_positions,
                sub[metric],
                width=bar_width,
                label=pipeline,
                color=colors[pipeline]
            )
        plt.xlabel("run")
        plt.ylabel(metric.capitalize())
        plt.title(f"{metric.capitalize()} per run")
        plt.ylim(0,1)
        plt.xticks(runs, run_labels, rotation=45, ha="right")
        plt.legend(
            fontsize=9,
            ncol=3,
            handlelength=1,
            handletextpad=0.4,
            borderpad=0.3,
            labelspacing=0.3
        )
        plt.tight_layout()
        outfile = f"{out_prefix}{metric}_by_run.png"
        plt.savefig(outfile, dpi=300)
        plt.close()
        print(f"Wrote plot: {outfile}")

def plot_precision_vs_recall(df: pd.DataFrame, out_prefix: str, colors: dict):
    plt.figure(figsize=(8,6))
    pipelines = sorted(df["pipeline"].unique())
    for pipeline in pipelines:
        if pipeline not in colors:
            continue
        sub = df[df["pipeline"] == pipeline]
        plt.scatter(
            sub["recall"],
            sub["precision"],
            label=pipeline,
            s=80,
            alpha=0.5,
            color=colors[pipeline]
        )
        for _, row in sub.iterrows():
            plt.annotate(
                row["run"],
                (row["recall"], row["precision"]),
                xytext=(4,4),
                textcoords="offset points",
                fontsize=9,
                alpha=0.9
            )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.legend(fontsize=9, ncol=1, borderpad=0.3, labelspacing=0.3)
    plt.tight_layout()
    outfile = f"{out_prefix}precision_vs_recall_scatter.png"
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Wrote plot: {outfile}")

# Heatmap Drug Class
TRUTH_CLASS_SET = [
    "aminoglycoside antibiotic", "nitroimidazole antibiotic", "peptide antibiotic",
    "pulvomycin", "rifamycin antibiotic", "carbapenem",
    "cephalosporin", "fluoroquinolone antibiotic", "fosfomycin",
    "macrolide antibiotic", "penicillin beta-lactam", "tetracycline antibiotic"
]

def plot_drugclass_detection_heatmap(all_df, out_prefix):
    pipelines = sorted(all_df["pipeline"].unique())
    # Inizialize matrix
    detection = {dc: {p: 0 for p in pipelines} for dc in TRUTH_CLASS_SET}
    for _, row in all_df.iterrows():
        pipe = row["pipeline"]
        if pd.isna(pipe):
            continue
        tp_str = str(row.get("TP_list", "")) if pd.notna(row.get("TP_list")) else ""
        tp_set = set([t.strip().lower() for t in tp_str.split(",") if t.strip()])
        for dc in TRUTH_CLASS_SET:
            if dc.lower() in tp_set:
                detection[dc][pipe] += 1
    heatmap_df = pd.DataFrame(detection).T[pipelines]
    plt.figure(figsize=(14, 10))
    im = plt.imshow(heatmap_df.values, aspect="auto", cmap="YlOrRd")
    plt.colorbar(im, label="Number of runs in which the class is detected")
    plt.yticks(np.arange(len(heatmap_df.index)), heatmap_df.index)
    plt.xticks(np.arange(len(heatmap_df.columns)), heatmap_df.columns, rotation=45, ha="right")
    plt.title("Heatmap: Drug-Class detection")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}heatmap_detection_total.png", dpi=300)
    plt.close()
    print(f"Wrote plot: {out_prefix}heatmap_detection_total.png")
# Barplot Drug Class frequency
def plot_drugclass_detection_frequency(all_df, out_prefix):

    pipelines = sorted(all_df["pipeline"].dropna().unique())

    # Count detection frequency across run and pipelines
    detection_counts = {
        dc: 0
        for dc in TRUTH_CLASS_SET
    }

    for _, row in all_df.iterrows():

        tp_str = (
            str(row.get("TP_list", ""))
            if pd.notna(row.get("TP_list"))
            else ""
        )

        tp_set = {
            t.strip().lower()
            for t in tp_str.split(",")
            if t.strip()
        }

        for dc in TRUTH_CLASS_SET:

            if dc.lower() in tp_set:
                detection_counts[dc] += 1

    frequency_df = (
        pd.DataFrame.from_dict(
            detection_counts,
            orient="index",
            columns=["detections"]
        )
        .sort_values(
            "detections",
            ascending=False
        )
    )

    plt.figure(figsize=(12, 6))

    plt.bar(
        frequency_df.index,
        frequency_df["detections"]
    )

    plt.xlabel("Drug class")
    plt.ylabel("Number of detections")
    plt.title("Drug-Class detection frequency")

    plt.xticks(
        rotation=45,
        ha="right"
    )

    plt.ylim(
        0,
        max(frequency_df["detections"]) + 1
    )

    plt.tight_layout()

    outfile = f"{out_prefix}drugclass_detection_frequency.png"

    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Wrote plot: {outfile}")
    frequency_df.to_csv(
        f"{out_prefix}drugclass_detection_frequency.tsv",
        sep="\t"
    )

# Main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--metrics_files",
        required=True,
        nargs="+",
        help="Lista dei file results_drug_metrics_drugclass_level_clean.tsv"
    )
    ap.add_argument("--out_prefix", default="all_runs_", help="Prefisso output plot")
    args = ap.parse_args()

    def extract_run_number(path):
        m = re.search(r"run(\d+)", path)
        return int(m.group(1)) if m else float('inf')

    all_metrics_files = sorted(args.metrics_files, key=extract_run_number)

    df_list = []
    for f in all_metrics_files:
        df = pd.read_csv(f, sep="\t", dtype=str)
        if "run" not in df.columns:
            run_id = os.path.basename(os.path.dirname(f))
            df["run"] = run_id
        for col in ["precision","recall","f1"]:
            df[col] = df[col].astype(float)
        df_list.append(df)

    all_df = pd.concat(df_list, ignore_index=True)

    # Build all run × pipeline combination
    all_pipelines = all_df["pipeline"].unique()
    all_runs = all_df["run"].unique()
    full_index = pd.MultiIndex.from_product([all_runs, all_pipelines], names=["run","pipeline"])
    all_df = all_df.set_index(["run","pipeline"]).reindex(full_index).reset_index()

    # Fill missing metrics
    for col in ["precision","recall","f1"]:
        all_df[col] = all_df[col].fillna(0)

    all_df["run_num"] = all_df["run"].str.extract(r'(\d+)')[0].fillna(0).astype(int)
    all_df = all_df.sort_values("run_num")
    pipeline_colors = {
        "Bacterial_genomes_AMRFinderPlus": "tab:blue",
        "Bacterial_genomes_RGI": "tab:orange",
        "Bacterial_genomes_resfinder": "tab:green",
        "MARTi": "tab:red",
        "Metagenomics_ResFinder": "tab:purple",
        "Metagenomics_RGI": "tab:brown",
        "Wf-metagenomics Abricate": "tab:pink",
        "Wf-metagenomics AMRFinderPlus": "tab:gray",
        "RUARGpore": "tab:olive",
        "NanoBsi_RGI": "tab:cyan"
    }

    #plot_metrics_by_run(all_df, args.out_prefix, pipeline_colors)
    #plot_precision_vs_recall(all_df, args.out_prefix, pipeline_colors)
    plot_drugclass_detection_heatmap(all_df, args.out_prefix)
    plot_drugclass_detection_frequency(
        all_df,
        args.out_prefix
    )
if __name__ == "__main__":
    main()
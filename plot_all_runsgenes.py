#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re
import os
import numpy as np

def plot_metrics_by_run(df: pd.DataFrame, out_prefix: str,colors:dict):

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
def plot_precision_vs_recall(df: pd.DataFrame, out_prefix: str, colors:dict):

    plt.figure(figsize=(8,6))

    pipelines = sorted(df["pipeline"].unique())

    for pipeline in pipelines:
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

    plt.legend(
        fontsize=9,
        ncol=1,
        borderpad=0.3,
        labelspacing=0.3
    )

    plt.tight_layout()

    outfile = f"{out_prefix}precision_vs_recall_scatter.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Wrote plot: {outfile}")
TRUTH_CLASS_SET = [
    "aminoglycoside antibiotic", "nitroimidazole antibiotic", "peptide antibiotic",
    "polymyxin", "pulvomycin", "rifamycin antibiotic", "carbapenem",
    "cephalosporin", "fluoroquinolone antibiotic", "fosfomycin",
    "macrolide antibiotic", "penicillin beta-lactam", "tetracycline antibiotic"
]
TRUTH_GENE_SET = [
    "abaq", "adef", "amph", "amrb", "arma", "bcii", "bes-1", 
    "carb-23", "crp", "ef-tu", "emrr", "fosb", "gyrb", "kpne", "kpnf", 
    "kpnh", "mcr-4.1", "mcr-4.5", "msba", "pbp3", "pmrf", "sgm", "soxr", 
    "uhpt", "tet(45)"
]
def plot_gene_detection_heatmap(all_df, out_prefix):
    """
    Crea la heatmap cumulativa (0-n_runs) basata sulla colonna TP_gene_list
    """
    import numpy as np
    pipelines = sorted(all_df["pipeline"].unique())
    detection = {}

    for gene in TRUTH_GENE_SET:
        detection[gene] = {}
        for pipeline in pipelines:
            detection[gene][pipeline] = 0
    # Fill counrs
    for _, row in all_df.iterrows():
        pipe = row["pipeline"]
        # TP gene list
        tp_str = str(row.get("TP_gene_list", "")) if pd.notna(row.get("TP_gene_list")) else ""
        tp_set = set([g.strip().lower() for g in tp_str.split(",") if g.strip()])

        for gene in TRUTH_GENE_SET:
            if gene.lower() in tp_set:
                detection[gene][pipe] += 1

    heatmap_df = pd.DataFrame(detection).T[pipelines]

    # Plotting
    plt.figure(figsize=(14, 12))
    im = plt.imshow(heatmap_df.values, aspect="auto", cmap="YlOrRd") 
    plt.colorbar(im, label="Number of runs in which the gene is detected")

    plt.yticks(np.arange(len(heatmap_df.index)), heatmap_df.index)
    plt.xticks(np.arange(len(heatmap_df.columns)), heatmap_df.columns, rotation=45, ha="right")

    plt.title("Heatmap: AMR gene detection")
    plt.tight_layout()
    
    outfile = f"{out_prefix}heatmap_gene_detection_total.png"
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Wrote gene heatmap: {outfile}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--metrics_files",
        required=True,
        nargs="+",  
        help="Lista dei file results_gene_metrics.tsv delle varie run"
    )
    ap.add_argument("--out_prefix", default="all_runs_", help="Prefisso output plot")
    args = ap.parse_args()


    def extract_run_number(path):
        m = re.search(r"run(\d+)", path)
        return int(m.group(1)) if m else float('inf')

    all_metrics_files = sorted(args.metrics_files, key=extract_run_number)

    # Read and concatenate files
    df_list = []
    for f in all_metrics_files:
        df = pd.read_csv(f, sep="\t", dtype=str)
        if "run" not in df.columns:
            run_id = os.path.basename(os.path.dirname(os.path.dirname(f)))
            df["run"] = run_id
        for col in ["precision","recall","f1"]:
            df[col] = df[col].astype(float)
        df_list.append(df)

    all_df = pd.concat(df_list, ignore_index=True)
    # build all run × pipeline combination
    all_pipelines = all_df["pipeline"].unique()
    all_runs = all_df["run"].unique()

    full_index = pd.MultiIndex.from_product(
    [all_runs, all_pipelines],
    names=["run", "pipeline"]
    )
    pipelines = sorted(all_df["pipeline"].unique())
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

    all_df = all_df.set_index(["run", "pipeline"]).reindex(full_index).reset_index()
    for col in ["precision", "recall", "f1"]:
        all_df[col] = all_df[col].fillna(0)
    all_df["run_num"] = all_df["run"].str.extract(r'(\d+)').astype(int)
    all_df = all_df.sort_values("run_num")

    plot_metrics_by_run(all_df, args.out_prefix, pipeline_colors)
    plot_precision_vs_recall(all_df, args.out_prefix, pipeline_colors)
    plot_gene_detection_heatmap(all_df, args.out_prefix)

if __name__ == "__main__":
    main()

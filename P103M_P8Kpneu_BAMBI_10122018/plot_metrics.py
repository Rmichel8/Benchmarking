#!/usr/bin/env python3

import re
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


ROOT = Path(
    "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018/benchmarking"
)


OUTDIR = ROOT / "results_plots"
OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

GENE_FILES = {

    "MARTi":
    ROOT / "MARTi/MARTi_genes_summary.txt",

    "Bacterial_genomes_ResFinder":
    ROOT / "bacterial_genomes_resfinder/summary_metrics.csv",

    "Bacterial_genomes_RGI":
    ROOT / "bacterial_genomes_rgi/summary_metrics.tsv",

    "Bacterial_genomes_AMRFinder":
    ROOT / "bacterial_genomes_amrfinder/summary_metrics.tsv",

    "Wf-metagenomic_Abricate":
    ROOT / "wf-metagenomic_output/abricate/P103M_P8Kpneu_all_10122018/summary_metrics_abricate.csv",

    "Wf-metagenomic_RGI":
    ROOT / "wf-metagenomic_output/rgi/P103M_P8Kpneu_all_10122018/summary_metrics.tsv",

    "Wf-metagenomic_AMRFinder":
    ROOT / "wf-metagenomic_output/amrfinder/P103M_P8Kpneu_all_10122018/summary_metrics.tsv",
    "Wf-metagenomic_ResFinder":
    ROOT / "wf-metagenomic_output/resfinder/P103M_P8Kpneu_all_10122018/metrics.tsv",
    "NanoBSI":
    ROOT / "nanobsi/summary_metrics.csv",
    "NanoBSI+RGI":
    ROOT / "nanobsi_rginudge/summary_metrics.csv",     
   "RUARGpore":
    ROOT/"RUARGpore/P103M_P8Kpneu_all_10122018/P103M_P8Kpneu_all_10122018.fa.uniq_ARGpore_2026-07-17.13:42:29/summary_metrics.txt"
}


DRUG_FILES = {

    "MARTi":
    ROOT / "MARTi/MARTi_drugclass_summary.txt",

    "Bacterial_genomes_ResFinder":
    ROOT / "bacterial_genomes_resfinder/summary_metrics.csv",

    "Bacterial_genomes_RGI":
    ROOT / "bacterial_genomes_rgi/summary_metrics.tsv",

    "Bacterial_genomes_AMRFinder":
    ROOT / "bacterial_genomes_amrfinder/summary_metrics.tsv",

    "Wf-metagenomic_Abricate":
    ROOT /"wf-metagenomic_output/abricate/P103M_P8Kpneu_all_10122018/summary_metrics_abricate.csv",
    "Wf-metagenomic_RGI":
    ROOT / "wf-metagenomic_output/rgi/P103M_P8Kpneu_all_10122018/summary_classes.tsv",
    "Wf-metagenomic_AMRFinder":
    ROOT / "wf-metagenomic_output/amrfinder/P103M_P8Kpneu_all_10122018/summary_metrics.tsv",
    "Wf-metagenomic_ResFinder":
    ROOT / "wf-metagenomic_output/resfinder/P103M_P8Kpneu_all_10122018/metrics.tsv",
    "NanoBSI":
    ROOT / "nanobsi/summary_metrics.csv",
    "NanoBSI+RGI":
    ROOT / "nanobsi_rginudge/summary_metrics.csv",  
   "RUARGpore":
    ROOT/"RUARGpore/P103M_P8Kpneu_all_10122018/P103M_P8Kpneu_all_10122018.fa.uniq_ARGpore_2026-07-17.13:42:29/summary_metrics.txt"
}


def calc_metrics(tp, fp, fn):

    precision = tp/(tp+fp) if tp+fp else 0
    recall = tp/(tp+fn) if tp+fn else 0

    f1 = (
        2*precision*recall/(precision+recall)
        if precision+recall
        else 0
    )

    return {
        "Precision": precision,
        "Recall": recall,
        "F1": f1
    }



def parse_marti_txt(path, mode):

    with open(path) as f:
        text = f.read()


    if mode == "gene":

        patterns = {
            "TP": r"TP:\s*(\d+)",
            "FP": r"FP:\s*(\d+)",
            "FN": r"FN:\s*(\d+)"
        }


    else:

        patterns = {
            "TP": r"TP:\s*(\d+)",
            "FP": r"FP:\s*(\d+)",
            "FN": r"FN:\s*(\d+)"
        }


    values = {}

    for k,p in patterns.items():

        m = re.search(p,text)

        values[k] = int(m.group(1)) if m else 0


    return calc_metrics(
        values["TP"],
        values["FP"],
        values["FN"]
    )
def parse_bacterial_resfinder(path, mode):

    df = pd.read_csv(path)

    target = (
        "genes"
        if mode == "gene"
        else "drug_classes"
    )

    row = df[
        df["level"] == target
    ].iloc[0]

    return {
        "Precision": float(row["precision"]),
        "Recall": float(row["recall"]),
        "F1": float(row["f1"])
    }
def parse_wf_amrfinder(path, mode):

    df = pd.read_csv(
        path,
        sep="\t"
    )

    row = df.iloc[0]

    if mode == "gene":

        return {
            "Precision": float(row["precision_genes"]),
            "Recall": float(row["recall_genes"]),
            "F1": float(row["f1_genes"])
        }

    else:

        return {
            "Precision": float(row["precision_classes"]),
            "Recall": float(row["recall_classes"]),
            "F1": float(row["f1_classes"])
        }
def parse_amrfinder(path, mode):

    df = pd.read_csv(
        path,
        sep="\t"
    )

    target = (
        "Genes"
        if mode == "gene"
        else "Drug_classes"
    )

    row = df[
        df["Category"] == target
    ].iloc[0]

    return {

        "Precision": float(row["Precision"]),

        "Recall": float(row["Recall"]),

        "F1": float(row["F1"])

    }
def parse_resfinder(path, mode):

    df = pd.read_csv(
        path,
        sep="\t"
    )

    row = df.iloc[0]

    if mode == "gene":

        return {
            "Precision": float(row["precision_genes"]),
            "Recall": float(row["recall_genes"]),
            "F1": float(row["f1_genes"])
        }

    else:

        return {
            "Precision": float(row["precision_classes"]),
            "Recall": float(row["recall_classes"]),
            "F1": float(row["f1_classes"])
        }
def parse_ruargpore(path, mode):

    with open(path) as f:
        text = f.read()

    if mode == "gene":

        block = text.split(
            "========================================"
        )[0]

    else:

        block = text.split(
            "RUARGpore DRUG CLASS BENCHMARK"
        )[1]

    tp = int(
        re.search(r"TP\s*:\s*(\d+)", block).group(1)
    )

    fp = int(
        re.search(r"FP\s*:\s*(\d+)", block).group(1)
    )

    fn = int(
        re.search(r"FN\s*:\s*(\d+)", block).group(1)
    )

    return calc_metrics(tp, fp, fn)
def parse_wfrgi(path):

    df = pd.read_csv(
        path,
        sep="\t",
        header=None
    )

    values = {}

    for _, row in df.iterrows():

        key = str(row.iloc[0]).strip()

        if len(row) < 2:
            continue

        value = str(row.iloc[1]).strip()

        try:
            value = float(value)
        except:
            continue

        values[key] = value

    precision = values.get("Precision", 0)
    recall = values.get("Recall", 0)

    f1 = (
        values.get("F1_score")
        if "F1_score" in values
        else values.get("F1", 0)
    )


    return {

        "Precision": precision,
        "Recall": recall,
        "F1": f1

    }
def parse_wfabricate(path, mode):

    df = pd.read_csv(path)

    row = df[
        df["type"] == (
            "genes" if mode == "gene"
            else "drug_classes"
        )
    ].iloc[0]

    return {
        "Precision": float(row["precision"]),
        "Recall": float(row["recall"]),
        "F1": float(row["f1"])
    }
def parse_nanobsi(path, mode):

    df = pd.read_csv(
        path
    )

    target = (
        "genes"
        if mode == "gene"
        else "drug_classes"
    )

    row = df[
        df["type"] == target
    ].iloc[0]

    return {
        "Precision": float(row["precision"]),
        "Recall": float(row["recall"]),
        "F1": float(row["f1"])
    }
def parse_rgi(path, mode):

    df = pd.read_csv(
        path,
        sep="\t"
    )

    target = "Genes" if mode == "gene" else "Drug_classes"

    row = df[
        df["Category"] == target
    ].iloc[0]

    return {
        "Precision": float(row["Precision"]),
        "Recall": float(row["Recall"]),
        "F1": float(row["F1"])
    }
def build_dataframe(files, mode):

    rows = []

    for tool, file in files.items():

        if not file.exists():
            print("Missing:", file)
            continue

        print("Reading:", tool)

        if tool == "MARTi":

            m = parse_marti_txt(file, mode)

        elif tool == "RUARGpore":

            m = parse_ruargpore(file, mode)
        elif tool in ["NanoBSI", "NanoBSI+RGI"]:

            m = parse_nanobsi(file, mode)
        elif tool == "Bacterial_genomes_ResFinder":

            m = parse_bacterial_resfinder(file, mode)
        elif tool == "Bacterial_genomes_AMRFinder":

            m = parse_amrfinder(file, mode)
        elif tool == "Wf-metagenomic_ResFinder":

            m = parse_resfinder(file, mode)
        elif tool == "Wf-metagenomic_AMRFinder":

            m = parse_wf_amrfinder(file, mode)

        elif tool == "Bacterial_genomes_RGI":       

            m = parse_rgi(file, mode)

        elif "Abricate" in tool:

            m = parse_wfabricate(file, mode)

        elif tool == "Wf-metagenomic_RGI":

            m = parse_wfrgi(file)

        else:

            continue

        m["Pipeline"] = tool
        rows.append(m)

    return pd.DataFrame(rows).set_index("Pipeline")


def make_barplot(df,name):

    ax=df.plot(
        kind="bar",
        figsize=(9,5),
        ylim=(0,1)
    )


    plt.ylabel("Score")
    plt.xlabel("")
    plt.xticks(
        rotation=30,
        ha="right"
    )

    plt.legend(
        loc="lower right"
    )

    plt.tight_layout()


    plt.savefig(
        OUTDIR/f"{name}_metrics_barplot.png",
        dpi=300,
        bbox_inches="tight"
    )


    plt.close()

from collections import defaultdict

def make_pr_plot(df,name):
    
    plt.figure(figsize=(7,6))

    offsets = defaultdict(int)

    for tool in df.index:

        x = df.loc[tool, "Recall"]
        y = df.loc[tool, "Precision"]

        key = (round(x, 4), round(y, 4))
        n = offsets[key]
        offsets[key] += 1
        x += 0.01 * n
        y += 0.01 * n

        plt.scatter(
            x,
            y,
            s=100,
            label=tool
        )

    plt.xlabel("Recall")
    plt.ylabel("Precision")

    plt.xlim(0,1)
    plt.ylim(0,1)

    plt.grid(alpha=0.3)

    plt.legend(fontsize=8)

    plt.tight_layout()


    plt.savefig(
        OUTDIR/f"{name}_precision_recall.png",
        dpi=300,
        bbox_inches="tight"
    )




gene_df=build_dataframe(
    GENE_FILES,
    "gene"
)


drug_df=build_dataframe(
    DRUG_FILES,
    "drug"
)


print("\nGENES")
print(gene_df)


print("\nDRUG CLASSES")
print(drug_df)



gene_df.to_csv(
    OUTDIR/"gene_metrics.csv"
)


drug_df.to_csv(
    OUTDIR/"drugclass_metrics.csv"
)



make_barplot(
    gene_df,
    "genes"
)


make_barplot(
    drug_df,
    "drug_classes"
)


make_pr_plot(
    gene_df,
    "genes"
)


make_pr_plot(
    drug_df,
    "drug_classes"
)


print("\nResults saved in:")
print(OUTDIR)
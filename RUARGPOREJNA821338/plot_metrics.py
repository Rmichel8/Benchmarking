#!/usr/bin/env python3

import re
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking")

FILES = {
    "Bacterial_genomes_ResFinder":
    ROOT/"bacterial_genomes_resfinder/RUARGPOREPRJNA821338_resfinder_results/split_by_sample/resfinder_summary_metrics.csv",


    "Bacterial_genomes_AMRFINDER":
    Path(
        "/nfsd/sysbiobig/Riccardom/"
        "RUARGPOREPRJNA869513/"
        "bacterial_genomes_amrfinder/"
        "global_comparison_results.txt"
    ),
    "Bacterial_genomes_RGI":
    ROOT/"bacterial_genomes_rgi/summary_metrics.csv",

    "wf-AMRFinder":
    ROOT/"wfmetagenomics_amrfinder/amrfinder_compare_results/summary_metrics.log",

    "wf-Abricate":
    ROOT/"wfmetagenomics_abricate2/abricate_compare_results/summary_metrics.log",

    "wf-RGI":
    ROOT/"wf-metagenomic_rgi/rgi_compare_results/summary_metrics.log",

    "wf-Resfinder":
    ROOT/"wfmetagenomics_resfinder/resfinder_compare_results/summary_metrics.log",

    "NanoBSI":
    ROOT/"nanobsi/global_metrics.csv",

    "NanoBSI+RGI":
    ROOT/"nanobsi_rgi/rgi_global_metrics.csv",

    "MARTi":
    ROOT/"MARTi/marti_summary.csv",

    "RUARGpore":
    ROOT/"RUARGpore/compare_results/RUARGpore_global_metrics.csv",
}

metric_names = [
    "Precision",
    "Recall",
    "F1",
    "Specificity",
    "Accuracy",
    "MCC"
]


def parse_summary_log(path):

    metrics = {}

    with open(path) as f:
        text = f.read()

    patterns = {

        "Precision":r"Precision\s*=\s*([0-9.]+)",
        "Recall":r"Recall\s*=\s*([0-9.]+)",
        "Specificity":r"Specificity\s*=\s*([0-9.]+)",
        "Accuracy":r"Accuracy\s*=\s*([0-9.]+)",
        "F1":r"F1-score\s*=\s*([0-9.]+)",
        "MCC":r"MCC\s*=\s*([0-9.]+)"
    }

    for m,p in patterns.items():

        x = re.search(p,text)

        metrics[m] = float(x.group(1)) if x else None

    return metrics
def parse_resfinder_global_csv(path):

    df = pd.read_csv(path)

    row = df.iloc[0]

    return {
        "Precision": row["Precision"],
        "Recall": row["Recall"],
        "Specificity": row["Specificity"],
        "Accuracy": row["Accuracy"],
        "F1": row["F1"],
        "MCC": row["MCC"]
    }

def parse_amrfinder_global(path):

    metrics = {}

    with open(path) as f:
        text = f.read()

    patterns = {

        "Precision": r"Precision\s*:\s*([0-9.]+)",
        "Recall": r"Recall\s*:\s*([0-9.]+)",
        "F1": r"F1\s*:\s*([0-9.]+)",
        "Specificity": r"Specificity\s*:\s*([0-9.]+)",
        "Accuracy":r"Accuracy\s*:\s*([0-9.]+)",
        "MCC": r"MCC\s*:\s*([0-9.]+)"

    }

    for metric, pattern in patterns.items():

        match = re.search(pattern, text)

        metrics[metric] = (
            float(match.group(1))
            if match
            else None
        )

    return metrics
def parse_rgi_global_csv(path):

    df = pd.read_csv(path)

    row = df.iloc[0]

    return {

        "Precision": row["precision"],
        "Recall": row["recall"],
        "Specificity": row["specificity"],
        "Accuracy": row["accuracy"],
        "F1": row["f1"],
        "MCC": row["MCC"]

    }
def parse_global_csv(path):

    df = pd.read_csv(path)

    row = df.iloc[0]

    return {

        "Precision":row["Precision"],
        "Recall":row["Recall"],
        "Specificity":row["Specificity"],
        "Accuracy":row["Accuracy"],
        "F1":row["F1-score"],
        "MCC":row["MCC"]

    }

def plot_precision_vs_recall(df, outdir):

    plt.figure(figsize=(8,6))

    for tool in df.index:

        recall = df.loc[tool, "Recall"]
        precision = df.loc[tool, "Precision"]

        plt.scatter(
            recall,
            precision,
            s=100,
            alpha=0.7,
            label=tool
        )


    plt.xlabel("Recall")
    plt.ylabel("Precision")

    plt.xlim(0,1)
    plt.ylim(0,1)

    plt.grid(True, alpha=0.3)

    plt.legend(
        fontsize=8,
        loc="best"
    )

    plt.tight_layout()

    plt.savefig(
        outdir / "precision_vs_recall.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
def plot_specificity_vs_recall(df, outdir):

    plt.figure(figsize=(8,6))

    for tool in df.index:

        recall = df.loc[tool, "Recall"]
        specificity = df.loc[tool, "Specificity"]

        plt.scatter(
            recall,
            specificity,
            s=100,
            alpha=0.7,
            label=tool
        )


    plt.xlabel("Recall")
    plt.ylabel("Specificity")

    plt.xlim(0,1)
    plt.ylim(0,1)

    plt.grid(True, alpha=0.3)

    plt.legend(
        fontsize=8,
        loc="best"
    )

    plt.tight_layout()

    plt.savefig(
        outdir / "specificity_vs_recall.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
def parse_marti(path):

    df = pd.read_csv(path)

    TP=df["TP"].sum()
    FP=df["FP"].sum()
    FN=df["FN"].sum()
    TN=df["TN"].sum()

    precision=TP/(TP+FP) if TP+FP else 0
    recall=TP/(TP+FN) if TP+FN else 0
    specificity=TN/(TN+FP) if TN+FP else 0
    tested=TP+TN+FP+FN
    accuracy = (TP+TN)/tested
    f1=2*precision*recall/(precision+recall) if precision+recall else 0

    import math
    den=math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

    mcc=((TP*TN-FP*FN)/den) if den else 0

    return {

        "Precision":precision,
        "Recall":recall,
        "Specificity":specificity,
        "Accuracy":accuracy,
        "F1":f1,
        "MCC":mcc

    }

rows=[]

for tool,file in FILES.items():

    if not file.exists():
        print("Missing:",file)
        continue

    if file.suffix==".log":

        m=parse_summary_log(file)
    elif tool == "Bacterial_genomes_ResFinder":


        m = parse_resfinder_global_csv(file)

    elif tool == "Bacterial_genomes_AMRFINDER":

        m = parse_amrfinder_global(file)

    elif tool=="MARTi":

        m=parse_marti(file)
    elif tool=="Bacterial_genomes_RGI":

        print("Parsing RGI CSV:", file)
        m=parse_rgi_global_csv(file)    

    else:
        print("Parsing CSV:", tool, file)
        m=parse_global_csv(file)

    m["Tool"]=tool

    rows.append(m)

metrics=pd.DataFrame(rows)

metrics=metrics.set_index("Tool")

print(metrics)

# OUTPUT DIRECTORY

OUTDIR = Path("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking")

OUTDIR.mkdir(parents=True, exist_ok=True)


metrics.to_csv(
    OUTDIR / "overall_metrics.csv"
)

ax = metrics[metric_names].plot(
    kind="bar",
    figsize=(10,6),
    ylim=(0,1)
)

plt.ylabel("Score")
plt.xlabel("")
plt.xticks(rotation=20, ha="right")
plt.legend(loc="lower right")
plt.tight_layout()
plot_precision_vs_recall(
    metrics,
    OUTDIR
)
plot_specificity_vs_recall(
    metrics,
    OUTDIR
)
plt.savefig(
    OUTDIR / "overall_metrics.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
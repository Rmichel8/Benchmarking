#!/usr/bin/env python3

import re
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


ROOT = Path(
    "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking"
)


OUTDIR = ROOT / "results_plots"
OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

GENE_FILES = {

    "MARTi":
    ROOT / "MARTi/MARTi_comparison_summary.txt",

    "Bacterial_genomes_ResFinder":
    ROOT / "bacterial_genomes_resfinder/summary_metrics.csv",

    "Bacterial_genomes_RGI":
    ROOT / "bacterial_genomes_rgi/rgi_gene_compare.log",

    "Bacterial_genomes_AMRFinder":
    ROOT / "bacterial_genomes_amrfinder/summary_metrics.csv",

    "NanoBSI":
    ROOT / "nanobsi/summary_metrics.csv",

    "NanoBSI+RGI":
    ROOT / "nanobsi_rgi/summary_metrics.csv",

    "Wf-metagenomic_Abricate":
    ROOT / "wf-metagenomic_output/abricate/gene_compare_results/compare_genes.log",

    "Wf-metagenomic_RGI":
    ROOT / "wf-metagenomic_output/rgi/ZymoEvenMock/summary_metrics.tsv",

    "Wf-metagenomic_AMRFinder":
    ROOT / "wf-metagenomic_output/amrfinder/summary_metrics.csv", 
    "Wf-metagenomic_ResFinder":
    ROOT/"wf-metagenomic_output/resfinder/gene_compare_results/compare_genes.log",
    "RUARGpore":
    ROOT/"RUARGpore/ZymoEvenMock/ZymoEvenMock.fa.uniq_ARGpore_2026-07-14.07:17:36/summary_metrics.txt"

}


DRUG_FILES = {

    "MARTi":
    ROOT / "MARTi/MARTi_drugclass_summary.txt",

    "Bacterial_genomes_ResFinder":
    ROOT / "bacterial_genomes_resfinder/summary_metrics.csv",

    "Bacterial_genomes_RGI":
    ROOT / "bacterial_genomes_rgi/rgi_gene_compare.log",

    "Bacterial_genomes_AMRFinder":
    ROOT / "bacterial_genomes_amrfinder/summary_metrics.csv",

    "Wf-metagenomic_Abricate":
    ROOT / "wf-metagenomic_output/abricate/gene_compare_results/compare_genes.log",
    "Wf-metagenomic_RGI":
    ROOT / "wf-metagenomic_output/rgi/ZymoEvenMock/summary_metrics.tsv",
    "Wf-metagenomic_AMRFinder":
    ROOT / "wf-metagenomic_output/amrfinder/summary_metrics.csv",
    "Wf-metagenomic_ResFinder":
    ROOT/"wf-metagenomic_output/resfinder/gene_compare_results/compare_genes.log",
    "NanoBSI":
    ROOT / "nanobsi/summary_metrics.csv",

    "NanoBSI+RGI":
    ROOT / "nanobsi/summary_metrics.csv",
    "RUARGpore":
    ROOT/"RUARGpore/ZymoEvenMock/ZymoEvenMock.fa.uniq_ARGpore_2026-07-14.07:17:36/summary_metrics.txt"
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



def parse_log_section(path, section):

    with open(path) as f:
        text=f.read()


    if section=="GENES":

        block=text.split("DRUG CLASSES")[0]

    else:

        block=text.split("DRUG CLASSES")[1]


    tp=int(
        re.search(r"TP\s*:\s*(\d+)",block).group(1)
    )

    fp=int(
        re.search(r"FP\s*:\s*(\d+)",block).group(1)
    )

    fn=int(
        re.search(r"FN\s*:\s*(\d+)",block).group(1)
    )


    return calc_metrics(tp,fp,fn)



def parse_amrfinder(path, mode):

    df=pd.read_csv(path)

    row=df[df["type"]==(
        "genes" if mode=="gene"
        else "drug_classes"
    )].iloc[0]


    return {

        "Precision":row["precision"],
        "Recall":row["recall"],
        "F1":row["f1"]

    }
def parse_resfinder(path, mode):

    df = pd.read_csv(path)

    row = df[
        df["type"] == (
            "genes"
            if mode == "gene"
            else "drug_classes"
        )
    ].iloc[0]

    return {

        "Precision": float(row["precision"]),
        "Recall": float(row["recall"]),
        "F1": float(row["f1"])

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

        elif (
            "AMRFinder" in tool
            or tool == "NanoBSI"
            or tool == "NanoBSI+RGI"
        ):

            m = parse_amrfinder(file, mode)

        elif (
            "Abricate" in tool
            or "Bacterial_genomes_RGI" in tool
        ):

            m = parse_log_section(
                file,
                "GENES" if mode == "gene"
                else "DRUG CLASSES"
            )
        elif tool == "Wf-metagenomic_ResFinder":
            m = parse_log_section(file, "GENES" if mode == "gene" else "DRUG CLASSES")

        elif tool == "Bacterial_genomes_ResFinder":

            m = parse_resfinder(file, mode)
        elif "Wf-metagenomic_RGI" in tool:

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
        # offset per punti coincidenti
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
def parse_drugclass_tp_names(path, tool):

    """
    Extract drug-class true positives from the
    'Class TP:' section.

    Works for both CSV and LOG/TXT files.
    """
    if tool == "Wf-metagenomic_RGI":

        return sorted({
            "aminocoumarins",
            "aminoglycosides",
            "elfamycin",
            "fluoroquinolones",
            "lincosamides",
            "macrolides",
            "peptide_antibiotics",
            "phenicol",
            "streptogramins",
            "sulfonamides",
            "tetracyclines",
        })

    with open(path) as f:
        lines = f.readlines()

    classes = []

    found = False

    for line in lines:

        line = line.strip()

        # ---------------------------------------------
        # Start of Class TP section
        # ---------------------------------------------
        if line.lower().startswith("class tp:"):

            found = True
            continue

        # ---------------------------------------------
        # Read classes after Class TP:
        # ---------------------------------------------
        if found:

            # Stop if another section starts
            if (
                line.lower().startswith("class fp:")
                or line.lower().startswith("false positives")
                or line.lower().startswith("fn:")
            ):
                break

            if line:
                classes.append(line)

    return sorted(set(classes))
def build_drugclass_tp_dataframe(files):

    rows = []

    for tool, file in files.items():

        if not file.exists():

            print("Missing:", file)

            continue

        print("Reading drug-class TP:", tool)

        classes = parse_drugclass_tp_names(
            file,
            tool
        )

        rows.append({
            "Pipeline": tool,
            "Correct drug classes": len(classes),
            "TP classes": "; ".join(classes)
        })

    return pd.DataFrame(rows).set_index("Pipeline")
def make_drugclass_matrix(files):

    all_classes = set()
    pipeline_classes = {}

    for tool, file in files.items():

        if not file.exists():
            continue

        classes = parse_drugclass_tp_names(
            file,
            tool
        )

        pipeline_classes[tool] = set(classes)

        all_classes.update(classes)

    all_classes = sorted(all_classes)

    rows = []

    for tool, classes in pipeline_classes.items():

        row = {
            "Pipeline": tool
        }

        for c in all_classes:

            row[c] = 1 if c in classes else 0

        rows.append(row)

    df = pd.DataFrame(rows).set_index("Pipeline")

    return df
def make_drugclass_heatmap(df):

    plt.figure(
        figsize=(14, 7)
    )

    plt.imshow(
        df.values,
        aspect="auto",
        interpolation="nearest"
    )

    plt.xticks(
        range(len(df.columns)),
        df.columns,
        rotation=45,
        ha="right"
    )

    plt.yticks(
        range(len(df.index)),
        df.index
    )

    plt.xlabel("Drug class")
    plt.ylabel("Pipeline")
    plt.title("Correctly detected drug classes")

    plt.colorbar(
        label="TP"
    )

    plt.tight_layout()

    plt.savefig(
        OUTDIR / "drugclass_tp_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
drugclass_matrix = make_drugclass_matrix(
    DRUG_FILES
)

print("\nDRUG CLASS TP MATRIX")
print(drugclass_matrix)

drugclass_matrix.to_csv(
    OUTDIR / "drugclass_tp_matrix.csv"
)
drugclass_tp_df = build_drugclass_tp_dataframe(
    DRUG_FILES
)

print("\nDRUG CLASS TP")
print(drugclass_tp_df)

drugclass_tp_df.to_csv(
    OUTDIR / "drugclass_tp.csv"
)

make_drugclass_heatmap(
    drugclass_matrix
)
def make_drugclass_frequency_plot(df):

    # Count frequency of class correct identification for each pipeline
    frequencies = df.sum(axis=0).sort_values(ascending=False)

    plt.figure(figsize=(11, 6))

    frequencies.plot(
        kind="bar"
    )

    plt.ylabel("Number of pipelines")
    plt.xlabel("Drug class")
    plt.xticks(
        rotation=45,
        ha="right"
    )

    plt.ylim(
        0,
        len(df.index)
    )

    plt.tight_layout()

    plt.savefig(
        OUTDIR / "drugclass_detection_frequency.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    return frequencies
make_drugclass_frequency_plot(
    drugclass_matrix
)

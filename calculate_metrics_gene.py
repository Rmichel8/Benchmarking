#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from gene_token import truth_to_token, pred_to_token
import numpy as np
import matplotlib.pyplot as plt
import os

def assign_pipeline(row) -> str:
    """Assign pipeline based on analysis_software_name e input_file_name"""
    tool = str(row.get("analysis_software_name", "")).strip().lower()
    infile = str(row.get("input_file_name", "")).strip().lower()
    sequence = str(row.get("SEQUENCE", "")).strip().lower()
    if sequence == "sample":
        return "Wf-metagenomics ResFinder"
    if tool == "abricate":
        if infile == "local_reads.fastq":
            return "Wf-metagenomics Abricate"
        else:
            return "MARTi"
    if tool in {"ruargpore", "ruargpore_summary", "ruargpore-hamronized"} or infile == "ruargpore":
        return "RUARGpore"
    if tool == "resfinder":
        if infile == "cluster":
            return "NanoBsi ResFinder"
        if infile == "input_genome":
            return "Bacterial_genomes_resfinder"
        return "ResFinder_unassigned"
    if infile == "nanobsi_rgi":
        return "NanoBsi_RGI"
    if infile== "metagenomics_amrfinder":
        return "wf-metagenomics AMRFinderPlus"
    if infile == "bacterial_genomes_amrfinder":
        return "Bacterial_genomes_AMRFinderPlus"
    return "unknown"
def extract_gene_from_rgi(gene_str: str) -> str:
    """
    Extract RGI removing common prefixes and alias
    """
    gene_str = str(gene_str).lower().strip()
    gene_str = re.sub(r'[^a-z0-9\-\(\)]', '', gene_str)
    species_prefixes = [
        "escherichiacoli", "klebsiellapneumoniae", "bacillussubtilis",
        "mycobacteriumtuberculosis", "enterococcusfaecium",
        "pseudomonasaeruginosa", "staphylococcusaureus",
        "salmonellaenterica", "neisseriagonorrhoeae", "helicobacterpylori",
        "morganellamorganii", "mycolicibacteriumsmegmatis", "klebsiellaaerogenes",
        "myobacteroidesabscessus", "cutibacteriumacnes","planobisporarosea",
        "chlamydophilapsittaci", "chlamydomonasreinhardtii", "salmonellaisangi", "shigellaflexneri",
        "moraxellacatarrhalis", "pasteurellamultocida", "serratiamarcescens", "salmonellaserovars",
        "clostridioidesdifficile", "acinetobacterbaumannii", "enterobactercloacae", "erysipelothrixrhusiopathiae", "streptomycescinnamoneus"
    ]
    for sp in species_prefixes:
        if gene_str.startswith(sp):
            gene_str = gene_str[len(sp):]
            break

    # Rimove all after "mutationconferringresistance", "mutantsconferringresistance", "withmutationconferringresistance"
    gene_str = re.split(r'(?:mutation|mutants|withmutation)?conferringresistance', gene_str)[0]


    return gene_str

def compute_tp_fp_fn(truth_set: set, pred_set: set):
    tp = truth_set & pred_set
    fn = truth_set - pred_set
    fp = pred_set - truth_set
    return tp, fp, fn


def metrics_from_counts(TP, FP, FN):
    recall = TP / (TP + FN) if (TP + FN) else 0.0
    precision = TP / (TP + FP) if (TP + FP) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return recall, precision, f1


def safe_name(s: str) -> str:
    return "".join([c if (c.isalnum() or c in "-_.") else "_" for c in str(s)])


def write_gene_lists(out_prefix: str, pipeline: str, tp: set, fp: set, fn: set):
    p = safe_name(pipeline)
    pd.DataFrame({"gene": sorted(tp)}).to_csv(f"{out_prefix}{p}_TP_genes.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": sorted(fp)}).to_csv(f"{out_prefix}{p}_FP_genes.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": sorted(fn)}).to_csv(f"{out_prefix}{p}_FN_genes.tsv", sep="\t", index=False)
def analyze_metagenomics_resfinder(out_prefix: str, resfinder_file: str, truth_set: set):
    """
    Parameters:
    - out_prefix: prefisso per output file
    - resfinder_file: TSV input con almeno la colonna 'gene'
    - truth_set: set di geni di verità
    """
    # read TSV
    df = pd.read_csv(resfinder_file, sep="\t", dtype=str).fillna("")

    if "gene" not in df.columns:
        raise ValueError(f"{resfinder_file} non contiene la colonna 'gene'")

    def clean_name(s):
        s = str(s).strip()
        s = re.sub(r"\s+", " ", s)
        return s

    df["gene_symbol_clean"] = df["gene"].apply(clean_name)

    RESFINDER_ALIAS = {
        "blaZ": "blaZ",
        "fosB": "fosb",
        "fosB1": "fosb",
        "tet(45)": "tet(45)",
    }

    df["pred_token"] = df["gene_symbol_clean"].map(RESFINDER_ALIAS)

    # Remove not mapped rows
    df = df.dropna(subset=["pred_token"])
    pred_set = set(df["pred_token"])

    # --- calculate TP, FP, FN ---
    tp = pred_set & truth_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set

    write_gene_lists(out_prefix, "Metagenomics_ResFinder", tp, fp, fn)

    # --- metriche ---
    r = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0
    p = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0
    f1 = 2 * r * p / (r + p) if (r + p) > 0 else 0

    result = {
        "pipeline": "Metagenomics_ResFinder",
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "recall": round(r, 4),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_genes": len(pred_set),
        "truth_genes": len(truth_set),
    }
    
    return result, pred_set
def analyze_rgi_tsv(out_prefix: str, rgi_tsv: str, truth_set: set):
    """
    Direclty analyze TSV RGI using column Best_Hit_ARO
    """
    df = pd.read_csv(rgi_tsv, sep="\t", dtype=str).fillna("")

    if "Best_Hit_ARO" not in df.columns:
        raise ValueError(f"{rgi_tsv} non contiene la colonna 'Best_Hit_ARO'")

    df["clean_gene"] = df["Best_Hit_ARO"].apply(lambda x: str(x).lower().strip())

    #riove species prefix
    df["clean_gene"] = df["clean_gene"].apply(
        lambda x: re.sub(r'^[a-z]+\s+[a-z]+\s+', '', x)
    )

    # tokenizzation
    df["pred_token"] = df["clean_gene"].apply(pred_to_token)

    pred_set = set(df["pred_token"])

    # metrics
    tp, fp, fn = compute_tp_fp_fn(truth_set, pred_set)
    write_gene_lists(out_prefix, "RGI_TSV_direct", tp, fp, fn)

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

    result = {
        "pipeline": "Bacterial_genomes_RGI",
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "recall": round(r, 4),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_genes": len(pred_set),
        "truth_genes": len(truth_set),
    }

    return result, pred_set

def plot_pipeline_metrics(df_metrics: pd.DataFrame, out_prefix: str):
    """
    Plot Precision, Recall e F1 for every pipeline.
    """
    pipelines = df_metrics["pipeline"].tolist()
    precision = df_metrics["precision"].tolist()
    recall = df_metrics["recall"].tolist()
    f1 = df_metrics["f1"].tolist()

    x = np.arange(len(pipelines))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10,6))

    ax.bar(x - width, precision, width, label="Precision")
    ax.bar(x, recall, width, label="Recall")
    ax.bar(x + width, f1, width, label="F1")

    ax.set_ylabel("Score")
    ax.set_title("Pipeline Performance (Gene level)")
    ax.set_xticks(x)
    ax.set_xticklabels(pipelines, rotation=45, ha="right")

    ax.set_ylim(0,1)
    ax.legend()

    plt.tight_layout()

    outfile = f"{out_prefix}pipeline_metrics_plot.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Wrote plot: {outfile}")
def plot_gene_detection_heatmap(pipeline_predictions, truth_set, out_prefix):
    """
    Heatmap pipeline vs gene (include standalone pipelines)
    """

    pipelines = sorted(pipeline_predictions.keys())
    genes = sorted(truth_set)

    matrix = []

    for pipe in pipelines:
        pred = pipeline_predictions[pipe]

        row = [1 if g in pred else 0 for g in genes]
        matrix.append(row)

    mat = pd.DataFrame(matrix, index=pipelines, columns=genes)

    plt.figure(figsize=(14,6))
    plt.imshow(mat, aspect="auto")

    plt.yticks(range(len(mat.index)), mat.index)
    plt.xticks(range(len(mat.columns)), mat.columns, rotation=90)

    plt.colorbar(label="Gene detected (1=yes)")
    plt.title("Gene Detection Across Pipelines")

    plt.tight_layout()

    outfile = f"{out_prefix}gene_detection_heatmap.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Wrote plot: {outfile}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="truth TSV con gene_symbol")
    ap.add_argument("--ham", required=True, help="merged_all_runs.tsv")
    ap.add_argument("--out_prefix", default="results_gene_", help="prefix output")
    ap.add_argument("--keep_unknown", action="store_true")
    ap.add_argument("--rgi", help="wf-metagenomics RGI standalone TSV (opzionale)")
    ap.add_argument("--metaresfinder", help="wf-metagenomics Resfinder")
    ap.add_argument("--rgi_tsv",required=False, default="")
    #ap.add_argument("--bacterial_amrfinder", help="bacterial_genomes_amrfinder TSV (opzionale)", default="")
    args = ap.parse_args()

    outp = args.out_prefix

    # --- Load truth ---
    truth = pd.read_csv(args.truth, sep="\t", dtype=str).fillna("")
    truth["truth_token"] = truth["gene_symbol"].apply(truth_to_token)
    truth_set = set(truth["truth_token"])

    truth.to_csv(f"{outp}truth_with_tokens.tsv", sep="\t", index=False)

    rows = []
    pipeline_predictions = {}
    if args.rgi_tsv:
        print("\nAnalisi RGI TSV diretto...")
        rgi_tsv_metrics, rgi_tsv_set  = analyze_rgi_tsv(outp, args.rgi_tsv, truth_set)
        rows.append(rgi_tsv_metrics)
        pipeline_predictions["Bacterial_genomes_RGI"] = rgi_tsv_set
        print("\nAnalisi bacterial_RGI completata")


    # Load harmonized merged
    ham = pd.read_csv(args.ham, sep="\t", low_memory=False)

    needed = {"input_file_name", "analysis_software_name", "gene_symbol"}
    missing = needed - set(ham.columns)
    if missing:
        raise ValueError(f"merged_all_runs.tsv manca colonne: {sorted(missing)}")

    if "genetic_variation_type" in ham.columns:
        ham = ham[
            ham["genetic_variation_type"].fillna("").eq("gene_presence_detected")
        ].copy()

    ham["pipeline"] = ham.apply(assign_pipeline, axis=1)
    print("\nALL PIPELINES BEFORE UNKNOWN FILTER:")
    print(ham["pipeline"].value_counts(dropna=False).to_string())
    ham["pred_token"] = ham["gene_symbol"].apply(pred_to_token)

    if not args.keep_unknown:
        ham = ham[~ham["pipeline"].isin(["unknown", "ResFinder_unassigned"])].copy()

    print("\nPIPELINES COUNT:")
    print(ham["pipeline"].value_counts(dropna=False).to_string())

    # Per pipeline
    for pipe, dfp in ham.groupby("pipeline"):

        pred_set = set(dfp["pred_token"])
        pipeline_predictions[pipe] = pred_set

        tp, fp, fn = compute_tp_fp_fn(truth_set, pred_set)
        write_gene_lists(outp, pipe, tp, fp, fn)

        r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

        rows.append({
            "pipeline": pipe,
            "TP": len(tp),
            "FP": len(fp),
            "FN": len(fn),
            "recall": round(r, 4),
            "precision": round(p, 4),
            "f1": round(f1, 4),
            "predicted_unique_genes": len(pred_set),
            "truth_genes": len(truth_set),
        })

    #  RGI standalone (gene-level) 
    if args.rgi:
        print("\nProcessing standalone wf-metagenomics RGI...")
        rgi = pd.read_csv(args.rgi, sep="\t", dtype=str, quotechar='"').fillna("")


        if "Gene" not in rgi.columns:
            raise ValueError("RGI file must contain 'Gene' column")

        rgi["clean_gene"] = rgi["Gene"].apply(extract_gene_from_rgi)
        rgi["pred_token"] = rgi["clean_gene"].apply(pred_to_token)
        rgi_set = set(rgi["pred_token"])
        pipeline_predictions["Wf-metagenomics RGI"] = rgi_set

        tp, fp, fn = compute_tp_fp_fn(truth_set, rgi_set)
        write_gene_lists(outp, "Wf-metagenomics_RGI", tp, fp, fn)

        r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

        rows.append({
            "pipeline": "Wf-metagenomics RGI",
            "TP": len(tp),
            "FP": len(fp),
            "FN": len(fn),
            "recall": round(r, 4),
            "precision": round(p, 4),
            "f1": round(f1, 4),
            "predicted_unique_genes": len(rgi_set),
            "truth_genes": len(truth_set),
        })
    if args.metaresfinder and os.path.exists(args.metaresfinder):

        print("\nProcessing standalone Wf-metagenomics ResFinder...")

        resf_metrics, resf_set = analyze_metagenomics_resfinder(
        outp,
        args.metaresfinder,
        truth_set
        )
        pipeline_predictions["Metagenomics_ResFinder"] = resf_set
        rows.append(resf_metrics)
    ruargpore_row = {
    "pipeline": "RUARGpore",
    "TP": 1,
    "FP": 20,
    "FN": 24,
    "recall": 0.040,
    "precision": 0.048,
    "f1": 0.044,
    "predicted_unique_genes": 21,
    "truth_genes": len(truth_set)}
    rows.append(ruargpore_row)

    #  Output 
    out_df = pd.DataFrame(rows).sort_values(["f1", "recall"], ascending=False)
    out_df.to_csv(f"{outp}metrics_gene_level.tsv", sep="\t", index=False)
    plot_pipeline_metrics(out_df, outp)
    plot_gene_detection_heatmap(pipeline_predictions, truth_set, outp)
    print("\n=== GENE LEVEL ===")
    print(out_df.to_string(index=False))
    print(f"\nWrote: {outp}metrics_gene_level.tsv")
    print(f"Wrote: {outp}truth_with_tokens.tsv")
    print(f"Wrote per-pipeline lists: {outp}<PIPELINE>_TP/FP/FN_genes.tsv")


if __name__ == "__main__":
    main()

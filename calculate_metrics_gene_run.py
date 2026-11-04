#!/usr/bin/env python3
import argparse
import pandas as pd
from gene_token import truth_to_token, pred_to_token
import matplotlib.pyplot as plt
import numpy as np
import re
import os
def assign_pipeline(row) -> str:
    """Assegna pipeline basandosi su analysis_software_name e input_file_name"""
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
        return "Wf-metagenomics AMRFinderPlus"
    if infile == "bacterial_genomes_amrfinder":
        return "Bacterial_genomes_AMRFinderPlus"
    return "unknown"

def extract_gene_from_rgi(gene_str: str) -> str:
    """
    Estrae il gene RGI pulito:
    - Rimuove prefisso specie (qualsiasi cosa fino a "_" o nome specie)
    - Rimuove tutto da mutanti/mutation/withmutation + conferringresistance in poi
    - Rimuove caratteri non alfanumerici finali
    """
    gene_str = str(gene_str).lower().strip()
    gene_str = re.sub(r'[^a-z0-9\-\(\)]', '', gene_str)
    # Remove common prefix
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

    # Remove all after "mutationconferringresistance", "mutantsconferringresistance", "withmutationconferringresistance"
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


def plot_pipeline_metrics(df_metrics: pd.DataFrame, out_prefix: str):
    """
    Plot Precision, Recall e F1 per ogni pipeline.
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
    ax.set_title("Pipeline Performance (Drug-class level)")
    ax.set_xticks(x)
    ax.set_xticklabels(pipelines, rotation=45, ha="right")

    ax.set_ylim(0,1)
    ax.legend()

    plt.tight_layout()

    outfile = f"{out_prefix}pipeline_best_metrics_plot.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Wrote plot: {outfile}")
def analyze_metagenomics_resfinder(out_prefix: str, resfinder_file: str, truth_set: set,run_id: str):
    """
    Analizza un file ResFinder e calcola metriche gene-level.
    
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
        "fosB": "fosb",
        "fosB1": "fosb",
        "tet(45)": "tet(45)",
       
    }

    df["pred_token"] = df["gene_symbol_clean"].map(RESFINDER_ALIAS)

    # Remove not mapped rows
    df = df.dropna(subset=["pred_token"])
    pred_set = set(df["pred_token"])

    # --- calcolate TP, FP, FN ---
    tp = pred_set & truth_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set

    write_gene_lists(out_prefix, "Metagenomics_ResFinder", tp, fp, fn)

    # --- metrics ---
    r = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0
    p = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0
    f1 = 2 * r * p / (r + p) if (r + p) > 0 else 0

    result = {
        "sample": run_id,
        "pipeline": "Metagenomics_ResFinder",
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "TP_set": list(tp),
        "recall": round(r, 4),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_genes": len(pred_set),
        "truth_genes": len(truth_set),
    }
    
    return result
def analyze_rgi_tsv(out_prefix: str, rgi_tsv: str, truth_set: set,run_id: str):
    """
    Analizza direttamente TSV RGI (tipo nanorgi_only.tsv o simili)
    usando la colonna Best_Hit_ARO
    """
    if not rgi_tsv or not os.path.exists(rgi_tsv):
        print(f"[WARNING] RGI file not found: {rgi_tsv} → skipping")

    df = pd.read_csv(rgi_tsv, sep="\t", dtype=str).fillna("")

    if "Best_Hit_ARO" not in df.columns:
        raise ValueError(f"{rgi_tsv} non contiene la colonna 'Best_Hit_ARO'")

    df["clean_gene"] = df["Best_Hit_ARO"].apply(lambda x: str(x).lower().strip())

    # remove prefix
    df["clean_gene"] = df["clean_gene"].apply(
        lambda x: re.sub(r'^[a-z]+\s+[a-z]+\s+', '', x)
    )
    print(df["clean_gene"].head())
    # --- tokenization ---
    df["pred_token"] = df["clean_gene"].apply(pred_to_token)

    pred_set = set(df["pred_token"])

    # --- metrics ---
    tp, fp, fn = compute_tp_fp_fn(truth_set, pred_set)
    write_gene_lists(out_prefix, "bacterial_genomes_RGI", tp, fp, fn)

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

    result = {
        "pipeline": "Bacterial_genomes_RGI",
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "recall": round(r, 4),
        "TP_set": list(tp),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_genes": len(pred_set),
        "truth_genes": len(truth_set),
        "sample": run_id,
    }

    return result
def analyze_ruargpore_genelevel(out_prefix, ruarg_file, truth_set):

    ru = pd.read_csv(ruarg_file, sep="\t", dtype=str, low_memory=False).fillna("")

    if "subtype" not in ru.columns:
        raise ValueError("RUARGpore file must contain 'subtype' column")

    ru["pred_token"] = ru["subtype"].apply(pred_to_token)

    ru = ru.drop_duplicates(subset=["pred_token"])

    pred_set = set(ru["pred_token"])

    tp, fp, fn = compute_tp_fp_fn(truth_set, pred_set)

    write_gene_lists(out_prefix, "RUARGpore", tp, fp, fn)

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

    print("\n=== RUARGpore (gene-level) ===")
    print("TP genes:", sorted(tp))
    print("FP genes:", sorted(fp))
    print("FN genes:", sorted(fn))

    return {
        "pipeline": "RUARGpore",
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "TP_set": list(tp),
        "recall": round(r,4),
        "precision": round(p,4),
        "f1": round(f1,4),
        "predicted_unique_genes": len(pred_set),
        "truth_genes": len(truth_set)
    }

def plot_metrics_by_run(df: pd.DataFrame, out_prefix: str):
    """
    Plot come cambiano precision, recall e F1 run-per-run per ogni pipeline.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = ["precision", "recall", "f1"]

    for metric in metrics:

        plt.figure(figsize=(10,6))

        for pipeline in sorted(df["pipeline"].unique()):

            sub = df[df["pipeline"] == pipeline]

            plt.plot(
                sub["run"],
                sub[metric],
                marker="o",
                label=pipeline
            )

        plt.xlabel("Run")
        plt.ylabel(metric.capitalize())
        plt.title(f"{metric.capitalize()} per run")
        plt.ylim(0,1)

        plt.legend()
        plt.tight_layout()

        outfile = f"{out_prefix}{metric}_by_run.png"
        plt.savefig(outfile, dpi=300)
        plt.close()

        print(f"Wrote plot: {outfile}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="truth TSV con gene_symbol")
    ap.add_argument("--ham", required=True, help="merged_run{i}.tsv")
    ap.add_argument("--out_prefix", default="results_gene_", help="prefix output")
    ap.add_argument("--keep_unknown", action="store_true")
    ap.add_argument("--rgi", help="wf-metagenomics RGI standalone TSV")
    ap.add_argument("--metaresfinder", help="wf-metagenomics Resfinder")
    ap.add_argument("--ruargpore", help="RUARGpore .w.taxa.tab file", default="")
    ap.add_argument("--run_id", required=True)
    ap.add_argument("--rgi_tsv", help="wf-bacterial genomes rgi", required=False, default="")
    args = ap.parse_args()
    run_id = args.run_id
    outp = args.out_prefix

    # --- Load truth ---
    truth = pd.read_csv(args.truth, sep="\t", dtype=str).fillna("")
    truth["truth_token"] = truth["gene_symbol"].apply(truth_to_token)

    truth_set = set(truth["truth_token"])
    truth.to_csv(f"{outp}truth_with_tokens.tsv", sep="\t", index=False)

    # --- Load harmonized merged ---
    ham = pd.read_csv(args.ham, sep="\t", low_memory=False)

    needed = {"input_file_name", "analysis_software_name", "gene_symbol"}
    missing = needed - set(ham.columns)
    if missing:
        raise ValueError(f"merged.tsv manca colonne: {sorted(missing)}")

    if "genetic_variation_type" in ham.columns:
        ham = ham[ham["genetic_variation_type"].fillna("").eq("gene_presence_detected")].copy()

    ham["pipeline"] = ham.apply(assign_pipeline, axis=1)
    ham["pred_token"] = ham["gene_symbol"].apply(pred_to_token)

    if not args.keep_unknown:
        ham = ham[~ham["pipeline"].isin(["unknown", "ResFinder_unassigned"])].copy()

    print("\nPIPELINES COUNT:")
    print(ham["pipeline"].value_counts(dropna=False).to_string())

       # --- Per pipeline ---
    rows = []

    for pipe, dfp in ham.groupby("pipeline"):

        pred_set = set(dfp["pred_token"])

        tp, fp, fn = compute_tp_fp_fn(truth_set, pred_set)
        write_gene_lists(outp, pipe, tp, fp, fn)

        r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

        rows.append({
            "sample": run_id,
            "pipeline": pipe,
            "TP": len(tp),
            "FP": len(fp),
            "FN": len(fn),
            "TP_set": list(tp),
            "recall": round(r, 4),
            "precision": round(p, 4),
            "f1": round(f1, 4),
            "predicted_unique_genes": len(pred_set),
            "truth_genes": len(truth_set)
        })
    if args.metaresfinder and os.path.exists(args.metaresfinder):

        print("\nProcessing standalone Wf-metagenomics ResFinder...")

        resf_metrics = analyze_metagenomics_resfinder(
        outp,
        args.metaresfinder,
        truth_set,
        run_id
        )

        rows.append(resf_metrics)

    else:
        print(f"\nResFinder results not found for {run_id}, skipping...")

    if args.rgi:
        print("\nProcessing standalone wf-metagenomics RGI...")

        rgi = pd.read_csv(args.rgi, sep="\t", dtype=str, quotechar='"').fillna("")


        if "Gene" not in rgi.columns:
            raise ValueError("RGI file must contain 'Gene' column")

        # --- Extract clean gene and compare with ground truth ---
        rgi["clean_gene"] = rgi["Gene"].apply(extract_gene_from_rgi)
        rgi["pred_token"] = rgi["clean_gene"].apply(pred_to_token)
        rgi_set = set(rgi["pred_token"])

        tp, fp, fn = compute_tp_fp_fn(truth_set, rgi_set)
        write_gene_lists(outp, "Wf-metagenomics_RGI", tp, fp, fn)

        r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

        rows.append({
        "sample": run_id,
        "pipeline": "Metagenomics_RGI",
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "TP_set": list(tp),
        "recall": round(r, 4),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_genes": len(rgi_set),
        "truth_genes": len(truth_set)})
    if args.rgi_tsv and os.path.exists(args.rgi_tsv):
        print("\nAnalisi RGI TSV diretto...")
        rgi_tsv_metrics = analyze_rgi_tsv(outp, args.rgi_tsv, truth_set, run_id)
        rows.append(rgi_tsv_metrics)
        print("\nAnalisi bacterial_RGI completata")
    else:
        print(f"\nRGI TSV not found for {run_id}, skipping...")
    if args.ruargpore:
        print("\nProcessing RUARGpore...")

        row = analyze_ruargpore_genelevel(
        outp,
        args.ruargpore,
        truth_set)

        row["sample"] = run_id
        rows.append(row)

    # --- Final metrics table and plots---
    out_df = pd.DataFrame(rows).sort_values(["f1", "recall"], ascending=False)
    if "TP_set" in out_df.columns:
        out_df["TP_gene_list"] = out_df["TP_set"].apply(lambda x: ",".join(x) if isinstance(x, list) else "")
    out_df.to_csv(f"{outp}gene_metrics.tsv", sep="\t", index=False)
    plot_pipeline_metrics(out_df, outp)
    #summarize_across_runs(out_df, outp)


    #plot_metrics_by_run(out_df, outp)
    print("\n=== GENE LEVEL METRICS ===")
    print(out_df.to_string(index=False))
    print(f"\nWrote: {outp}gene_metrics.tsv")
    print(f"Wrote: {outp}truth_with_tokens.tsv")
    print(f"Wrote per-pipeline lists: {outp}<PIPELINE>_TP/FP/FN_genes.tsv")


if __name__ == "__main__":
    main()

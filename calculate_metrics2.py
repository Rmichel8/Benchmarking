#!/usr/bin/env python3
import argparse
import math
import pandas as pd

from gene_token import truth_to_token, pred_to_token
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import gene_token
print("USING gene_token.py from:", gene_token.__file__)
print("TEST pred_to_token(Klebsiella pneumoniae KpnE) =", pred_to_token("Klebsiella pneumoniae KpnE"))


def assign_pipeline(row) -> str:
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
    if infile == "nanorgi":
        return "NanoBsi_RGI"
    if infile== "amrfinderplus":
        return "wf-metagenomics AMRFinderPlus"
    if infile == "bacterial_genomes_amrfinder":
        return "Bacterial_genomes_AMRFinderPlus"
    return "unknown"

def compute_confusion(truth_set: set, pred_set: set):
    tp = truth_set & pred_set
    fn = truth_set - pred_set
    fp = pred_set - truth_set
    tn = 0
    return tp, fp, fn, tn

def metrics_from_counts(TP, FP, FN, TN):

    # Recall (Sensitivity)
    if (TP + FN) != 0:
        recall = TP / (TP + FN)
    else:
        recall = 0.0

    # Precision
    if (TP + FP) != 0:
        precision = TP / (TP + FP)
    else:
        precision = 0.0

    # F1 score
    if (precision + recall) != 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = 0.0

    return recall, precision, f1


def safe_name(s: str) -> str:
    # filename-safe
    return "".join([c if (c.isalnum() or c in "-_.") else "_" for c in str(s)])


def write_gene_lists(out_prefix: str, pipeline: str, tp: set, fp: set, fn: set):
    p = safe_name(pipeline)
    pd.DataFrame({"gene_token": sorted(tp)}).to_csv(f"{out_prefix}{p}_TP.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_token": sorted(fp)}).to_csv(f"{out_prefix}{p}_FP.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_token": sorted(fn)}).to_csv(f"{out_prefix}{p}_FN.tsv", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="FPES_mock_gene_truth.tsv")
    ap.add_argument("--ham", required=True, help="merged_all_runs.tsv (hAMRonized)")
    ap.add_argument("--card_aro", required=True, help="CARD aro.tsv")
    ap.add_argument("--out_prefix", default="results_", help="prefix output")
    ap.add_argument("--keep_unknown", action="store_true")
    args = ap.parse_args()

    outp = args.out_prefix  

    # Load truth 
    truth = pd.read_csv(args.truth, sep="\t", dtype=str).fillna("")
    if "gene_symbol" not in truth.columns:
        raise ValueError("Truth deve avere colonna 'gene_symbol'.")
    truth["truth_token"] = truth["gene_symbol"].apply(truth_to_token)
    truth_gene_set = set(truth["truth_token"].dropna().unique())

    # Load ham 
    ham = pd.read_csv(args.ham, sep="\t", low_memory=False)
    needed = {"input_file_name", "analysis_software_name", "gene_symbol"}
    missing = needed - set(ham.columns)
    if missing:
        raise ValueError(f"merged_all_runs.tsv manca colonne: {sorted(missing)}")

    if "genetic_variation_type" in ham.columns:
        ham = ham[ham["genetic_variation_type"].fillna("").eq("gene_presence_detected")].copy()

    ham["pipeline"] = ham.apply(assign_pipeline, axis=1)
    ham["pred_token"] = ham["gene_symbol"].apply(pred_to_token)

    if not args.keep_unknown:
        ham = ham[~ham["pipeline"].isin(["unknown", "ResFinder_unassigned"])].copy()

    print("\nPIPELINES COUNT:")
    print(ham["pipeline"].value_counts(dropna=False).to_string())


    # Compute per pipeline + write TP/FP/FN lists 
    rows_gene = []
    for pipe, dfp in ham.groupby("pipeline"):
        pred_gene_set = set(dfp["pred_token"].dropna().unique())
        tp, fp, fn, tn = compute_confusion(truth_gene_set, pred_gene_set)
        write_gene_lists(outp, pipe, tp, fp, fn)

        recall, precision, f1 = metrics_from_counts(len(tp), len(fp), len(fn), tn)
        rows_gene.append({
            "pipeline": pipe,
            "TP": len(tp), "FP": len(fp), "FN": len(fn), "TN": tn,
            "recall": round(recall, 4),
            "precision": round(precision, 4),
            "f1": round(f1, 4),
            "predicted_unique_genes": len(pred_gene_set),
            "truth_genes": len(truth_gene_set),
        })

    gene_df = pd.DataFrame(rows_gene).sort_values([ "f1", "recall"], ascending=False)
    gene_df.to_csv(f"{outp}metrics_gene_level.tsv", sep="\t", index=False)

    print("\n=== GENE LEVEL ===")
    print(gene_df.to_string(index=False))
    print(f"\nWrote: {outp}metrics_gene_level.tsv")
    print(f"Wrote per-pipeline lists: {outp}<PIPELINE>_TP.tsv, {outp}<PIPELINE>_FP.tsv, {outp}<PIPELINE>_FN.tsv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from gene_token import truth_to_token, pred_to_token
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


def analyze_nanobsi_rgi_drugclass(out_prefix: str, rgi_file: str, truth_class_set: set):
    """
    Analizza NanoBsi RGI standalone a livello drug-class.
    """
    df = pd.read_csv(rgi_file, sep="\t", dtype=str).fillna("")
    if "drug_class" not in df.columns:
        raise ValueError(f"{rgi_file} non contiene la colonna 'drug_class'")

    pred_class_set = explode_drug_classes(df["drug_class"])
    
    tp, fp, fn = compute_tp_fp_fn(truth_class_set, pred_class_set)
    write_class_lists(out_prefix, "NanoBsi_RGI", tp, fp, fn)

    print("\n=== NanoBsi RGI (drug-class) ===")
    print("TP drug classes:", sorted(tp))
    print("FP drug classes:", sorted(fp))

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))
    return {
        "pipeline": "NanoBsi_RGI",
        "TP": len(tp), "FP": len(fp), "FN": len(fn),
        "recall": round(r, 4), "precision": round(p, 4), "f1": round(f1, 4),
        "predicted_unique_drug_classes": len(pred_class_set),
        "truth_drug_classes": len(truth_class_set),
    }
def analyze_bacterial_genomes_rgi_drugclass(out_prefix: str, rgi_file: str, truth_class_set: set):
    """
    Analizza RGI (wf-bacterial genomes) a livello drug-class
    con gestione robusta (file missing, skip).
    """

    # --- CHECK FILE ---
    if not rgi_file or not os.path.exists(rgi_file):
        print(f"[WARNING] Bacterial RGI file not found: {rgi_file} → skipping")
        return {
            "pipeline": "Bacterial_genomes_RGI",
            "TP": 0, "FP": 0, "FN": 0,
            "recall": 0.0,
            "precision": 0.0,
            "f1": 0.0,
            "predicted_unique_drug_classes": 0,
            "truth_drug_classes": len(truth_class_set)
        }

    df = pd.read_csv(rgi_file, sep="\t", dtype=str).fillna("")

    # CHECK COLUMNS
    if "Drug Class" not in df.columns:
        raise ValueError(f"{rgi_file} non contiene la colonna 'drug_class'")

    # CLASS EXTRACTION
    pred_class_set = explode_drug_classes(df["Drug Class"])

    # METRICS
    tp, fp, fn = compute_tp_fp_fn(truth_class_set, pred_class_set)

    write_class_lists(out_prefix, "Bacterial_genomes_RGI", tp, fp, fn)

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))

    print("\n=== Bacterial Genomes RGI (drug-class) ===")
    print("TP drug classes:", sorted(tp))
    print("FP drug classes:", sorted(fp))
    print("FN drug classes:", sorted(fn))

    return {
        "pipeline": "Bacterial_genomes_RGI",
        "TP_set": tp,
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "recall": round(r, 4),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_drug_classes": len(pred_class_set),
        "truth_drug_classes": len(truth_class_set),
    }
def analyze_metagenomics_resfinder_drugclass(out_prefix: str, resfinder_file: str, truth_class_set: set):
    """
    Analizza Wf-metagenomics ResFinder standalone a livello drug-class.
    """
    df = pd.read_csv(resfinder_file, sep="\t", dtype=str).fillna("")
    if "Phenotype" not in df.columns:
        raise ValueError(f"{resfinder_file} non contiene la colonna 'Phenotype'")

    pred_class_set = explode_drug_classes(df["Phenotype"])
    
    tp, fp, fn = compute_tp_fp_fn(truth_class_set, pred_class_set)
    write_class_lists(out_prefix, "Metagenomics_ResFinder", tp, fp, fn)

    print("\n=== Wf-metagenomics ResFinder (drug-class) ===")
    print("TP drug classes:", sorted(tp))
    print("FP drug classes:", sorted(fp))

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))
    return {
        "pipeline": "Metagenomics_ResFinder",
        "TP_set": tp,
        "TP": len(tp), "FP": len(fp), "FN": len(fn),
        "recall": round(r, 4), "precision": round(p, 4), "f1": round(f1, 4),
        "predicted_unique_drug_classes": len(pred_class_set),
        "truth_drug_classes": len(truth_class_set),
    }

def analyze_metagenomics_rgi_drugclass(out_prefix: str, rgi_file: str, truth_class_set: set):
    """
    Analizza wf-metagenomics RGI standalone a livello drug-class.
    """
    df = pd.read_csv(rgi_file, sep="\t", dtype=str).fillna("")
    if "drug_class" not in df.columns:
        raise ValueError(f"{rgi_file} non contiene la colonna 'drug_class'")

    pred_class_set = explode_drug_classes(df["drug_class"])
    
    tp, fp, fn = compute_tp_fp_fn(truth_class_set, pred_class_set)
    write_class_lists(out_prefix, "Metagenomics_RGI", tp, fp, fn)

    print("\n=== Wf-metagenomics RGI (drug-class) ===")
    print("TP drug classes:", sorted(tp))
    print("FP drug classes:", sorted(fp))

    r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))
    return {
        "pipeline": "Metagenomics_RGI",
        "TP_set": tp,
        "TP": len(tp), "FP": len(fp), "FN": len(fn),
        "recall": round(r, 4),
        "precision": round(p, 4),
        "f1": round(f1, 4),
        "predicted_unique_drug_classes": len(pred_class_set),
        "truth_drug_classes": len(truth_class_set),
    }

def canon_dc(s: str) -> str:
    s = str(s).strip().lower()
    greek_to_ascii = {
        "β": "beta",
        "α": "alpha",
        "γ": "gamma",
        "δ": "delta",
    }
    for k, v in greek_to_ascii.items():
        s = s.replace(k, v)
    s = s.replace(",", ";")
    s = re.sub(r"\s+", " ", s).strip()
    return s
EXCLUDE_DRUG_CLASSES = {
    canon_dc("disinfecting agents and antiseptics")} 
    
def normalize_dc(dc: str) -> str:
    """
    Normalizza vocabolari diversi (tool) verso un set comparabile con la truth.
    Non inventa classi: quando non c'è equivalenza sensata, lascia invariato (quindi resta FP).
    """
    dc = canon_dc(dc)
    mapping = {
        # short form ->truth
        "aminoglycoside": "aminoglycoside antibiotic",
        "fluoroquinolone": "fluoroquinolone antibiotic",
        "quinolone": "fluoroquinolone antibiotic",  # NanoBsi ResFinder uses 'quinolone'
        "macrolide": "macrolide antibiotic",
        "peptide": "peptide antibiotic",
        "polymyxin": "peptide antibiotic", #subclass
        "glycopeptide antibiotic": "peptide antibiotic",
        "phosphonic acid antibiotic": "fosfomycin",

        # beta-lactams
        "cephamycin": "cephalosporin",
        "monobactam": "penicillin beta-lactam",
        "penam": "penicillin beta-lactam",
        "penem": "carbapenem",

        "elfamycin antibiotic": "pulvomycin",

        # tetraciclines
        "glycylcycline": "tetracycline antibiotic",
        # streptogramin family
        "streptogramin b antibiotic": "streptogramin antibiotic",
        "streptogramin a antibiotic": "streptogramin antibiotic",
        # phenicol/amphenicol
        "amphenicol": "phenicol antibiotic",
        #beta-lactam to penicilin
        "beta-lactam": "penicillin beta-lactam",
        # common variants
        "beta lactam": "penicillin beta-lactam",
        "penicillin": "penicillin beta-lactam",
        "tetracycline": "tetracycline antibiotic",
    }

    return mapping.get(dc, dc)



def explode_drug_classes(series: pd.Series) -> set:
    out = set()
    for x in series.fillna("").astype(str):
        x = canon_dc(x)
        if not x:
            continue
        for p in x.split(";"):
            p = p.strip()
            if not p:
                continue
            norm = normalize_dc(p)
            if norm in EXCLUDE_DRUG_CLASSES:
                continue
            out.add(norm)

    return out


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


#Manual mapping: truth gene token -> drug classes 
TRUTH_GENE_TO_DRUGCLASS = {
    "bcii": {"penicillin beta-lactam", "cephalosporin"},
    "fosb": {"fosfomycin"},
    "mcr-4.5": {"peptide antibiotic"},
    "tet(45)": {"tetracycline antibiotic", "tetracycline"},
    "crp": {"fluoroquinolone antibiotic", "macrolide antibiotic", "penicillin beta-lactam"},
    "ef-tu": {"pulvomycin"},
    "pbp3": {"cephalosporin", "penicillin beta-lactam"},
    "kpnf": {"tetracycline antibiotic","aminoglycoside antibiotic", "macrolide antibiotic","cephalosporin","rifamycin antibiotic"},
    "kpnh": {"fluoroquinolone antibiotic", "tetracycline antibiotic", "penicillin beta-lactam", "carbapenem","aminoglycoside antibiotic", "macrolide antibiotic","cephalosporin"},
    "adef": {"fluoroquinolone antibiotic", "tetracycline antibiotic"},
    "emrr": {"fluoroquinolone antibiotic"},
    "msba": {"nitroimidazole antibiotic"},
    "abaq": {"fluoroquinolone antibiotic"},
    "soxr": {"fluoroquinolone antibiotic","tetracycline antibiotic", "penicillin beta-lactam"},
    "arma": {"aminoglycoside antibiotic"},
    "mcr-4.1": {"peptide antibiotic"},
    "sgm": {"aminoglycoside antibiotic"},
    "carb-23": {"penicillin beta-lactam"},
    "amph": {"penicillin beta-lactam"},
    "gyrb": {"fluoroquinolone antibiotic"},
    "pmrf": {"peptide antibiotic"},
    "bes-1": {"penicillin beta-lactam"},
    "uhpt": {"fosfomycin"},
    "kpne": { "peptide antibiotic", "tetracycline antibiotic","aminoglycoside antibiotic", "macrolide antibiotic","cephalosporin","rifamycin antibiotic"},
    "amrb": {"aminoglycoside antibiotic", "macrolide antibiotic"}
}


def truth_token_to_drugclass(token: str) -> str:
    token = str(token).strip().lower()
    classes = TRUTH_GENE_TO_DRUGCLASS.get(token, set())
    classes = {normalize_dc(c) for c in classes}
    return ";".join(sorted(classes))


def safe_name(s: str) -> str:
    return "".join([c if (c.isalnum() or c in "-_.") else "_" for c in str(s)])


def write_class_lists(out_prefix: str, pipeline: str, tp: set, fp: set, fn: set):
    p = safe_name(pipeline)
    pd.DataFrame({"drug_class": sorted(tp)}).to_csv(f"{out_prefix}{p}_TP_drugclass.tsv", sep="\t", index=False)
    pd.DataFrame({"drug_class": sorted(fp)}).to_csv(f"{out_prefix}{p}_FP_drugclass.tsv", sep="\t", index=False)
    pd.DataFrame({"drug_class": sorted(fn)}).to_csv(f"{out_prefix}{p}_FN_drugclass.tsv", sep="\t", index=False)
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

    outfile = f"{out_prefix}pipeline_metrics_plot.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Wrote plot: {outfile}")
def plot_drugclass_detection_heatmap(results, truth_class_set, out_prefix):
    """
    Heatmap drug-class detection:
        1 = drug class della truth rilevata (TP)
        0 = non rilevata
    """

    pipelines = []
    detection = {}

    for res in results:
        pipe = res["pipeline"]
        pipelines.append(pipe)

        tp_set = set(res.get("TP_set", []))

        for dc in truth_class_set:
            detection.setdefault(dc, {})
            detection[dc][pipe] = 1 if dc in tp_set else 0

    heatmap_df = pd.DataFrame(detection).T.fillna(0)

    # order columns
    heatmap_df = heatmap_df[pipelines]

    plt.figure(figsize=(12, 8))

    plt.imshow(heatmap_df.values, aspect="auto")

    plt.colorbar(label="Detection (1=TP, 0=missed)")

    plt.yticks(
        np.arange(len(heatmap_df.index)),
        heatmap_df.index
    )

    plt.xticks(
        np.arange(len(heatmap_df.columns)),
        heatmap_df.columns,
        rotation=45,
        ha="right"
    )

    plt.title("Drug-class Detection Heatmap")

    plt.tight_layout()

    outfile = f"{out_prefix}drugclass_detection_heatmap.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

    print(f"Wrote heatmap: {outfile}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="truth TSV con gene_symbol/present/copy_number")
    ap.add_argument("--ham", required=True, help="merged_all_runs.tsv (hAMRonized)")
    ap.add_argument("--out_prefix", default="results_drug_", help="prefix output")
    ap.add_argument("--keep_unknown", action="store_true")
    ap.add_argument("--metargi", help="wf-metagenomics RGI standalone TSV (opzionale)")
    ap.add_argument("--metaresfinder", help="Wf-metagenomics ResFinder TSV (opzionale)", default="")
    ap.add_argument("--rgi_tsv", help="wf-bacterial genomes rgi", required=False, default="")
    args = ap.parse_args()

    outp = args.out_prefix

    # Load truth 
    truth = pd.read_csv(args.truth, sep="\t", dtype=str).fillna("")
    if "gene_symbol" not in truth.columns:
        raise ValueError("Truth deve avere colonna 'gene_symbol'.")

    truth["truth_token"] = truth["gene_symbol"].apply(truth_to_token)
    truth["truth_drug_class"] = truth["truth_token"].apply(truth_token_to_drugclass)
    truth.to_csv(f"{outp}truth_with_drugclass.tsv", sep="\t", index=False)

    truth_class_set = explode_drug_classes(truth["truth_drug_class"])

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

    if "drug_class" in ham.columns:
        ham["drug_class_norm_set"] = ham["drug_class"].fillna("").astype(str).apply(
            lambda x: ";".join(sorted(explode_drug_classes(pd.Series([x]))))
        )
    else:
        ham["drug_class_norm_set"] = ""

    if not args.keep_unknown:
        ham = ham[~ham["pipeline"].isin(["unknown", "ResFinder_unassigned"])].copy()

    print("\nPIPELINES COUNT:")
    print(ham["pipeline"].value_counts(dropna=False).to_string())

    print("\nTruth drug classes (n={}):".format(len(truth_class_set)))
    print(sorted(truth_class_set))

    # --- Per pipeline ---
    rows = []
    # pipeline standalone
    if args.metaresfinder and os.path.exists(args.metaresfinder):
        rows.append(analyze_metagenomics_resfinder_drugclass(outp, args.metaresfinder, truth_class_set))
    if args.metargi:
        rows.append(analyze_metagenomics_rgi_drugclass(outp, args.metargi, truth_class_set))
    if args.rgi_tsv and os.path.exists(args.rgi_tsv):
        rows.append(
        analyze_bacterial_genomes_rgi_drugclass(
            outp,
            args.rgi_tsv,
            truth_class_set
        )
        )
    for pipe, dfp in ham.groupby("pipeline"):
        pred_class_set = explode_drug_classes(dfp.get("drug_class", pd.Series([], dtype=str)))

        tp, fp, fn = compute_tp_fp_fn(truth_class_set, pred_class_set)
        write_class_lists(outp, pipe, tp, fp, fn)

        r, p, f1 = metrics_from_counts(len(tp), len(fp), len(fn))
        rows.append({
            "pipeline": pipe,
            "TP": len(tp), "FP": len(fp), "FN": len(fn),
            "recall": round(r, 4),
            "precision": round(p, 4),
            "TP_set": tp,
            "f1": round(f1, 4),
            "predicted_unique_drug_classes": len(pred_class_set),
            "truth_drug_classes": len(truth_class_set),
        })
    ruargpore_tp = {
    normalize_dc("aminoglycoside antibiotic"),
    normalize_dc("carbapenem"),
    normalize_dc("cephalosporin"),
    normalize_dc("fluoroquinolone antibiotic"),
    normalize_dc("fosfomycin"),
    normalize_dc("macrolide antibiotic"),
    normalize_dc("penicillin beta-lactam"),
    normalize_dc("peptide antibiotic"),
    normalize_dc("rifamycin antibiotic"),
    normalize_dc("tetracycline antibiotic"),}
    ruargpore_tp = ruargpore_tp & truth_class_set
    ruargpore_fp = set()
    ruargpore_row = {
    "pipeline": "RUARGpore", "TP_set": ruargpore_tp,
    "TP": 9, "FP": 4, "FN": 3, "recall": 0.7500, "precision": 0.6923, "f1": 0.7200, "predicted_unique_drug_classes": len(pred_class_set),
    "truth_drug_classes": len(truth_class_set),
    }
    rows.append(ruargpore_row)
    out_df = pd.DataFrame(rows).sort_values(["f1", "recall"], ascending=False)
    out_df.to_csv(f"{outp}metrics_drugclass_level.tsv", sep="\t", index=False)
    plot_pipeline_metrics(out_df, outp)
    plot_drugclass_detection_heatmap(rows, truth_class_set, outp)
    print("\n=== DRUG-CLASS LEVEL ===")
    print(out_df.to_string(index=False))
    print(f"\nWrote: {outp}metrics_drugclass_level.tsv")
    print(f"Wrote: {outp}truth_with_drugclass.tsv")
    print(f"Wrote per-pipeline lists: {outp}<PIPELINE>_TP_drugclass.tsv / FP / FN")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from pathlib import Path
from gene_token import truth_to_token

def canon_dc(s: str) -> str:
    s = str(s).strip().lower()
    s = s.replace(",", ";")
    s = re.sub(r"\s+", " ", s).strip()
    return s

def normalize_dc(dc: str) -> str:
    dc = canon_dc(dc)
    mapping = {
        "aminoglycoside": "aminoglycoside antibiotic",
        "fluoroquinolone": "fluoroquinolone antibiotic",
        "quinolone": "fluoroquinolone antibiotic",
        "macrolide": "macrolide antibiotic",
        "peptide": "peptide antibiotic",
        "polymyxin": "polymyxin",
        "cephamycin": "cephalosporin",
        "monobactam": "penicillin beta-lactam",
        "penam": "penicillin beta-lactam",
        "penem": "carbapenem",
        "elfamycin antibiotic": "pulvomycin",
        "glycylcycline": "tetracycline antibiotic",
        "tetracycline": "tetracycline antibiotic",
        "streptogramin b antibiotic": "streptogramin antibiotic",
        "streptogramin a antibiotic": "streptogramin antibiotic",
        "amphenicol": "phenicol antibiotic",
        "beta-lactam": "penicillin beta-lactam",
        "beta lactam": "penicillin beta-lactam",
        "penicillin": "penicillin beta-lactam",
        "fosfomycin": "fosfomycin",
        "phosphonic acid antibiotic": "phosphonic acid antibiotic",
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
            if p == "macrolide-lincosamide-streptogramin":
                out.update(["macrolide antibiotic", "lincosamide antibiotic", "streptogramin antibiotic"])
                continue
            out.add(normalize_dc(p))
    return out

# Truth mapping

TRUTH_GENE_TO_DRUGCLASS = {
    "bcii": {"penicillin beta-lactam", "cephalosporin"},
    "fosb": {"phosphonic acid antibiotic", "fosfomycin"},
    "mcr-4.5": {"peptide antibiotic", "polymyxin"},
    "tet(45)": {"tetracycline antibiotic"},
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
    "pmrf": {"peptide antibiotic", "polymyxin"},
    "bes-1": {"penicillin beta-lactam"},
    "uhpt": {"phosphonic acid antibiotic", "fosfomycin"},
    "kpne": {"peptide antibiotic", "tetracycline antibiotic","aminoglycoside antibiotic", "macrolide antibiotic","cephalosporin","rifamycin antibiotic"},
    "amrb": {"aminoglycoside antibiotic", "macrolide antibiotic"}
}

CARD_FP_GENE_TO_CLASSES = {
    "acra": {"tetracycline antibiotic", "penicillin beta-lactam", "cephalosporin", "disinfecting agents and antiseptics", "phenicol antibiotic", "rifamycin antibiotic", "glycylcycline", "fluoroquinolone antibiotic"},
    "acrb": {"fluoroquinolone antibiotic", "penicillin beta-lactam", "cephalosporin","disinfecting agents and antiseptics", "tetracycline antibiotic", "rifamycin antibiotic","phenicol antibiotic"},
    "arna": {"peptide antibiotic"},
    "baca": {"peptide antibiotic"},
    "bcra": {"peptide antibiotic"},
    "bicyclomycin-multidrug_efflux_protein_bcr": {"bicyclomycin-like antibiotic"},
    "ceob": {"aminoglycoside antibiotic", "fluoroquinolone antibiotic"},
    "classabeta-lactamase": {"penicillin beta-lactam"},
    "classbbeta-lactamase": {"penicillin beta-lactam"},
    "classcbeta-lactamase": {"penicillin beta-lactam"},
    "emra": {"fluoroquinolone antibiotic"},
    "emrb": {"fluoroquinolone antibiotic"},
    "emrb-qacafamilymajorfacilitatortransporter": {"fluoroquinolone antibiotic"},
    "emrd": {"phenicol antibiotic", "disinfecting agents and antiseptics"},
    "kasugamycinresistanceproteinksga": {"aminoglycoside antibiotic"},
    "maca": {"macrolide antibiotic"},
    "macb": {"macrolide antibiotic"},
    "mdfa": {"disinfecting agents and antiseptics", "phenicol antibiotic", "tetracycline antibiotic"},
    "mdta": {"aminocoumarin antibiotic"},
    "mdtb": {"aminocoumarin antibiotic"},
    "mdtc": {"aminocoumarin antibiotic"},
    "mdtg": {"phosphonic acid antibiotic"},
    "mdth": {"fluoroquinolone antibiotic"},
    "mdtk": {"fluoroquinolone antibiotic"},
    "mexd": {"phenicol antibiotic","diaminopyrimidine antibiotic","aminocoumarin antibiotic", "tetracycline antibiotic", "aminoglycoside antibiotic", "penicillin beta-lactam", "cephalosporin","fluoroquinolone antibiotic", "macrolide antibiotic"},
    "mexe": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "mexf": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "mext": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "mexw": {"disinfecting agents and antiseptics", "phenicol antibiotic", "tetracycline antibiotic", "macrolide antibiotic", "fluoroquinolone antibiotic"},
    "omp36": {"cephalosporin", "penicillin beta-lactam", "carbapenem", "monobactam"},
    "ompf": {"cephalosporin", "penicillin beta-lactam", "carbapenem", "monobactam"},
    "ompr": {"cephalosporin", "penicillin beta-lactam", "carbapenem", "monobactam"},
    "oprn": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "qnrb": {"fluoroquinolone antibiotic"},
    "rosa": {"peptide antibiotic"},
    "rosb": {"peptide antibiotic"},
    "rpsd(ramaorsud2)": {"aminoglycoside antibiotic"},
    "tet34": {"tetracycline antibiotic"},
    "tetl": {"tetracycline antibiotic"},
    "tolc": {"peptide antibiotic", "tetracycline antibiotic", "penicillin beta-lactam", "cephalosporin", "aminoglycoside antibiotic", "disinfecting agents and antiseptics", "phenicol antibiotic", "rifamycin antibiotic", "aminocoumarin antibiotic", "glycylcycline", "carbapenem", "fluoroquinolone antibiotic", "macrolide antibiotic"},
    "transcriptionalregulatoryproteincpxrcpxr": {"peptide antibiotic", "sulfonamide antibiotic", "diaminopyrimidine antibiotic", "aminocoumarin antibiotic", "penicillin beta-lactam", "cephalosporin", "aminoglycoside antibiotic", "phenicol antibiotic", "tetracycline antibiotic", "carbapenem", "monobactam", "fluoroquinolone antibiotic", "macrolide antibiotic"},
    "vanr": {"vancomycin"},
    "vans": {"vancomycin"},
    "vany": {"vancomycin"},
    "vatg": {"vancomycin"}
}

# Normalizes CARD keys
CARD_FP_GENE_TO_CLASSES_TOKEN = {
    truth_to_token(k): {normalize_dc(c) for c in v} for k, v in CARD_FP_GENE_TO_CLASSES.items()
}

# Metrics

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

# Load RUARG raw

def load_ruarg_raw(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    needed = {"subtype", "type"}
    if not needed.issubset(set(df.columns)):
        raise ValueError(f"File RUARG raw deve avere colonne {sorted(needed)}")
    out = pd.DataFrame({
        "gene_symbol": df["subtype"].astype(str)
    })
    out["gene_symbol"] = out["gene_symbol"].str.strip()
    out = out[out["gene_symbol"] != ""].copy()
    return out

# Main

def main():
    ap = argparse.ArgumentParser(description="RUARGpore metrics (gene + drug-class from truth mapping)")
    ap.add_argument("--truth", required=True, help="FPES_moch_gene_truth.tsv (gene_symbol,present,copy_number)")
    ap.add_argument("--ruarg_raw", required=True, help="RUARGpore raw file")
    ap.add_argument("--out_dir", default=".", help="Output directory")
    ap.add_argument("--out_prefix", default="ruargpore", help="Prefix output files")
    args = ap.parse_args()

    out_dir = args.out_dir.rstrip("/")
    outp = args.out_prefix.rstrip("_") + "_"
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    #Load truth
    truth = pd.read_csv(args.truth, sep="\t", dtype=str).fillna("")
    truth = truth[truth["present"].astype(str).eq("1")].copy()
    truth["truth_token"] = truth["gene_symbol"].apply(truth_to_token)

    # map drug classes from truth
    def token_to_drugclass(token):
        classes = TRUTH_GENE_TO_DRUGCLASS.get(token, set())
        return ";".join(sorted(normalize_dc(c) for c in classes))
    truth["truth_drug_class"] = truth["truth_token"].apply(token_to_drugclass)
    truth_with_file = f"{out_dir}/{outp}truth_with_tokens_and_drugclass.tsv"
    truth.to_csv(truth_with_file, sep="\t", index=False)

    truth_gene_set = set(truth["truth_token"].tolist())
    truth_class_set = explode_drug_classes(truth["truth_drug_class"])

    # Load RUARG 
    ru = load_ruarg_raw(args.ruarg_raw)
    ru["pred_gene_token"] = ru["gene_symbol"].apply(truth_to_token)
    ru = ru.drop_duplicates(subset=["pred_gene_token"])
    pred_gene_set = set(ru["pred_gene_token"].tolist())

    # Gene-level metrics 
    tp_g, fp_g, fn_g = compute_tp_fp_fn(truth_gene_set, pred_gene_set)
    r_g, p_g, f1_g = metrics_from_counts(len(tp_g), len(fp_g), len(fn_g))

    #Drug-class-level metrics 
    tp_class_set = set()
    pred_class_set = set()
    for gene in ru["pred_gene_token"]:
        if gene in TRUTH_GENE_TO_DRUGCLASS:
            classes = {normalize_dc(c) for c in TRUTH_GENE_TO_DRUGCLASS[gene]}
            tp_class_set.update(classes)
            pred_class_set.update(classes)
    # Integrate CARD classes for FP genes
    for gene in fp_g:
        if gene in CARD_FP_GENE_TO_CLASSES_TOKEN:
            pred_class_set.update(CARD_FP_GENE_TO_CLASSES_TOKEN[gene])

    fn_class_set = truth_class_set - tp_class_set
    fp_class_set = pred_class_set - truth_class_set
    r_c, p_c, f1_c = metrics_from_counts(len(tp_class_set), len(fp_class_set), len(fn_class_set))

    # Write outputs 
    gene_metrics = pd.DataFrame([{
        "pipeline": "RUARGpore_raw",
        "TP": len(tp_g),
        "FP": len(fp_g),
        "FN": len(fn_g),
        "recall": round(r_g,4),
        "precision": round(p_g,4),
        "f1": round(f1_g,4),
        "truth_unique_genes": len(truth_gene_set),
        "pred_unique_genes": len(pred_gene_set),
        "n_rows_ruargpore": len(ru),
        "ruarg_raw_file": args.ruarg_raw,
    }])

    drug_metrics = pd.DataFrame([{
        "pipeline": "RUARGpore_raw",
        "TP": len(tp_class_set),
        "FP": len(fp_class_set),
        "FN": len(fn_class_set),
        "recall": round(r_c,4),
        "precision": round(p_c,4),
        "f1": round(f1_c,4),
        "truth_unique_drug_classes": len(truth_class_set),
        "pred_unique_drug_classes": len(pred_class_set),
        "n_rows_ruargpore": len(ru),
        "ruarg_raw_file": args.ruarg_raw,
    }])

    # lists
    pd.DataFrame({"gene_token": sorted(tp_g)}).to_csv(f"{out_dir}/{outp}gene_TP.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_token": sorted(fp_g)}).to_csv(f"{out_dir}/{outp}gene_FP.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_token": sorted(fn_g)}).to_csv(f"{out_dir}/{outp}gene_FN.tsv", sep="\t", index=False)

    pd.DataFrame({"drug_class": sorted(tp_class_set)}).to_csv(f"{out_dir}/{outp}drugclass_TP.tsv", sep="\t", index=False)
    pd.DataFrame({"drug_class": sorted(fp_class_set)}).to_csv(f"{out_dir}/{outp}drugclass_FP.tsv", sep="\t", index=False)
    pd.DataFrame({"drug_class": sorted(fn_class_set)}).to_csv(f"{out_dir}/{outp}drugclass_FN.tsv", sep="\t", index=False)

    print("\n=== GENE-LEVEL METRICS ===")
    print(gene_metrics.to_string(index=False))
    print("\n=== DRUG-CLASS-LEVEL METRICS ===")
    print(drug_metrics.to_string(index=False))
    print(f"\nWrote: {truth_with_file}")

if __name__ == "__main__":
    main()

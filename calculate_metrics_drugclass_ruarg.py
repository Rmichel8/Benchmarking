#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from pathlib import Path

def canon_dc(s: str) -> str:
    s = str(s).strip().lower()
    s = s.replace(",", ";")
    s = re.sub(r"\s+", " ", s).strip()
    return s

def normalize_dc(dc: str) -> str:
    mapping = {
        "aminoglycoside": "aminoglycoside antibiotic",
        "fluoroquinolone": "fluoroquinolone antibiotic",
        "quinolone": "fluoroquinolone antibiotic",
        "macrolide": "macrolide antibiotic",
        "peptide": "peptide antibiotic",
        "polymyxin": "peptide antibiotic",
        "cephamycin": "cephalosporin",
        "monobactam": "penicillin beta-lactam",
        "penam": "penicillin beta-lactam",
        "penem": "carbapenem",
        "elfamycin antibiotic": "pulvomycin",
        "glycylcycline": "tetracycline antibiotic",
        "streptogramin b antibiotic": "streptogramin antibiotic",
        "streptogramin a antibiotic": "streptogramin antibiotic",
        "amphenicol": "phenicol antibiotic",
        "beta-lactam": "penicillin beta-lactam",
        "beta lactam": "penicillin beta-lactam",
        "penicillin": "penicillin beta-lactam",
        "tetracycline": "tetracycline antibiotic",
        "fosfomycin": "fosfomycin",
        "phosphonic acid antibiotic": "fosfomycin",
    }
    return mapping.get(dc.lower(), dc.lower())

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
                out.add("macrolide antibiotic")
                out.add("lincosamide antibiotic")
                out.add("streptogramin antibiotic")
                continue
            out.add(normalize_dc(p))
    return out

def truth_to_token(gene: str) -> str:
    g = str(gene).strip().lower()
    g = re.sub(r"[\s_]+", "", g)
    g = re.sub(r"[^a-z0-9\-\.\(\)]+", "", g)
    return g
EXCLUDE_DRUG_CLASSES = {
    canon_dc("disinfecting agents and antiseptics"),
}

# Truth mapping
TRUTH_GENE_TO_DRUGCLASS = {
    "bcii": {"penicillin beta-lactam", "cephalosporin"},
    "fosb": { "fosfomycin"},
    "mcr-4.5": {"peptide antibiotic"},
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
    "pmrf": {"peptide antibiotic"},
    "bes-1": {"penicillin beta-lactam"},
    "uhpt": { "fosfomycin"},
    "kpne": { "peptide antibiotic", "tetracycline antibiotic","aminoglycoside antibiotic", "macrolide antibiotic","cephalosporin","rifamycin antibiotic"},
    "amrb": {"aminoglycoside antibiotic", "macrolide antibiotic"}
}
CARD_FP_GENE_TO_CLASSES = {
    "acra": {"tetracycline antibiotic", "penicillin beta-lactam", "cephalosporin", "phenicol antibiotic", "rifamycin antibiotic", "fluoroquinolone antibiotic"},
    "acrb": {"fluoroquinolone antibiotic", "penicillin beta-lactam", "cephalosporin", "tetracycline antibiotic", "rifamycin antibiotic","phenicol antibiotic"},
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
    "emrd": {"phenicol antibiotic"},
    "kasugamycinresistanceproteinksga": {"aminoglycoside antibiotic"},
    "maca": {"macrolide antibiotic"},
    "macb": {"macrolide antibiotic"},
    "mdfa": {"phenicol antibiotic", "tetracycline antibiotic"},
    "mdta": {"aminocoumarin antibiotic"},
    "mdtb": {"aminocoumarin antibiotic"},
    "mdtc": {"aminocoumarin antibiotic"},
    "mdtg": { "fosfomycin"},
    "mdth": {"fluoroquinolone antibiotic"},
    "mdtk": {"fluoroquinolone antibiotic"},
    "mexd": {"phenicol antibiotic","diaminopyrimidine antibiotic","aminocoumarin antibiotic", "tetracycline antibiotic", "aminoglycoside antibiotic", "penicillin beta-lactam", "cephalosporin","fluoroquinolone antibiotic", "macrolide antibiotic"},
    "mexe": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "mexf": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "mext": {"phenicol antibiotic", "diaminopyrimidine antibiotic", "fluoroquinolone antibiotic"},
    "mexw": {"phenicol antibiotic", "tetracycline antibiotic", "macrolide antibiotic", "fluoroquinolone antibiotic"},
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
    "tolc": {"peptide antibiotic", "tetracycline antibiotic", "penicillin beta-lactam", "cephalosporin", "aminoglycoside antibiotic", "phenicol antibiotic", "rifamycin antibiotic", "aminocoumarin antibiotic", "carbapenem", "fluoroquinolone antibiotic", "macrolide antibiotic"},
    "transcriptionalregulatoryproteincpxrcpxr": {"peptide antibiotic", "sulfonamide antibiotic", "diaminopyrimidine antibiotic", "aminocoumarin antibiotic", "penicillin beta-lactam", "cephalosporin", "aminoglycoside antibiotic", "phenicol antibiotic", "tetracycline antibiotic", "carbapenem", "monobactam", "fluoroquinolone antibiotic", "macrolide antibiotic"},
    "vanr": {"vancomycin"},
    "vans": {"vancomycin"},
    "vany": {"vancomycin"},
    "vatg": {"vancomycin"}
}
def truth_token_to_drugclass(token: str) -> set:
    token = str(token).strip().lower()
    classes = TRUTH_GENE_TO_DRUGCLASS.get(token, set())
    return {normalize_dc(c) for c in classes}

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

# Main

def main():
    ap = argparse.ArgumentParser(description="RUARGpore drug-class metrics basate sul mapping truth")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--ruarg_ham", required=True, nargs="+")
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--out_prefix", default="ruargpore")
    ap.add_argument("--keep_only_gene_presence", action="store_true")
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    outp = args.out_prefix.rstrip("_") + "_"

    #Truth 
    truth = pd.read_csv(args.truth, sep="\t", dtype=str).fillna("")
    truth = truth[truth["present"].astype(str).eq("1")].copy()
    truth["truth_token"] = truth["gene_symbol"].apply(truth_to_token)
    truth_gene_set = set(truth["truth_token"])
    truth_class_set = {
    "fosfomycin",
    "nitroimidazole antibiotic",
    "aminoglycoside antibiotic",
    "carbapenem",
    "cephalosporin",
    "fluoroquinolone antibiotic",
    "macrolide antibiotic",
    "penicillin beta-lactam",
    "peptide antibiotic",
    "pulvomycin",
    "rifamycin antibiotic",
    "tetracycline antibiotic"}


    #RUARG harmonized
    dfs = []
    for f in args.ruarg_ham:
        df = pd.read_csv(f, sep="\t", low_memory=False).fillna("")
        dfs.append(df)
    ru = pd.concat(dfs, ignore_index=True)

    #Gene-level
    ru["pred_gene_token"] = ru["subtype"].apply(truth_to_token)
    ru = ru.drop_duplicates(subset=["pred_gene_token"])
    pred_gene_set = set(ru["pred_gene_token"])

    tp_g, fp_g, fn_g = compute_tp_fp_fn(truth_gene_set, pred_gene_set)
    r_g, p_g, f1_g = metrics_from_counts(len(tp_g), len(fp_g), len(fn_g))

    # Drug-class-level 
    pred_class_set = set()
    for gene in ru["pred_gene_token"]:
        pred_class_set.update(truth_token_to_drugclass(gene))
    pred_class_set = {
        canon_dc(c)
        for c in pred_class_set
        if canon_dc(c) not in EXCLUDE_DRUG_CLASSES
    }
    # Integrate CARD classes for FP genes
    for gene in fp_g:
        if gene in CARD_FP_GENE_TO_CLASSES:
            pred_class_set.update(CARD_FP_GENE_TO_CLASSES[gene])
    tp_c, fp_c, fn_c = compute_tp_fp_fn(truth_class_set, pred_class_set)
    r_c, p_c, f1_c = metrics_from_counts(len(tp_c), len(fp_c), len(fn_c))

    # ---- Write ----
    pd.DataFrame([{
        "pipeline": "RUARGpore",
        "TP": len(tp_g),
        "FP": len(fp_g),
        "FN": len(fn_g),
        "recall": round(r_g,4),
        "precision": round(p_g,4),
        "f1": round(f1_g,4),
        "truth_unique_genes": len(truth_gene_set),
        "pred_unique_genes": len(pred_gene_set),
        "n_rows_ruargpore": len(ru),
    }]).to_csv(f"{args.out_dir}/{outp}gene_metrics.tsv", sep="\t", index=False)

    pd.DataFrame([{
        "pipeline": "RUARGpore",
        "TP": len(tp_c),
        "FP": len(fp_c),
        "FN": len(fn_c),
        "recall": round(r_c,4),
        "precision": round(p_c,4),
        "f1": round(f1_c,4),
        "truth_unique_drug_classes": len(truth_class_set),
        "pred_unique_drug_classes": len(pred_class_set),
        "n_rows_ruargpore": len(ru),
    }]).to_csv(f"{args.out_dir}/{outp}drugclass_metrics.tsv", sep="\t", index=False)

    # lists
    pd.DataFrame({"gene_token": sorted(tp_g)}).to_csv(f"{args.out_dir}/{outp}gene_TP.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_token": sorted(fp_g)}).to_csv(f"{args.out_dir}/{outp}gene_FP.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_token": sorted(fn_g)}).to_csv(f"{args.out_dir}/{outp}gene_FN.tsv", sep="\t", index=False)

    pd.DataFrame({"drug_class": sorted(tp_c)}).to_csv(f"{args.out_dir}/{outp}drugclass_TP.tsv", sep="\t", index=False)
    pd.DataFrame({"drug_class": sorted(fp_c)}).to_csv(f"{args.out_dir}/{outp}drugclass_FP.tsv", sep="\t", index=False)
    pd.DataFrame({"drug_class": sorted(fn_c)}).to_csv(f"{args.out_dir}/{outp}drugclass_FN.tsv", sep="\t", index=False)

    print("\n=== Gene-level metrics ===")
    print(f"TP: {len(tp_g)}, FP: {len(fp_g)}, FN: {len(fn_g)}, recall: {r_g:.4f}, precision: {p_g:.4f}, f1: {f1_g:.4f}")
    print("\n=== Drug-class-level metrics (mapping truth) ===")
    print(f"TP: {len(tp_c)}, FP: {len(fp_c)}, FN: {len(fn_c)}, recall: {r_c:.4f}, precision: {p_c:.4f}, f1: {f1_c:.4f}")
    print("\n=== Gene-level ===")
    print("TP genes:", sorted(tp_g))
    print("FP genes:", sorted(fp_g))
    print("FN genes:", sorted(fn_g))
   
    print("TP drug classes:", sorted(tp_c))
    print("FP drug classes:", sorted(fp_c))
    print("FN drug classes:", sorted(fn_c))

if __name__ == "__main__":
    main()

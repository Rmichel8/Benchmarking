#!/usr/bin/env python3

import csv
import re
import os
import sys

# Ground truth

amr_genes = sorted({
    "mexI",
    "mexB",
    "mexF",
    "mexD",
    "mexQ",
    "ceoB",
    "amrB",
    "carA",
    "OXA-114a",
    "sul2",
    "srmB",
    "mdtF",
    "mdtC",
    "mdtD",
    "mdtB",
    "acrB",
    "acrA",
    "acrD",
    "mfd",
    "PBP1a",
    "tetH",
    "cpxR",
    "CRP",
    "emrD",
    "emrB",
    "H-NS",
    "arnA",
    "bacA",
    "smeE",
    "vgaC",
    "FosA5",
    "SHV-152",
    "SHV-187",
    "SHV-185",
    "SHV-16",
    "SHV-43",
    "oqxB",
    "oleB"
})

amr_genes = list(amr_genes)

print(f"Total unique AMR-related genes: {len(amr_genes)}")
print(amr_genes)
drug_classes = [
    "penicillin beta-lactam",
    "tetracyclines",
    "sulfonamides",
    "fluoroquinolones",
    "fosfomycin",
    "polymyxins",
    "aminoglycosides",
    "macrolides",
    "lincosamides",
    "streptogramins",
    "aminocoumarins",
    "peptide_antibiotics",
    "phenicol",
    "cephalosporin",
    "elfamycin"]

GT_GENES = {
    g.lower().replace("-", "").replace("_", "").replace(" ", "")
    for g in amr_genes
}

GT_CLASSES = {
    c.lower().strip()
    for c in drug_classes
}

# Gene normalization


def normalize_gene(gene):

    gene = gene.lower().strip()
    gene = re.sub(r"^bla", "", gene)
    
    # remove parentesis keep names
    gene = gene.replace("(", "")
    gene = gene.replace(")", "")
    
    gene = gene.replace("_", "")
    gene = gene.replace("-", "")
    gene = gene.replace(" ", "")
    
    gene = re.sub(r"[^a-z0-9]", "", gene)

    return gene

# Drug-class mapping

CLASS_MAP = {

    "amoxicillin": "penicillin beta-lactam",
    "ampicillin": "penicillin beta-lactam",
    "piperacillin": "penicillin beta-lactam",
    "ticarcillin": "penicillin beta-lactam",
    "amoxicillin+clavulanic acid": "penicillin beta-lactam",
    "ampicillin+clavulanic acid": "penicillin beta-lactam",
    "aminocoumarin antibiotic": "aminocoumarins",
    "cephalothin": "cephalosporin",
    "lincosamide antibiotic": "lincosamides",
    "ciprofloxacin": "fluoroquinolones",
    "nalidixic acid": "fluoroquinolones",
    "streptogramin antibiotic": "streptogramins",
    "streptogramin a antibiotic": "streptogramins",
    "chloramphenicol": "phenicol",
    "glycylcycline": "glycylcycline",
    "fosfomycin": "fosfomycin",

    "sulfamethoxazole": "sulfonamides",

    "doxycycline": "tetracyclines",
    "tetracycline": "tetracyclines",
    "ceftriaxone": "cephalosporin",
    "cefepime": "cephalosporin",
    "cefotaxime": "cephalosporin",
    "ceftazidime": "cephalosporin",
    "norfloxacin": "fluoroquinolones",
    "levofloxacin": "fluoroquinolones",
    "ofloxacin": "fluoroquinolones",
    "moxifloxacin": "fluoroquinolones",
    "oxytetracycline": "tetracyclines",
    "minocycline": "tetracyclines",
}

def precision(tp, fp):
    return tp / (tp + fp) if tp + fp else 0

def recall(tp, fn):
    return tp / (tp + fn) if tp + fn else 0

def f1(tp, fp, fn):
    p = precision(tp, fp)
    r = recall(tp, fn)
    return 2 * p * r / (p + r) if p + r else 0


def extract_predictions(tsv):

    genes = set()
    classes = set()

    with open(tsv) as f:

        reader = csv.DictReader(f)

        for row in reader:

            # genes

            gene = normalize_gene(row["Resistance gene"])
            genes.add(gene)

            # phenotypes

            for drug in row["Phenotype"].split(","):

                drug = drug.lower().strip()

                classes.add(CLASS_MAP.get(drug, drug))

    return genes, classes


def evaluate(tsv):

    pred_genes, pred_classes = extract_predictions(tsv)

    TPg = pred_genes & GT_GENES
    FPg = pred_genes - GT_GENES
    FNg = GT_GENES - pred_genes

    TPc = pred_classes & GT_CLASSES
    FPc = pred_classes - GT_CLASSES
    FNc = GT_CLASSES - pred_classes

    print("Gene TP:", sorted(TPg))
    print("Gene FP:", sorted(FPg))
    print("Class FP:", sorted(FPc))
    print("Class TP:", sorted(TPc))
    print("Class TP:")
    for cls in sorted(TPc):
        print(f"  - {cls}")



    out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "summary_metrics.csv"
    )

    with open(out, "w", newline="") as f:

        w = csv.writer(f)

        w.writerow([
            "type",
            "TP",
            "FP",
            "FN",
            "precision",
            "recall",
            "f1"
        ])

        w.writerow([
            "genes",
            len(TPg),
            len(FPg),
            len(FNg),
            round(precision(len(TPg), len(FPg)),3),
            round(recall(len(TPg), len(FNg)),3),
            round(f1(len(TPg), len(FPg), len(FNg)),3)
        ])

        w.writerow([
            "drug_classes",
            len(TPc),
            len(FPc),
            len(FNc),
            round(precision(len(TPc), len(FPc)),3),
            round(recall(len(TPc), len(FNc)),3),
            round(f1(len(TPc), len(FPc), len(FNc)),3)
        ])

    print(f"\nResults written to {out}")


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python compare.py Resfinder_resistance_report.csv")
        sys.exit(1)

    evaluate(sys.argv[1])
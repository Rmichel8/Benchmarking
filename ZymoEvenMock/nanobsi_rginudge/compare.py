#!/usr/bin/env python3

import csv
import re
import os
import sys

# Ground truth


amr_genes = sorted({
    "mexI","mexB","mexF","mexD","mexQ","ceoB","amrB","carA",
    "OXA-114a","sul2","srmB","mdtF","mdtC","mdtD","mdtB",
    "acrB","acrA","acrD","mfd","PBP1a","tetH","cpxR","CRP",
    "emrD","emrB","H-NS","arnA","bacA","smeE","vgaC",
    "FosA5","SHV-152","SHV-187","SHV-185","SHV-16","SHV-43",
    "oqxB","oleB"
})

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
    "elfamycin"
]

GT_GENES = {
    g.lower().replace("-", "").replace("_", "").replace(" ", "")
    for g in amr_genes
}

GT_CLASSES = {c.lower() for c in drug_classes}

# Normalization genes

def normalize_gene(g):

    g = g.lower()

    g = g.replace("bla", "")
    g = g.replace("(", "")
    g = g.replace(")", "")
    
    g = g.replace("-", "")
    g = g.replace("_", "")
    g = g.replace(" ", "")
    
    g = re.sub(r"[^a-z0-9]", "", g)

    return g

# Mapping classes

CLASS_MAP = {

    # beta-lactams
    "amoxicillin": "penicillin beta-lactam",
    "ampicillin": "penicillin beta-lactam",
    "piperacillin": "penicillin beta-lactam",
    "ticarcillin": "penicillin beta-lactam",
    "amoxicillin+clavulanic acid": "penicillin beta-lactam",
    "ampicillin+clavulanic acid": "penicillin beta-lactam",

    # cephalosporins
    "cephalothin": "cephalosporin",
    "ceftriaxone": "cephalosporin",
    "cefepime": "cephalosporin",
    "cefotaxime": "cephalosporin",
    "ceftazidime": "cephalosporin",

    # tetracyclines
    "tetracycline antibiotic": "tetracyclines",
    "tetracycline": "tetracyclines",
    "doxycycline": "tetracyclines",
    "oxytetracycline": "tetracyclines",
    "minocycline": "tetracyclines",
    "glycylcycline": "glycylcycline",

    # sulfonamides
    "sulfonamide antibiotic": "sulfonamides",
    "sulfamethoxazole": "sulfonamides",

    # fluoroquinolones
    "fluoroquinolone antibiotic": "fluoroquinolones",
    "ciprofloxacin": "fluoroquinolones",
    "nalidixic acid": "fluoroquinolones",
    "norfloxacin": "fluoroquinolones",
    "levofloxacin": "fluoroquinolones",
    "ofloxacin": "fluoroquinolones",
    "moxifloxacin": "fluoroquinolones",

    # fosfomycin
    "fosfomycin": "fosfomycin",
    "phosphonic acid antibiotic": "fosfomycin",

    # aminoglycosides
    "aminoglycoside antibiotic": "aminoglycosides",

    # macrolides
    "macrolide antibiotic": "macrolides",

    # lincosamides
    "lincosamide antibiotic": "lincosamides",

    # streptogramins
    "streptogramin antibiotic": "streptogramins",
    "streptogramin a antibiotic": "streptogramins",

    # aminocoumarins
    "aminocoumarin antibiotic": "aminocoumarins",

    # peptide antibiotics
    "peptide antibiotic": "peptide_antibiotics",
    "glycopeptide antibiotic": "peptide_antibiotics",

    # phenicol
    "phenicol antibiotic": "phenicol",
    "chloramphenicol": "phenicol",

    # elfamycin
    "elfamycin antibiotic": "elfamycin",
}


def precision(tp, fp):
    return tp/(tp+fp) if tp+fp else 0

def recall(tp, fn):
    return tp/(tp+fn) if tp+fn else 0

def f1(tp, fp, fn):

    p = precision(tp,fp)
    r = recall(tp,fn)

    return 2*p*r/(p+r) if p+r else 0


def extract_predictions(tsv):

    genes = set()
    classes = set()

    with open(tsv) as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            
            # gene

            gene = normalize_gene(row["Best_Hit_ARO"])

            genes.add(gene)
            
            # class

            drug = row["Drug Class"].lower()

            for c in re.split(r"[;,]", drug):

                c = c.strip()

                classes.add(CLASS_MAP.get(c, c))
            for c in re.split(r"[;,]", row["Drug Class"]):

                c = c.lower().strip()

                mapped = CLASS_MAP.get(c)
                if c in GT_CLASSES:
                    classes.add(c)
                    continue
                if mapped is None:
                    print("UNMAPPED CLASS:", c)
                    classes.add(c)
                else:
                    classes.add(mapped)

    return genes, classes

def evaluate(tsv):

    pred_genes, pred_classes = extract_predictions(tsv)

    TPg = pred_genes & GT_GENES
    FPg = pred_genes - GT_GENES
    FNg = GT_GENES - pred_genes

    TPc = pred_classes & GT_CLASSES
    FPc = pred_classes - GT_CLASSES
    FNc = GT_CLASSES - pred_classes

    print("\n========== GENE TP ==========")
    for g in sorted(TPg):
        print(g)

    print("\n========== GENE FP ==========")
    for g in sorted(FPg):
        print(g)

    print("\n========== CLASS TP ==========")
    for c in sorted(TPc):
        print(c)

    print("\n========== CLASS FP ==========")
    for c in sorted(FPc):
        print(c)

    out = os.path.join(
        os.path.dirname(__file__),
        "summary_metrics.csv"
    )

    with open(out,"w",newline="") as f:

        w = csv.writer(f)

        w.writerow([
            "type","TP","FP","FN",
            "precision","recall","f1"
        ])

        w.writerow([
            "genes",
            len(TPg),
            len(FPg),
            len(FNg),
            round(precision(len(TPg),len(FPg)),3),
            round(recall(len(TPg),len(FNg)),3),
            round(f1(len(TPg),len(FPg),len(FNg)),3)
        ])

        w.writerow([
            "drug_classes",
            len(TPc),
            len(FPc),
            len(FNc),
            round(precision(len(TPc),len(FPc)),3),
            round(recall(len(TPc),len(FNc)),3),
            round(f1(len(TPc),len(FPc),len(FNc)),3)
        ])
    for c in sorted(TPc):
        print(" ", c)   
    print(f"\nResults written to {out}")

if __name__ == "__main__":

    if len(sys.argv)!=2:
        print("Usage: python compare.py barcode01.txt")
        sys.exit(1)

    evaluate(sys.argv[1])
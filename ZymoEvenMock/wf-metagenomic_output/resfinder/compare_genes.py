#!/usr/bin/env python3

import os
import re
import sys
import csv

from ground_truth import amr_genes, drug_classes

GT_GENES = {
    normalize.lower().replace("-", "").replace("_", "").replace(" ", "")
    for normalize in amr_genes
}

GT_CLASSES = {
    c.lower().strip()
    for c in drug_classes
}

OUTDIR = os.path.dirname(os.path.abspath(__file__))

SUMMARY_CSV = os.path.join(
    OUTDIR,
    "summary_gene_compare.csv"
)


def normalize_gene(gene):

    gene = gene.lower().strip()

    gene = re.sub(r"\(.*?\)", "", gene)

    gene = gene.replace("_", "")
    gene = gene.replace("-", "")
    gene = gene.replace(" ", "")

    gene = re.sub(r"[^a-z0-9]", "", gene)
    species = [
        "klebsiellapneumoniae",
        "escherichiacoli",
        "pseudomonasaeruginosa",
        "agrobacteriumfabrum",
        "achromobacterxylosoxidans",
        "acinetobacterbaumannii",
        "enterobactercloacae",
        "staphylococcusaureus",
        "salmonellaenterica",
    ]

    for s in species:
        if gene.startswith(s):
            gene = gene[len(s):]
            break
    return gene

CLASS_MAP = {

    
    # Penicillin beta-lactam
    

    "amoxicillin": "penicillin beta-lactam",
    "amoxicillin+clavulanic acid": "penicillin beta-lactam",
    "ampicillin": "penicillin beta-lactam",
    "ampicillin+clavulanic acid": "penicillin beta-lactam",
    "benzylpenicillin": "penicillin beta-lactam",
    "penicillin": "penicillin beta-lactam",
    "oxacillin": "penicillin beta-lactam",
    "methicillin": "penicillin beta-lactam",
    "nafcillin": "penicillin beta-lactam",
    "dicloxacillin": "penicillin beta-lactam",
    "cloxacillin": "penicillin beta-lactam",
    "piperacillin": "penicillin beta-lactam",
    "piperacillin+tazobactam": "penicillin beta-lactam",
    "ticarcillin": "penicillin beta-lactam",
    "ticarcillin+clavulanic acid": "penicillin beta-lactam",


    # Cephalosporins


    "cefazolin": "cephalosporin",
    "cefalexin": "cephalosporin",
    "cephalexin": "cephalosporin",
    "cefuroxime": "cephalosporin",
    "cefotaxime": "cephalosporin",
    "ceftriaxone": "cephalosporin",
    "ceftazidime": "cephalosporin",
    "cefepime": "cephalosporin",
    "cefoxitin": "cephalosporin",
    "cefpodoxime": "cephalosporin",
    "cefixime": "cephalosporin",
    "cephalothin": "cephalosporin",
    "cephalosporin": "cephalosporin",


    # Tetracyclines

    "tetracycline": "tetracyclines",
    "doxycycline": "tetracyclines",
    "minocycline": "tetracyclines",
    "tigecycline": "tetracyclines",

    # Aminoglycosides

    "amikacin": "aminoglycosides",
    "gentamicin": "aminoglycosides",
    "tobramycin": "aminoglycosides",
    "kanamycin": "aminoglycosides",
    "neomycin": "aminoglycosides",
    "streptomycin": "aminoglycosides",
    "spectinomycin": "aminoglycosides",
    "apramycin": "aminoglycosides",
    "hygromycin": "aminoglycosides",


    # Fluoroquinolones

    "ciprofloxacin": "fluoroquinolones",
    "levofloxacin": "fluoroquinolones",
    "norfloxacin": "fluoroquinolones",
    "ofloxacin": "fluoroquinolones",
    "moxifloxacin": "fluoroquinolones",
    "nalidixic acid": "fluoroquinolones",
    "fluoroquinolone": "fluoroquinolones",

    # Sulfonamides

    "sulfamethoxazole": "sulfonamides",
    "sulfisoxazole": "sulfonamides",
    "sulfadiazine": "sulfonamides",
    "trimethoprim+sulfamethoxazole": "sulfonamides",
    "sulfonamide": "sulfonamides",

    # Macrolides

    "erythromycin": "macrolides",
    "azithromycin": "macrolides",
    "clarithromycin": "macrolides",
    "telithromycin": "macrolides",
    "macrolide": "macrolides",

    # Lincosamides

    "clindamycin": "lincosamides",
    "lincomycin": "lincosamides",
    "lincosamide": "lincosamides",
    # Streptogramins

    "quinupristin": "streptogramins",
    "dalfopristin": "streptogramins",
    "quinupristin/dalfopristin": "streptogramins",
    "streptogramin": "streptogramins",
    "pristinamycin ia": "streptogramins",
    "virginiamycin s": "streptogramins",

    # Aminocoumarins

    "novobiocin": "aminocoumarins",
    "coumermycin": "aminocoumarins",
    "aminocoumarin": "aminocoumarins",

    # Fosfomycin

    "fosfomycin": "fosfomycin",

    # Polymyxins

    "colistin": "polymyxins",
    "polymyxin b": "polymyxins",
    "polymyxin": "polymyxins",

    # Phenicols

    "chloramphenicol": "phenicol",
    "florfenicol": "phenicol",
    "phenicol": "phenicol",

    # Peptide antibiotics

    "vancomycin": "peptide_antibiotics",
    "teicoplanin": "peptide_antibiotics",
    "daptomycin": "peptide_antibiotics",
    "bacitracin": "peptide_antibiotics",
    "peptide": "peptide_antibiotics",

    # Elfamycin

    "kirromycin": "elfamycin",
    "pulvomycin": "elfamycin",
    "elfamycin": "elfamycin",
}
def lookup_gene(gene):

    g = normalize_gene(gene)

    if g in GT_GENES:
        return g

    # substring
    for truth in GT_GENES:

        if g in truth or truth in g:
            return truth

def precision(tp, fp):
    return tp / (tp + fp) if tp + fp else 0


def recall(tp, fn):
    return tp / (tp + fn) if tp + fn else 0


def f1(tp, fp, fn):

    p = precision(tp, fp)
    r = recall(tp, fn)

    return 2 * p * r / (p + r) if p + r else 0
def extract_resfinder(tsv):

    predicted_genes = set()
    predicted_classes = set()

    with open(tsv) as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            # gene

            gene = row["Resistance gene"]

            gene = normalize_gene(gene)
            if gene.startswith("bla"):
                gene = gene[3:]
            predicted_genes.add(gene)

            # drug classes

            resistance = row["Phenotype"].strip()

            if resistance:

                for c in re.split(r"[;,]", resistance):

                    c = c.lower().strip()

                    if c not in CLASS_MAP:
                        print("UNMAPPED PHENOTYPE:", c)

                    predicted_classes.add(CLASS_MAP.get(c, c))

    return predicted_genes, predicted_classes
def evaluate(tsv):

    pred_genes, pred_classes = extract_resfinder(tsv)

    # GENES

    TPg = pred_genes & GT_GENES
    FPg = pred_genes - GT_GENES
    FNg = GT_GENES - pred_genes

    # DRUG CLASSES

    TPc = pred_classes & GT_CLASSES
    FPc = pred_classes - GT_CLASSES
    FNc = GT_CLASSES - pred_classes
    print("\nTrue positives")
    print("----------------")
    for g in sorted(TPg):
        print(" ", g)
    print("\nFalse positives")
    print("----------------")
    for g in sorted(FPg):
        print(" ", g)
    print("\n==============================")
    print(os.path.basename(tsv))
    print("==============================")
    print("\nFalse positives")
    print("----------------")
    for c in sorted(FPc):
        print(" ", c)
    print("\nGENES")
    print("Predicted :", len(pred_genes))
    print("TP :", len(TPg))
    print("FP :", len(FPg))
    print("FN :", len(FNg))
    print("Precision :", round(precision(len(TPg), len(FPg)),3))
    print("Recall    :", round(recall(len(TPg), len(FNg)),3))
    print("F1         :", round(f1(len(TPg), len(FPg), len(FNg)),3))

    print("\nDRUG CLASSES")
    print("Predicted :", len(pred_classes))
    print("TP :", len(TPc))
    print("FP :", len(FPc))
    print("FN :", len(FNc))
    print("Precision :", round(precision(len(TPc), len(FPc)),3))
    print("Recall    :", round(recall(len(TPc), len(FNc)),3))
    print("F1         :", round(f1(len(TPc), len(FPc), len(FNc)),3))
    print(TPc)
if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: compare_genes.py <resfinder.tsv>")
        sys.exit(1)

    evaluate(sys.argv[1])
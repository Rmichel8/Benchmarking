#!/usr/bin/env python3
import csv
import re
import os
import sys
# Ground truth

genes = [
    "patA",
    "emrB",
    "acrA",
    "acrB",
    "marA",
    "emrR",
    "CRP",
    "H-NS",
    "mdtB",
    "mdtC",
    "ramA",
    "baeR",
    "cpxA",
    "cpxR",
    "vgaC",
    "smeB",
    "SHV-187",
    "FosA5",
    "alaS",
    "z",
    "mfd",
    "oqxA",
    "oqxB",
    "lsaa",
    "tetO",
    "ileS",
    "emeA",
    "dfrE",
    "PC1",
    "tet32",
    "norA",
    "tetM",
    "qacA",
    "aac6aph2",
    "ErmC"
]
drug_class= [
    "aminocoumarin",
    "aminoglycoside",
    "penicillin beta-lactam",
    "cephalosporin",
    "carbapenem",
    "fluoroquinolone",
    "tetracycline",
    "glycylcycline",
    "phenicol",
    "macrolide",
    "lincosamide",
    "streptogramin",
    "pleuromutilin",
    "rifamycin",
    "fosfomycin",
    "diaminopyrimidine",
    "mupirocin",
    "disinfecting agents and antiseptics",
    "nitrofuran"
]


GT_GENES = {
    g.lower()
    .replace("-", "")
    .replace("_", "")
    .replace(" ", "")
    for g in genes
}


GT_CLASSES = {
    c.lower()
    for c in drug_class
}
def normalize_gene(gene):

    gene = gene.lower().strip()

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
# Penicillins

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
"unknown beta-lactam": "penicillin beta-lactam",
# Monobactams
"aztreonam": "monobactam",
# Aminoglycosides

"netilmicin": "aminoglycoside",
"dibekacin": "aminoglycoside",
"isepamicin": "aminoglycoside",
"fortimicin": "aminoglycoside",
# Fosfomycin
"fosfomycin": "fosfomycin",
# Streptogramins

"pristinamycin ia": "streptogramin",
"pristinamycin iia": "streptogramin",
"virginiamycin s": "streptogramin",
"virginiamycin m": "streptogramin",
# Macrolides

"erythromycin": "macrolide",
"azithromycin": "macrolide",
"telithromycin": "macrolide",
# Lincosamides


"lincomycin": "lincosamide",
"clindamycin": "lincosamide",
# Fluoroquinolones / Quinolones

"ciprofloxacin": "fluoroquinolone",
"nalidixic acid": "fluoroquinolone",
# Diaminopyrimidines

"trimethoprim": "diaminopyrimidine",

# Phenicols

"chloramphenicol": "phenicol",
# Tetracyclines

"doxycycline": "tetracycline",
"minocycline": "tetracycline",
"tetracycline": "tetracycline"}

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
                        continue

                    predicted_classes.add(CLASS_MAP[c])

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
    return {
        "sample": os.path.basename(tsv),

        "predicted_genes": len(pred_genes),
        "tp_genes": len(TPg),
        "fp_genes": len(FPg),
        "fn_genes": len(FNg),
        "precision_genes": precision(len(TPg), len(FPg)),
        "recall_genes": recall(len(TPg), len(FNg)),
        "f1_genes": f1(len(TPg), len(FPg), len(FNg)),

        "predicted_classes": len(pred_classes),
        "tp_classes": len(TPc),
        "fp_classes": len(FPc),
        "fn_classes": len(FNc),
        "precision_classes": precision(len(TPc), len(FPc)),
        "recall_classes": recall(len(TPc), len(FNc)),
        "f1_classes": f1(len(TPc), len(FPc), len(FNc)),
    }
if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: compare_genes.py <resfinder.tsv>")
        sys.exit(1)

    results = evaluate(sys.argv[1])

    outfile = "metrics.tsv"

    with open(outfile, "w", newline="") as f:

        writer = csv.DictWriter(f, fieldnames=results.keys(), delimiter="\t")

        writer.writeheader()
        writer.writerow(results)

    print(f"\nMetrics saved to {outfile}")
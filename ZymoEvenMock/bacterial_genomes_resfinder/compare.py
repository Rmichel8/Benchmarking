#!/usr/bin/env python3

import csv
import re
import os
import sys


############################################################
# Ground truth
############################################################

amr_genes = {
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
}


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


############################################################
# Normalizzazione ground truth
############################################################
def normalize_gene(gene):

    gene = gene.lower().strip()

    # Rimuove il prefisso bla
    gene = re.sub(r"^bla", "", gene)

    # Rimuove solo parentesi e non il loro contenuto
    gene = gene.replace("(", "")
    gene = gene.replace(")", "")

    # Rimuove _, -, spazi e altri caratteri speciali
    gene = re.sub(r"[^a-z0-9]", "", gene)

    return gene


GT_GENES = {
    normalize_gene(gene)
    for gene in amr_genes
}


# CORREZIONE: costruzione corretta del set delle classi
GT_CLASSES = {
    c.lower().strip()
    for c in drug_classes
}


############################################################
# Drug-class mapping
############################################################

CLASS_MAP = {

    # Beta-lactams / penicillins
    "amoxicillin": "penicillin beta-lactam",
    "ampicillin": "penicillin beta-lactam",
    "piperacillin": "penicillin beta-lactam",
    "ticarcillin": "penicillin beta-lactam",
    "amoxicillin+clavulanic acid": "penicillin beta-lactam",
    "ampicillin+clavulanic acid": "penicillin beta-lactam",

    # Aminocoumarins
    "aminocoumarin antibiotic": "aminocoumarins",

    # Cephalosporins
    "cephalothin": "cephalosporin",
    "ceftriaxone": "cephalosporin",
    "cefepime": "cephalosporin",
    "cefotaxime": "cephalosporin",
    "ceftazidime": "cephalosporin",

    # Lincosamides
    "lincosamide antibiotic": "lincosamides",

    # Fluoroquinolones
    "ciprofloxacin": "fluoroquinolones",
    "nalidixic acid": "fluoroquinolones",
    "norfloxacin": "fluoroquinolones",
    "levofloxacin": "fluoroquinolones",
    "ofloxacin": "fluoroquinolones",
    "moxifloxacin": "fluoroquinolones",

    # Streptogramins
    "streptogramin antibiotic": "streptogramins",
    "streptogramin a antibiotic": "streptogramins",

    # Phenicols
    "chloramphenicol": "phenicol",

    # Tetracyclines
    "doxycycline": "tetracyclines",
    "tetracycline": "tetracyclines",
    "oxytetracycline": "tetracyclines",
    "minocycline": "tetracyclines",

    # Sulfonamides
    "sulfamethoxazole": "sulfonamides",

    # Fosfomycin
    "fosfomycin": "fosfomycin",

    # Other classes
    "trimethoprim": "trimethoprim",
}


############################################################
# Metriche
############################################################

def precision(tp, fp):

    return tp / (tp + fp) if tp + fp else 0


def recall(tp, fn):

    return tp / (tp + fn) if tp + fn else 0


def f1(tp, fp, fn):

    p = precision(tp, fp)
    r = recall(tp, fn)

    return 2 * p * r / (p + r) if p + r else 0


############################################################
# Estrazione predizioni ResFinder
############################################################

def extract_predictions(tsv):

    genes = set()
    classes = set()

    with open(tsv, newline="", encoding="utf-8") as f:

        # Il file ResFinder è TAB-delimited
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            ################################################
            # GENE
            ################################################

            raw_gene = row["Resistance gene"]

            gene = normalize_gene(raw_gene)

            genes.add(gene)


            ################################################
            # PHENOTYPE
            ################################################

            phenotype = row["Phenotype"].strip()

            # Ignora fenotipi non specifici
            if phenotype.lower() in {
                "unknown beta-lactam",
                "unknown"
            }:
                continue


            # ResFinder separa i fenotipi con virgole
            for drug in phenotype.split(","):

                drug = drug.lower().strip()

                if not drug:
                    continue

                drug_class = CLASS_MAP.get(drug, drug)

                classes.add(drug_class)


    return genes, classes


############################################################
# Valutazione
############################################################

def evaluate(tsv):

    pred_genes, pred_classes = extract_predictions(tsv)


    ########################################################
    # GENI
    ########################################################

    TPg = pred_genes & GT_GENES
    FPg = pred_genes - GT_GENES
    FNg = GT_GENES - pred_genes


    ########################################################
    # CLASSI
    ########################################################

    TPc = pred_classes & GT_CLASSES
    FPc = pred_classes - GT_CLASSES
    FNc = GT_CLASSES - pred_classes


    ########################################################
    # OUTPUT DETTAGLIATO
    ########################################################

    print("\n==============================================")
    print("GENE-LEVEL RESULTS")
    print("==============================================")

    print("\nGround truth genes:", len(GT_GENES))
    print("Detected genes:", len(pred_genes))

    print("\nTP genes:")
    print(sorted(TPg))

    print("\nFP genes:")
    print(sorted(FPg))

    print("\nFN genes:")
    print(sorted(FNg))


    print("\n==============================================")
    print("DRUG-CLASS RESULTS")
    print("==============================================")

    print("\nGround truth classes:", len(GT_CLASSES))
    print("Detected classes:", len(pred_classes))

    print("\nTP classes:")
    print(sorted(TPc))

    print("\nFP classes:")
    print(sorted(FPc))

    print("\nFN classes:")
    print(sorted(FNc))


    ########################################################
    # CSV OUTPUT
    ########################################################

    out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "summary_metrics.csv"
    )


    with open(out, "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow([
            "type",
            "TP",
            "FP",
            "FN",
            "precision",
            "recall",
            "f1"
        ])


        writer.writerow([
            "genes",
            len(TPg),
            len(FPg),
            len(FNg),
            round(precision(len(TPg), len(FPg)), 3),
            round(recall(len(TPg), len(FNg)), 3),
            round(f1(len(TPg), len(FPg), len(FNg)), 3)
        ])


        writer.writerow([
            "drug_classes",
            len(TPc),
            len(FPc),
            len(FNc),
            round(precision(len(TPc), len(FPc)), 3),
            round(recall(len(TPc), len(FNc)), 3),
            round(f1(len(TPc), len(FPc), len(FNc)), 3)
        ])


    print(f"\nResults written to: {out}")


############################################################
# MAIN
############################################################

if __name__ == "__main__":

    if len(sys.argv) != 2:

        print(
            "Usage: python compare.py "
            "ResFinder_results_tab.txt"
        )

        sys.exit(1)


    evaluate(sys.argv[1])
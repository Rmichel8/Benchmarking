#!/usr/bin/env python3
import os
import sys
import csv
import json
import re

B1_TP = ["Amikacin" ,"Amoxicillin/clavulanic acid", "Tetracycline"]
B1_TN = ["Aztreonam", "Ceftazidime", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin", "Ceftriaxone",
         "Cefuroxime", "Colistin", "Ertapenem", "Imipenem",
         "Meropenem", "Minocycline", "Nitrofurantoin", "Norfloxacin",
         "Tegacyclin", "BTobramycin"]


B2_TP = []
B2_TN = ["Oxacillin", "Penicillin", "Chloramphenicol", "Ciprofloxacin",
         "Fusidic acid", "Mupirocin", "Clindamycin", "Daptomycin",
         "Erythromycin", "Gentamicin", "Linezolid", "Minocycline",
         "Tetracycline", "Nitrofurantoin", "Rifampicin", "Tegacyclin",
         "Teicoplanin", "Vancomycin"]


B3_TP = []
B3_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
         "Cefoperazone/sulbactam", "Aztreonam", "Ceftazidime",
         "Cefatriaxone", "Cefepime", "Chloramphenicol",
         "Ciprofloxacin", "Levofloxacin", "Colistin", "Imipenem"]


B4_TP = ["Oxacillin", "Penicillin", "Ciprofloxacin",
         "Fusidic acid", "Gentamicin"]
B4_TN = ["Chloramphenicol", "Mupirocin", "Clindamycin", "Daptomycin",
         "Erythromycin", "Linezolid", "Minocycline", "Tetracycline",
         "Nitrofurantoin", "Rifampicin", "Tegacyclin", "Teicoplanin",
         "Vancomycin"]


B5_TP = ["Ceftazidime", "Cefepime", "Chloramphenicol"]
B5_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
         "Cefoperazone/sulbactam", "Aztreonam", "Cefatriaxone",
         "Ciprofloxacin", "Levofloxacin", "Colistin", "Imipenem"]


B6_TP = []
B6_TN = [
    "Amikacin",
    "Gentamicin",
    "Ampicillin/sulbactam",
    "Cefoperazone/sulbactam",
    "Aztreonam",
    "Ceftazidime",
    "Cefatriaxone",
    "Cefepime",
    "Chloramphenicol",
    "Ciprofloxacin",
    "Levofloxacin",
    "Colistin",
    "Imipenem",
    "Oxacillin",
    "Penicillin",
    "Fusidic acid",
    "Mupirocin",
    "Clindamycin",
    "Daptomycin",
    "Erythromycin",
    "Linezolid",
    "Minocycline",
    "Tetracycline",
    "Nitrofurantoin",
    "Rifampicin",
    "Tegacyclin",
    "Teicoplanin",
    "Vancomycin"
]


B7_TP = []
B7_TN = ["Oxacillin", "Penicillin", "Chloramphenicol", "Ciprofloxacin",
         "Fusidic acid", "Mupirocin", "Clindamycin", "Daptomycin",
         "Erythromycin", "Gentamicin", "Linezolid", "Minocycline",
         "Tetracycline", "Nitrofurantoin", "Rifampicin",
         "Tegacyclin", "Teicoplanin", "Vancomycin"]


B8_TP = ["Amoxicillin/clavulanic acid", "Tetracycline", "Tobramycin", "Ceftazidime"] #"Ceftazidime" Intermedio
B8_TN = ["Amikacin", "Aztreonam", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin",
         "Ceftriaxone", "Cefuroxime", "Colistin", "Ertapenem",
         "Imipenem", "Meropenem", "Minocycline",
         "Nitrofurantoin", "Norfloxacin", "Tegacyclin"]


B9_TP = ["Amoxicillin/clavulanic acid", "Tetracycline"]
B9_TN = ["Amikacin", "Aztreonam", "Ceftazidime", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin",
         "Ceftriaxone", "Cefuroxime", "Colistin", "Ertapenem",
         "Imipenem", "Meropenem", "Minocycline",
         "Nitrofurantoin", "Norfloxacin", "Tegacyclin", "Tobramycin"]


B10_TP = []
B10_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
          "Cefoperazone/sulbactam", "Aztreonam", "Ceftazidime",
          "Cefatriaxone", "Cefepime", "Chloramphenicol",
          "Ciprofloxacin", "Levofloxacin", "Colistin", "Imipenem",
          "Amoxicillin/clavulanic acid", "Tetracycline",
          "Cefuroxime", "Ertapenem", "Meropenem", "Minocycline",
          "Nitrofurantoin", "Norfloxacin", "Tegacyclin", "Tobramycin"]


B11_TP = ["Ceftazidime", "Cefepime", "Clindamycin", "Erythromycin"]
B11_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
          "Cefoperazone/sulbactam", "Aztreonam", "Cefatriaxone",
          "Chloramphenicol", "Ciprofloxacin", "Levofloxacin",
          "Colistin", "Imipenem", "Oxacillin", "Penicillin",
          "Fusidic acid", "Mupirocin", "Daptomycin", "Linezolid",
          "Minocycline", "Tetracycline", "Nitrofurantoin",
          "Rifampicin", "Tegacyclin", "Teicoplanin", "Vancomycin"]

GT = {
    "B1": (set(B1_TP), set(B1_TN)),
    "B2": (set(B2_TP), set(B2_TN)),
    "B3": (set(B3_TP), set(B3_TN)),
    "B4": (set(B4_TP), set(B4_TN)),
    "B5": (set(B5_TP), set(B5_TN)),
    "B6": (set(B6_TP), set(B6_TN)),
    "B7": (set(B7_TP), set(B7_TN)),
    "B8": (set(B8_TP), set(B8_TN)),
    "B9": (set(B9_TP), set(B9_TN)),
    "B10": (set(B10_TP), set(B10_TN)),
    "B11": (set(B11_TP), set(B11_TN)),
}



CARD_JSON = "/nfsd/sysbiobig/Riccardom/database_CARD/card-data/card.json"
RESULTS_DIR = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "RUARGpore90/"
    "compare_results"
)

SUMMARY_CSV = os.path.join(
    RESULTS_DIR,
    "summary_CARD_compare.csv"
)
from difflib import get_close_matches

CARD_TO_GT = {
    # beta-lattamici
    "amoxicillin":"Amoxicillin/clavulanic acid",
    "ampicillin":"Amoxicillin/clavulanic acid",
    "benzylpenicillin":"Amoxicillin/clavulanic acid",
    "cefepime":"Cefepime",
    "ceftriaxone":"Ceftriaxone",
    "cefuroxime":"Cefuroxime",
    "ceftazidime":"Ceftazidime",
    "aztreonam":"Aztreonam",
    "imipenem":"Imipenem",
    "meropenem":"Meropenem",
    "ertapenem":"Ertapenem",
    # quinoloni
    "ciprofloxacin":
        "Ciprofloxacin",

    "levofloxacin":
        "Levofloxacin",

    "norfloxacin":
        "Norfloxacin",

    # aminoglicosidi

    "amikacin":
        "Amikacin",

    "gentamicin":
        "Gentamicin",

    "tobramycin":
        "Tobramycin",

    # tetracicline

    "tetracycline":
        "Tetracycline",

    "minocycline":
        "Minocycline",

    "tigecycline":
        "Tigecycline",
    # altri

    "chloramphenicol":
        "Chloramphenicol",

    "rifampin":
        "Rifampicin",

    "rifampicin":
        "Rifampicin",

    "vancomycin":
        "Vancomycin",

    "teicoplanin":
        "Teicoplanin",

    "colistin A":
        "Colistin",

    "colistin B":
        "Colistin",

    "polymyxin B":
        "Colistin",

    "nitrofurantoin":
        "Nitrofurantoin",

    "fosfomycin":
        "Fosfomycin",

    "linezolid":
        "Linezolid",

    "daptomycin":
        "Daptomycin"
}

# ============================================================
# COSTRUZIONE INDICI CARD
# ============================================================

with open(CARD_JSON) as f:
    CARD = json.load(f)


GENE_TO_ABX = {}
ARO_TO_ABX = {}


def normalize_gene(g):

    g = g.lower()

    g = re.sub(r"[^a-z0-9]", "", g)

    return g


for obj in CARD.values():

    if "ARO_name" not in obj:
        continue

    aro = obj.get("ARO_accession")

    gene = obj["ARO_name"]

    antibiotics = set()

    if "ARO_category" not in obj:
        continue

    for cat in obj["ARO_category"].values():

        if cat["category_aro_class_name"] != "Antibiotic":
            continue

        ab = cat["category_aro_name"].lower()

        if ab in CARD_TO_GT:

            antibiotics.add(CARD_TO_GT[ab])

    if not antibiotics:
        continue

    ##########################################################
    # costruzione di tutte le forme del gene
    ##########################################################

    gene_norm = normalize_gene(gene)

    variants = set()

    variants.add(gene_norm)

    # senza numeri
    variants.add(
        re.sub(r"[0-9]", "", gene_norm)
    )

    # solo lettere iniziali
    m = re.match(r"[a-z]+", gene_norm)

    if m:
        variants.add(m.group())

    # rimuove eventuale "bla"
    if gene_norm.startswith("bla"):

        variants.add(
            gene_norm[3:]
        )

    ##########################################################

    for g in variants:

        if not g:
            continue

        GENE_TO_ABX.setdefault(g, set()).update(antibiotics)

    if aro:

        ARO_TO_ABX[aro] = antibiotics


print()

print("CARD entries :", len(CARD))

print("Indexed genes:", len(GENE_TO_ABX))

print("Indexed AROs :", len(ARO_TO_ABX))

def lookup_gene(gene):

    gene = normalize_gene(gene)

    if gene in GENE_TO_ABX:
        return GENE_TO_ABX[gene]

    matches = get_close_matches(
        gene,
        GENE_TO_ABX.keys(),
        n=1,
        cutoff=0.90
    )

    if matches:
        return GENE_TO_ABX[matches[0]]

    return set()
def get_benchmark(path):

    for part in os.path.normpath(path).split(os.sep):

        m = re.search(r"_B(\d+)", part)

        if m:
            return "B" + m.group(1)

    return None

def extract_ruargpore(tsv):

    predicted = set()
    found_genes = set()
    unmapped_genes = set()

    with open(tsv) as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            gene = row["subtype"]

            gene = normalize_gene(gene)

            found_genes.add(gene)

            antibiotics = lookup_gene(gene)

            if antibiotics:
                predicted.update(antibiotics)
            else:
                unmapped_genes.add(gene)
        print("\nMapped genes:", len(found_genes)-len(unmapped_genes))
        print("Unmapped genes:", len(unmapped_genes))
        if unmapped_genes:
            print("\nUnmapped:")
            for g in sorted(unmapped_genes):
                print(g)

    return predicted, found_genes, unmapped_genes

def evaluate(tsv):

    bench = get_benchmark(tsv)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    if bench not in GT:
        print("Cannot identify benchmark:", tsv)
        return

    tp_truth, tn_truth = GT[bench]

    predicted, genes, unmapped = extract_ruargpore(tsv)
    tool = "RUARGPORE"

    ############################################################
    # confronto con la ground truth
    ############################################################

    TP = predicted & tp_truth
    FP = predicted & tn_truth
    FN = tp_truth - predicted
    TN = tn_truth - predicted

    ############################################################
    # report a video
    ############################################################

    sample = os.path.basename(tsv)

    print()
    print("=" * 70)
    print(sample)
    print("=" * 70)

    print("Benchmark :", bench)
    print("Tool      :", tool)

    print()

    print("Predicted antibiotics")
    print("---------------------")
    for a in sorted(predicted):
        print("  ", a)

    print()

    print("True positives")
    print("---------------------")
    for a in sorted(TP):
        print("  ", a)

    print()

    print("False positives")
    print("---------------------")
    for a in sorted(FP):
        print("  ", a)

    print()

    print("False negatives")
    print("---------------------")
    for a in sorted(FN):
        print("  ", a)

    print()

    print("True negatives")
    print("---------------------")
    for a in sorted(TN):
        print("  ", a)

    ############################################################
    # geni trovati
    ############################################################
    gene_file = os.path.join(
    RESULTS_DIR,
    sample + "_predicted_genes.txt"
)
    with open(gene_file, "w") as out:
        for g in sorted(genes):
            out.write(g + "\n")

    ############################################################
    # summary csv
    ############################################################

    os.makedirs(
        os.path.dirname(SUMMARY_CSV),
        exist_ok=True
    )

    write_header = not os.path.exists(SUMMARY_CSV)

    with open(
        SUMMARY_CSV,
        "a",
        newline=""
    ) as out:

        writer = csv.writer(out)

        if write_header:
            writer.writerow([
                "sample",
                "benchmark",
                "tool",
                "TP",
                "FP",
                "TN",
                "FN",
                "predicted_antibiotics",
                "tp_antibiotics",
                "fp_antibiotics",
                "fn_antibiotics",
                "genes_found",
                "genes_unmapped"
            ])

        writer.writerow([
            sample,
            bench,
            tool,
            len(TP),
            len(FP),
            len(TN),
            len(FN),
            ";".join(sorted(predicted)),
            ";".join(sorted(TP)),
            ";".join(sorted(FP)),
            ";".join(sorted(FN)),
            ";".join(sorted(genes)),
            ";".join(sorted(unmapped))
        ])

    print()
    print("Summary updated:")
    print(SUMMARY_CSV)

    print("Mapped genes   :", len(genes) - len(unmapped))
    print("Unmapped genes :", len(unmapped))
    
if __name__=="__main__":


    if len(sys.argv)!=2:

        print(
        "Usage: compare.py <tsv>"
        )

        sys.exit(1)

    evaluate(sys.argv[1])
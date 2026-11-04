#!/usr/bin/env python3

import os
import sys
import csv
import json
from collections import defaultdict
import re
B1_TP = ["Amikacin" ,"Amoxicillin/clavulanic acid", "Tetracycline"]
B1_TN = ["Aztreonam", "Ceftazidime", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin", "Ceftriaxone",
         "Cefuroxime", "Colistin", "Ertapenem", "Imipenem",
         "Meropenem", "Minocycline", "Nitrofurantoin", "Norfloxacin",
         "Tegacyclin", "Tobramycin"]


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


RGI_FILE = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_rgi/"
    "rgi_results/"
    "RUARGPOREPRJNA821338_rgi.txt"
)

INPUT_DIR = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_rgi/"
    "amr_fastq_assignment/"
    "split_by_fastq_annotated"
)

OUTDIR = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_rgi"
)

OUTFILE = os.path.join(
    OUTDIR,
    "per_sample_results.tsv"
)

SUMMARY = os.path.join(
    OUTDIR,
    "summary_metrics.csv"
)
CARD_JSON = "/nfsd/sysbiobig/Riccardom/database_CARD/card-data/card.json"

# NORMALIZZAZIONE CARD -> Ground truth



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
        "Tegacycline",


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
        "Daptomycin",

}


with open(CARD_JSON) as f:
    CARD = json.load(f)


GENE_TO_ABX = {}
ARO_TO_ABX = {}
CARD_GENES = []


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

    gene_norm = normalize_gene(gene)

    variants = set()

    variants.add(gene_norm)

    variants.add(
        re.sub(r"[0-9]", "", gene_norm)
    )


    m = re.match(r"[a-z]+", gene_norm)

    if m:
        variants.add(m.group())


    if gene_norm.startswith("bla"):

        variants.add(
            gene_norm[3:]
        )


    for g in variants:

        if not g:
            continue

        GENE_TO_ABX[g] = antibiotics

        CARD_GENES.append(
            (g, antibiotics)
        )

    if aro:

        ARO_TO_ABX[aro] = antibiotics


print()

print("CARD entries :", len(CARD))

print("Indexed genes:", len(GENE_TO_ABX))

print("Indexed AROs :", len(ARO_TO_ABX))


def get_benchmark(path):

    m = re.search(r"_B(\d+)", path)

    if m:
        return "B" + m.group(1)

    return None
def predicted_antibiotics(tsv):

    predicted = set()

    with open(tsv) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            gene = normalize_gene(
                row["GENE"]
            )


            if gene in GENE_TO_ABX:

                predicted.update(
                    GENE_TO_ABX[gene]
                )

            else:

                for card_gene, abx in CARD_GENES:

                    if gene.startswith(card_gene) or card_gene.startswith(gene):

                        predicted.update(abx)


    return predicted
def evaluate(sample, predicted):

    gt_tp, gt_tn = GT[sample]

    tp = len(predicted & gt_tp)

    fp = len(predicted & gt_tn)

    fn = len(gt_tp - predicted)

    tn = len(gt_tn - predicted)

    return tp, fp, tn, fn
results = []

TOTAL_TP = 0
TOTAL_FP = 0
TOTAL_TN = 0
TOTAL_FN = 0

for file in sorted(os.listdir(INPUT_DIR)):

    if not file.endswith(".tsv"):
        continue

    sample = get_benchmark(file)

    if sample is None:
        continue

    preds = predicted_antibiotics(
        os.path.join(INPUT_DIR, file)
    )

    tp, fp, tn, fn = evaluate(sample, preds)

    TOTAL_TP += tp
    TOTAL_FP += fp
    TOTAL_TN += tn
    TOTAL_FN += fn

    results.append([
        sample,
        tp,
        fp,
        tn,
        fn,
        ";".join(sorted(preds))
    ])
with open(OUTFILE, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "sample",
        "TP",
        "FP",
        "TN",
        "FN",
        "predicted_antibiotics"
    ])

    writer.writerows(results)
precision = (
    TOTAL_TP /
    (TOTAL_TP + TOTAL_FP)
    if TOTAL_TP + TOTAL_FP else 0
)

recall = (
    TOTAL_TP /
    (TOTAL_TP + TOTAL_FN)
    if TOTAL_TP + TOTAL_FN else 0
)

specificity = (
    TOTAL_TN /
    (TOTAL_TN + TOTAL_FP)
    if TOTAL_TN + TOTAL_FP else 0
)

accuracy = (
    (TOTAL_TP + TOTAL_TN) /
    (TOTAL_TP + TOTAL_TN + TOTAL_FP + TOTAL_FN)
    if TOTAL_TP + TOTAL_TN + TOTAL_FP + TOTAL_FN else 0
)

f1 = (
    2 * precision * recall /
    (precision + recall)
    if precision + recall else 0
)
denominator = (
    (TOTAL_TP + TOTAL_FP) *
    (TOTAL_TP + TOTAL_FN) *
    (TOTAL_TN + TOTAL_FP) *
    (TOTAL_TN + TOTAL_FN)
)

if denominator:

    mcc = (
        (TOTAL_TP * TOTAL_TN - TOTAL_FP * TOTAL_FN)
        /
        (denominator ** 0.5)
    )

else:

    mcc = 0
with open(SUMMARY, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "TP",
        "FP",
        "TN",
        "FN",
        "precision",
        "recall",
        "f1",
        "specificity",
        "accuracy",
        "MCC"
    ])

    writer.writerow([
        TOTAL_TP,
        TOTAL_FP,
        TOTAL_TN,
        TOTAL_FN,
        round(precision,4),
        round(recall,4),
        round(f1,4),
        round(specificity,4),
        round(accuracy,4),
        round(mcc,4)
    ])
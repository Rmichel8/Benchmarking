#!/usr/bin/env python3

import csv
import sys

BASE_DIR = "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018"

NANOBSI_FILE = (
    BASE_DIR
    + "/benchmarking/bacterial_genomes_rgi/"
      "barcode01_rgi_results3/barcode01.txt"
)

sys.path.append(BASE_DIR)
from ground_truth_AST import antibiotics

GT = set(antibiotics)


# Broad mapping from RGI drug classes to the AST antibiotics
DRUG_CLASS_TO_ABX = {
    "aminoglycoside": {"Gentamicin"},
    "aminoglycoside antibiotic": {"Gentamicin"},

    "beta-lactam": {
        "Benzylpenicillin",
        "Amoxicillin",
        "Meropenem"
    },
    "penicillin beta-lactam": {
        "Benzylpenicillin",
        "Amoxicillin",
        "Meropenem"
    },

    "carbapenem": {"Meropenem"},

    "glycopeptide": {"Vancomycin"},
    "glycopeptide antibiotic": {"Vancomycin"},
    "vancomycin": {"Vancomycin"},

    "mupirocin": {"Mupirocin"},

    "nitroimidazole": {"Metronidazole"},
}


def predict_antibiotics(row):
    predicted = set()

    # ---------------------------------------------------------
    # 1. Direct matching using the "Antibiotic" column
    # ---------------------------------------------------------
    antibiotic_field = row.get("Antibiotic", "")

    if antibiotic_field:
        antibiotic_text = antibiotic_field.lower()

        for abx in GT:
            if abx.lower() in antibiotic_text:
                predicted.add(abx)

    # ---------------------------------------------------------
    # 2. Broad mapping using the "Drug Class" column
    # ---------------------------------------------------------
    drug_class_field = row.get("Drug Class", "")

    if drug_class_field:
        classes = [
            x.strip().lower()
            for x in drug_class_field.split(";")
            if x.strip()
        ]

        for drug_class in classes:
            if drug_class in DRUG_CLASS_TO_ABX:
                predicted.update(
                    DRUG_CLASS_TO_ABX[drug_class]
                )

    # Keep only antibiotics belonging to the AST ground truth
    return predicted & GT


# -------------------------------------------------------------
# Process NanoBSI/RGI results
# -------------------------------------------------------------

predicted_antibiotics = set()
num_hits = 0

with open(NANOBSI_FILE, "r", newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        num_hits += 1

        predicted_antibiotics.update(
            predict_antibiotics(row)
        )


# -------------------------------------------------------------
# Calculate TP, FN and recall
# -------------------------------------------------------------

TP = predicted_antibiotics & GT
FN = GT - predicted_antibiotics

recall = len(TP) / len(GT) if GT else 0.0


# -------------------------------------------------------------
# Print results
# -------------------------------------------------------------


print(f"\nNumber of RGI hits: {num_hits}")

print("\nGround truth antibiotics:")
for abx in sorted(GT):
    print(f"  - {abx}")

print("\nPredicted antibiotics:")
for abx in sorted(predicted_antibiotics):
    print(f"  - {abx}")

print("\nTrue positives:")
for abx in sorted(TP):
    print(f"  - {abx}")

print("\nFalse negatives:")
for abx in sorted(FN):
    print(f"  - {abx}")

print("\nResults:")
print(f"TP = {len(TP)}")
print(f"FN = {len(FN)}")
print(f"Recall = {recall:.4f}")
print(f"Recall = {recall * 100:.2f}%")

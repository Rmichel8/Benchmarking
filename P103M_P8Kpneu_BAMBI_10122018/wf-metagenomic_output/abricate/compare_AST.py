#!/usr/bin/env python3

import sys
from pathlib import Path
import csv


# ============================================================
# Ground truth
# ============================================================

GT_DIR = Path(
    "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_19122018"
)

sys.path.insert(0, str(GT_DIR))

from ground_truth_AST import antibiotics


# ============================================================
# Input / output
# ============================================================

INPUT_FILE = Path(
    "/nfsd/sysbiobig/Riccardom/"
    "P103M_P8Kpneu_BAMBI_10122018/"
    "benchmarking/wf-metagenomic_output/abricate/"
    "P103M_P8Kpneu_all_10122018/amr/tsv/"
    "P103M_P8Kpneu_all_10122018_amr_results.tsv"
)

OUTPUT_FILE = INPUT_FILE.parent / (
    "P103M_P8Kpneu_all_10122018_amr_AST_comparison.txt"
)


# ============================================================
# Ground truth
# ============================================================

GT = set(antibiotics)


# ============================================================
# Resistance class -> antibiotics
# ============================================================

RESISTANCE_TO_ABX = {

    # Aminoglycosides
    "aminoglycoside": {
        "Gentamicin"
    },

    # Beta-lactams
    "carbapenem": {
        "Meropenem"
    },

    "penam": {
        "Benzylpenicillin",
        "Amoxicillin"
    },

    "penicillin": {
        "Benzylpenicillin",
        "Amoxicillin"
    },

    "penem": {
        "Meropenem"
    },

    # Glycopeptides
    "glycopeptide": {
        "Vancomycin"
    },

    "vancomycin": {
        "Vancomycin"
    },

    # Mupirocin
    "mupirocin": {
        "Mupirocin"
    },

    # Nitroimidazoles
    "nitroimidazole": {
        "Metronidazole"
    },

    "nitroimidazole derivative": {
        "Metronidazole"
    },

    # Direct antibiotic names
    "gentamicin": {
        "Gentamicin"
    },

    "benzylpenicillin": {
        "Benzylpenicillin"
    },

    "amoxicillin": {
        "Amoxicillin"
    },

    "metronidazole": {
        "Metronidazole"
    },

    "meropenem": {
        "Meropenem"
    },
}


# ============================================================
# Predict antibiotics
# ============================================================

def predict_antibiotics(resistance, product):
    """
    Predict AST antibiotics from the CARD resistance classes
    and product description.

    If a resistance class contains an antibiotic represented
    in the AST ground truth, that antibiotic is considered
    predicted.
    """

    predicted = set()

    # --------------------------------------------------------
    # Resistance classes
    # --------------------------------------------------------

    categories = [
        x.strip().lower()
        for x in resistance.split(";")
    ]

    for category in categories:

        if category in RESISTANCE_TO_ABX:
            predicted.update(
                RESISTANCE_TO_ABX[category]
            )

    # --------------------------------------------------------
    # Product description
    # --------------------------------------------------------

    product_lower = product.lower()

    for antibiotic in GT:

        if antibiotic.lower() in product_lower:
            predicted.add(antibiotic)

    return predicted


# ============================================================
# Read Abricate results
# ============================================================

predicted_antibiotics = set()
hit_details = []

with open(INPUT_FILE, "r") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        gene = row["GENE"]
        resistance = row["RESISTANCE"]
        product = row["PRODUCT"]

        predicted = predict_antibiotics(
            resistance,
            product
        )

        predicted_antibiotics.update(predicted)

        if predicted:

            hit_details.append(
                {
                    "gene": gene,
                    "resistance": resistance,
                    "product": product,
                    "predicted": sorted(predicted),
                }
            )


# ============================================================
# Calculate metrics
# ============================================================

TP = sorted(
    GT.intersection(predicted_antibiotics)
)

FN = sorted(
    GT.difference(predicted_antibiotics)
)

tp = len(TP)
fn = len(FN)

if tp + fn > 0:
    recall = tp / (tp + fn)
else:
    recall = 0.0


# ============================================================
# Save results
# ============================================================

with open(OUTPUT_FILE, "w") as out:

    out.write(
        "Abricate / wf-metagenomics AST evaluation\n"
    )
    out.write(
        "========================================\n\n"
    )

    out.write(
        f"Input file: {INPUT_FILE}\n\n"
    )

    out.write("Ground truth antibiotics:\n")
    for abx in sorted(GT):
        out.write(f"  - {abx}\n")

    out.write("\nPredicted antibiotics:\n")
    for abx in sorted(predicted_antibiotics):
        out.write(f"  - {abx}\n")

    out.write("\nTrue positives (TP):\n")
    for abx in TP:
        out.write(f"  - {abx}\n")

    out.write("\nFalse negatives (FN):\n")
    for abx in FN:
        out.write(f"  - {abx}\n")

    out.write("\nMetrics:\n")
    out.write(f"TP = {tp}\n")
    out.write(f"FN = {fn}\n")
    out.write(f"Recall = {recall:.4f}\n")

    out.write(
        "\nAbricate hits contributing to predictions:\n"
    )

    for hit in hit_details:

        out.write(
            f"\nGene: {hit['gene']}\n"
        )

        out.write(
            f"Resistance: {hit['resistance']}\n"
        )

        out.write(
            f"Product: {hit['product']}\n"
        )

        out.write(
            "Predicted AST antibiotics: "
            + ", ".join(hit["predicted"])
            + "\n"
        )


# ============================================================
# Print summary
# ============================================================

print("Evaluation completed.")
print(f"Ground truth: {sorted(GT)}")
print(f"Predicted:    {sorted(predicted_antibiotics)}")
print(f"TP:           {TP}")
print(f"FN:           {FN}")
print(f"Recall:       {recall:.3f}")
print(f"\nResults saved to:\n{OUTPUT_FILE}")
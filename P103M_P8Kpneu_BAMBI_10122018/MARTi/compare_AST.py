
#!/usr/bin/env python3

import csv
import sys
import os

# ============================================================
# PATHS
# ============================================================

BASE_DIR = "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018"

MARTI_FILE = os.path.join(
    BASE_DIR,
    "benchmarking/MARTi/MARTi_AMR_results.tsv"
)

sys.path.append(BASE_DIR)

from ground_truth_AST import antibiotics


# ============================================================
# GROUND TRUTH
# ============================================================

GT = set(antibiotics)


# ============================================================
# MARTi DRUG CLASS -> ANTIBIOTIC MAPPING
# ============================================================

CLASS_TO_ABX = {
    # Aminoglycosides
    "aminoglycoside": {
        "Gentamicin"
    },

    # Beta-lactams
    # Broad mapping, as used for the other workflows
    "penicillin beta-lactam": {
        "Benzylpenicillin",
        "Amoxicillin",
        "Meropenem"
    },

    "beta-lactam": {
        "Benzylpenicillin",
        "Amoxicillin",
        "Meropenem"
    },

    # Carbapenems
    "carbapenem": {
        "Meropenem"
    },

    # Glycopeptides
    "glycopeptide": {
        "Vancomycin"
    },

    # Mupirocin
    "mupirocin": {
        "Mupirocin"
    },

    "mupirocin-like antibiotic": {
        "Mupirocin"
    },

    # Nitroimidazoles
    "nitroimidazole": {
        "Metronidazole"
    },

    "nitroimidazole antibiotic": {
        "Metronidazole"
    },
}


# ============================================================
# MARTi -> ANTIBIOTICS
# ============================================================

def predict_antibiotics(row):
    """
    Convert one MARTi result into one or more antibiotics
    represented in the AST ground truth.

    The mapping follows the broad class-based strategy used
    for the other workflows.
    """

    predicted = set()

    name = row.get("name", "").strip().lower()
    short_name = row.get("shortName", "").strip().lower()
    drug_class = row.get("drugClass", "").strip().lower()
    gene_family = row.get("geneFamily", "").strip().lower()
    resistance_mechanism = row.get(
        "resistanceMechanism", ""
    ).strip().lower()
    description = row.get(
        "description", ""
    ).strip().lower()

    # --------------------------------------------------------
    # Combine relevant MARTi fields
    # --------------------------------------------------------

    text = " ".join([
        name,
        short_name,
        drug_class,
        gene_family,
        resistance_mechanism,
        description
    ])

    # --------------------------------------------------------
    # Match drug classes
    # --------------------------------------------------------

    # drugClass may contain multiple classes separated by ";"
    classes = [
        x.strip()
        for x in drug_class.split(";")
        if x.strip()
    ]

    for cls in classes:

        if cls in CLASS_TO_ABX:
            predicted.update(CLASS_TO_ABX[cls])

    # --------------------------------------------------------
    # Explicit antibiotic names
    # --------------------------------------------------------

    if "gentamicin" in text:
        predicted.add("Gentamicin")

    if "benzylpenicillin" in text:
        predicted.add("Benzylpenicillin")

    if "amoxicillin" in text:
        predicted.add("Amoxicillin")

    if "meropenem" in text:
        predicted.add("Meropenem")

    if "vancomycin" in text:
        predicted.add("Vancomycin")

    if "mupirocin" in text:
        predicted.add("Mupirocin")

    if "metronidazole" in text:
        predicted.add("Metronidazole")

    # --------------------------------------------------------
    # Return only antibiotics present in the GT
    # --------------------------------------------------------

    return predicted & GT


# ============================================================
# READ MARTi RESULTS
# ============================================================

predicted_antibiotics = set()

n_rows = 0

with open(MARTI_FILE, "r", encoding="utf-8") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        n_rows += 1

        predicted = predict_antibiotics(row)

        predicted_antibiotics.update(predicted)


# ============================================================
# EVALUATION
# ============================================================

TP = predicted_antibiotics & GT
FN = GT - predicted_antibiotics

tp = len(TP)
fn = len(FN)

recall = tp / len(GT) if GT else 0.0


# ============================================================
# OUTPUT
# ============================================================

print("=" * 60)
print("MARTi AST-LEVEL EVALUATION")
print("=" * 60)

print(f"\nInput file:")
print(MARTI_FILE)

print(f"\nMARTi records analyzed: {n_rows}")

print("\nGround truth antibiotics:")
for abx in sorted(GT):
    print(f"  - {abx}")

print("\nPredicted antibiotics:")
for abx in sorted(predicted_antibiotics):
    print(f"  - {abx}")

print("\nTrue positives:")
for abx in sorted(TP):
    print(f"  - {abx}")

print("\n" + "-" * 60)

print(f"TP      = {tp}")
print(f"FN      = {fn}")
print(f"Recall  = {recall:.4f}")
print(f"Recall  = {recall * 100:.2f}%")

print("=" * 60)

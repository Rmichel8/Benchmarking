
#!/usr/bin/env python3

import csv
import os
import sys

# ============================================================
# PATHS
# ============================================================

BASE_DIR = "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018"

RUARGPORE_FILE = (
    "/nfsd/sysbiobig/Riccardom/"
    "P103M_P8Kpneu_BAMBI_10122018/"
    "benchmarking/RUARGpore/"
    "P103M_P8Kpneu_all_10122018/"
    "P103M_P8Kpneu_all_10122018.fa.uniq_ARGpore_2026-07-17.13:42:29/"
    "P103M_P8Kpneu_all_10122018.fa.uniq_arg.w.taxa.tab"
)

sys.path.append(BASE_DIR)

from ground_truth_AST import antibiotics


# ============================================================
# GROUND TRUTH
# ============================================================

GT = set(antibiotics)


# ============================================================
# RUARGPORE SUBTYPE -> ANTIBIOTIC MAPPING
# ============================================================

SUBTYPE_TO_ABX = {

    # --------------------------------------------------------
    # Aminoglycosides
    # --------------------------------------------------------
    "aminoglycoside": {
        "Gentamicin"
    },

    # --------------------------------------------------------
    # Beta-lactams
    # Broad mapping, consistent with the other workflows
    # --------------------------------------------------------
    "beta-lactam": {
        "Benzylpenicillin",
        "Amoxicillin",
        "Meropenem"
    },

    # --------------------------------------------------------
    # Carbapenems
    # --------------------------------------------------------
    "carbapenem": {
        "Meropenem"
    },

    # --------------------------------------------------------
    # Glycopeptides
    # --------------------------------------------------------
    "glycopeptide": {
        "Vancomycin"
    },

    "vancomycin": {
        "Vancomycin"
    },

    # --------------------------------------------------------
    # Mupirocin
    # --------------------------------------------------------
    "mupirocin": {
        "Mupirocin"
    },

    "mupirocin-like antibiotic": {
        "Mupirocin"
    },

    # --------------------------------------------------------
    # Nitroimidazoles
    # --------------------------------------------------------
    "nitroimidazole": {
        "Metronidazole"
    },

    "nitroimidazole antibiotic": {
        "Metronidazole"
    },

    # --------------------------------------------------------
    # Explicit antibiotic names
    # --------------------------------------------------------
    "gentamicin": {
        "Gentamicin"
    },

    "benzylpenicillin": {
        "Benzylpenicillin"
    },

    "amoxicillin": {
        "Amoxicillin"
    },

    "meropenem": {
        "Meropenem"
    },

    "metronidazole": {
        "Metronidazole"
    },

    "mupirocin": {
        "Mupirocin"
    },

    "vancomycin": {
        "Vancomycin"
    },
}


# ============================================================
# PREDICT ANTIBIOTICS FROM ONE RUARGPORE RECORD
# ============================================================

def predict_antibiotics(row):
    """
    Convert one RUARGpore result into one or more antibiotics
    represented in the AST ground truth.
    """

    predicted = set()

    subtype = row.get("subtype", "").strip().lower()
    type_field = row.get("type", "").strip().lower()

    # --------------------------------------------------------
    # Use subtype
    # --------------------------------------------------------

    if subtype in SUBTYPE_TO_ABX:
        predicted.update(SUBTYPE_TO_ABX[subtype])

    # --------------------------------------------------------
    # Also inspect the type column
    # --------------------------------------------------------

    if type_field in SUBTYPE_TO_ABX:
        predicted.update(SUBTYPE_TO_ABX[type_field])

    # --------------------------------------------------------
    # Explicit antibiotic names in subtype/type
    # --------------------------------------------------------

    text = f"{subtype} {type_field}"

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
    # Keep only antibiotics in the AST ground truth
    # --------------------------------------------------------

    return predicted & GT


# ============================================================
# READ RUARGPORE RESULTS
# ============================================================

predicted_antibiotics = set()

n_rows = 0

with open(RUARGPORE_FILE, "r", encoding="utf-8") as f:

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
print("RUARGpore AST-LEVEL EVALUATION")
print("=" * 60)

print(f"\nInput file:")
print(RUARGPORE_FILE)

print(f"\nRUARGpore records analyzed: {n_rows}")

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

print("\n" + "-" * 60)

print(f"TP      = {tp}")
print(f"FN      = {fn}")
print(f"Recall  = {recall:.4f}")
print(f"Recall  = {recall * 100:.2f}%")

print("=" * 60)

#!/usr/bin/env python3

import csv
import os
import sys
import re
import sys

sys.path.append(
    "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018"
)
from ground_truth_AST import antibiotics


# ============================================================
# GROUND TRUTH
# ============================================================

GT = set(antibiotics)


# ============================================================
# AMRFINDER -> ANTIBIOTIC MAPPING
# ============================================================

def predict_antibiotics(row):
    """
    Convert one AMRFinder result into one or more antibiotics
    represented in the AST ground truth.

    The mapping is intentionally restricted to the antibiotics
    present in ground_truth_AST.py.
    """

    predicted = set()

    element_symbol = row.get("Element symbol", "").strip().lower()
    element_name = row.get("Element name", "").strip().lower()
    subclass = row.get("Subclass", "").strip().lower()
    amr_class = row.get("Class", "").strip().lower()
    closest_name = row.get("Closest reference name", "").strip().lower()

    text = " ".join([
        element_symbol,
        element_name,
        subclass,
        amr_class,
        closest_name
    ])

    # --------------------------------------------------------
    # Beta-lactams
    # --------------------------------------------------------
    #
    # The AST ground truth contains three beta-lactams:
    # Benzylpenicillin, Amoxicillin and Meropenem.
    #
    # A beta-lactam resistance determinant is therefore
    # considered evidence for the beta-lactam antibiotics
    # represented in the ground truth.
    #
    if amr_class == "beta-lactam":
        predicted.update({
            "Benzylpenicillin",
            "Amoxicillin",
            "Meropenem"
        })

    # --------------------------------------------------------
    # Gentamicin
    # --------------------------------------------------------
    #
    # Avoid mapping every aminoglycoside determinant to
    # Gentamicin. We require explicit evidence that the
    # determinant/reference is associated with gentamicin.
    #
    elif amr_class == "aminoglycoside":
        if "gentamicin" in text:
            predicted.add("Gentamicin")

    # --------------------------------------------------------
    # Vancomycin
    # --------------------------------------------------------

    if (
        "vancomycin" in text
        or subclass == "vancomycin"
        or amr_class == "glycopeptide"
    ):
        predicted.add("Vancomycin")

    # --------------------------------------------------------
    # Mupirocin
    # --------------------------------------------------------

    if (
        "mupirocin" in text
        or "mupA" in element_symbol
        or "mupB" in element_symbol
    ):
        predicted.add("Mupirocin")

    # --------------------------------------------------------
    # Metronidazole
    # --------------------------------------------------------

    if (
        "metronidazole" in text
        or "nitroimidazole" in text
    ):
        predicted.add("Metronidazole")

    # --------------------------------------------------------
    # Keep only antibiotics present in the GT
    # --------------------------------------------------------

    return predicted & GT


# ============================================================
# READ AMRFINDER RESULTS
# ============================================================

def read_amrfinder(tsv_file):

    predicted = set()
    detected_rows = []

    with open(tsv_file, encoding="utf-8") as f:

        reader = csv.DictReader(f, delimiter="\t")

        required_columns = [
            "Element symbol",
            "Element name",
            "Class",
            "Subclass",
            "Closest reference name"
        ]

        missing = [
            col for col in required_columns
            if col not in reader.fieldnames
        ]

        if missing:
            raise RuntimeError(
                "Missing columns in AMRFinder file: "
                + ", ".join(missing)
            )

        for row in reader:

            row_predicted = predict_antibiotics(row)

            if row_predicted:
                predicted.update(row_predicted)

                detected_rows.append({
                    "Element symbol": row["Element symbol"],
                    "Element name": row["Element name"],
                    "Class": row["Class"],
                    "Subclass": row["Subclass"],
                    "Closest reference name":
                        row["Closest reference name"],
                    "Predicted antibiotics":
                        ";".join(sorted(row_predicted))
                })

    return predicted, detected_rows


# ============================================================
# EVALUATION
# ============================================================

def evaluate(tsv_file):

    predicted, detected_rows = read_amrfinder(tsv_file)

    tp = predicted & GT
    fn = GT - predicted

    recall = (
        len(tp) / len(GT)
        if GT
        else 0
    )

    sample = os.path.basename(tsv_file)

    if sample.endswith("_amrfinder_results.tsv"):
        sample = sample.replace(
            "_amrfinder_results.tsv", ""
        )

    print()
    print("=" * 70)
    print("AMRFinder AST evaluation")
    print("=" * 70)

    print("Sample:", sample)

    print()
    print("Ground truth:")
    print(sorted(GT))

    print()
    print("Predicted antibiotics:")
    print(sorted(predicted))

    print()
    print("TP:")
    print(sorted(tp))

    print()
    print("FN:")
    print(sorted(fn))

    print()
    print("TP:", len(tp))
    print("FN:", len(fn))
    print(f"Recall: {recall:.4f}")

    # --------------------------------------------------------
    # Detailed evidence
    # --------------------------------------------------------

    print()
    print("=" * 70)
    print("AMRFinder hits contributing to predictions")
    print("=" * 70)

    if detected_rows:

        for hit in detected_rows:

            print()
            print("Element:", hit["Element symbol"])
            print("Name:", hit["Element name"])
            print("Class:", hit["Class"])
            print("Subclass:", hit["Subclass"])
            print(
                "Closest reference:",
                hit["Closest reference name"]
            )
            print(
                "Predicted:",
                hit["Predicted antibiotics"]
            )

    else:
        print("No AMRFinder hits mapped to the AST ground truth.")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 2:

        print(
            "Usage: python compare_AST.py "
            "amrfinder_results.tsv"
        )

        sys.exit(1)

    evaluate(sys.argv[1])


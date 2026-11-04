#!/usr/bin/env python3
import os
import sys

from ground_truth import (
    B1_TP, B1_TN,
    B2_TP, B2_TN,
    B3_TP, B3_TN,
    B4_TP, B4_TN,
    B5_TP, B5_TN,
    B6_TP, B6_TN,
    B7_TP, B7_TN,
    B8_TP, B8_TN,
    B9_TP, B9_TN,
    B10_TP, B10_TN,
    B11_TP, B11_TN,
)

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

CARD_TO_ABX = {
    "tetracycline": ["Tetracycline", "Minocycline", "Tigecycline"],
    "fluoroquinolone": ["Ciprofloxacin", "Levofloxacin", "Norfloxacin"],
    "penam": ["Penicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
              "Oxacillin", "Cefuroxime"],
    "cephalosporin": ["Ceftazidime", "Cefepime", "Ceftriaxone", "Cefuroxime"],
    "aminoglycoside": ["Gentamicin", "Amikacin", "Tobramycin"],
    "glycylcycline": ["Tigecycline"],
    "colistin": ["Colistin"],
    "peptide": ["Colistin"],
    "acridine_dye": [],
}


def extract_card_classes(tsv_path):
    classes = set()

    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            for c in cols[-1].split(";"):
                c = c.strip()
                if c:
                    classes.add(c)
    return classes


def map_to_antibiotics(card_classes):
    abx = set()

    for c in card_classes:
        if c in CARD_TO_ABX:
            abx.update(CARD_TO_ABX[c])
        else:
            # se non mappato, lo ignoriamo ma potresti loggarlo
            pass

    return abx

import re

def get_benchmark(path):
    base = os.path.basename(path)

    match = re.search(r"_B(\d+)_", base)
    if match:
        return f"B{match.group(1)}"

    match = re.search(r"(B10|B11|B\d)", base)
    if match:
        return match.group(1)

    parent = os.path.basename(os.path.dirname(path))
    match = re.search(r"(B10|B11|B\d)", parent)
    if match:
        return match.group(1)

    return None
import csv

def evaluate(tsv):

    bench = get_benchmark(tsv)

    if bench is None:
        print("No benchmark found")
        return

    tp_truth, tn_truth = GT[bench]

    # CARD classes -> antibiotici
    predicted = map_to_antibiotics(extract_card_classes(tsv))

    tp_set = predicted & tp_truth
    fp_set = predicted & tn_truth
    fn_set = tp_truth - predicted
    tn_set = tn_truth - (tn_truth & predicted)

    TP = len(tp_set)
    FP = len(fp_set)
    FN = len(fn_set)
    TN = len(tn_set)

    sample = os.path.basename(tsv).replace("_amr_results.tsv", "")

    print("Sample:", sample)
    print("Benchmark:", bench)
    print("Predicted:", sorted(predicted))
    print("TP:", TP)
    print("FP:", FP)
    print("FN:", FN)
    print("TN:", TN)

    out_csv = tsv.replace(".tsv", "_compare.csv")

    with open(out_csv, "w", newline="") as out:

        writer = csv.writer(out)

        writer.writerow([
            "sample",
            "benchmark",
            "TP",
            "FP",
            "TN",
            "FN"
        ])

        writer.writerow([
            sample,
            bench,
            TP,
            FP,
            TN,
            FN
        ])

    print(f"Wrote {out_csv}")
if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: compare.py <abricate_results.tsv>")
        sys.exit(1)

    evaluate(sys.argv[1])
#!/usr/bin/env python3

import os
import re
import sys
import csv
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

RESFINDER_TO_GT = {

    # β-lattamici
    "Amoxicillin": "Amoxicillin/clavulanic acid",
    "Ampicillin": "Amoxicillin/clavulanic acid",
    "Piperacillin": "Amoxicillin/clavulanic acid",

    "Ceftriaxone": "Ceftriaxone",
    "Cefuroxime": "Cefuroxime",
    "Ceftazidime": "Ceftazidime",
    "Cefepime": "Cefepime",
    "Cefoperazone": "Cefoperazone/sulbactam",
    "Aztreonam": "Aztreonam",

    "Imipenem": "Imipenem",
    "Meropenem": "Meropenem",
    "Ertapenem": "Ertapenem",

    # aminoglicosidi
    "Gentamicin": "Gentamicin",
    "Amikacin": "Amikacin",
    "Tobramycin": "Tobramycin",

    # fluorochinoloni
    "Ciprofloxacin": "Ciprofloxacin",
    "Levofloxacin": "Levofloxacin",
    "Norfloxacin": "Norfloxacin",

    # altri
    "Colistin": "Colistin",
    "Chloramphenicol": "Chloramphenicol",

    # tetracicline
    "Tetracycline": "Tetracycline",
    "Minocycline": "Minocycline",
    "Doxycycline": "Tetracycline",
    "Tigecycline": "Tegacyclin",   

    # macrolidi / altri
    "Erythromycin": "Erythromycin",
    "Clindamycin": "Clindamycin",
    "Linezolid": "Linezolid",
    "Vancomycin": "Vancomycin",
    "Teicoplanin": "Teicoplanin",

    # antifosforici ecc
    "Rifampicin": "Rifampicin",
    "Nitrofurantoin": "Nitrofurantoin",
    "Fusidic acid": "Fusidic acid",
    "Mupirocin": "Mupirocin",
    "Daptomycin": "Daptomycin",
}

def get_benchmark(path):

    base = os.path.basename(path)

    m = re.search(r"(B10|B11|B\d+)", base)
    if m:
        return m.group(1)

    m = re.search(r"(B10|B11|B\d+)", path)
    if m:
        return m.group(1)

    return None


def extract_resfinder_antibiotics(tsv_path):

    antibiotics = set()

    with open(tsv_path) as f:

        header = f.readline().rstrip().split("\t")

        try:
            pheno_idx = header.index("Phenotype")
        except ValueError:
            raise RuntimeError("Phenotype column not found")

        for line in f:

            if not line.strip():
                continue

            cols = line.rstrip().split("\t")

            if len(cols) <= pheno_idx:
                continue

            phenos = cols[pheno_idx]

            for ab in phenos.split(","):

                ab = ab.strip()

                if not ab:
                    continue

                if ab in RESFINDER_TO_GT:
                    antibiotics.add(RESFINDER_TO_GT[ab])

    return antibiotics

def evaluate(tsv):

    bench = get_benchmark(tsv)

    if bench is None:
        print("No benchmark found")
        return

    tp_truth, tn_truth = GT[bench]

    predicted_all = extract_resfinder_antibiotics(tsv)
    predicted = predicted_all & (tp_truth | tn_truth)

    tp_set = predicted & tp_truth
    fp_set = predicted & tn_truth
    fn_set = tp_truth - predicted
    tn_set = tn_truth - predicted

    TP = len(tp_set)
    FP = len(fp_set)
    FN = len(fn_set)
    TN = len(tn_set)

    sample = os.path.basename(tsv).replace("_amr_results.tsv", "")

    print("Sample:", sample)
    print("Benchmark:", bench)
    print("Total:", sorted(predicted_all))
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

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: compare_resfinder.py results.tsv")
        sys.exit(1)

    evaluate(sys.argv[1])
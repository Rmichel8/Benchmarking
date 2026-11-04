#!/usr/bin/env python3

import os
import re
import csv
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




INPUT_DIR = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_resfinder/"
    "RUARGPOREPRJNA821338_resfinder_results/"
    "split_by_sample"
)

OUTDIR = INPUT_DIR

OUTFILE = os.path.join(
    OUTDIR,
    "resfinder_sample_results.tsv"
)

SUMMARY = os.path.join(
    OUTDIR,
    "resfinder_summary_metrics.csv"
)

# GROUND TRUTH

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

# RESFINDER -> GROUND TRUTH MAPPING

RESFINDER_TO_GT = {

    # Beta-lactams

    "Amoxicillin": "Amoxicillin/clavulanic acid",
    "Amoxicillin+Clavulanic acid":
        "Amoxicillin/clavulanic acid",

    "Ampicillin": "Amoxicillin/clavulanic acid",
    "Ampicillin+Clavulanic acid":
        "Amoxicillin/clavulanic acid",

    "Piperacillin": "Amoxicillin/clavulanic acid",
    "Piperacillin+Tazobactam":
        "Amoxicillin/clavulanic acid",

    "Ceftriaxone": "Ceftriaxone",
    "Cefuroxime": "Cefuroxime",
    "Ceftazidime": "Ceftazidime",
    "Cefepime": "Cefepime",

    "Cefoperazone": "Cefoperazone/sulbactam",
    "Cefoperazone+Sulbactam":
        "Cefoperazone/sulbactam",

    "Aztreonam": "Aztreonam",

    "Imipenem": "Imipenem",
    "Meropenem": "Meropenem",
    "Ertapenem": "Ertapenem",

    "Cefoxitin": "Cefoxitin",
    "Cefotaxime": "Cefotaxime",

    "Cephalothin": "Cephalothin",
    "Cephalotin": "Cephalothin",

    # Aminoglycosides

    "Gentamicin": "Gentamicin",
    "Amikacin": "Amikacin",
    "Tobramycin": "Tobramycin",
    "Kanamycin": "Kanamycin",
    "Neomycin": "Neomycin",
    "Paromomycin": "Paromomycin",
    "Ribostamycin": "Ribostamycin",

    # Fluoroquinolones

    "Ciprofloxacin": "Ciprofloxacin",
    "Levofloxacin": "Levofloxacin",
    "Norfloxacin": "Norfloxacin",

    "Fluoroquinolone": None,

    # Tetracyclines

    "Tetracycline": "Tetracycline",
    "Minocycline": "Minocycline",
    "Doxycycline": "Tetracycline",
    "Tigecycline": "Tegacyclin",

    # Macrolides / lincosamides

    "Erythromycin": "Erythromycin",
    "Azithromycin": "Azithromycin",
    "Clindamycin": "Clindamycin",
    "Lincomycin": "Lincomycin",

    # Glycopeptides

    "Vancomycin": "Vancomycin",
    "Teicoplanin": "Teicoplanin",

    # Other antibiotics

    "Colistin": "Colistin",
    "Chloramphenicol": "Chloramphenicol",
    "Florfenicol": "Florfenicol",

    "Rifampicin": "Rifampicin",
    "Nitrofurantoin": "Nitrofurantoin",
    "Fusidic acid": "Fusidic acid",
    "Mupirocin": "Mupirocin",
    "Daptomycin": "Daptomycin",

    "Fosfomycin": "Fosfomycin",
    "Sulfamethoxazole": "Sulfamethoxazole",
    "Trimethoprim": "Trimethoprim",

    # Other drugs

    "Spectinomycin": "Spectinomycin",
    "Streptomycin": "Streptomycin",

    "Nalidixic acid": "Nalidixic acid",

    "Penicillin": "Penicillin",

    "Tiamulin": "Tiamulin",

    "Quinupristin": "Quinupristin",
    "Pristinamycin IA": "Pristinamycin IA",
    "Pristinamycin IIA": "Pristinamycin IIA",
    "Virginiamycin M": "Virginiamycin M",
    "Virginiamycin S": "Virginiamycin S",

    # Ignored

    "Unknown Aminoglycoside": None,
    "Unknown Beta-lactam": None,
}

# GET BENCHMARK

def get_benchmark(filename):

    m = re.search(
        r"(B10|B11|B[1-9])",
        filename
    )

    if m:
        return m.group(1)

    return None

# NORMALIZATION

def normalize_antibiotic(name):

    name = name.strip()

    if not name:
        return None

    name = re.sub(r"\s+", " ", name)

    replacements = {

        "Amoxicillin + Clavulanic acid":
            "Amoxicillin+Clavulanic acid",

        "Ampicillin + Clavulanic acid":
            "Ampicillin+Clavulanic acid",

        "Piperacillin + Tazobactam":
            "Piperacillin+Tazobactam",

        "Cefoperazone + Sulbactam":
            "Cefoperazone+Sulbactam",

        "Cephalotin":
            "Cephalothin",
    }

    return replacements.get(name, name)

# EXTRACT RESFINDER ANTIBIOTICS

def predicted_antibiotics(tsv_path):

    antibiotics = set()

    with open(
        tsv_path,
        encoding="utf-8"
    ) as f:

        header = (
            f.readline()
            .rstrip("\n")
            .split("\t")
        )

        try:
            pheno_idx = header.index("Phenotype")

        except ValueError:

            raise RuntimeError(
                "Phenotype column not found in: "
                + tsv_path
            )

        for line in f:

            if not line.strip():
                continue

            cols = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(cols) <= pheno_idx:
                continue

            phenotypes = cols[pheno_idx]

            for antibiotic in phenotypes.split(","):

                antibiotic = normalize_antibiotic(
                    antibiotic
                )

                if not antibiotic:
                    continue

                if antibiotic not in RESFINDER_TO_GT:
                    continue

                mapped = RESFINDER_TO_GT[
                    antibiotic
                ]

                if mapped is not None:
                    antibiotics.add(mapped)

    return antibiotics

# EVALUATE ONE SAMPLE

def evaluate(sample, predicted):

    gt_tp, gt_tn = GT[sample]

    tp = len(predicted & gt_tp)

    fp = len(predicted & gt_tn)

    fn = len(gt_tp - predicted)

    tn = len(gt_tn - predicted)

    return tp, fp, tn, fn

# MAIN

def main():

    if not os.path.isdir(INPUT_DIR):

        print(
            "ERROR: Input directory not found:"
        )

        print(INPUT_DIR)

        sys.exit(1)


    results = []

    TOTAL_TP = 0
    TOTAL_FP = 0
    TOTAL_TN = 0
    TOTAL_FN = 0

    # PROCESS ALL SAMPLES

    print()
    print("=" * 60)
    print("ResFinder benchmark")
    print("=" * 60)
    print()
    print("Input directory:")
    print(INPUT_DIR)
    print()


    for file in sorted(os.listdir(INPUT_DIR)):

        if not file.endswith(".tsv"):
            continue

        sample = get_benchmark(file)

        if sample is None:
            print(
                "WARNING: benchmark not found:",
                file
            )
            continue

        path = os.path.join(
            INPUT_DIR,
            file
        )

        preds = predicted_antibiotics(path)

        tp, fp, tn, fn = evaluate(
            sample,
            preds
        )


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


        print(
            f"{sample:>3} | "
            f"TP={tp:<3} "
            f"FP={fp:<3} "
            f"TN={tn:<3} "
            f"FN={fn:<3}"
        )

    # SAVE PER-SAMPLE RESULTS
    with open(
        OUTFILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.writer(
            f,
            delimiter="\t"
        )

        writer.writerow([
            "sample",
            "TP",
            "FP",
            "TN",
            "FN",
            "predicted_antibiotics"
        ])

        writer.writerows(results)

    # CALCULATE GLOBAL METRICS

    Precision = (
        TOTAL_TP /
        (TOTAL_TP + TOTAL_FP)
        if TOTAL_TP + TOTAL_FP
        else 0
    )

    Recall = (
        TOTAL_TP /
        (TOTAL_TP + TOTAL_FN)
        if TOTAL_TP + TOTAL_FN
        else 0
    )

    Specificity = (
        TOTAL_TN /
        (TOTAL_TN + TOTAL_FP)
        if TOTAL_TN + TOTAL_FP
        else 0
    )

    Accuracy = (
        (TOTAL_TP + TOTAL_TN) /
        (
            TOTAL_TP +
            TOTAL_TN +
            TOTAL_FP +
            TOTAL_FN
        )
        if (
            TOTAL_TP +
            TOTAL_TN +
            TOTAL_FP +
            TOTAL_FN
        )
        else 0
    )

    F1 = (
        2 * Precision * Recall /
        (Precision + Recall)
        if Precision + Recall
        else 0
    )

    # MCC

    denominator = (
        (TOTAL_TP + TOTAL_FP) *
        (TOTAL_TP + TOTAL_FN) *
        (TOTAL_TN + TOTAL_FP) *
        (TOTAL_TN + TOTAL_FN)
    )

    if denominator:

        MCC = (
            (
                TOTAL_TP * TOTAL_TN
                -
                TOTAL_FP * TOTAL_FN
            )
            /
            (denominator ** 0.5)
        )

    else:

        MCC = 0


    print()
    print("=" * 60)
    print("GLOBAL RESULTS")
    print("=" * 60)

    print()
    print("TP:", TOTAL_TP)
    print("FP:", TOTAL_FP)
    print("TN:", TOTAL_TN)
    print("FN:", TOTAL_FN)

    print()
    print("Precision :", round(Precision, 4))
    print("Recall    :", round(Recall, 4))
    print("F1        :", round(F1, 4))
    print("Specificity:", round(Specificity, 4))
    print("Accuracy  :", round(Accuracy, 4))
    print("MCC       :", round(MCC, 4))

    # SAVE SUMMARY

    with open(
        SUMMARY,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.writer(f)

        writer.writerow([
            "TP",
            "FP",
            "TN",
            "FN",
            "Precision",
            "Recall",
            "F1",
            "Specificity",
            "Accuracy",
            "MCC"
        ])

        writer.writerow([
            TOTAL_TP,
            TOTAL_FP,
            TOTAL_TN,
            TOTAL_FN,
            round(Precision, 4),
            round(Recall, 4),
            round(F1, 4),
            round(Specificity, 4),
            round(Accuracy, 4),
            round(MCC, 4)
        ])


    print()
    print("Per-sample results:")
    print(OUTFILE)

    print()
    print("Summary metrics:")
    print(SUMMARY)

    print()


if __name__ == "__main__":
    main()
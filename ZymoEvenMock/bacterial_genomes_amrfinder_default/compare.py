TPg=2
FPg=32
FNg=36
drug_classes = [
    "penicillin beta-lactam",
    "tetracyclines",
    "sulfonamides",
    "fluoroquinolones",
    "fosfomycin",
    "polymyxins",
    "aminoglycosides",
    "macrolides",
    "lincosamides",
    "streptogramins",
    "aminocoumarins",
    "peptide_antibiotics",
    "phenicol",
    "cephalosporin",
    "elfamycin"]
import csv
import re
import sys
GT_CLASSES = {
    c.lower().strip()
    for c in drug_classes
}

# Drug class mapping


CLASS_MAP = {

    "beta-lactam": "penicillin beta-lactam",
    "penam": "penicillin beta-lactam",

    "cephalosporin": "cephalosporin",

    "tetracycline": "tetracyclines",

    "aminoglycoside": "aminoglycosides",

    "fluoroquinolone": "fluoroquinolones",

    "sulfonamide": "sulfonamides",

    "macrolide": "macrolides",

    "lincosamide": "lincosamides",

    "streptogramin": "streptogramins",

    "aminocoumarin": "aminocoumarins",

    "glycopeptide": "peptide_antibiotics",

    "phenicol": "phenicol",

    "phosphonic acid": "fosfomycin",

    "polymyxin": "polymyxins",

    "elfamycin": "elfamycin",
}

def precision(tp, fp):
    return tp / (tp + fp) if tp + fp else 0


def recall(tp, fn):
    return tp / (tp + fn) if tp + fn else 0


def f1(tp, fp, fn):

    p = precision(tp, fp)
    r = recall(tp, fn)

    return 2 * p * r / (p + r) if p + r else 0


def extract_amrfinder(tsv):
    pred_classes = set()

    with open(tsv) as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            cls = row["Class"].lower().strip()

            pred_classes.add(CLASS_MAP.get(cls, cls))

    return pred_classes


import os

OUT_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "summary_metrics.csv"
)
def evaluate(tsv):

    pred_classes = extract_amrfinder(tsv)

    TPc = pred_classes & GT_CLASSES
    FPc = pred_classes - GT_CLASSES
    FNc = GT_CLASSES - pred_classes

    print("\nGENES")
    print("TP =", TPg)
    print("FP =", FPg)
    print("FN =", FNg)
    print("Precision =", round(precision(TPg, FPg), 3))
    print("Recall    =", round(recall(TPg, FNg), 3))
    print("F1         =", round(f1(TPg,FPg, FNg), 3))

    print("\nDRUG CLASSES")
    print("TP =", len(TPc))
    print("FP =", len(FPc))
    print("FN =", len(FNc))
    print("Precision =", round(precision(len(TPc), len(FPc)), 3))
    print("Recall    =", round(recall(len(TPc), len(FNc)), 3))
    print("F1         =", round(f1(len(TPc), len(FPc), len(FNc)), 3))
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
        "type",
        "TP",
        "FP",
        "FN",
        "precision",
        "recall",
        "f1"
        ])

        writer.writerow([
        "genes",
        TPg,
        FPg,
        FNg,
        round(precision(TPg, FPg), 3),
        round(recall(TPg, FNg), 3),
        round(f1(TPg, FPg, FNg), 3)
        ])

        writer.writerow([
        "drug_classes",
        len(TPc),
        len(FPc),
        len(FNc),
        round(precision(len(TPc), len(FPc)), 3),
        round(recall(len(TPc), len(FNc)), 3),
        round(f1(len(TPc), len(FPc), len(FNc)), 3)
        ])
    print(f"\nResults written to {OUT_CSV}")



if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: compare_genes.py <amrfinder.tsv>")
        sys.exit(1)

    evaluate(sys.argv[1])


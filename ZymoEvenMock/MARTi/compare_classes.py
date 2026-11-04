#!/usr/bin/env python3

import pandas as pd

# Ground truth drug classes


ground_truth = {
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
    "elfamycin"
}


# CARD -> Ground truth mapping


CLASS_MAP = {

    "penicillin beta-lactam": "penicillin beta-lactam",

    "cephalosporin": "cephalosporin",

    "fluoroquinolone antibiotic": "fluoroquinolones",

    "aminoglycoside antibiotic": "aminoglycosides",

    "macrolide antibiotic": "macrolides",

    "lincosamide antibiotic": "lincosamides",

    "tetracycline antibiotic": "tetracyclines",

    "phenicol antibiotic": "phenicol",

    "phosphonic acid antibiotic": "fosfomycin",

    "aminocoumarin antibiotic": "aminocoumarins",

    "peptide antibiotic": "peptide_antibiotics",

}


marti_file = "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/MARTi/MARTi_AMR_results.tsv"

output_prefix = "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/MARTi/MARTi_drugclass"


# Load

df = pd.read_csv(marti_file, sep="\t")

found = set()

for value in df["drugClass"].fillna(""):

    if value == "":
        continue

    for cls in value.split(";"):

        cls = cls.strip().lower()

        if cls in CLASS_MAP:
            found.add(CLASS_MAP[cls])

        else:
            # mantiene le classi non riconosciute
            found.add(cls)

truth = {x.lower() for x in ground_truth}

# ============================
# Comparison
# ============================

TP = sorted(found & truth)
FP = sorted(found - truth)
FN = sorted(truth - found)

precision = len(TP) / (len(TP) + len(FP)) if TP or FP else 0
recall = len(TP) / (len(TP) + len(FN)) if TP or FN else 0

f1 = (
    2 * precision * recall / (precision + recall)
    if precision + recall > 0
    else 0
)

# Results


print("\n===== MARTi Drug Class benchmark =====\n")

print(f"Ground truth classes : {len(truth)}")
print(f"Detected classes     : {len(found)}")

print(f"\nTP: {len(TP)}")
print(f"FP: {len(FP)}")
print(f"FN: {len(FN)}")

print(f"\nPrecision : {precision:.4f}")
print(f"Recall    : {recall:.4f}")
print(f"F1-score  : {f1:.4f}")

print("\nTrue positives")
for x in TP:
    print(" ", x)

print("\nFalse positives")
for x in FP:
    print(" ", x)

print("\nFalse negatives")
for x in FN:
    print(" ", x)

with open(output_prefix + "_TP.txt", "w") as f:
    f.write("\n".join(TP))

with open(output_prefix + "_FP.txt", "w") as f:
    f.write("\n".join(FP))

with open(output_prefix + "_FN.txt", "w") as f:
    f.write("\n".join(FN))

with open(output_prefix + "_summary.txt", "w") as f:

    f.write("MARTi Drug Class benchmark\n\n")

    f.write(f"Ground truth: {len(truth)}\n")
    f.write(f"Detected: {len(found)}\n\n")

    f.write(f"TP: {len(TP)}\n")
    f.write(f"FP: {len(FP)}\n")
    f.write(f"FN: {len(FN)}\n\n")

    f.write(f"Precision: {precision:.3f}\n")
    f.write(f"Recall: {recall:.3f}\n")
    f.write(f"F1-score: {f1:.3f}\n")
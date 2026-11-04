#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


INPUT = Path(
    "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/RUARGpore65/ZymoEvenMock/ZymoEvenMock.fa.uniq_ARGpore_2026-10-09.16:00:31/ZymoEvenMock.fa.uniq_arg.w.taxa.tab"
)

OUTDIR = INPUT.parent

SUMMARY_FILE = OUTDIR / "summary_metrics.txt"


drug_classes = {
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


# Classi troppo generiche da ignorare
EXCLUDED_CLASSES = {
    "multidrug",
    "unclassified"
}


CLASS_ALIAS = {

    "beta-lactam":
        {"penicillin beta-lactam"},

    "beta lactam":
        {"penicillin beta-lactam"},

    "penicillin":
        {"penicillin beta-lactam"},

    "cephalosporins":
        {"cephalosporin"},

    "aminoglycoside":
        {"aminoglycosides"},

    "macrolide":
        {"macrolides"},

    "tetracycline":
        {"tetracyclines"},

    "fluoroquinolone":
        {"fluoroquinolones"},

    "quinolone":
        {"fluoroquinolones"},

    "polymyxin":
        {"polymyxins"},

    "sulfonamide":
        {"sulfonamides"},

    "chloramphenicol":
        {"phenicol"},
    "macrolide-lincosamide-streptogramin":
        {
            "macrolides",
            "lincosamides",
            "streptogramins"
        },
    "kasugamycin":
    {"aminoglycosides"},
     "bacitracin":
    {"peptide_antibiotics"}

}


def normalize_class(x):

    x = str(x).strip().lower()

    # Ignora classi generiche
    if x in EXCLUDED_CLASSES:
        return set()

    # Classe già presente nella ground truth
    if x in drug_classes:
        return {x}

    # Alias
    if x in CLASS_ALIAS:
        return CLASS_ALIAS[x]

    # Classe non riconosciuta
    return {x}


df = pd.read_csv(
    INPUT,
    sep="\t"
)


predicted = set()

for x in df["type"].dropna():

    normalized_classes = normalize_class(x)

    predicted.update(normalized_classes)


tp = sorted(
    predicted & drug_classes
)

fp = sorted(
    predicted - drug_classes
)

fn = sorted(
    drug_classes - predicted
)


precision = (
    len(tp) / (len(tp) + len(fp))
    if tp or fp
    else 0
)


recall = (
    len(tp) / (len(tp) + len(fn))
    if tp or fn
    else 0
)


f1 = (
    2 * precision * recall / (precision + recall)
    if precision + recall
    else 0
)



rows = []


for c in sorted(drug_classes | predicted):

    rows.append({

        "drug_class": c,

        "ground_truth":
        c in drug_classes,

        "predicted":
        c in predicted,

        "status":
        "TP" if c in tp else
        "FP" if c in fp else
        "FN"

    })


summary = pd.DataFrame(rows)


summary.to_csv(
    OUTDIR / "RUARGpore_drugclass_compare.csv",
    index=False
)


with open(
    SUMMARY_FILE,
    "a"
) as f:

    f.write(
        "\n\n"
        "========================================\n"
        "RUARGpore DRUG CLASS BENCHMARK\n"
        "========================================\n\n"
    )

    f.write(
        f"Ground truth : {len(drug_classes)}\n"
    )

    f.write(
        f"Predicted    : {len(predicted)}\n\n"
    )

    f.write(
        f"TP : {len(tp)}\n"
    )

    f.write(
        f"FP : {len(fp)}\n"
    )

    f.write(
        f"FN : {len(fn)}\n\n"
    )

    f.write(
        f"Precision : {precision:.3f}\n"
    )

    f.write(
        f"Recall    : {recall:.3f}\n"
    )

    f.write(
        f"F1        : {f1:.3f}\n"
    )

print(summary)

print("\nTP:", tp)
print("FP:", fp)
print("FN:", fn)

print()

print(
    f"Precision = {precision:.3f}"
)

print(
    f"Recall    = {recall:.3f}"
)

print(
    f"F1        = {f1:.3f}"
)

print(
    f"\nDrug-class metrics appended to:\n{SUMMARY_FILE}"
)
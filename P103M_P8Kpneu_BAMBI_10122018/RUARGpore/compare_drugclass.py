#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


INPUT = Path(
    "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018/benchmarking/RUARGpore/P103M_P8Kpneu_all_10122018/P103M_P8Kpneu_all_10122018.fa.uniq_ARGpore_2026-07-17.13:42:29/P103M_P8Kpneu_all_10122018.fa.uniq_arg.w.taxa.tab"
)

OUTDIR = INPUT.parent

SUMMARY_FILE = OUTDIR / "summary_metrics.txt"



drug_classes = {
    "aminocoumarin",
    "aminoglycoside",
    "penicillin beta-lactam",
    "cephalosporin",
    "carbapenem",
    "fluoroquinolone",
    "tetracycline",
    "glycylcycline",
    "phenicol",
    "macrolide",
    "lincosamide",
    "streptogramin",
    "pleuromutilin",
    "rifamycin",
    "fosfomycin",
    "diaminopyrimidine",
    "nitrofuran",
    "mupirocin",
    "disinfecting agents and antiseptics"
}



EXCLUDED_CLASSES = {
    "multidrug",
    "unclassified",
    "unknown",
    "unknown antibiotic",
    "unknown beta-lactam"
}



CLASS_ALIAS = {

    # Beta-lactams
    "beta-lactam":
        {"penicillin beta-lactam"},

    "beta lactam":
        {"penicillin beta-lactam"},

    "penicillin":
        {"penicillin beta-lactam"},

    "penicillins":
        {"penicillin beta-lactam"},

    "penicillin beta-lactam":
        {"penicillin beta-lactam"},

    # Cephalosporins
    "cephalosporin":
        {"cephalosporin"},

    "cephalosporins":
        {"cephalosporin"},

    # Carbapenems
    "carbapenem":
        {"carbapenem"},

    "carbapenems":
        {"carbapenem"},

    # Aminoglycosides
    "aminoglycoside":
        {"aminoglycoside"},

    "aminoglycosides":
        {"aminoglycoside"},

    # Macrolides
    "macrolide":
        {"macrolide"},

    "macrolides":
        {"macrolide"},

    # Tetracyclines
    "tetracycline":
        {"tetracycline"},

    "tetracyclines":
        {"tetracycline"},

    # Glycylcyclines
    "glycylcycline":
        {"glycylcycline"},

    "glycylcyclines":
        {"glycylcycline"},

    # Fluoroquinolones
    "fluoroquinolone":
        {"fluoroquinolone"},

    "fluoroquinolones":
        {"fluoroquinolone"},

    "quinolone":
        {"fluoroquinolone"},

    "quinolones":
        {"fluoroquinolone"},


    # Lincosamides
    "lincosamide":
        {"lincosamide"},

    "lincosamides":
        {"lincosamide"},

    # Streptogramins
    "streptogramin":
        {"streptogramin"},

    "streptogramins":
        {"streptogramin"},

    # Aminocoumarins
    "aminocoumarin":
        {"aminocoumarin"},

    "aminocoumarins":
        {"aminocoumarin"},

    # Phenicols
    "phenicol":
        {"phenicol"},

    "phenicols":
        {"phenicol"},

    "chloramphenicol":
        {"phenicol"},

    # Fosfomycin
    "fosfomycin":
        {"fosfomycin"},

    # Rifamycins
    "rifamycin":
        {"rifamycin"},

    "rifamycins":
        {"rifamycin"},

    # Pleuromutilins
    "pleuromutilin":
        {"pleuromutilin"},

    "pleuromutilins":
        {"pleuromutilin"},

    # Diaminopyrimidines
    "diaminopyrimidine":
        {"diaminopyrimidine"},

    "diaminopyrimidines":
        {"diaminopyrimidine"},

    # Nitrofurans
    "nitrofuran":
        {"nitrofuran"},

    "nitrofurans":
        {"nitrofuran"},

    # Combinazioni
    "macrolide-lincosamide-streptogramin":
        {
            "macrolide",
            "lincosamide",
            "streptogramin"
        },

    # Alcune classi ARGPORE specifiche
    "kasugamycin":
        {"aminoglycoside"},

    "bacitracin":
        {"peptide antibiotic"}
}



def normalize_class(x):

    x = str(x).strip().lower()

    if x in EXCLUDED_CLASSES:
        return set()

    if x in drug_classes:
        return {x}

    # Alias
    if x in CLASS_ALIAS:
        return CLASS_ALIAS[x]

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
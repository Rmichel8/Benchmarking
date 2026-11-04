#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re




INPUT = Path(
    "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018/benchmarking/RUARGpore/P103M_P8Kpneu_all_10122018/P103M_P8Kpneu_all_10122018.fa.uniq_ARGpore_2026-07-17.13:42:29/P103M_P8Kpneu_all_10122018.fa.uniq_arg.w.taxa.tab"
)

OUTDIR = INPUT.parent



genes = [
    
        "emrB",
        "acrA",
        "acrB",
        "marA",
        "emrR",
        "CRP",
        "H-NS",
        "mdtB",
        "mdtC",
        "ramA",
        "baeR",
        "cpxA",
        "cpxR",
        "vgaC",
        "smeB",
        "SHV-187",
        "FosA5",
        "alaS",
        "blaZ",
        "mfd",
        "oqxA",
        "oqxB",
        "lsa",
        "tetO",
        "ileS",
        "emeA",
        "dfrE",
        "class A beta-lactamase",
        "tet32",
        "norA",
        "tetM",
        "qacA",
        "aac6-ie-aph2-ia",
        "ErmC"
]

drug_class = [
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
    "mupirocin",
    "disinfecting agents and antiseptics",
    "nitrofuran",
]


GROUND_TRUTH = set(genes)



ALIASES = {

    # Regulatory proteins
    "transcriptional regulatory protein cpxr cpxr": "cpxR",
    "cpxr": "cpxR",

    "camp-regulatory protein": "CRP",
    "camp regulatory protein": "CRP",
    "crp": "CRP",

    "h-ns": "H-NS",
    "hns": "H-NS",
    "DNA-binding_protein_H-NS": "H-NS",

    # EmrB
    "EmrB-QacA family major facilitator transporter": "emrB",
    "emrb": "emrB",

    # Case-insensitive variants
    "fosa5": "FosA5",
    "shv-187": "SHV-187",
    "oqxa": "oqxA",
    "oqxb": "oqxB",

}


# GENE NORMALIZATION

def normalize(name):

    name = str(name).strip()

    # Match exact ground-truth gene
    for gene in GROUND_TRUTH:

        if name.lower() == gene.lower():

            return gene

    low = name.lower()

    # Explicit aliases
    if low in {
        alias.lower(): target
        for alias, target in ALIASES.items()
    }:

        for alias, target in ALIASES.items():

            if low == alias.lower():

                return target

    # FosA family
    #
    # IMPORTANT:
    # FosA5 is the only FosA allele in the ground truth.
    #
    # If RUARGpore reports FosA5 specifically,
    # it is a TP.


    m = re.search(
        r"\b(FosA\d+)\b",
        name,
        re.I
    )

    if m:

        candidate = m.group(1)

        for gene in GROUND_TRUTH:

            if candidate.lower() == gene.lower():

                return gene

        return candidate

    # SHV alleles
    m = re.search(
        r"\b(SHV-\d+)\b",
        name,
        re.I
    )

    if m:

        candidate = m.group(1).upper()

        for gene in GROUND_TRUTH:

            if candidate.lower() == gene.lower():

                return gene

        return candidate

    # Search gene symbols inside longer descriptions
    #
    # Longest names first avoids partial matches.

    for gene in sorted(
        GROUND_TRUTH,
        key=len,
        reverse=True
    ):

        if gene.lower() in low:

            return gene

    return name




df = pd.read_csv(
    INPUT,
    sep="\t"
)

# PREDICTED GENES


predicted = {

    normalize(g)

    for g in df["subtype"]

}


# METRICS


tp = sorted(
    predicted & GROUND_TRUTH
)

fp = sorted(
    predicted - GROUND_TRUTH
)

fn = sorted(
    GROUND_TRUTH - predicted
)


precision = (

    len(tp) /
    (len(tp) + len(fp))

    if len(tp) + len(fp)

    else 0

)


recall = (

    len(tp) /
    (len(tp) + len(fn))

    if len(tp) + len(fn)

    else 0

)


f1 = (

    2 * precision * recall /
    (precision + recall)

    if precision + recall

    else 0

)

rows = []


for gene in sorted(
    GROUND_TRUTH | predicted
):

    if gene in tp:

        status = "TP"

    elif gene in fp:

        status = "FP"

    else:

        status = "FN"


    rows.append({

        "gene": gene,

        "ground_truth":
        gene in GROUND_TRUTH,

        "predicted":
        gene in predicted,

        "status":
        status

    })


summary = pd.DataFrame(
    rows
)


summary.to_csv(

    OUTDIR /
    "summary_compare.csv",

    index=False

)


with open(

    OUTDIR /
    "summary_metrics.txt",

    "w"

) as f:

    f.write(
        "RUARGpore Gene Benchmark\n\n"
    )

    f.write(
        f"Ground truth : "
        f"{len(GROUND_TRUTH)}\n"
    )

    f.write(
        f"Predicted    : "
        f"{len(predicted)}\n\n"
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
        f"Precision : "
        f"{precision:.3f}\n"
    )

    f.write(
        f"Recall    : "
        f"{recall:.3f}\n"
    )

    f.write(
        f"F1        : "
        f"{f1:.3f}\n"
    )


print(summary)

print()

print(
    "TP:",
    len(tp)
)

print(
    "FP:",
    len(fp)
)

print(
    "FN:",
    len(fn)
)

print(
    "Precision:",
    round(precision, 3)
)

print(
    "Recall:",
    round(recall, 3)
)

print(
    "F1:",
    round(f1, 3)
)
#!/usr/bin/env python3

import os
import re
import sys
import csv

from ground_truth import amr_genes, drug_classes

def normalize_gene(gene):

    gene = gene.lower().strip()
    gene = gene.replace("(", "")
    gene = gene.replace(")", "")

    gene = gene.replace("_", "")
    gene = gene.replace("-", "")
    gene = gene.replace(" ", "")

    gene = re.sub(r"[^a-z0-9]", "", gene)

    species = [
        "escherichiacoli",
        "klebsiellapneumoniae",
        "pseudomonasaeruginosa",
        "agrobacteriumfabrum",
        "achromobacterxylosoxidans",
        "acinetobacterbaumannii",
        "enterobactercloacae",
        "haemophilusinfluenzae",
        "morganellamorganii",
        "streptomycesrimosus",
        "bifidobacteriumadolescentis",
        "bifidobacteriumbifidum",
    ]

    for s in species:
        if gene.startswith(s):
            gene = gene[len(s):]
            break

    gene = re.sub(r"conferring.*", "", gene)
    gene = re.sub(r"mutants?.*", "", gene)
    gene = re.sub(r"mutation.*", "", gene)
    gene = re.sub(r"with.*", "", gene)

    gene = re.sub(r"genein.*", "", gene)

    gene = gene.removesuffix("gene")

    if gene.startswith("bla"):
        gene = gene[3:]

    return gene


GT_GENES = {
    normalize_gene(gene)
    for gene in amr_genes
}

GT_CLASSES = {
    c.lower().strip()
    for c in drug_classes
}



CLASS_MAP = {

    "penam":
        "penicillin beta-lactam",

    "penicillin":
        "penicillin beta-lactam",

    "penicillin beta-lactam":
        "penicillin beta-lactam",

    "tetracycline":
        "tetracycline",

    "aminoglycoside":
        "aminoglycoside",

    "fluoroquinolone":
        "fluoroquinolone",

    "macrolide":
        "macrolide",

    "lincosamide":
        "lincosamide",

    "streptogramin":
        "streptogramin",

    "polymyxin":
        "polymyxin",

    "aminocoumarin":
        "aminocoumarin",

    "phenicol":
        "phenicol",

    "glycylcycline":
        "glycylcycline",

    "fosfomycin":
        "fosfomycin",

    "phosphonic acid":
        "fosfomycin",

    "phosphonic acid antibiotic":
        "fosfomycin",

    "cephalosporin":
        "cephalosporin",

    "carbapenem":
        "carbapenem",

    "diaminopyrimidine":
        "diaminopyrimidine",

    "nitrofuran":
        "nitrofuran",

    "pleuromutilin":
        "pleuromutilin",

    "rifamycin":
        "rifamycin",

    "rifamycin antibiotic":
        "rifamycin",

    "aminocoumarin antibiotic":
        "aminocoumarin",

    "aminoglycoside antibiotic":
        "aminoglycoside",

    "fluoroquinolone antibiotic":
        "fluoroquinolone",

    "tetracycline antibiotic":
        "tetracycline",

    "glycylcycline antibiotic":
        "glycylcycline",

    "phenicol antibiotic":
        "phenicol",

    "macrolide antibiotic":
        "macrolide",

    "lincosamide antibiotic":
        "lincosamide",

    "streptogramin antibiotic":
        "streptogramin",

    "diaminopyrimidine antibiotic":
        "diaminopyrimidine",

    "nitrofuran antibiotic":
        "nitrofuran",

    "carbapenem antibiotic":
        "carbapenem",

    "cephalosporin antibiotic":
        "cephalosporin",

    "pleuromutilin antibiotic":
        "pleuromutilin",

    "peptide antibiotic":
        "peptide antibiotic",

    "disinfecting agents and antiseptics":
        "disinfecting agents and antiseptics",
}



OUTDIR = os.path.dirname(os.path.abspath(__file__))

SUMMARY_CSV = os.path.join(
    OUTDIR,
    "summary_gene_compare.csv"
)



def lookup_gene(gene):

    g = normalize_gene(gene)

    if g in GT_GENES:

        return g

    for truth in GT_GENES:

        if g in truth or truth in g:

            return truth

    return None



def precision(tp, fp):

    return tp / (tp + fp) if tp + fp else 0


def recall(tp, fn):

    return tp / (tp + fn) if tp + fn else 0


def f1(tp, fp, fn):

    p = precision(tp, fp)

    r = recall(tp, fn)

    return (
        2 * p * r / (p + r)
        if p + r
        else 0
    )


def extract_rgi(tsv):

    predicted_genes = set()

    predicted_classes = set()

    with open(tsv) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            # ==================================================
            # GENE
            # ==================================================

            raw_gene = row["Best_Hit_ARO"]

            gene = normalize_gene(raw_gene)

            predicted_genes.add(gene)

            # ==================================================
            # DRUG CLASSES
            # ==================================================

            resistance = row["Drug Class"].strip()

            if resistance:

                for c in re.split(
                    r"[;,]",
                    resistance
                ):

                    c = c.lower().strip()

                    # Rimuove la parola antibiotic/antibiotics
                    c = c.replace(
                        " antibiotic",
                        ""
                    )

                    c = c.replace(
                        " antibiotics",
                        ""
                    )

                    c = c.strip()

                    predicted_classes.add(
                        CLASS_MAP.get(
                            c,
                            c
                        )
                    )

    return predicted_genes, predicted_classes

# ============================================================
# SAVE METRICS
# ============================================================

def save_metrics(
    pred_genes,
    TPg,
    FPg,
    FNg,
    pred_classes,
    TPc,
    FPc,
    FNc
):

    metrics_file = os.path.join(
        OUTDIR,
        "summary_metrics.tsv"
    )

    with open(
        metrics_file,
        "w",
        newline=""
    ) as f:

        writer = csv.writer(
            f,
            delimiter="\t"
        )

        writer.writerow([
            "Category",
            "Predicted",
            "TP",
            "FP",
            "FN",
            "Precision",
            "Recall",
            "F1"
        ])

        writer.writerow([
            "Genes",
            len(pred_genes),
            len(TPg),
            len(FPg),
            len(FNg),
            round(
                precision(
                    len(TPg),
                    len(FPg)
                ),
                4
            ),
            round(
                recall(
                    len(TPg),
                    len(FNg)
                ),
                4
            ),
            round(
                f1(
                    len(TPg),
                    len(FPg),
                    len(FNg)
                ),
                4
            )
        ])

        writer.writerow([
            "Drug_classes",
            len(pred_classes),
            len(TPc),
            len(FPc),
            len(FNc),
            round(
                precision(
                    len(TPc),
                    len(FPc)
                ),
                4
            ),
            round(
                recall(
                    len(TPc),
                    len(FNc)
                ),
                4
            ),
            round(
                f1(
                    len(TPc),
                    len(FPc),
                    len(FNc)
                ),
                4
            )
        ])

    print(
        f"\nMetrics saved to:\n{metrics_file}"
    )

def evaluate(tsv):

    pred_genes, pred_classes = extract_rgi(tsv)


    TPg = pred_genes & GT_GENES

    FPg = pred_genes - GT_GENES

    FNg = GT_GENES - pred_genes



    TPc = pred_classes & GT_CLASSES

    FPc = pred_classes - GT_CLASSES

    FNc = GT_CLASSES - pred_classes



    print("\n")
    print("=" * 60)
    print(os.path.basename(tsv))
    print("=" * 60)


    print("\nTrue positive genes")
    print("-------------------")

    for g in sorted(TPg):

        print(" ", g)



    print("\nFalse positive genes")
    print("--------------------")

    for g in sorted(FPg):

        print(" ", g)



    print("\nFalse negative genes")
    print("--------------------")

    for g in sorted(FNg):

        print(" ", g)

    print("\nGENES")

    print("Predicted :", len(pred_genes))

    print("TP :", len(TPg))

    print("FP :", len(FPg))

    print("FN :", len(FNg))

    print(
        "Precision :",
        round(
            precision(
                len(TPg),
                len(FPg)
            ),
            3
        )
    )

    print(
        "Recall    :",
        round(
            recall(
                len(TPg),
                len(FNg)
            ),
            3
        )
    )

    print(
        "F1        :",
        round(
            f1(
                len(TPg),
                len(FPg),
                len(FNg)
            ),
            3
        )
    )


    print("\nTrue positive drug classes")
    print("--------------------------")

    for c in sorted(TPc):

        print(" ", c)



    print("\nFalse positive drug classes")
    print("---------------------------")

    for c in sorted(FPc):

        print(" ", c)


    print("\nFalse negative drug classes")
    print("---------------------------")

    for c in sorted(FNc):

        print(" ", c)


    print("\nDRUG CLASSES")

    print("Predicted :", len(pred_classes))

    print("TP :", len(TPc))

    print("FP :", len(FPc))

    print("FN :", len(FNc))

    print(
        "Precision :",
        round(
            precision(
                len(TPc),
                len(FPc)
            ),
            3
        )
    )

    print(
        "Recall    :",
        round(
            recall(
                len(TPc),
                len(FNc)
            ),
            3
        )
    )

    print(
        "F1        :",
        round(
            f1(
                len(TPc),
                len(FPc),
                len(FNc)
            ),
            3
        )
    )
    save_metrics(
        pred_genes,
        TPg,
        FPg,
        FNg,
        pred_classes,
        TPc,
        FPc,
        FNc
    )

if __name__ == "__main__":

    if len(sys.argv) != 2:

        print(
            "Usage: compare_genes.py <rgi.tsv>"
        )

        sys.exit(1)

    evaluate(
        sys.argv[1]
    )
#!/usr/bin/env python3

import os
import re
import sys
import csv

from ground_truth import amr_genes, drug_classes


# ============================================================
# GROUND TRUTH
# ============================================================

GT_GENES = {
    gene.lower()
    .replace("-", "")
    .replace("_", "")
    .replace(" ", "")
    for gene in amr_genes
}
GENE_ALIASES = {
    "shv": "shv187",
}

GT_CLASSES = {
    c.lower().strip()
    for c in drug_classes
}


# ============================================================
# OUTPUT
# ============================================================

OUTDIR = os.path.dirname(os.path.abspath(__file__))

SUMMARY_CSV = os.path.join(
    OUTDIR,
    "summary_amrfinder_compare.csv"
)


# ============================================================
# AMRFINDER CLASS -> GROUND TRUTH CLASS
# ============================================================

CLASS_MAP = {

    "beta-lactam":
        "penicillin beta-lactam",

    "tetracycline":
        "tetracycline",

    "aminoglycoside":
        "aminoglycoside",

    "fluoroquinolone":
        "fluoroquinolone",
    "quinolone":
        "fluoroquinolone",
    "efflux":
        "fluoroquinolone",   

    "macrolide":
        "macrolide",

    "lincosamide":
        "lincosamide",

    "streptogramin":
        "streptogramin",

    "pleuromutilin":
        "pleuromutilin",

    "rifamycin":
        "rifamycin",

    "aminocoumarin":
        "aminocoumarin",

    "phenicol":
        "phenicol",

    "fosfomycin":
        "fosfomycin",

    "diaminopyrimidine":
        "diaminopyrimidine",

    "nitrofuran":
        "nitrofuran",

    "glycylcycline":
        "glycylcycline",

    "carbapenem":
        "carbapenem",

    "cephalosporin":
        "cephalosporin",

    "glycopeptide":
        "glycopeptide",

    "peptide":
        "peptide",

    "nucleoside":
        "nucleoside",

    "disinfecting agents and antiseptics":
        "disinfecting agents and antiseptics",

}


# ============================================================
# GENE NORMALIZATION
# ============================================================
def normalize_gene(gene):

    gene = gene.strip().lower()

    # remove bla prefix
    gene = re.sub(r"^bla", "", gene)

    # keep text inside parentheses
    gene = gene.replace("(", "")
    gene = gene.replace(")", "")

    # remove punctuation
    gene = re.sub(r"[-_ ]", "", gene)

    # remove apostrophes
    gene = gene.replace("'", "")

    return gene



# ============================================================
# GENE LOOKUP
# ============================================================

def lookup_gene(gene):

    g = normalize_gene(gene)

    # Match esatto
    if g in GT_GENES:

        return g

    # Match substring
    for truth in GT_GENES:

        if g in truth or truth in g:

            return truth

    return None


# ============================================================
# METRICS
# ============================================================

def precision(tp, fp):

    return (
        tp / (tp + fp)
        if tp + fp
        else 0
    )


def recall(tp, fn):

    return (
        tp / (tp + fn)
        if tp + fn
        else 0
    )


def f1(tp, fp, fn):

    p = precision(tp, fp)
    r = recall(tp, fn)

    return (
        2 * p * r / (p + r)
        if p + r
        else 0
    )


# ============================================================
# EXTRACT AMRFINDER RESULTS
# ============================================================

def extract_amrfinder(tsv):

    predicted_genes = set()
    predicted_classes = set()

    with open(tsv, newline="") as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            # ==================================================
            # Consideriamo solo le righe AMR
            # ==================================================

            if row["Type"].strip().upper() != "AMR":

                continue


            # ==================================================
            # GENE
            # ==================================================

            gene_name = row["Element symbol"].strip()

            if gene_name:

                gene = normalize_gene(gene_name)
                gene = GENE_ALIASES.get(gene, gene)
                predicted_genes.add(gene)


            # ==================================================
            # DRUG CLASS
            # ==================================================

            drug_class = row["Class"].strip().lower()

            if drug_class:

                # Alcune classi possono essere separate da /
                for c in re.split(
                    r"[;,/]",
                    drug_class
                ):

                    c = c.strip()

                    if not c:

                        continue

                    predicted_classes.add(
                        CLASS_MAP.get(c, c)
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

# ============================================================
# EVALUATION
# ============================================================

def evaluate(tsv):

    pred_genes, pred_classes = extract_amrfinder(tsv)


    # ========================================================
    # GENES
    # ========================================================

    TPg = pred_genes & GT_GENES

    FPg = pred_genes - GT_GENES

    FNg = GT_GENES - pred_genes


    # ========================================================
    # DRUG CLASSES
    # ========================================================

    TPc = pred_classes & GT_CLASSES

    FPc = pred_classes - GT_CLASSES

    FNc = GT_CLASSES - pred_classes


    # ========================================================
    # GENES
    # ========================================================

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


    print("\n==============================")
    print(os.path.basename(tsv))
    print("==============================")


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


    # ========================================================
    # DRUG CLASSES
    # ========================================================

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


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 2:

        print(
            "Usage: compare_amrfinder.py <amrfinder.tsv>"
        )

        sys.exit(1)


    evaluate(sys.argv[1])
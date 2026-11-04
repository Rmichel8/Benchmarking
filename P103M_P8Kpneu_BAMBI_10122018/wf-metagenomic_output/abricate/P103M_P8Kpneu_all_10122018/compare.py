#!/usr/bin/env python3
import csv
import re
import os
import sys

# Ground truth
genes = [
    "patA",
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
    "lsaa",
    "tetO",
    "ileS",
    "emeA",
    "dfrE",
    "pc1betalactamaseblaz",
    "tet32",
    "norA",
    "tetM",
    "qacA",
    "aac6-ie-aph2-ia",
    "ErmC"
]
drug_class= [
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
    "nitrofuran"
]



GT_GENES = {
    g.lower()
    .replace("-", "")
    .replace("_", "")
    .replace(" ", "")
    for g in genes
}


GT_CLASSES = {
    c.lower()
    for c in drug_class
}


def normalize_gene(g):

    g = g.lower().strip()

    g = re.sub(r"^[a-z]+_[a-z]+_", "", g)

    g = g.replace("-", "")
    g = g.replace("_", "")
    g = g.replace(" ", "")

    g = re.sub(r"[^a-z0-9]", "", g)

    return g

# Mapping classi ABRicate -> truth

CLASS_MAP = {
    "aminoglycoside": "aminoglycoside",

    "fluoroquinolone": "fluoroquinolone",

    "tetracycline": "tetracycline",

    "glycylcycline": "glycylcycline",

    "phenicol": "phenicol",

    "macrolide": "macrolide",

    "lincosamide": "lincosamide",

    "streptogramin": "streptogramin",
    "cephamycin":"cephalosporin",
    
    "triclosan": "disinfecting agents and antiseptics",

    "pleuromutilin": "pleuromutilin",

    "rifamycin": "rifamycin",

    "fosfomycin": "fosfomycin",

    "diaminopyrimidine": "diaminopyrimidine",

    "nitrofuran": "nitrofuran",

    "aminocoumarin": "aminocoumarin",

    "penam": "penicillin beta-lactam",

    "penem": "carbapenem",

    "penicillin": "penicillin beta-lactam",

    "cephalosporin": "cephalosporin",

    "cephamycin": "cephalosporin",

    "carbapenem": "carbapenem",

}

def precision(tp, fp):

    return tp/(tp+fp) if tp+fp else 0



def recall(tp, fn):

    return tp/(tp+fn) if tp+fn else 0



def f1(tp, fp, fn):

    p = precision(tp, fp)
    r = recall(tp, fn)

    return 2*p*r/(p+r) if p+r else 0




def extract_predictions(tsv):

    pred_genes = set()
    pred_classes = set()


    with open(tsv) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )


        for row in reader:
            # genes

            gene = normalize_gene(
                row["GENE"]
            )

            pred_genes.add(gene)

            # classes

            resistance = row["RESISTANCE"].lower()


            for c in resistance.split(";"):

                c = c.strip()


                mapped = CLASS_MAP.get(c)


                if mapped:

                    pred_classes.add(mapped)

                else:

                    #FP
                    pred_classes.add(c)


    return pred_genes, pred_classes


def evaluate(tsv):


    pred_genes, pred_classes = extract_predictions(tsv)

    # genes


    TPg = pred_genes & GT_GENES

    FPg = pred_genes - GT_GENES

    FNg = GT_GENES - pred_genes


    # classes

    TPc = pred_classes & GT_CLASSES

    FPc = pred_classes - GT_CLASSES

    FNc = GT_CLASSES - pred_classes



    print("\n========== GENE TP ==========")

    for x in sorted(TPg):
        print(x)



    print("\n========== GENE FP ==========")

    for x in sorted(FPg):
        print(x)



    print("\n========== CLASS TP ==========")

    for x in sorted(TPc):
        print(x)



    print("\n========== CLASS FP ==========")

    for x in sorted(FPc):
        print(x)



    print("\n========== METRICS ==========")


    print(
        "GENES:",
        "precision=",round(precision(len(TPg),len(FPg)),3),
        "recall=",round(recall(len(TPg),len(FNg)),3),
        "f1=",round(f1(len(TPg),len(FPg),len(FNg)),3)
    )


    print(
        "CLASSES:",
        "precision=",round(precision(len(TPc),len(FPc)),3),
        "recall=",round(recall(len(TPc),len(FNc)),3),
        "f1=",round(f1(len(TPc),len(FPc),len(FNc)),3)
    )



    out = os.path.join(
        os.path.dirname(__file__),
        "summary_metrics_abricate.csv"
    )


    with open(out,"w",newline="") as f:

        w = csv.writer(f)

        w.writerow([
            "type",
            "TP",
            "FP",
            "FN",
            "precision",
            "recall",
            "f1"
        ])


        w.writerow([
            "genes",
            len(TPg),
            len(FPg),
            len(FNg),
            round(precision(len(TPg),len(FPg)),3),
            round(recall(len(TPg),len(FNg)),3),
            round(f1(len(TPg),len(FPg),len(FNg)),3)
        ])


        w.writerow([
            "drug_classes",
            len(TPc),
            len(FPc),
            len(FNc),
            round(precision(len(TPc),len(FPc)),3),
            round(recall(len(TPc),len(FNc)),3),
            round(f1(len(TPc),len(FPc),len(FNc)),3)
        ])


    print("\nResults written to:", out)

if __name__ == "__main__":


    if len(sys.argv)!=2:

        print(
            "Usage: python compare.py abricate_merged.tsv"
        )

        sys.exit(1)


    evaluate(sys.argv[1])
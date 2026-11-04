#!/usr/bin/env python3

import csv
import math
import re
from pathlib import Path

##############################################################################
# INPUT / OUTPUT
##############################################################################

ROOT = Path(
    "/nfsd/sysbiobig/Riccardom/"
    "P103M_P8Kpneu_BAMBI_10122018/"
    "benchmarking/"
    "bacterial_genomes_resfinder"
)

INPUT = ROOT / "barcode01_resfinder_results" / "ResFinder_results_tab.txt"

OUT_RESULTS = ROOT / "per_sample_results.tsv"

OUT_SUMMARY = ROOT / "summary_metrics.csv"

##############################################################################
# GROUND TRUTH
##############################################################################

GT_GENES = {
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
    "PC1",
    "tet32",
    "norA",
    "tetM",
    "qacA",
    "aac6-ie-aph2-ia",
    "ErmC"
}

GT_CLASSES = {

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

}

##############################################################################
# GENE NORMALIZATION
##############################################################################

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


##############################################################################
# NORMALIZED GROUND TRUTH
##############################################################################

GT_GENES = {

    normalize_gene(g)

    for g in GT_GENES

}

print()

print("Ground truth genes :", len(GT_GENES))

print("Ground truth drug classes :", len(GT_CLASSES))
##############################################################################
# ANTIBIOTIC -> DRUG CLASS
##############################################################################

ABX_TO_CLASS = {

    # ----------------------------------------------------------
    # Aminoglycosides
    # ----------------------------------------------------------

    "amikacin": "aminoglycoside",
    "gentamicin": "aminoglycoside",
    "kanamycin": "aminoglycoside",
    "tobramycin": "aminoglycoside",
    "streptomycin": "aminoglycoside",
    "netilmicin": "aminoglycoside",
    "dibekacin": "aminoglycoside",
    "fortimicin": "aminoglycoside",
    "isepamicin": "aminoglycoside",

    # ----------------------------------------------------------
    # Penicillins
    # ----------------------------------------------------------

    "penicillin": "penicillin beta-lactam",
    "amoxicillin": "penicillin beta-lactam",
    "ampicillin": "penicillin beta-lactam",
    "piperacillin": "penicillin beta-lactam",

    "amoxicillin+clavulanic acid": "penicillin beta-lactam",
    "ampicillin+clavulanic acid": "penicillin beta-lactam",
    "piperacillin+tazobactam": "penicillin beta-lactam",

    # ----------------------------------------------------------
    # Cephalosporins
    # ----------------------------------------------------------

    "cefepime": "cephalosporin",
    "ceftazidime": "cephalosporin",
    "cefotaxime": "cephalosporin",
    "cefixime": "cephalosporin",
    "cefoxitin": "cephalosporin",

    # ----------------------------------------------------------
    # Carbapenems
    # ----------------------------------------------------------

    "imipenem": "carbapenem",
    "meropenem": "carbapenem",
    "ertapenem": "carbapenem",

    # ----------------------------------------------------------
    # Fluoroquinolones
    # ----------------------------------------------------------

    "ciprofloxacin": "fluoroquinolone",
    "nalidixic acid": "fluoroquinolone",

    # ----------------------------------------------------------
    # Tetracyclines
    # ----------------------------------------------------------

    "tetracycline": "tetracycline",
    "doxycycline": "tetracycline",
    "minocycline": "tetracycline",

    # ----------------------------------------------------------
    # Glycylcyclines
    # ----------------------------------------------------------

    "tigecycline": "glycylcycline",

    # ----------------------------------------------------------
    # Phenicols
    # ----------------------------------------------------------

    "chloramphenicol": "phenicol",

    # ----------------------------------------------------------
    # Macrolides
    # ----------------------------------------------------------

    "erythromycin": "macrolide",
    "azithromycin": "macrolide",
    "spiramycin": "macrolide",
    "telithromycin": "macrolide",

    # ----------------------------------------------------------
    # Lincosamides
    # ----------------------------------------------------------

    "clindamycin": "lincosamide",
    "lincomycin": "lincosamide",

    # ----------------------------------------------------------
    # Streptogramins
    # ----------------------------------------------------------

    "quinupristin": "streptogramin",
    "pristinamycin ia": "streptogramin",
    "pristinamycin iia": "streptogramin",
    "virginiamycin s": "streptogramin",
    "virginiamycin m": "streptogramin",
    "dalfopristin": "streptogramin",

    # ----------------------------------------------------------
    # Pleuromutilins
    # ----------------------------------------------------------

    "tiamulin": "pleuromutilin",

    # ----------------------------------------------------------
    # Rifamycins
    # ----------------------------------------------------------

    "rifampicin": "rifamycin",
    "rifampin": "rifamycin",

    # ----------------------------------------------------------
    # Fosfomycin
    # ----------------------------------------------------------

    "fosfomycin": "fosfomycin",

    # ----------------------------------------------------------
    # Diaminopyrimidines
    # ----------------------------------------------------------

    "trimethoprim": "diaminopyrimidine",

    # ----------------------------------------------------------
    # Nitrofurans
    # ----------------------------------------------------------

    "nitrofurantoin": "nitrofuran"

}

##############################################################################
# PARSER
##############################################################################

def parse_resfinder(path):

    predicted_genes = set()

    predicted_classes = set()

    with open(path) as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            ####################################################
            # Gene
            ####################################################

            gene = normalize_gene(
                row["Resistance gene"]
            )

            predicted_genes.add(gene)

            ####################################################
            # Drug class
            ####################################################

            phenotype = row["Phenotype"]

            for ab in phenotype.split(","):

                ab = ab.strip().lower()

                if ab in ABX_TO_CLASS:

                    predicted_classes.add(
                        ABX_TO_CLASS[ab]
                    )

    return predicted_genes, predicted_classes


##############################################################################
# LOAD PREDICTIONS
##############################################################################

predicted_genes, predicted_classes = parse_resfinder(INPUT)

print()

print("Predicted genes:", len(predicted_genes))

print("Predicted drug classes:", len(predicted_classes))

##############################################################################
# GENE MATCHING
##############################################################################

def gene_matches_truth(pred_gene, truth_gene):

    #
    # exact match
    #

    if pred_gene == truth_gene:

        return True

    return False


##############################################################################
# GENE BENCHMARK
##############################################################################

matched_truth = set()

TP_gene = 0
FP_gene = 0

for pred in predicted_genes:

    match = False

    for gt in GT_GENES:

        if gene_matches_truth(pred, gt):

            TP_gene += 1

            matched_truth.add(gt)
             
            print(f"TP: {pred}  ->  {gt}")

            match = True

            break

    if not match:

        FP_gene += 1


FN_gene = len(GT_GENES - matched_truth)


##############################################################################
# DRUG CLASS BENCHMARK
##############################################################################

TP_class = len(predicted_classes & GT_CLASSES)

FP_class = len(predicted_classes - GT_CLASSES)

FN_class = len(GT_CLASSES - predicted_classes)


##############################################################################
# METRICS
##############################################################################

def calculate_metrics(tp, fp, fn):

    precision = (
        tp / (tp + fp)
        if tp + fp else 0
    )

    recall = (
        tp / (tp + fn)
        if tp + fn else 0
    )

    f1 = (
        2 * precision * recall /
        (precision + recall)
        if precision + recall else 0
    )

    return {

        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1

    }


gene_metrics = calculate_metrics(
    TP_gene,
    FP_gene,
    FN_gene
)

class_metrics = calculate_metrics(
    TP_class,
    FP_class,
    FN_class
)


##############################################################################
# PRINT
##############################################################################

print("\n==============================")
print("GENE BENCHMARK")
print("==============================")

for k, v in gene_metrics.items():

    print(f"{k:10s} {v}")


print("\n==============================")
print("DRUG CLASS BENCHMARK")
print("==============================")

for k, v in class_metrics.items():

    print(f"{k:10s} {v}")
##############################################################################
# SAVE RESULTS
##############################################################################

##############################################################################
# SUMMARY CSV
##############################################################################

summary_file = OUT_SUMMARY


with open(summary_file, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "level",
        "TP",
        "FP",
        "FN",
        "precision",
        "recall",
        "f1"
    ])


    writer.writerow([

        "genes",
        gene_metrics["TP"],
        gene_metrics["FP"],
        gene_metrics["FN"],
        round(gene_metrics["precision"],4),
        round(gene_metrics["recall"],4),
        round(gene_metrics["f1"],4)

    ])


    writer.writerow([

        "drug_classes",
        class_metrics["TP"],
        class_metrics["FP"],
        class_metrics["FN"],
        round(class_metrics["precision"],4),
        round(class_metrics["recall"],4),
        round(class_metrics["f1"],4)

    ])


print()
print(predicted_classes & GT_CLASSES)
print("Saved summary:")
print(summary_file)



##############################################################################
# DETAILED TSV
##############################################################################

with open(OUT_RESULTS, "w", newline="") as f:

    writer = csv.writer(
        f,
        delimiter="\t"
    )

    writer.writerow([

        "type",
        "predicted"

    ])


    for g in sorted(predicted_genes):

        writer.writerow([

            "gene",
            g

        ])


    for c in sorted(predicted_classes):

        writer.writerow([

            "drug_class",
            c

        ])


print()

print("Saved detailed results:")
print(OUT_RESULTS)
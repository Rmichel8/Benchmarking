#!/usr/bin/env python3

import os
import sys
import csv
import json
import re
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

CARD_JSON = "/nfsd/sysbiobig/Riccardom/database_CARD/card-data/card.json"

#AMRFINDER
#SUMMARY_CSV = ( "/nfsd/sysbiobig/Riccardom/" "RUARGPOREPRJNA821338/benchmarking/" "wfmetagenomics_amrfinder/" "amrfinder_compare_results/""summary_CARD_compare.csv")
#RGI

SUMMARY_CSV = (
"/nfsd/sysbiobig/Riccardom/"
"RUARGPOREPRJNA821338/benchmarking/"
"wf-metagenomic_rgi/"
"rgi_compare_results/"
"summary_CARD_compare.csv"
)


# ============================================================
# NORMALIZATION SYNONIMS CARD -> Ground truth
# ============================================================


CARD_TO_GT = {


    # beta-lattamici
    "amoxicillin":"Amoxicillin/clavulanic acid",
    "ampicillin":"Amoxicillin/clavulanic acid",
    "benzylpenicillin":"Amoxicillin/clavulanic acid",
    "cefepime":"Cefepime",
    "ceftriaxone":"Ceftriaxone",
    "cefuroxime":"Cefuroxime",
    "ceftazidime":"Ceftazidime",
    "aztreonam":"Aztreonam",
    "imipenem":"Imipenem",
    "meropenem":"Meropenem",
    "ertapenem":"Ertapenem",


    # quinoloni

    "ciprofloxacin":
        "Ciprofloxacin",

    "levofloxacin":
        "Levofloxacin",

    "norfloxacin":
        "Norfloxacin",


    # aminoglicosidi

    "amikacin":
        "Amikacin",

    "gentamicin":
        "Gentamicin",

    "tobramycin":
        "Tobramycin",


    # tetracicline

    "tetracycline":
        "Tetracycline",

    "minocycline":
        "Minocycline",

    "tigecycline":
        "Tigecycline",


    # altri

    "chloramphenicol":
        "Chloramphenicol",

    "rifampin":
        "Rifampicin",

    "rifampicin":
        "Rifampicin",

    "vancomycin":
        "Vancomycin",

    "teicoplanin":
        "Teicoplanin",

    "colistin A":
        "Colistin",

    "colistin B":
        "Colistin",

    "polymyxin B":
        "Colistin",

    "nitrofurantoin":
        "Nitrofurantoin",

    "fosfomycin":
        "Fosfomycin",

    "linezolid":
        "Linezolid",

    "daptomycin":
        "Daptomycin",

}



with open(CARD_JSON) as f:
    CARD = json.load(f)


GENE_TO_ABX = {}
ARO_TO_ABX = {}
CARD_GENES = []


def normalize_gene(g):

    g = g.lower()

    g = re.sub(r"[^a-z0-9]", "", g)

    return g


for obj in CARD.values():

    if "ARO_name" not in obj:
        continue

    aro = obj.get("ARO_accession")

    gene = obj["ARO_name"]

    antibiotics = set()

    if "ARO_category" not in obj:
        continue

    for cat in obj["ARO_category"].values():

        if cat["category_aro_class_name"] != "Antibiotic":
            continue

        ab = cat["category_aro_name"].lower()

        if ab in CARD_TO_GT:

            antibiotics.add(CARD_TO_GT[ab])

    if not antibiotics:
        continue


    gene_norm = normalize_gene(gene)

    variants = set()

    variants.add(gene_norm)
    variants.add(
        re.sub(r"[0-9]", "", gene_norm)
    )
    m = re.match(r"[a-z]+", gene_norm)

    if m:
        variants.add(m.group())
    if gene_norm.startswith("bla"):

        variants.add(
            gene_norm[3:]
        )


    for g in variants:

        if not g:
            continue

        GENE_TO_ABX[g] = antibiotics

        CARD_GENES.append(
            (g, antibiotics)
        )

    if aro:

        ARO_TO_ABX[aro] = antibiotics


print()

print("CARD entries :", len(CARD))

print("Indexed genes:", len(GENE_TO_ABX))

print("Indexed AROs :", len(ARO_TO_ABX))

def lookup_gene(gene):

    gene = normalize_gene(gene)


    # 1) exact match


    if gene in GENE_TO_ABX:
        return GENE_TO_ABX[gene]

    # 2) substring


    for card_gene, abx in CARD_GENES:

        if gene in card_gene:

            return abx

        if card_gene in gene:

            return abx


    # 3) same family

    m = re.match(r"[a-z]+", gene)

    if m:

        prefix = m.group()

        for card_gene, abx in CARD_GENES:

            if card_gene.startswith(prefix):

                return abx

    # 4) fuzzy matching


    matches = get_close_matches(

        gene,

        GENE_TO_ABX.keys(),

        n=1,

        cutoff=0.90

    )

    if matches:

        return GENE_TO_ABX[matches[0]]

    return set()


def get_benchmark(path):

    parts=os.path.normpath(path).split(os.sep)

    for p in parts:

        if p.startswith("B") and p[1:].isdigit():

            return p

    return None

# DETECTION INPUT TYPE

def detect_file_type(tsv):

    with open(tsv) as f:

        header=f.readline()

    if "Element symbol" in header:
        return "AMRFINDER"

    if "ARO" in header:
        return "RGI"

    return None

# PARSING AMRFINDER

def extract_amrfinder(tsv):

    predicted = set()

    found_genes = set()

    unmapped_genes = set()

    with open(tsv) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            gene = row["Element symbol"]

            gene = normalize_gene(gene)

            found_genes.add(gene)

            antibiotics = lookup_gene(gene)

            if antibiotics:

                predicted.update(antibiotics)

            else:

                unmapped_genes.add(gene)

    return predicted, found_genes, unmapped_genes

# PARSING RGI

from difflib import get_close_matches, SequenceMatcher
def extract_rgi(tsv):

    predicted = set()
    genes = set()
    unmapped = set()

    with open(tsv) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            cols = line.split("\t")

            first = cols[0]

            # ARO accession
            aro = first.split("|")[0].replace("ARO:", "")

            if aro in ARO_TO_ABX:
                predicted.update(ARO_TO_ABX[aro])

            # Gene
            gene = None

            for field in first.split("|"):

                if field.startswith("Name:"):
                    gene = field.replace("Name:", "")
                    break

            if gene:

                gene = normalize_gene(gene)

                genes.add(gene)

                antibiotics = lookup_gene(gene)

                if antibiotics:
                    predicted.update(antibiotics)
                else:
                    unmapped.add(gene)

    return predicted, genes, unmapped

# EVALUATION

def evaluate(tsv):

    bench = get_benchmark(tsv)

    if bench not in GT:

        print("Cannot identify benchmark:", tsv)
        return

    tp_truth, tn_truth = GT[bench]

    tool = detect_file_type(tsv)

    if tool == "AMRFINDER":

        predicted, genes, unmapped = extract_amrfinder(tsv)

    elif tool == "RGI":

        predicted, genes, unmapped = extract_rgi(tsv)

    else:

        print("Unknown file type")
        return


    TP = predicted & tp_truth

    FP = predicted & tn_truth

    FN = tp_truth - predicted

    TN = tn_truth - predicted

    sample = os.path.basename(tsv)

    print()
    print("=" * 70)
    print(sample)
    print("=" * 70)

    print("Benchmark :", bench)
    print("Tool      :", tool)

    print()

    print("Predicted antibiotics")
    print("---------------------")

    for a in sorted(predicted):

        print("  ", a)

    print()

    print("True positives")
    print("---------------------")

    for a in sorted(TP):

        print("  ", a)

    print()

    print("False positives")
    print("---------------------")

    for a in sorted(FP):

        print("  ", a)

    print()

    print("False negatives")
    print("---------------------")

    for a in sorted(FN):

        print("  ", a)

    print()

    print("True negatives")
    print("---------------------")

    for a in sorted(TN):

        print("  ", a)


    gene_file = tsv.replace(
        ".tsv",
        "_predicted_genes.txt"
    )

    with open(gene_file, "w") as out:

        for g in sorted(genes):

            out.write(g + "\n")

    unmapped_file = tsv.replace(
        ".tsv",
        "_unmapped_genes.txt"
    )

    with open(unmapped_file, "w") as out:

        for g in sorted(unmapped):

            out.write(g + "\n")

    # summary csv


    os.makedirs(
        os.path.dirname(SUMMARY_CSV),
        exist_ok=True
    )

    write_header = not os.path.exists(
        SUMMARY_CSV
    )

    with open(
        SUMMARY_CSV,
        "a",
        newline=""
    ) as out:

        writer = csv.writer(out)

        if write_header:

            writer.writerow([
                "sample",
                "benchmark",
                "tool",
                "TP",
                "FP",
                "TN",
                "FN",
                "predicted_antibiotics",
                "tp_antibiotics",
                "fp_antibiotics",
                "fn_antibiotics",
                "genes_found",
                "genes_unmapped"
            ])


        writer.writerow([

            sample,

            bench,

            tool,

            len(TP),

            len(FP),

            len(TN),

            len(FN),

            ";".join(sorted(predicted)),

            ";".join(sorted(TP)),

            ";".join(sorted(FP)),

            ";".join(sorted(FN)),

            ";".join(sorted(genes)),

            ";".join(sorted(unmapped))

        ])


    print()
    print("Summary updated:")
    print(SUMMARY_CSV)

    print("Mapped genes   :", len(genes) - len(unmapped))
    print("Unmapped genes :", len(unmapped))
if __name__=="__main__":


    if len(sys.argv)!=2:

        print(
        "Usage: compare.py <tsv>"
        )

        sys.exit(1)


    evaluate(sys.argv[1])
#!/usr/bin/env python3
import pandas as pd

# Ground truth AMR genes

ground_truth = sorted({
    "mexI",
    "mexB",
    "mexF",
    "mexD",
    "mexQ",
    "ceoB",
    "amrB",
    "carA",
    "OXA-114a",
    "sul2",
    "srmB",
    "mdtF",
    "mdtC",
    "mdtD",
    "mdtB",
    "acrB",
    "acrA",
    "acrD",
    "mfd",
    "PBP1a",
    "tetH",
    "cpxR",
    "CRP",
    "emrD",
    "emrB",
    "H-NS",
    "arnA",
    "bacA",
    "smeE",
    "vgaC",
    "FosA5",
    "SHV-152",
    "SHV-187",
    "SHV-185",
    "SHV-16",
    "SHV-43",
    "oqxB",
    "oleB"
})

# Input MARTi TSV


marti_file = "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/MARTi/MARTi_AMR_results.tsv"

output_prefix = "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/MARTi/MARTi_comparison"
import re
ALIASES = {
    "klebsiellapneumoniaeacra": "acra",
    "pseudomonasaeruginosacpxr": "cpxr"
}
#Normalize genes
def normalize_gene_name(gene):
    """
    Normalizza i nomi AMR per confrontarli indipendentemente
    da maiuscole/minuscole e descrizioni aggiuntive.
    """

    gene = str(gene).strip()

    gene = gene.lower()

    # rimove spaces and other symbols
    gene = re.sub(r"[^a-z0-9\-]", "", gene)
    gene = ALIASES.get(gene, gene)
    return gene

# Load MARTi results


df = pd.read_csv(
    marti_file,
    sep="\t"
)


# AMR names founded by MARTi
found_genes_raw = set(
    df["name"]
    .dropna()
    .astype(str)
)

found_genes = {
    normalize_gene_name(x)
    for x in found_genes_raw
}


truth = {
    normalize_gene_name(x)
    for x in ground_truth
}

# Comparison


TP = sorted(found_genes & truth)
FP = sorted(found_genes - truth)
FN = sorted(truth - found_genes)


precision = len(TP) / (len(TP) + len(FP)) if (len(TP)+len(FP)) > 0 else 0
recall = len(TP) / (len(TP) + len(FN)) if (len(TP)+len(FN)) > 0 else 0

f1 = (
    2 * precision * recall / (precision + recall)
    if (precision + recall) > 0
    else 0
)


print("\n===== MARTi AMR benchmark =====")

print(f"Ground truth genes : {len(truth)}")
print(f"MARTi genes found  : {len(found_genes)}")

print("\nTP:", len(TP))
print("FP:", len(FP))
print("FN:", len(FN))

print("\nPrecision :", round(precision, 4))
print("Recall    :", round(recall, 4))
print("F1-score  :", round(f1, 4))


# Save results


with open(output_prefix + "_TP.txt", "w") as f:
    f.write("\n".join(TP))

with open(output_prefix + "_FP.txt", "w") as f:
    f.write("\n".join(FP))

with open(output_prefix + "_FN.txt", "w") as f:
    f.write("\n".join(FN))


# Save complete summary

with open(output_prefix + "_summary.txt", "w") as f:
    f.write("MARTi AMR benchmark\n\n")
    f.write(f"Ground truth: {len(truth)}\n")
    f.write(f"Detected: {len(found_genes)}\n\n")
    f.write(f"TP: {len(TP)}\n")
    f.write(f"FP: {len(FP)}\n")
    f.write(f"FN: {len(FN)}\n\n")
    f.write(f"Precision: {precision:.4f}\n")
    f.write(f"Recall: {recall:.4f}\n")
    f.write(f"F1: {f1:.4f}\n")

print("\nResults saved:")
print(output_prefix + "_TP.txt")
print(output_prefix + "_FP.txt")
print(output_prefix + "_FN.txt")
print(output_prefix + "_summary.txt")
mapping = {}

for g in found_genes_raw:
    ng = normalize_gene_name(g)
    mapping.setdefault(ng, []).append(g)


for g in sorted(TP):
    print(
        g,
        "->",
        mapping.get(g, [])
    )
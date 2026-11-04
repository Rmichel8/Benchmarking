import re

# Ground truth

genes = sorted({
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
    "oleB",
})
drug_classes = [
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
    "elfamycin"]
rgi_file = "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/wf-metagenomic_output/rgi/ZymoEvenMock/rgi_genes.txt"

output_file = "/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/wf-metagenomic_output/rgi/ZymoEvenMock/summary_metrics.tsv"


# Normalization


def normalize_gene(gene):

    gene = gene.strip().lower()

    # rimuove eventuale prefisso bla
    gene = re.sub(r"^bla", "", gene)

    # elimina caratteri speciali
    gene = re.sub(r"[^a-z0-9]", "", gene)

    return gene

#known gene normalization
known_genes = {normalize_gene(g) for g in genes}



def extract_gene_name(gene):

    gene = normalize_gene(gene)

    # remove CARD description
    gene = re.sub(r"conferring.*$", "", gene)
    gene = re.sub(r"withmutation.*$", "", gene)
    gene = re.sub(r"mutation.*$", "", gene)
    gene = re.sub(r"mutants.*$", "", gene)
    gene = re.sub(r"genein.*$", "", gene)

    # search a ground tuth gene at the end of the string 
    for g in sorted(known_genes, key=len, reverse=True):
        if gene.endswith(g):
            return g

    return gene



rgi_genes = set()

with open(rgi_file) as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        rgi_genes.add(extract_gene_name(line))

truth_genes = {normalize_gene(g) for g in genes}
print("\n===== Ground truth =====")
for g in genes:
    print(f"{g:10s} -> {normalize_gene(g)}")

print("\nGround truth normalizzata:")
for g in sorted(truth_genes):
    print(" ", g)

TP = truth_genes & rgi_genes
FP = rgi_genes - truth_genes
FN = truth_genes - rgi_genes


# Metrics

precision = len(TP) / (len(TP) + len(FP)) if TP or FP else 0

recall = len(TP) / (len(TP) + len(FN)) if TP or FN else 0

f1 = (
    2 * precision * recall / (precision + recall)
    if precision + recall
    else 0
)



print("\n===== AMR Benchmark =====\n")

print(f"Ground truth genes : {len(truth_genes)}")
print(f"RGI detected genes : {len(rgi_genes)}\n")

print("TP")
for g in sorted(TP):
    print(" ", g)

print("\nFP")
for g in sorted(FP):
    print(" ", g)

print("\nFN")
for g in sorted(FN):
    print(" ", g)

print("\n===== Metrics =====")

print(f"TP = {len(TP)}")
print(f"FP = {len(FP)}")
print(f"FN = {len(FN)}")

print(f"\nPrecision = {precision:.3f}")
print(f"Recall    = {recall:.3f}")
print(f"F1-score  = {f1:.3f}")


# Save TSV

with open(output_file, "w") as out:

    out.write("Metric\tValue\n")

    out.write(f"Ground_truth_genes\t{len(truth_genes)}\n")
    out.write(f"RGI_detected_genes\t{len(rgi_genes)}\n")

    out.write(f"TP\t{len(TP)}\n")
    out.write(f"FP\t{len(FP)}\n")
    out.write(f"FN\t{len(FN)}\n")

    out.write(f"Precision\t{precision:.4f}\n")
    out.write(f"Recall\t{recall:.4f}\n")
    out.write(f"F1_score\t{f1:.4f}\n")

    out.write("\nCategory\tGene\n")

    for g in sorted(TP):
        out.write(f"TP\t{g}\n")

    for g in sorted(FP):
        out.write(f"FP\t{g}\n")

    for g in sorted(FN):
        out.write(f"FN\t{g}\n")

print(f"\nResults saved in: {output_file}")
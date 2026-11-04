import csv
import re
from collections import defaultdict
# Files
rgi_file = "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018/benchmarking/wf-metagenomic_output/rgi/P103M_P8Kpneu_all_10122018/amr/tsv/P103M_P8Kpneu_all_10122018_amr_results_filtered.tsv"
card_file = "/nfsd/sysbiobig/Riccardom/database_CARD/card-data/aro_index.tsv"
output = "/nfsd/sysbiobig/Riccardom/P103M_P8Kpneu_BAMBI_10122018/benchmarking/wf-metagenomic_output/rgi/P103M_P8Kpneu_all_10122018/summary_classes.tsv"
# Ground truth drug classes
truth_classes = {
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
    "nitrofuran"
}

# Normalization

def normalize_name(x):

    x = x.lower().strip()

    x = re.sub(
        r"[^a-z0-9]",
        "",
        x
    )

    return x



def normalize_class(x):

    x = x.lower().strip()


    mapping = {

    # beta-lactam
    "beta-lactam antibiotic":
        "penicillin beta-lactam",

    "penicillin antibiotic":
        "penicillin beta-lactam",

    "cephalosporin antibiotic":
        "cephalosporin",

    "carbapenem antibiotic":
        "carbapenem",


    # tetracycline
    "tetracycline antibiotic":
        "tetracycline",

    "glycylcycline antibiotic":
        "glycylcycline",


    # quinolone
    "fluoroquinolone antibiotic":
        "fluoroquinolone",


    # folate pathway
    "diaminopyrimidine antibiotic":
        "diaminopyrimidine",

    "sulfonamide antibiotic":
        "sulfonamide",


    # protein synthesis
    "aminoglycoside antibiotic":
        "aminoglycoside",

    "macrolide antibiotic":
        "macrolide",

    "lincosamide antibiotic":
        "lincosamide",

    "streptogramin antibiotic":
        "streptogramin",

    "phenicol antibiotic":
        "phenicol",

    "pleuromutilin antibiotic":
        "pleuromutilin",


    # other
    "aminocoumarin antibiotic":
        "aminocoumarin",

    "rifamycin antibiotic":
        "rifamycin",

    "fosfomycin antibiotic":
        "fosfomycin",

    "nitrofuran antibiotic":
        "nitrofuran",


    # eventuali CARD names
    "elfamycin antibiotic":
        "rifamycin",

    "peptide antibiotic":
        "peptide_antibiotics",

    "polymyxin antibiotic":
        "polymyxin"}

    return mapping.get(x)

# Load CARD

print("Loading CARD...")


aro_to_class = defaultdict(set)

name_to_class = defaultdict(set)



with open(card_file) as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )


    for row in reader:


        drug = row["Drug Class"].strip()


        if not drug:
            continue


        classes = set()


        for d in drug.split(";"):

            c = normalize_class(d)

            if c:
                classes.add(c)



        if not classes:
            continue



        # ARO mapping

        aro = row["ARO Accession"].strip()

        aro_to_class[aro].update(classes)



        # Name mapping

        for field in [
            "Model Name",
            "ARO Name",
            "AMR Gene Family"
        ]:


            value = row[field].strip()


            if value:

                value = normalize_name(value)


                name_to_class[value].update(
                    classes
                )



print(
    "CARD ARO loaded:",
    len(aro_to_class)
)


print(
    "CARD names loaded:",
    len(name_to_class)
)

# Parse RGI

detected_classes = set()

detected_aros = set()

missing = []

aro_report = []


with open(rgi_file) as f:


    for line in f:


        if not line.strip():
            continue

        # ARO

        aro = line.split("|")[0].strip()

        detected_aros.add(aro)

        # NAME

        match = re.search(
            r"Name:([^|]+)",
            line
        )

        gene_name = None

        if match:

            gene_name = match.group(1)

        found = False

        # --------------------
        # 1) ARO lookup
        # --------------------
        if aro in aro_to_class:

            detected_classes.update(
                aro_to_class[aro]
                )

            for c in aro_to_class[aro]:
                aro_report.append(
                (aro, c, "ARO")
                )

            found = True

      

        if not found and gene_name:


            name_norm = normalize_name(
                gene_name
            )


            if name_norm in name_to_class:
                detected_classes.update(
                name_to_class[name_norm]
                )

                for c in name_to_class[name_norm]:
                    aro_report.append(
                    (aro, c, "NAME")
                    )

                found = True


        if not found and gene_name:


            name_norm = normalize_name(
                gene_name
            )


            for card_name, classes in name_to_class.items():


                if (
                    name_norm in card_name
                    or card_name in name_norm
                ):

                    detected_classes.update(classes)

                    for c in classes:
                        aro_report.append(
                        (aro, c, "PARTIAL")
                        )
                    found=True
                    break


        if not found:


            missing.append(
                (
                    aro,
                    gene_name
                )
            )




print("\nDetected ARO:")
print(len(detected_aros))


print("\nDrug classes detected:")

for x in sorted(detected_classes):
    print(x)



print("\n======================")
print("Missing CARD mapping")
print("======================")


for aro,name in missing:

    print(
        aro,
        name
    )



TP = truth_classes & detected_classes

FP = detected_classes - truth_classes

FN = truth_classes - detected_classes



precision = (
    len(TP)/(len(TP)+len(FP))
    if len(TP)+len(FP)
    else 0
)


recall = (
    len(TP)/(len(TP)+len(FN))
    if len(TP)+len(FN)
    else 0
)


f1 = (
    2*precision*recall/(precision+recall)
    if precision+recall
    else 0
)




print("\n======================")
print("Drug class benchmark")
print("======================")


print("\nTP:")

for x in sorted(TP):
    print(x)



print("\nFP:")

for x in sorted(FP):
    print(x)



print("\nFN:")

for x in sorted(FN):
    print(x)



print("\nMetrics")
print("----------------")

print("TP:",len(TP))
print("FP:",len(FP))
print("FN:",len(FN))

print(
    "\nPrecision:",
    round(precision,3)
)

print(
    "Recall:",
    round(recall,3)
)

print(
    "F1:",
    round(f1,3)
)

with open(output,"w") as out:


    out.write(
        "Section\tARO_or_Class\tDrug_class\tSource\n"
    )


    for aro,c,s in aro_report:

        out.write(
            f"ARO\t{aro}\t{c}\t{s}\n"
        )


    out.write("\n")


    for x in sorted(TP):

        out.write(
            f"TP\t{x}\t\t\n"
        )


    for x in sorted(FP):

        out.write(
            f"FP\t{x}\t\t\n"
        )


    for x in sorted(FN):

        out.write(
            f"FN\t{x}\t\t\n"
        )


    out.write("\n")

    out.write(
        f"Precision\t{precision}\n"
    )

    out.write(
        f"Recall\t{recall}\n"
    )

    out.write(
        f"F1\t{f1}\n"
    )


print(
    "\nSaved:",
    output
)
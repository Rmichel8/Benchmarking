#!/usr/bin/env python3

import csv
import os
import re
from collections import defaultdict

############################################################
# INPUT FILES
############################################################

ASSIGNMENT_FILE = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_rgi/"
    "amr_fastq_assignment/"
    "rgi_fastq_assignment_final.tsv"
)

RGI_FILE = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_rgi/"
    "rgi_results/"
    "RUARGPOREPRJNA821338_rgi.txt"
)

############################################################
# OUTPUT DIRECTORY
############################################################

OUTDIR = os.path.join(
    os.path.dirname(ASSIGNMENT_FILE),
    "split_by_fastq_annotated"
)

os.makedirs(OUTDIR, exist_ok=True)


############################################################
# READ RGI RESULTS
############################################################

def read_rgi():

    contig_hits = defaultdict(list)

    with open(
        RGI_FILE,
        newline="",
        encoding="utf-8"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            # ORF_ID:
            #
            # contig_103 # 32767 # 33870 # +
            #
            contig = (
                row["ORF_ID"]
                .split("#")[0]
                .strip()
            )

            contig_hits[contig].append({

                "RGI_ORF_ID":
                    row["ORF_ID"],

                "RGI_BEST_HIT":
                    row["Best_Hit_ARO"],

                "RGI_ARO":
                    row["ARO"],

                "RGI_DRUG_CLASS":
                    row["Drug Class"],

                "RGI_RESISTANCE_MECHANISM":
                    row["Resistance Mechanism"],

                "RGI_AMR_GENE_FAMILY":
                    row["AMR Gene Family"],

                "RGI_PERCENT_IDENTITY":
                    row["Best_Identities"],

                "RGI_PERCENT_LENGTH":
                    row["Percentage Length of Reference Sequence"],

                "RGI_CUTOFF":
                    row["Cut_Off"],

                "RGI_MODEL_TYPE":
                    row["Model_type"],

                "RGI_SNP":
                    row["SNPs_in_Best_Hit_ARO"]

            })

    return contig_hits


############################################################
# SPLIT FILES
############################################################

def split_results():

    rgi = read_rgi()

    writers = {}
    handles = {}

    total_assignment = 0
    matched = 0
    unmatched = 0

    with open(
        ASSIGNMENT_FILE,
        newline="",
        encoding="utf-8"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            total_assignment += 1

            ####################################################
            # FASTQ
            ####################################################

            fastq = row["FASTQ"]

            sample = re.sub(
                r"\.fastq(\.gz)?$",
                "",
                fastq,
                flags=re.IGNORECASE
            )

            ####################################################
            # CONTIG
            ####################################################

            contig = row["CONTIG"].strip()

            hits = rgi.get(contig)

            ####################################################
            # NO MATCH
            ####################################################

            if not hits:

                unmatched += 1

                hits = [{

                    "RGI_ORF_ID": "NA",
                    "RGI_BEST_HIT": "NA",
                    "RGI_ARO": "NA",
                    "RGI_DRUG_CLASS": "NA",
                    "RGI_RESISTANCE_MECHANISM": "NA",
                    "RGI_AMR_GENE_FAMILY": "NA",
                    "RGI_PERCENT_IDENTITY": "NA",
                    "RGI_PERCENT_LENGTH": "NA",
                    "RGI_CUTOFF": "NA",
                    "RGI_MODEL_TYPE": "NA",
                    "RGI_SNP": "NA"

                }]

            else:

                matched += 1

            ####################################################
            # CREATE WRITER
            ####################################################

            if sample not in writers:

                outfile = os.path.join(
                    OUTDIR,
                    f"{sample}.tsv"
                )

                handles[sample] = open(
                    outfile,
                    "w",
                    newline="",
                    encoding="utf-8"
                )

                fieldnames = (
                    list(row.keys()) +
                    list(hits[0].keys())
                )

                writers[sample] = csv.DictWriter(
                    handles[sample],
                    fieldnames=fieldnames,
                    delimiter="\t"
                )

                writers[sample].writeheader()

            ####################################################
            # WRITE ALL RGI HITS FOR THIS CONTIG
            ####################################################

            for hit in hits:

                output = dict(row)
                output.update(hit)

                writers[sample].writerow(output)

    ############################################################
    # CLOSE FILES
    ############################################################

    for h in handles.values():

        h.close()

    ############################################################
    # SUMMARY
    ############################################################

    print("\n=======================================")
    print("SUMMARY")
    print("=======================================")

    print(f"Assignment rows : {total_assignment}")
    print(f"Contigs matched : {matched}")
    print(f"Contigs missing : {unmatched}")
    print(f"FASTQ files     : {len(handles)}")

    print("\nOutput:")
    print(OUTDIR)


############################################################
# MAIN
############################################################

if __name__ == "__main__":

    split_results()
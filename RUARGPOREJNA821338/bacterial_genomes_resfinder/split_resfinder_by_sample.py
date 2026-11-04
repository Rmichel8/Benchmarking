#!/usr/bin/env python3

import csv
import os
from collections import defaultdict

# INPUT

BASE = (
    "/nfsd/sysbiobig/Riccardom/"
    "RUARGPOREPRJNA821338/"
    "benchmarking/"
    "bacterial_genomes_resfinder/"
)

RESFINDER_FILE = (
    BASE +
    "RUARGPOREPRJNA821338_resfinder_results/"
    "ResFinder_results_tab.txt"
)

ASSIGNMENT_FILE = (
    BASE +
    "contig_sample_assignment/"
    "contig_to_sample.tsv"
)

# OUTPUT

OUTDIR = (
    BASE +
    "RUARGPOREPRJNA821338_resfinder_results/"
    "split_by_sample/"
)

os.makedirs(OUTDIR, exist_ok=True)

# READ CONTIG -> SAMPLE

def read_assignment():

    contig_to_sample = {}

    with open(
        ASSIGNMENT_FILE,
        encoding="utf-8"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            contig = row["contig"].strip()
            sample = row["sample"].strip()

            contig_to_sample[contig] = sample

    return contig_to_sample

# SPLIT RESFINDER

def split_results():

    contig_to_sample = read_assignment()

    writers = {}
    handles = {}

    total = 0
    assigned = 0
    unassigned = 0

    unassigned_contigs = defaultdict(int)

    with open(
        RESFINDER_FILE,
        encoding="utf-8"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        fieldnames = reader.fieldnames

        if fieldnames is None:
            raise RuntimeError(
                "Could not read ResFinder header"
            )

        if "Contig" not in fieldnames:
            raise RuntimeError(
                "'Contig' column not found in ResFinder file"
            )

        for row in reader:

            total += 1

            contig = row["Contig"].strip()

            sample = contig_to_sample.get(contig)

            # CONTIG NOT ASSIGNED

            if sample is None:

                unassigned += 1
                unassigned_contigs[contig] += 1

                continue

            assigned += 1
            # CREATE OUTPUT FILE

            if sample not in writers:

                outfile = os.path.join(
                    OUTDIR,
                    f"{sample}_amr_results.tsv"
                )

                handles[sample] = open(
                    outfile,
                    "w",
                    newline="",
                    encoding="utf-8"
                )

                writers[sample] = csv.DictWriter(
                    handles[sample],
                    fieldnames=fieldnames,
                    delimiter="\t"
                )

                writers[sample].writeheader()


            writers[sample].writerow(row)

    for handle in handles.values():
        handle.close()

    print()
    print("=======================================")
    print("RESFINDER SPLITTING SUMMARY")
    print("=======================================")

    print(f"ResFinder rows       : {total}")
    print(f"Rows assigned        : {assigned}")
    print(f"Rows not assigned    : {unassigned}")
    print(f"Samples found        : {len(writers)}")

    print()
    print("Rows per sample:")

    for sample in sorted(writers):

        outfile = os.path.join(
            OUTDIR,
            f"{sample}_amr_results.tsv"
        )

        with open(
            outfile,
            encoding="utf-8"
        ) as f:

            n = sum(1 for _ in f) - 1

        print(f"  {sample:20s} {n}")

    if unassigned_contigs:

        print()
        print("WARNING: unassigned contigs:")

        for contig, count in sorted(
            unassigned_contigs.items()
        ):

            print(
                f"  {contig}: {count} rows"
            )

    print()
    print("Output:")
    print(OUTDIR)


if __name__ == "__main__":
    split_results()
#!/usr/bin/env python3
import pandas as pd
import glob
import os

BASE = "/nfsd/sysbiobig/Riccardom/benchmarking"
OUTPUT_FILE = os.path.join(BASE, "rgi_best_merged.tsv")

all_dfs = []

# Search files rgi_parsed.tsv nelle run1..run10
for i in range(1, 11):
    pattern = os.path.join(BASE, f"run{i}",f"wf-metagenomics", f"rgi", f"rgi_bwt_parsed.tsv")
    files = glob.glob(pattern)

if not files:
    print("Nessun file trovato. Esco.")
    exit(1)

for f in sorted(files):
    print(f"Leggo {f}")
    df = pd.read_csv(f, sep="\t", dtype=str).fillna("")
    all_dfs.append(df)

# Concat dataframe
merged_df = pd.concat(all_dfs, ignore_index=True)

# Remove duplicates
merged_df = merged_df.drop_duplicates()

# Save
merged_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"Merged RGI scritto in: {OUTPUT_FILE}")
print(f"Totale righe dopo dedup: {len(merged_df)}")

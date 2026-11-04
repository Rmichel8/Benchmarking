#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re
import math

INPUT = Path(
"/nfsd/sysbiobig/Riccardom/ZymoEvenMock/benchmarking/RUARGpore65/ZymoEvenMock/ZymoEvenMock.fa.uniq_ARGpore_2026-10-09.16:00:31/ZymoEvenMock.fa.uniq_arg.w.taxa.tab"
)

OUTDIR = INPUT.parent

GROUND_TRUTH = {
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
}
ALIASES = {

    "transcriptional regulatory protein cpxr cpxr":"cpxR",
    "camp-regulatory protein":"CRP",
    "camp regulatory protein":"CRP",
    "h-ns":"H-NS",
    "DNA-binding_protein_H-NS": "H-NS",
    "EmrB-QacA family major facilitator transporter":"emrB"

}

def normalize(name):

    name = str(name).strip()

    if name in GROUND_TRUTH:
        return name

    low = name.lower()
    if name.startswith("OXA-114"):
        return "OXA-114a"
    if name.lower().startswith("fosa"):
        return "FosA5"
    if low in ALIASES:
        return ALIASES[low]

    m = re.search(r"(SHV-\d+)", name, re.I)
    if m:
        return m.group(1).upper()

    m = re.search(r"(OXA-\d+[A-Za-z]?)", name, re.I)
    if m:
        return m.group(1)

    tokens = re.split(r"[ ;,/()]+", name)

    for t in tokens:

        if t in GROUND_TRUTH:
            return t

    return name

df = pd.read_csv(INPUT, sep="\t")

predicted = {
    normalize(g)
    for g in df["subtype"]
}

tp = sorted(predicted & GROUND_TRUTH)
fp = sorted(predicted - GROUND_TRUTH)
fn = sorted(GROUND_TRUTH - predicted)

precision = len(tp)/(len(tp)+len(fp)) if tp or fp else 0
recall = len(tp)/(len(tp)+len(fn)) if tp or fn else 0

f1 = (
    2*precision*recall/(precision+recall)
    if precision+recall
    else 0
)

rows=[]

for g in sorted(GROUND_TRUTH | predicted):

    rows.append({

        "gene":g,
        "ground_truth":g in GROUND_TRUTH,
        "predicted":g in predicted,

        "status":
        "TP" if g in tp else
        "FP" if g in fp else
        "FN"

    })

summary = pd.DataFrame(rows)

summary.to_csv(
    OUTDIR/"summary_compare.csv",
    index=False
)

with open(OUTDIR/"summary_metrics.txt","w") as f:

    f.write("RUARGpore Gene Benchmark\n\n")

    f.write(f"Ground truth : {len(GROUND_TRUTH)}\n")
    f.write(f"Predicted    : {len(predicted)}\n\n")

    f.write(f"TP : {len(tp)}\n")
    f.write(f"FP : {len(fp)}\n")
    f.write(f"FN : {len(fn)}\n\n")

    f.write(f"Precision : {precision:.3f}\n")
    f.write(f"Recall    : {recall:.3f}\n")
    f.write(f"F1         : {f1:.3f}\n")

print(summary)

print()
print("TP:",len(tp))
print("FP:",len(fp))
print("FN:",len(fn))
print("Precision:",round(precision,3))
print("Recall   :",round(recall,3))
print("F1       :",round(f1,3))
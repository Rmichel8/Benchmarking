#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path
if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} summary_CARD_compare.csv")
    sys.exit(1)

summary = sys.argv[1]

df = pd.read_csv(summary)

TP = int(df["TP"].sum())
FP = int(df["FP"].sum())
TN = int(df["TN"].sum())
FN = int(df["FN"].sum())

accuracy = (TP + TN) / (TP + TN + FP + FN) if (TP + TN + FP + FN) else 0
precision = TP / (TP + FP) if (TP + FP) else 0
recall = TP / (TP + FN) if (TP + FN) else 0
specificity = TN / (TN + FP) if (TN + FP) else 0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0

den = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
mcc = (TP * TN - FP * FN) / den if den else 0

print("\n===== PER SAMPLE =====")
print(df)

print("\n===== OVERALL RESULTS =====")
print(f"TP = {TP}")
print(f"FP = {FP}")
print(f"TN = {TN}")
print(f"FN = {FN}")
print(f"Accuracy    = {accuracy:.3f}")
print(f"Precision   = {precision:.3f}")
print(f"Recall      = {recall:.3f}")
print(f"Specificity = {specificity:.3f}")
print(f"F1-score    = {f1:.3f}")
print(f"MCC         = {mcc:.3f}")
metrics = pd.DataFrame([{

    "benchmark": "RUARGpore",

    "TP": TP,
    "FP": FP,
    "TN": TN,
    "FN": FN,

    "Accuracy": accuracy,
    "Precision": precision,
    "Recall": recall,
    "Specificity": specificity,
    "F1-score": f1,
    "MCC": mcc

}])


outdir = Path(
    "/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/RUARGpore/compare_results"
)

outdir.mkdir(
    parents=True,
    exist_ok=True
)


outfile = outdir / "RUARGpore_global_metrics.csv"


metrics.to_csv(
    outfile,
    index=False
)


print("\nSaved:")
print(outfile)
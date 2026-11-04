#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} summary_CARD_compare.csv")
    sys.exit(1)

summary = sys.argv[1]

df = pd.read_csv(summary)

TP = df["TP"].sum()
FP = df["FP"].sum()
TN = df["TN"].sum()
FN = df["FN"].sum()
tested=TP+TN+FP+FN
precision = TP / (TP + FP) if (TP + FP) else 0
recall = TP / (TP + FN) if (TP + FN) else 0
specificity = TN / (TN + FP) if (TN + FP) else 0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0
accuracy = (TP+TN)/tested
den = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
mcc = (TP * TN - FP * FN) / den if den else 0

print(df)

print("\n===== OVERALL RESULTS =====")
print(f"TP = {TP}")
print(f"FP = {FP}")
print(f"TN = {TN}")
print(f"FN = {FN}")
print(f"Precision  = {precision:.3f}")
print(f"Recall     = {recall:.3f}")
print(f"Specificity= {specificity:.3f}")
print(f"Accuracy= {accuracy:.3f}")
print(f"F1-score   = {f1:.3f}")
print(f"MCC        = {mcc:.3f}")
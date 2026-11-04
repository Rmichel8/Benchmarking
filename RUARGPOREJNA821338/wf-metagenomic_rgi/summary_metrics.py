#!/usr/bin/env python3

import sys
import math
import pandas as pd

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} summary_CARD_compare.csv")
    sys.exit(1)

summary = sys.argv[1]

df = pd.read_csv(summary)

TP = int(df["TP"].sum())
FP = int(df["FP"].sum())
TN = int(df["TN"].sum())
FN = int(df["FN"].sum())

precision = TP/(TP+FP) if TP+FP else 0
recall = TP/(TP+FN) if TP+FN else 0
specificity = TN/(TN+FP) if TN+FP else 0
accuracy = (TP+TN)/(TP+TN+FP+FN)

balanced_accuracy = (recall + specificity)/2

f1 = (
    2*precision*recall/(precision+recall)
    if precision+recall else 0
)

den = math.sqrt(
    (TP+FP)*
    (TP+FN)*
    (TN+FP)*
    (TN+FN)
)

mcc = (TP*TN-FP*FN)/den if den else 0


print("\n========== GLOBAL RESULTS ==========\n")

print(f"TP = {TP}")
print(f"FP = {FP}")
print(f"TN = {TN}")
print(f"FN = {FN}")

print()

print(f"Precision          {precision:.4f}")
print(f"Recall             {recall:.4f}")
print(f"Specificity        {specificity:.4f}")
print(f"Accuracy           {accuracy:.4f}")
print(f"BalancedAccuracy   {balanced_accuracy:.4f}")
print(f"F1-score           {f1:.4f}")
print(f"MCC                {mcc:.4f}")

metrics = pd.DataFrame({

    "Metric":[
        "TP",
        "FP",
        "TN",
        "FN",
        "Precision",
        "Recall",
        "Specificity",
        "Accuracy",
        "BalancedAccuracy",
        "F1",
        "MCC"
    ],

    "Value":[
        TP,
        FP,
        TN,
        FN,
        precision,
        recall,
        specificity,
        accuracy,
        balanced_accuracy,
        f1,
        mcc
    ]

})

outfile = summary.replace(
    ".csv",
    "_overall_metrics.csv"
)

metrics.to_csv(
    outfile,
    index=False
)

print()
print("Metrics written to")
print(outfile)
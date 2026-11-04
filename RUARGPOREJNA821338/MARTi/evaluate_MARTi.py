#!/usr/bin/env python3

import json
import csv
import importlib.util
from pathlib import Path
from collections import defaultdict

# PATHS

ROOT = Path("/nfsd/sysbiobig/Riccardom/runKraken_mock")
CARD_JSON = Path("/nfsd/sysbiobig/Riccardom/database_CARD/card-data/card.json")
GROUND_TRUTH = Path("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/ground_truth.py")
CLASS_TO_AB = {

    "aminoglycoside antibiotic": {
        "Amikacin",
        "Gentamicin",
        "Tobramycin"
    },

    "cephalosporin": {
        "Ceftazidime",
        "Cefepime",
        "Ceftriaxone",
        "Cefuroxime"
    },

    "carbapenem": {
        "Imipenem",
        "Meropenem",
        "Ertapenem"
    },

    "penam": {
        "Amoxicillin/clavulanic acid",
        "Ampicillin/sulbactam",
        "Oxacillin",
        "Penicillin"
    },

    "monobactam": {
        "Aztreonam"
    },

    "phenicol antibiotic": {
        "Chloramphenicol"
    },

    "fluoroquinolone antibiotic": {
        "Ciprofloxacin",
        "Levofloxacin",
        "Norfloxacin"
    },

    "glycylcycline": {
        "Tegacyclin"
    },

    "tetracycline antibiotic": {
        "Tetracycline",
        "Minocycline",
        "Tegacyclin"
    },

    "polymyxin": {
        "Colistin"
    },

    "rifamycin antibiotic": {
        "Rifampicin"
    },

    "lincosamide antibiotic": {
        "Clindamycin"
    },

    "macrolide antibiotic": {
        "Erythromycin"
    },

    "glycopeptide antibiotic": {
        "Vancomycin",
        "Teicoplanin"
    },

    "oxazolidinone antibiotic": {
        "Linezolid"
    },

    "lipopeptide antibiotic": {
        "Daptomycin"
    },

    "mupirocin": {
        "Mupirocin"
    },

    "fusidane antibiotic": {
        "Fusidic acid"
    },

    "nitrofuran antibiotic": {
        "Nitrofurantoin"
    }

}

# load ground truth
spec = importlib.util.spec_from_file_location("gt", GROUND_TRUTH)
gt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gt)

# load CARD

print("Loading CARD...")

with open(CARD_JSON) as f:
    card = json.load(f)

aro2drugclass = {}

for k, v in card.items():

    if not isinstance(v, dict):
        continue

    aro = v.get("ARO_accession")
    if aro is None:
        continue

    aro = f"ARO:{aro}"

    drug_classes = set()

    if "ARO_category" in v:

        for c in v["ARO_category"].values():

            if c.get("category_aro_class_name") == "Drug Class":

                drug_classes.add(c["category_aro_name"])

    aro2drugclass[aro] = drug_classes

print("Loaded", len(aro2drugclass), "ARO entries")
# helper
def aro_to_antibiotics(aro):

    abs_found = set()

    for dc in aro2drugclass.get(aro, []):

        if dc in CLASS_TO_AB:
            abs_found.update(CLASS_TO_AB[dc])

    return abs_found

# output

summary = []
details = []

# samples

for i in range(1,12):

    sample = f"B{i}"

    tp_gt = set(getattr(gt, f"{sample}_TP"))
    tn_gt = set(getattr(gt, f"{sample}_TN"))

    tested = tp_gt | tn_gt

    amr = ROOT / f"B{i:02d}_results" / "marti" / "RUARG"/ "amr.json"

    with open(amr) as f:
        j = json.load(f)

    predicted = set()

    for gene in j["amr"]["geneList"]:

        aro = gene["cardId"]

        antibiotics = aro_to_antibiotics(aro)

        antibiotics &= tested

        predicted |= antibiotics

        for ab in antibiotics:

            details.append({
                "sample":sample,
                "ARO":aro,
                "gene":gene["name"],
                "antibiotic":ab
            })

    TP = predicted & tp_gt
    FP = predicted & tn_gt
    FN = tp_gt - predicted
    TN = tn_gt - predicted

    precision = len(TP)/(len(TP)+len(FP)) if len(TP)+len(FP)>0 else 0
    recall = len(TP)/(len(TP)+len(FN)) if len(TP)+len(FN)>0 else 0

    specificity = len(TN)/(len(TN)+len(FP)) if len(TN)+len(FP)>0 else 0

    accuracy = (len(TP)+len(TN))/len(tested)

    f1 = 0
    if precision+recall>0:
        f1 = 2*precision*recall/(precision+recall)

    summary.append({

        "sample":sample,

        "TP":len(TP),
        "FP":len(FP),
        "FN":len(FN),
        "TN":len(TN),

        "Precision":round(precision,3),
        "Recall":round(recall,3),
        "Specificity":round(specificity,3),
        "Accuracy":round(accuracy,3),
        "F1":round(f1,3),

        "Predicted":",".join(sorted(predicted))
    })


with open("marti_summary.csv","w",newline="") as f:

    w = csv.DictWriter(f,fieldnames=summary[0].keys())
    w.writeheader()
    w.writerows(summary)

with open("marti_predictions.csv","w",newline="") as f:

    w = csv.DictWriter(f,fieldnames=details[0].keys())
    w.writeheader()
    w.writerows(details)

print("Done.")
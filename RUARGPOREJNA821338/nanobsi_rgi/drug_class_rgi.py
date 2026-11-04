# RGI drug classes per barcode

drug_classes_rgi = {
    "barcode01": ['aminoglycoside antibiotic', 'cephalosporin', 'diaminopyrimidine antibiotic', 'disinfecting agents and antiseptics', 'fluoroquinolone antibiotic', 'macrolide antibiotic', 'monobactam', 'nitroimidazole antibiotic', 'penicillin beta-lactam', 'peptide antibiotic', 'phenicol antibiotic', 'phosphonic acid antibiotic', 'rifamycin antibiotic', 'tetracycline antibiotic'],
    "barcode02": ['glycopeptide antibiotic', 'nitroimidazole antibiotic'],
    "barcode03": [],
    "barcode04": [],
    "barcode05": ['aminoglycoside antibiotic', 'glycopeptide antibiotic', 'peptide antibiotic'],
    "barcode06": ['fluoroquinolone antibiotic', 'glycopeptide antibiotic'],
    "barcode07": [],
    "barcode08": [],
    "barcode09": [],
    "barcode10": [],
    "barcode11": [],
}
B1_TP = ["Amikacin" ,"Amoxicillin/clavulanic acid", "Tetracycline"]
B1_TN = ["Aztreonam", "Ceftazidime", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin", "Ceftriaxone",
         "Cefuroxime", "Colistin", "Ertapenem", "Imipenem",
         "Meropenem", "Minocycline", "Nitrofurantoin", "Norfloxacin",
         "Tegacyclin", "Tobramycin"]


B2_TP = []
B2_TN = ["Oxacillin", "Penicillin", "Chloramphenicol", "Ciprofloxacin",
         "Fusidic acid", "Mupirocin", "Clindamycin", "Daptomycin",
         "Erythromycin", "Gentamicin", "Linezolid", "Minocycline",
         "Tetracycline", "Nitrofurantoin", "Rifampicin", "Tegacyclin",
         "Teicoplanin", "Vancomycin"]


B3_TP = []
B3_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
         "Cefoperazone/sulbactam", "Aztreonam", "Ceftazidime",
         "Cefatriaxone", "Cefepime", "Chloramphenicol",
         "Ciprofloxacin", "Levofloxacin", "Colistin", "Imipenem"]


B4_TP = ["Oxacillin", "Penicillin", "Ciprofloxacin",
         "Fusidic acid", "Gentamicin"]
B4_TN = ["Chloramphenicol", "Mupirocin", "Clindamycin", "Daptomycin",
         "Erythromycin", "Linezolid", "Minocycline", "Tetracycline",
         "Nitrofurantoin", "Rifampicin", "Tegacyclin", "Teicoplanin",
         "Vancomycin"]


B5_TP = ["Ceftazidime", "Cefepime", "Chloramphenicol"]
B5_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
         "Cefoperazone/sulbactam", "Aztreonam", "Cefatriaxone",
         "Ciprofloxacin", "Levofloxacin", "Colistin", "Imipenem"]


B6_TP = []
B6_TN = [
    "Amikacin",
    "Gentamicin",
    "Ampicillin/sulbactam",
    "Cefoperazone/sulbactam",
    "Aztreonam",
    "Ceftazidime",
    "Cefatriaxone",
    "Cefepime",
    "Chloramphenicol",
    "Ciprofloxacin",
    "Levofloxacin",
    "Colistin",
    "Imipenem",
    "Oxacillin",
    "Penicillin",
    "Fusidic acid",
    "Mupirocin",
    "Clindamycin",
    "Daptomycin",
    "Erythromycin",
    "Linezolid",
    "Minocycline",
    "Tetracycline",
    "Nitrofurantoin",
    "Rifampicin",
    "Tegacyclin",
    "Teicoplanin",
    "Vancomycin"
]


B7_TP = []
B7_TN = ["Oxacillin", "Penicillin", "Chloramphenicol", "Ciprofloxacin",
         "Fusidic acid", "Mupirocin", "Clindamycin", "Daptomycin",
         "Erythromycin", "Gentamicin", "Linezolid", "Minocycline",
         "Tetracycline", "Nitrofurantoin", "Rifampicin",
         "Tegacyclin", "Teicoplanin", "Vancomycin"]


B8_TP = ["Amoxicillin/clavulanic acid", "Tetracycline", "Tobramycin", "Ceftazidime"] #"Ceftazidime" Intermedio
B8_TN = ["Amikacin", "Aztreonam", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin",
         "Ceftriaxone", "Cefuroxime", "Colistin", "Ertapenem",
         "Imipenem", "Meropenem", "Minocycline",
         "Nitrofurantoin", "Norfloxacin", "Tegacyclin"]


B9_TP = ["Amoxicillin/clavulanic acid", "Tetracycline"]
B9_TN = ["Amikacin", "Aztreonam", "Ceftazidime", "Chloramphenicol",
         "Ciprofloxacin", "Gentamicin", "Levofloxacin",
         "Ceftriaxone", "Cefuroxime", "Colistin", "Ertapenem",
         "Imipenem", "Meropenem", "Minocycline",
         "Nitrofurantoin", "Norfloxacin", "Tegacyclin", "Tobramycin"]


B10_TP = []
B10_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
          "Cefoperazone/sulbactam", "Aztreonam", "Ceftazidime",
          "Cefatriaxone", "Cefepime", "Chloramphenicol",
          "Ciprofloxacin", "Levofloxacin", "Colistin", "Imipenem",
          "Amoxicillin/clavulanic acid", "Tetracycline",
          "Cefuroxime", "Ertapenem", "Meropenem", "Minocycline",
          "Nitrofurantoin", "Norfloxacin", "Tegacyclin", "Tobramycin"]


B11_TP = ["Ceftazidime", "Cefepime", "Clindamycin", "Erythromycin"]
B11_TN = ["Amikacin", "Gentamicin", "Ampicillin/sulbactam",
          "Cefoperazone/sulbactam", "Aztreonam", "Cefatriaxone",
          "Chloramphenicol", "Ciprofloxacin", "Levofloxacin",
          "Colistin", "Imipenem", "Oxacillin", "Penicillin",
          "Fusidic acid", "Mupirocin", "Daptomycin", "Linezolid",
          "Minocycline", "Tetracycline", "Nitrofurantoin",
          "Rifampicin", "Tegacyclin", "Teicoplanin", "Vancomycin"]
CLASS_TO_ANTIBIOTICS = {
    "aminoglycoside antibiotic": {
        "Amikacin", "Gentamicin", "Tobramycin"
    },

    "penicillin beta-lactam": {
        "Penicillin", "Oxacillin", "Amoxicillin/clavulanic acid"
    },

    "cephalosporin": {
        "Cefepime", "Ceftazidime", "Ceftriaxone", "Cefuroxime"
    },

    "carbapenem": {
        "Imipenem", "Meropenem", "Ertapenem"
    },

    "monobactam": {
        "Aztreonam"
    },

    "phenicol antibiotic": {
        "Chloramphenicol"
    },

    "fluoroquinolone antibiotic": {
        "Ciprofloxacin", "Levofloxacin", "Norfloxacin"
    },

    "tetracycline antibiotic": {
        "Tetracycline", "Minocycline", "Tegacyclin"
    },

    "glycylcycline": {
        "Tegacyclin"
    },

    "macrolide antibiotic": {
        "Erythromycin"
    },

    "lincosamide antibiotic": {
        "Clindamycin"
    },

    "glycopeptide antibiotic": {
        "Vancomycin", "Teicoplanin"
    },

    "rifamycin antibiotic": {
        "Rifampicin"
    },

    "nitrofuran antibiotic": {
        "Nitrofurantoin"
    },

    "diaminopyrimidine antibiotic": {
        "Trimethoprim"
    },

    "sulfonamide antibiotic": {
        "Sulfamethoxazole"
    }
}

from math import sqrt
import pandas as pd

TP = FP = FN = TN = 0

rows = []

for barcode in range(1, 12):

    sample = f"barcode{barcode:02d}"

    predicted_classes = drug_classes_rgi[sample]

    predicted_all = set()
    for c in predicted_classes:
        predicted_all |= CLASS_TO_ANTIBIOTICS.get(c, set())

    tp_truth = set(globals()[f"B{barcode}_TP"])
    tn_truth = set(globals()[f"B{barcode}_TN"])

    predicted = predicted_all & (tp_truth | tn_truth)

    tp_set = predicted & tp_truth
    fp_set = predicted & tn_truth
    fn_set = tp_truth - predicted
    tn_set = tn_truth - predicted

    tp = len(tp_set)
    fp = len(fp_set)
    fn = len(fn_set)
    tn = len(tn_set)

    TP += tp
    FP += fp
    FN += fn
    TN += tn

    rows.append({
        "Sample": sample,
        "Predicted classes": ", ".join(sorted(predicted_classes)),
        "Predicted antibiotics": ", ".join(sorted(predicted)),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "TP antibiotics": ", ".join(sorted(tp_set)),
        "FP antibiotics": ", ".join(sorted(fp_set)),
        "FN antibiotics": ", ".join(sorted(fn_set)),
    })
tested=TP+TN+FP+FN
precision = TP / (TP + FP) if TP + FP else 0
recall = TP / (TP + FN) if TP + FN else 0
specificity = TN / (TN + FP) if TN + FP else 0
f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
accuracy = (TP+TN)/tested
den = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
mcc = ((TP * TN) - (FP * FN)) / den if den else 0

print("\n===== GLOBAL RESULTS =====")
print(f"TP = {TP}")
print(f"FP = {FP}")
print(f"FN = {FN}")
print(f"TN = {TN}")
print(f"Precision   = {precision:.3f}")
print(f"Recall      = {recall:.3f}")
print(f"F1-score    = {f1:.3f}")
print(f"Specificity = {specificity:.3f}")
print(f"Accuracy= {accuracy:.3f}")
print(f"MCC         = {mcc:.3f}")

pd.DataFrame(rows).to_csv("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/nanobsi_rgi/rgi_barcode_results.csv", index=False)

pd.DataFrame([{
    "TP": TP,
    "FP": FP,
    "FN": FN,
    "TN": TN,
    "Precision": precision,
    "Recall": recall,
    "F1-score": f1,
    "Specificity": specificity,
    "Accuracy":accuracy,
    "MCC": mcc,
}]).to_csv("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/nanobsi_rgi/rgi_global_metrics.csv", index=False)
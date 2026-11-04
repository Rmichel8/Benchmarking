drug_classes_by_barcode = {
    "barcode01_spe": ['Amikacin', 'Amoxicillin', 'Amoxicillin/clavulanic acid', 'Ampicillin', 'Ampicillin/clavulanic acid', 'Azithromycin', 'Aztreonam', 'Cefepime', 'Cefotaxime', 'Cefoxitin', 'Ceftazidime', 'Ceftriaxone', 'Cephalothin', 'Chloramphenicol', 'Ciprofloxacin', 'Dibekacin', 'Doxycycline', 'Ertapenem', 'Erythromycin', 'Florfenicol', 'Fluoroquinolone', 'Fosfomycin', 'Imipenem', 'Meropenem', 'Minocycline', 'Nalidixic acid', 'Netilmicin', 'Piperacillin', 'Rifampicin', 'Sisomicin', 'Spectinomycin', 'Spiramycin', 'Streptomycin', 'Sulfamethoxazole', 'Telithromycin', 'Tetracycline', 'Ticarcillin', 'Tobramycin', 'Trimethoprim'],
    "barcode02_spe": ['Amoxicillin', 'Amoxicillin/clavulanic acid', 'Ampicillin', 'Ampicillin/clavulanic acid', 'Cefepime', 'Cefixime', 'Cefotaxime', 'Cefoxitin', 'Ceftazidime', 'Clindamycin', 'Doxycycline', 'Ertapenem', 'Erythromycin', 'Imipenem', 'Lincomycin', 'Meropenem', 'Minocycline', 'Piperacillin', 'Pristinamycin IA', 'Quinupristin', 'Tetracycline', 'Tegacyclin', 'Virginiamycin S'],
    "barcode03_spe": [],
    "barcode04_spe": ['Amoxicillin', 'Ampicillin', 'Penicillin', 'Piperacillin'],
    "barcode05_spe": ['Amoxicillin', 'Ampicillin', 'Butiromycin', 'Cefepime', 'Ceftazidime', 'Chloramphenicol', 'Gentamicin', 'Kanamycin', 'Neomycin', 'Paromomycin', 'Ribostamycin', 'Unknown Aminoglycoside'],
    "barcode06_spe": ['Doxycycline', 'Minocycline', 'Tetracycline'],
    "barcode07_spe": ['Clindamycin', 'Doxycycline', 'Erythromycin', 'Lincomycin', 'Minocycline', 'Pristinamycin IA', 'Quinupristin', 'Tetracycline', 'Virginiamycin S'],
    "barcode08_spe": [],
    "barcode09_spe": ['Amoxicillin', 'Ampicillin', 'Cephalothin', 'Piperacillin', 'Ticarcillin'],
    "barcode10_spe": [],
    "barcode11_spe": [],
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


from math import sqrt
import pandas as pd

TP = FP = FN = TN = 0

rows = []

for barcode in range(1, 12):

    sample = f"barcode{barcode:02d}_spe"

    predicted_all = set(drug_classes_by_barcode[sample])

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
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "Predicted": ", ".join(sorted(predicted)),
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

df = pd.DataFrame(rows)
df.to_csv("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/nanobsi/barcode_results.csv", index=False)

metrics = pd.DataFrame([{
    "TP": TP,
    "FP": FP,
    "FN": FN,
    "TN": TN,
    "Precision": precision,
    "Recall": recall,
    "F1-score": f1,
    "Specificity": specificity,
    "Accuracy": accuracy,
    "MCC": mcc,
}])

metrics.to_csv("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/nanobsi/global_metrics.csv", index=False)

print("\nSalvati:")
print(" - barcode_results.csv")
print(" - global_metrics.csv")
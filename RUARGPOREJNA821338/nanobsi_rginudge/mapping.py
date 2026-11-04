from pathlib import Path
import csv

base_dir = Path("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/nanobsi_rginudge")

output = {}

for i in range(1, 12):
    barcode = f"barcode{i:02d}"
    file_path = base_dir / barcode / f"{barcode}.txt"

    if not file_path.exists():
        print(f"[WARN] missing {file_path}")
        continue

    drug_classes = set()

    with open(file_path, newline="") as f:
        # RGI è TSV quasi sempre
        reader = csv.DictReader(f, delimiter="\t")

        if not reader.fieldnames:
            print(f"[ERROR] no header in {file_path}")
            continue

        if "Drug Class" not in reader.fieldnames:
            print(f"[ERROR] 'Drug Class' not found in {file_path}")
            print("Columns:", reader.fieldnames)
            continue

        for row in reader:
            val = row.get("Drug Class")
            if val and val.strip():
                # RGI può avere più classi separate da ;
                for v in val.split(";"):
                    v = v.strip()
                    if v:
                        drug_classes.add(v)

    output[barcode] = sorted(drug_classes)

# save result
out_file = base_dir / "drug_class_rgi.py"

with open(out_file, "w") as f:
    f.write("# RGI drug classes per barcode\n\n")
    f.write("drug_classes_rgi = {\n")
    for k, v in output.items():
        f.write(f'    "{k}": {v},\n')
    f.write("}\n")

print("Saved:", out_file)
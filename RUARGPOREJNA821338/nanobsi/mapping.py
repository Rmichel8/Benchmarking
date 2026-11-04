from pathlib import Path

base_dir = Path("/nfsd/sysbiobig/Riccardom/RUARGPOREPRJNA821338/benchmarking/nanobsi")

output = {}

for bdir in sorted(base_dir.glob("barcode*_spe")):
    barcode = bdir.name
    file_path = bdir / "resfinderc" / "ResFinder_results_tab.txt"

    if not file_path.exists():
        print(f"[WARN] missing {file_path}")
        continue

    phenotypes = set()

    with open(file_path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # skip header
            if line.startswith("Resistance gene"):
                continue

            cols = line.split("\t")

            # safety check
            if len(cols) < 8:
                continue

            phenotype_field = cols[7].strip()

            if phenotype_field and phenotype_field.lower() != "phenotype":
                # split multiple classes separated by comma
                for ph in phenotype_field.split(","):
                    ph = ph.strip()
                    if ph:
                        phenotypes.add(ph)

    output[barcode] = sorted(phenotypes)

# save
out_file = base_dir / "drug_classes.py"

with open(out_file, "w") as f:
    f.write("drug_classes_by_barcode = {\n")
    for k, v in output.items():
        f.write(f'    "{k}": {v},\n')
    f.write("}\n")

print("Saved:", out_file)
#!/usr/bin/env python3
import os, json, csv, statistics
from pathlib import Path

BASE = "/data/habib/EcoSurfBGC"
AS_DIR = os.path.join(BASE, "antismash_output")
RESULTS = os.path.join(BASE, "results")
OUT_FILE = os.path.join(RESULTS, "bgc_matrix_raw.csv")
BGC_DETAIL_FILE = os.path.join(RESULTS, "bgc_detailed_list.csv")
LOG_FILE = os.path.join(BASE, "logs", "parse_antismash.log")
os.makedirs(RESULTS, exist_ok=True)

BGC_CLASSES = [
    "NRPS","PKS","PKS-NRP-Hybrid","T1PKS","T2PKS","T3PKS",
    "transAT-PKS","NRPS-like","PKS-like",
    "lanthipeptide","lanthipeptide-class-i","lanthipeptide-class-ii",
    "lanthipeptide-class-iii","lanthipeptide-class-iv","lanthipeptide-class-v",
    "lassopeptide","sactipeptide","thiopeptide","RiPP-like",
    "head_to_tail","cyclic-lactone-autoinducer","proteusin",
    "terpene","betalactone","phosphonate","ectoine",
    "hserlactone","fatty_acid","arylpolyene","resorcinol",
    "ladderane","PUFA","melanin","oligosaccharide",
    "furan","nucleoside","siderophore","other"
]

def parse_genome(genome_id, as_genome_dir):
    json_files = list(Path(as_genome_dir).glob("*.json"))
    if not json_files:
        return None, None
    with open(json_files[0]) as f:
        try:
            data = json.load(f)
        except:
            return None, None
    record = {"genome_id": genome_id, "total_bgcs": 0}
    for cls in BGC_CLASSES:
        record[cls] = 0
    bgc_list = []
    for seq_record in data.get("records", []):
        contig_id = seq_record.get("id", "")
        for area in seq_record.get("areas", []):
            record["total_bgcs"] += 1
            products = area.get("products", [])
            if isinstance(products, str):
                products = [products]
            bgc_entry = {
                "genome_id": genome_id,
                "contig": contig_id,
                "bgc_number": record["total_bgcs"],
                "products": ";".join(products),
                "start": area.get("start",""),
                "end": area.get("end","")
            }
            bgc_list.append(bgc_entry)
            for prod in products:
                prod_lower = prod.lower().strip()
                matched = False
                for cls in BGC_CLASSES:
                    if cls.lower() == prod_lower:
                        record[cls] += 1
                        matched = True
                        break
                if not matched:
                    for cls in BGC_CLASSES:
                        if cls.lower() in prod_lower:
                            record[cls] += 1
                            break
    return record, bgc_list

print("Parsing antiSMASH v7 output (areas-based)...")
all_records = []
all_bgcs = []
failed = []
genome_dirs = [d for d in Path(AS_DIR).iterdir() if d.is_dir()]
total = len(genome_dirs)
print(f"Found {total} genome directories")

for i, gdir in enumerate(genome_dirs):
    genome_id = gdir.name
    if i % 200 == 0:
        print(f"  Processing {i}/{total}...")
    record, bgc_list = parse_genome(genome_id, gdir)
    if record is None:
        failed.append(genome_id)
        continue
    all_records.append(record)
    all_bgcs.extend(bgc_list)

print("Saving outputs...")
keys = ["genome_id","total_bgcs"] + BGC_CLASSES
with open(OUT_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=keys)
    writer.writeheader()
    writer.writerows(all_records)
print(f"BGC matrix: {OUT_FILE}")

if all_bgcs:
    with open(BGC_DETAIL_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_bgcs[0].keys())
        writer.writeheader()
        writer.writerows(all_bgcs)
    print(f"BGC detail: {BGC_DETAIL_FILE}")

with open(LOG_FILE, "w") as f:
    f.write("\n".join(failed))

print(f"\n=== RESULTS ===")
print(f"Genomes parsed:   {len(all_records)}")
print(f"Failed:           {len(failed)}")
print(f"Total BGCs:       {len(all_bgcs)}")
if all_records:
    bgc_counts = [r["total_bgcs"] for r in all_records]
    print(f"Mean BGCs/genome: {statistics.mean(bgc_counts):.2f}")
    print(f"Median:           {statistics.median(bgc_counts):.1f}")
    print(f"Max:              {max(bgc_counts)}")
    print(f"Genomes 0 BGCs:   {sum(1 for x in bgc_counts if x==0)}")
    classes_to_show = ["NRPS","NRPS-like","T1PKS","PKS-NRP-Hybrid",
        "terpene","siderophore","lanthipeptide","lassopeptide",
        "sactipeptide","RiPP-like","fatty_acid","betalactone"]
    for cls in classes_to_show:
        count = sum(1 for r in all_records if r.get(cls,0)>0)
        print(f"  {cls:<25} {count}")

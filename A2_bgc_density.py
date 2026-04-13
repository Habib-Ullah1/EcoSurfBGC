#!/usr/bin/env python3
"""
Analysis 2: Normalize BGC counts by genome size (BGCs per Mb)
This removes the confound that larger genomes have more BGCs
"""
import os, csv
from pathlib import Path
from collections import defaultdict

BASE = '/data/habib/EcoSurfBGC'
FASTA_DIR = os.path.join(BASE, 'genomes/fasta')
BGC_MATRIX = os.path.join(BASE, 'results/bgc_matrix_final.csv')
OUT_DENSITY = os.path.join(BASE, 'results/bgc_density.csv')

def get_genome_size_mb(fasta_file):
    """Calculate genome size in Mb from FASTA file."""
    total_bp = 0
    n_contigs = 0
    contig_lengths = []
    try:
        with open(fasta_file) as f:
            current_len = 0
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_len > 0:
                        contig_lengths.append(current_len)
                        total_bp += current_len
                    current_len = 0
                    n_contigs += 1
                else:
                    current_len += len(line)
            if current_len > 0:
                contig_lengths.append(current_len)
                total_bp += current_len
    except:
        return None, None, None, None

    size_mb = total_bp / 1_000_000
    # Calculate N50
    if contig_lengths:
        sorted_lens = sorted(contig_lengths, reverse=True)
        cumsum = 0
        n50 = 0
        for l in sorted_lens:
            cumsum += l
            if cumsum >= total_bp / 2:
                n50 = l
                break
    else:
        n50 = 0

    return size_mb, n_contigs, total_bp, n50

print("Calculating genome sizes from FASTA files...")

# Load BGC matrix
bgc_data = {}
with open(BGC_MATRIX) as f:
    reader = csv.DictReader(f)
    for row in reader:
        bgc_data[row['genome_id']] = int(row['total_bgcs'])

print(f"Genomes in BGC matrix: {len(bgc_data)}")

results = []
missing_fasta = 0

fasta_files = list(Path(FASTA_DIR).glob('*.fasta'))
total = len(fasta_files)
print(f"FASTA files found: {total}")

for i, fasta in enumerate(fasta_files):
    if i % 200 == 0:
        print(f"  Processing {i}/{total}...")
    genome_id = fasta.stem
    size_mb, n_contigs, total_bp, n50 = get_genome_size_mb(fasta)
    if size_mb is None:
        missing_fasta += 1
        continue

    total_bgcs = bgc_data.get(genome_id, 0)
    bgc_per_mb = round(total_bgcs / size_mb, 4) if size_mb > 0 else 0

    results.append({
        'genome_id': genome_id,
        'genome_size_mb': round(size_mb, 4),
        'genome_size_bp': total_bp,
        'n_contigs': n_contigs,
        'n50_bp': n50,
        'total_bgcs': total_bgcs,
        'bgc_per_mb': bgc_per_mb
    })

with open(OUT_DENSITY, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print(f"\n=== GENOME SIZE AND BGC DENSITY ===")
sizes = [r['genome_size_mb'] for r in results]
densities = [r['bgc_per_mb'] for r in results]
n_contigs_list = [r['n_contigs'] for r in results]
import statistics
print(f"Genomes processed: {len(results)}")
print(f"Mean genome size: {statistics.mean(sizes):.2f} Mb")
print(f"Median genome size: {statistics.median(sizes):.2f} Mb")
print(f"Size range: {min(sizes):.2f} - {max(sizes):.2f} Mb")
print(f"Mean BGC density: {statistics.mean(densities):.2f} BGCs/Mb")
print(f"Median BGC density: {statistics.median(densities):.2f} BGCs/Mb")
print(f"Mean contigs: {statistics.mean(n_contigs_list):.0f}")
print(f"Saved: {OUT_DENSITY}")

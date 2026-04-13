#!/usr/bin/env python3
"""
Analysis 1: Extract MIBiG similarity from knownclusterblast output
For each BGC, parse the txt file to get:
- Best MIBiG hit accession and description
- Number of protein hits (proxy for similarity)
- Cumulative BLAST score
- Novelty classification
"""
import os, re, csv
from pathlib import Path
from collections import defaultdict

BASE = '/data/habib/EcoSurfBGC'
AS_DIR = os.path.join(BASE, 'antismash_output')
BGC_LIST = os.path.join(BASE, 'results', 'bgc_detailed_list.csv')
OUT_NOVELTY = os.path.join(BASE, 'results', 'bgc_mibig_novelty.csv')
OUT_GENOME = os.path.join(BASE, 'results', 'genome_novelty_summary.csv')

def parse_knownclusterblast(txt_file):
    """Parse one knownclusterblast txt file.
    Returns dict with best hit info."""
    try:
        with open(txt_file) as f:
            content = f.read()
    except:
        return None

    result = {
        'has_hit': False,
        'n_hits': 0,
        'best_bgc_acc': '',
        'best_bgc_desc': '',
        'best_bgc_type': '',
        'best_n_proteins': 0,
        'best_cumulative_score': 0.0,
        'all_hits': []
    }

    # Find significant hits section
    sig_match = re.search(r'Significant hits:\s*\n(.*?)\n\nDetails:', 
                           content, re.DOTALL)
    if not sig_match:
        return result

    hits_text = sig_match.group(1).strip()
    if not hits_text or 'No hits' in hits_text:
        return result

    # Parse hit list
    hit_lines = [l.strip() for l in hits_text.split('\n') if l.strip()]
    bgc_hits = []
    for line in hit_lines:
        m = re.match(r'\d+\.\s+(BGC\d+)\s+(.*)', line)
        if m:
            bgc_hits.append({'acc': m.group(1), 'desc': m.group(2).strip()})

    if not bgc_hits:
        return result

    result['has_hit'] = True
    result['n_hits'] = len(bgc_hits)

    # Parse details for each hit
    details_blocks = re.findall(
        r'>>\n(\d+)\.\s+(BGC\d+)\nSource:\s*(.*?)\nType:\s*(.*?)\n'
        r'Number of proteins with BLAST hits to this cluster:\s*(\d+)\n'
        r'Cumulative BLAST score:\s*([\d.]+)',
        content, re.DOTALL
    )

    best_score = 0
    for block in details_blocks:
        rank, acc, source, bgc_type, n_proteins, cum_score = block
        score = float(cum_score)
        n_prot = int(n_proteins)
        result['all_hits'].append({
            'rank': int(rank),
            'acc': acc,
            'source': source.strip(),
            'type': bgc_type.strip(),
            'n_proteins': n_prot,
            'cumulative_score': score
        })
        if score > best_score:
            best_score = score
            result['best_bgc_acc'] = acc
            result['best_bgc_desc'] = source.strip()
            result['best_bgc_type'] = bgc_type.strip()
            result['best_n_proteins'] = n_prot
            result['best_cumulative_score'] = score

    # If no details parsed but hits listed, use first hit
    if not result['best_bgc_acc'] and bgc_hits:
        result['best_bgc_acc'] = bgc_hits[0]['acc']
        result['best_bgc_desc'] = bgc_hits[0]['desc']

    return result

# ── Load BGC list to map contig → genome + bgc_number ────
print("Loading BGC list...")
bgc_map = {}  # (genome_id, bgc_number) → row
with open(BGC_LIST) as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = (row['genome_id'], int(row['bgc_number']))
        bgc_map[key] = row

print(f"Total BGCs in list: {len(bgc_map)}")

# ── Process all knownclusterblast files ──────────────────
print("Parsing knownclusterblast files...")

genome_dirs = [d for d in Path(AS_DIR).iterdir() if d.is_dir()]
total = len(genome_dirs)

all_bgc_rows = []
genome_stats = defaultdict(lambda: {
    'total_bgcs': 0, 'bgcs_with_hit': 0,
    'bgcs_novel': 0, 'bgcs_highly_similar': 0,
    'best_scores': []
})

for i, gdir in enumerate(genome_dirs):
    genome_id = gdir.name
    if i % 200 == 0:
        print(f"  {i}/{total} genomes processed...")

    kb_dir = gdir / 'knownclusterblast'
    if not kb_dir.exists():
        continue

    # Get all BGC txt files (not region subdirs)
    txt_files = sorted(kb_dir.glob('*.txt'))

    for bgc_idx, txt_file in enumerate(txt_files):
        result = parse_knownclusterblast(txt_file)
        if result is None:
            continue

        genome_stats[genome_id]['total_bgcs'] += 1

        # Novelty classification based on n_proteins hit
        # Novel = 0 hits; Low similarity = 1-2 proteins hit
        # Moderate = 3-4; High = 5+
        n_prot = result['best_n_proteins']
        if not result['has_hit']:
            novelty_class = 'novel'
        elif n_prot <= 2:
            novelty_class = 'low_similarity'
        elif n_prot <= 4:
            novelty_class = 'moderate_similarity'
        else:
            novelty_class = 'high_similarity'

        if result['has_hit']:
            genome_stats[genome_id]['bgcs_with_hit'] += 1
            genome_stats[genome_id]['best_scores'].append(
                result['best_cumulative_score'])
            if n_prot >= 5:
                genome_stats[genome_id]['bgcs_highly_similar'] += 1
        else:
            genome_stats[genome_id]['bgcs_novel'] += 1

        all_bgc_rows.append({
            'genome_id': genome_id,
            'bgc_file': txt_file.stem,
            'bgc_idx': bgc_idx + 1,
            'has_mibig_hit': 1 if result['has_hit'] else 0,
            'n_mibig_hits': result['n_hits'],
            'best_bgc_accession': result['best_bgc_acc'],
            'best_bgc_description': result['best_bgc_desc'],
            'best_bgc_type': result['best_bgc_type'],
            'best_n_proteins': result['best_n_proteins'],
            'best_cumulative_score': result['best_cumulative_score'],
            'novelty_class': novelty_class
        })

print(f"\nTotal BGC records processed: {len(all_bgc_rows)}")

# ── Save BGC-level novelty table ─────────────────────────
with open(OUT_NOVELTY, 'w', newline='') as f:
    if all_bgc_rows:
        writer = csv.DictWriter(f, fieldnames=all_bgc_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_bgc_rows)
print(f"BGC novelty saved: {OUT_NOVELTY}")

# ── Save genome-level novelty summary ────────────────────
import statistics
genome_rows = []
for gid, stats in genome_stats.items():
    total_b = stats['total_bgcs']
    with_hit = stats['bgcs_with_hit']
    novel = stats['bgcs_novel']
    high_sim = stats['bgcs_highly_similar']
    scores = stats['best_scores']

    genome_rows.append({
        'genome_id': gid,
        'total_bgcs': total_b,
        'bgcs_with_mibig_hit': with_hit,
        'bgcs_novel': novel,
        'bgcs_highly_similar': high_sim,
        'pct_novel': round(novel/total_b*100, 1) if total_b > 0 else 0,
        'pct_with_hit': round(with_hit/total_b*100, 1) if total_b > 0 else 0,
        'mean_best_score': round(statistics.mean(scores), 1) if scores else 0,
        'max_best_score': round(max(scores), 1) if scores else 0
    })

with open(OUT_GENOME, 'w', newline='') as f:
    if genome_rows:
        writer = csv.DictWriter(f, fieldnames=genome_rows[0].keys())
        writer.writeheader()
        writer.writerows(genome_rows)
print(f"Genome novelty saved: {OUT_GENOME}")

# ── Print summary statistics ─────────────────────────────
total_bgcs = len(all_bgc_rows)
novel_bgcs = sum(1 for r in all_bgc_rows if r['novelty_class'] == 'novel')
low_sim = sum(1 for r in all_bgc_rows if r['novelty_class'] == 'low_similarity')
mod_sim = sum(1 for r in all_bgc_rows if r['novelty_class'] == 'moderate_similarity')
high_sim = sum(1 for r in all_bgc_rows if r['novelty_class'] == 'high_similarity')

print(f"\n=== MIBiG NOVELTY SUMMARY ===")
print(f"Total BGCs analyzed:          {total_bgcs}")
print(f"Novel (no MIBiG hit):         {novel_bgcs} ({round(novel_bgcs/total_bgcs*100,1)}%)")
print(f"Low similarity (1-2 proteins):{low_sim} ({round(low_sim/total_bgcs*100,1)}%)")
print(f"Moderate (3-4 proteins):      {mod_sim} ({round(mod_sim/total_bgcs*100,1)}%)")
print(f"High similarity (5+ proteins):{high_sim} ({round(high_sim/total_bgcs*100,1)}%)")

print(f"\nTop 10 most common MIBiG hits:")
from collections import Counter
hit_types = Counter(r['best_bgc_description'] for r in all_bgc_rows 
                    if r['best_bgc_description'])
for desc, count in hit_types.most_common(10):
    print(f"  {count:5d}  {desc}")

print(f"\nFiles saved:")
print(f"  {OUT_NOVELTY}")
print(f"  {OUT_GENOME}")

import pandas as pd
import networkx as nx
from collections import defaultdict

BASE = '/data/habib/EcoSurfBGC'
NET_DIR = f'{BASE}/bigscape_output/network_files/2026-03-16_10-57-49_hybrids_auto_EcoSurfBGC_2057genomes'

print("Loading files...")
net    = pd.read_csv(f'{NET_DIR}/mix/mix_c0.30.network', sep='\t')
annot  = pd.read_csv(f'{NET_DIR}/Network_Annotations_Full.tsv', sep='\t')
master = pd.read_csv(f'{BASE}/results/integrated_master_table_v2.csv')
master['IMG.Genome.ID'] = master['IMG.Genome.ID'].astype(str)

our_bgcs = annot[~annot['BGC'].str.startswith('BGC')].copy()
our_bgcs['genome_id'] = our_bgcs['BGC'].str.split('_').str[0]
our_bgc_set = set(our_bgcs['BGC'])
print(f"Our BGCs: {len(our_bgcs)} | MIBiG refs: {len(annot)-len(our_bgcs)}")

print("Building network...")
G = nx.Graph()
G.add_nodes_from(our_bgc_set)

mibig_connected = set()
edges_added = 0
for _, row in net.iterrows():
    n1 = row['Clustername 1']
    n2 = row['Clustername 2']
    d  = row['Raw distance']
    if d <= 0.30:
        if n1 in our_bgc_set and n2 in our_bgc_set:
            G.add_edge(n1, n2)
            edges_added += 1
        if n1.startswith('BGC') and n2 in our_bgc_set:
            mibig_connected.add(n2)
        elif n2.startswith('BGC') and n1 in our_bgc_set:
            mibig_connected.add(n1)

components = list(nx.connected_components(G))
singletons = [c for c in components if len(c)==1]
families   = [c for c in components if len(c)>1]

print(f"\n=== GCF SUMMARY (cutoff 0.30) ===")
print(f"Total GCFs:                    {len(components)}")
print(f"Singleton BGCs (unique):       {len(singletons)}")
print(f"Multi-member GCFs:             {len(families)}")
print(f"BGCs in multi-member GCFs:     {sum(len(c) for c in families)}")
print(f"Largest GCF:                   {max(len(c) for c in components)} BGCs")
mean_fam = sum(len(c) for c in families) / max(len(families),1)
print(f"Mean size (non-singleton):     {mean_fam:.1f} BGCs")
print(f"Edges used:                    {edges_added}")

print(f"\n=== MIBiG CONNECTIVITY ===")
n_novel = len(our_bgc_set) - len(mibig_connected)
print(f"BGCs connected to MIBiG:       {len(mibig_connected)}")
print(f"BGCs with no MIBiG match:      {n_novel}")
print(f"Pct novel (no MIBiG):          {n_novel/len(our_bgc_set)*100:.1f}%")

components_sorted = sorted(components, key=len, reverse=True)
bgc_to_gcf  = {}
gcf_to_size = {}
for gcf_id, comp in enumerate(components_sorted):
    gid = f"GCF_{gcf_id:05d}"
    gcf_to_size[gid] = len(comp)
    for bgc in comp:
        bgc_to_gcf[bgc] = gid

our_bgcs['GCF_id']          = our_bgcs['BGC'].map(bgc_to_gcf)
our_bgcs['GCF_size']        = our_bgcs['GCF_id'].map(gcf_to_size)
our_bgcs['is_singleton']    = our_bgcs['GCF_size'] == 1
our_bgcs['mibig_connected'] = our_bgcs['BGC'].isin(mibig_connected)

def top_class(x):
    m = x.mode()
    return m.iloc[0] if len(m) > 0 else 'Unknown'

genome_gcf = our_bgcs.groupby('genome_id').agg(
    n_bgcs=('BGC','count'),
    n_gcfs=('GCF_id','nunique'),
    n_singletons=('is_singleton','sum'),
    n_mibig_conn=('mibig_connected','sum'),
    bigscape_class=('BiG-SCAPE class', top_class)
).reset_index()

genome_gcf['pct_singleton']  = round(genome_gcf['n_singletons']  / genome_gcf['n_bgcs'] * 100, 1)
genome_gcf['pct_mibig_conn'] = round(genome_gcf['n_mibig_conn']  / genome_gcf['n_bgcs'] * 100, 1)

print(f"\n=== GENOME-LEVEL GCF DIVERSITY ===")
print(f"Genomes with GCF data:         {len(genome_gcf)}")
print(f"Mean BGCs per genome:          {genome_gcf['n_bgcs'].mean():.1f}")
print(f"Mean GCFs per genome:          {genome_gcf['n_gcfs'].mean():.1f}")
print(f"Mean pct singleton BGCs:       {genome_gcf['pct_singleton'].mean():.1f}%")

genome_gcf_eco = genome_gcf.merge(
    master[['IMG.Genome.ID','Ecosystem_Category','C_Phylum',
            'bgc_per_mb','Culture_Status']],
    left_on='genome_id',
    right_on='IMG.Genome.ID',
    how='left'
)

print(f"\n=== GCF DIVERSITY BY ECOSYSTEM ===")
eco_gcf = genome_gcf_eco.groupby('Ecosystem_Category').agg(
    n_genomes=('genome_id','count'),
    mean_bgcs=('n_bgcs','mean'),
    mean_gcfs=('n_gcfs','mean'),
    mean_pct_singleton=('pct_singleton','mean'),
    mean_pct_mibig=('pct_mibig_conn','mean')
).round(2).reset_index().sort_values('mean_gcfs', ascending=False)
print(eco_gcf.to_string(index=False))

print(f"\n=== BiG-SCAPE CLASS DISTRIBUTION ===")
print(our_bgcs['BiG-SCAPE class'].value_counts().head(10))

print(f"\n=== GCF SIZE DISTRIBUTION ===")
for size in [1,2,3,4,5,10,20,50]:
    n = sum(1 for c in components if len(c)==size)
    if n > 0:
        print(f"  Size {size:3d}: {n:5d} GCFs")
print(f"  Size 50+: {sum(1 for c in components if len(c)>=50):5d} GCFs")

our_bgcs.to_csv(f'{BASE}/results/bgc_gcf_assignments.csv', index=False)
genome_gcf.to_csv(f'{BASE}/results/genome_gcf_summary.csv', index=False)
eco_gcf.to_csv(f'{BASE}/results/ecosystem_gcf_summary.csv', index=False)
print(f"\nAll files saved.")

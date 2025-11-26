### Figure 2F
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load RNA-seq data
rnaseq = pd.read_csv('results/GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt', sep='\t')
rnaseq_clean = rnaseq.dropna(subset=['padj'])

# Apply thresholds from paper
sig_genes = rnaseq_clean[(rnaseq_clean['padj'] < 0.01) & 
                         (abs(rnaseq_clean['log2FoldChange']) > 1)]
up_genes = sig_genes[sig_genes['log2FoldChange'] > 1]
down_genes = sig_genes[sig_genes['log2FoldChange'] < -1]

up_gene_set = set(up_genes['genename'])
down_gene_set = set(down_genes['genename'])

# Load annotated peaks
peaks = pd.read_csv('results/homer/annotations/annotated_peaks.txt', 
                    sep='\t', comment='#')

# Check what column contains gene names
#print("Peak annotation columns:")
#print(peaks.columns.tolist())
#print("\nFirst few rows:")
#print(peaks.head())

# Adjust this column name based on your actual file
# Common HOMER column names: 'Gene Name', 'Nearest PromoterID', 'Gene Symbol'
gene_col = 'Gene Name'  # CHANGE THIS if needed

# Function to check if gene has peak within distance
def genes_with_peaks_in_region(peaks_df, distance_col, max_distance):
    """Get genes with peaks within specified distance from TSS"""
    nearby_peaks = peaks_df[abs(peaks_df[distance_col]) <= max_distance]
    return set(nearby_peaks[gene_col].dropna().unique())

def genes_with_peaks_in_gene_body(peaks_df, max_distance=20000):
    """Get genes with peaks within gene body +/- 20kb"""
    # Peaks annotated as promoter, exon, intron, TTS, or within 20kb
    gene_body_annotations = ['promoter', 'exon', 'intron', '5', '3', 'UTR']
    gene_body_peaks = peaks_df[
        (peaks_df['Annotation'].str.contains('|'.join(gene_body_annotations), case=False, na=False)) |
        (abs(peaks_df['Distance to TSS']) <= max_distance)
    ]
    return set(gene_body_peaks[gene_col].dropna().unique())

# Get genes with peaks near TSS (+/- 5kb)
genes_peaks_tss = genes_with_peaks_in_region(peaks, 'Distance to TSS', 5000)

# Get genes with peaks in gene body (+/- 20kb)
genes_peaks_gene = genes_with_peaks_in_gene_body(peaks, 20000)

# Calculate overlaps for each category
# TSS +/- 5kb
up_tss_with = up_gene_set & genes_peaks_tss
up_tss_without = up_gene_set - genes_peaks_tss
down_tss_with = down_gene_set & genes_peaks_tss
down_tss_without = down_gene_set - genes_peaks_tss

# Gene body +/- 20kb
up_gene_with = up_gene_set & genes_peaks_gene
up_gene_without = up_gene_set - genes_peaks_gene
down_gene_with = down_gene_set & genes_peaks_gene
down_gene_without = down_gene_set - genes_peaks_gene

# Print statistics
print(f"\n========== Figure 2F Data ==========")
print(f"\nTSS +/- 5kb:")
print(f"  Up-regulated with RUNX1: {len(up_tss_with)}/{len(up_gene_set)} ({100*len(up_tss_with)/len(up_gene_set):.1f}%)")
print(f"  Down-regulated with RUNX1: {len(down_tss_with)}/{len(down_gene_set)} ({100*len(down_tss_with)/len(down_gene_set):.1f}%)")

print(f"\nGene body +/- 20kb:")
print(f"  Up-regulated with RUNX1: {len(up_gene_with)}/{len(up_gene_set)} ({100*len(up_gene_with)/len(up_gene_set):.1f}%)")
print(f"  Down-regulated with RUNX1: {len(down_gene_with)}/{len(down_gene_set)} ({100*len(down_gene_with)/len(down_gene_set):.1f}%)")

# Create Figure 2F
fig, ax = plt.subplots(figsize=(10, 7))

categories = ['Up-regulated\nTSS ±5kb', 'Down-regulated\nTSS ±5kb',
              'Up-regulated\nGene ±20kb', 'Down-regulated\nGene ±20kb']

# Data for each bar
bound_counts = [len(up_tss_with), len(down_tss_with), 
                len(up_gene_with), len(down_gene_with)]
not_bound_counts = [len(up_tss_without), len(down_tss_without),
                    len(up_gene_without), len(down_gene_without)]
totals = [len(up_gene_set), len(down_gene_set), 
          len(up_gene_set), len(down_gene_set)]

# Calculate percentages
bound_pct = [100 * b / t for b, t in zip(bound_counts, totals)]
not_bound_pct = [100 * nb / t for nb, t in zip(not_bound_counts, totals)]

x = np.arange(len(categories))
width = 0.6

# Create stacked bars
p1 = ax.bar(x, bound_pct, width, label='RUNX1 Bound', 
            color='#d62728', edgecolor='black', linewidth=1.5)
p2 = ax.bar(x, not_bound_pct, width, bottom=bound_pct, label='Not Bound',
            color='#808080', edgecolor='black', linewidth=1.5)

# Add count labels on bars
for i, (b_pct, b_count, nb_pct, nb_count) in enumerate(zip(bound_pct, bound_counts, 
                                                             not_bound_pct, not_bound_counts)):
    # Label for bound (red) section
    ax.text(i, b_pct/2, str(b_count), 
            ha='center', va='center', fontweight='bold', fontsize=12, color='black')
    
    # Label for not bound (grey) section
    ax.text(i, b_pct + nb_pct/2, str(nb_count),
            ha='center', va='center', fontweight='bold', fontsize=12, color='black')

# Formatting
ax.set_ylabel('Percentage of Genes', fontsize=13, fontweight='bold')
ax.set_xlabel('')
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=11)
ax.set_ylim(0, 105)
ax.legend(loc='upper right', fontsize=11, frameon=True, fancybox=True)
ax.set_title('RUNX1 Peak Binding at Differentially Expressed Genes', 
             fontsize=14, fontweight='bold', pad=20)

# Add grid for easier reading
ax.yaxis.grid(True, linestyle='--', alpha=0.3)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig('results/figure_2F_accurate.png', dpi=300, bbox_inches='tight')
plt.show()

# Summary table
summary_data = {
    'Category': categories,
    'RUNX1 Bound': bound_counts,
    'Not Bound': not_bound_counts,
    'Total': totals,
    '% Bound': [f"{p:.1f}%" for p in bound_pct]
}
summary_df = pd.DataFrame(summary_data)
print("\n========== Summary Table ==========")
print(summary_df.to_string(index=False))

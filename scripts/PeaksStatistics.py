import pandas as pd
import numpy as np

# Load data
rnaseq = pd.read_csv('results/GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt', sep='\t')
rnaseq_clean = rnaseq.dropna(subset=['padj']).copy()  # Use .copy() to avoid warning

peaks = pd.read_csv('results/homer/annotations/annotated_peaks.txt', 
                    sep='\t', comment='#')

print("PEAKS STATISTICS")
print("="*80)

print(f"\nTotal peaks: {len(peaks)}")
print(f"Total genes in RNA-seq: {len(rnaseq_clean)}")
print(f"\nPeak Score statistics:")
print(peaks['Peak Score'].describe())
print(f"\nDistance to TSS statistics:")
print(peaks['Distance to TSS'].describe())

# Check a few peak examples
print("\nFirst 5 peaks:")
print(peaks[['Gene Name', 'Peak Score', 'Distance to TSS', 'Annotation']].head())

# Get significant DE genes - use .copy() to avoid warning
sig_genes = rnaseq_clean[(rnaseq_clean['padj'] < 0.01) & 
                         (abs(rnaseq_clean['log2FoldChange']) > 1)].copy()

print(f"\nSignificant DE genes: {len(sig_genes)}")
print(f"  Up-regulated: {len(sig_genes[sig_genes['log2FoldChange'] > 1])}")
print(f"  Down-regulated: {len(sig_genes[sig_genes['log2FoldChange'] < -1])}")

# Lower the peak score threshold since your peaks seem to have low scores
print("\n" + "="*80)
print("FINDING CANDIDATES (using lower thresholds)")
print("="*80)

# Use lower peak score threshold
promoter_peaks = peaks[
    (abs(peaks['Distance to TSS']) <= 5000) &
    (peaks['Peak Score'] > 0)  # Changed from 10 to 0 - accept all peaks
].copy()

print(f"\nPromoter peaks (within 5kb of TSS): {len(promoter_peaks)}")

# Clean gene names
promoter_peaks['Gene Name Clean'] = promoter_peaks['Gene Name'].str.strip()
sig_genes['genename_clean'] = sig_genes['genename'].str.strip()

# Find overlap
overlap = promoter_peaks.merge(
    sig_genes, 
    left_on='Gene Name Clean', 
    right_on='genename_clean',
    how='inner'
)

print(f"Genes with both peaks and significant DE: {len(overlap)}")

if len(overlap) > 0:
    # Sort by significance and fold change
    overlap['abs_log2FC'] = abs(overlap['log2FoldChange'])
    overlap_sorted = overlap.sort_values(['abs_log2FC', 'Peak Score'], ascending=[False, False])
    
    print("\n" + "="*80)
    print("TOP 10 CANDIDATES FOR FIGURES 2D/2E (Down-regulated)")
    print("="*80)
    
    down_candidates = overlap_sorted[overlap_sorted['log2FoldChange'] < -1].head(10)
    
    if len(down_candidates) == 0:
        print("No down-regulated candidates found. Showing all candidates:")
        down_candidates = overlap_sorted.head(10)
    
    for i, (idx, row) in enumerate(down_candidates.iterrows(), 1):
        print(f"\n{i}. Gene: {row['Gene Name Clean']}")
        print(f"   Location: chr{row['Chr']}:{row['Start']:,}-{row['End']:,}")
        print(f"   Peak Score: {row['Peak Score']:.2f}")
        print(f"   Distance to TSS: {row['Distance to TSS']:,} bp")
        print(f"   log2 Fold Change: {row['log2FoldChange']:.2f}")
        print(f"   Adjusted p-value: {row['padj']:.2e}")
        print(f"   Annotation: {row['Annotation']}")
        if 'Gene Description' in row:
            desc = row['Gene Description']
            if isinstance(desc, str) and len(desc) > 0:
                print(f"   Description: {desc[:100]}...")
else:
    print("\nNO OVERLAP FOUND. Checking why...")
    
    # Check some gene name examples from each dataset
    print("\nSample gene names from peaks:")
    print(promoter_peaks['Gene Name Clean'].head(10).tolist())
    
    print("\nSample gene names from RNA-seq:")
    print(sig_genes['genename_clean'].head(10).tolist())

# Show MALAT1 and NEAT1 info
print("\n" + "="*80)
print("PAPER'S GENES (MALAT1, NEAT1, FN1, FBN2, BMP2)")
print("="*80)

paper_genes = ['MALAT1', 'NEAT1', 'FN1', 'FBN2', 'BMP2']
for gene in paper_genes:
    gene_peaks = peaks[peaks['Gene Name'].str.strip() == gene]
    gene_rnaseq = rnaseq_clean[rnaseq_clean['genename'].str.strip() == gene]
    
    print(f"\n{gene}:")
    
    if len(gene_peaks) > 0:
        best_peak = gene_peaks.sort_values('Peak Score', ascending=False).iloc[0]
        print(f"  ✓ Found {len(gene_peaks)} peak(s)")
        print(f"  Location: chr{best_peak['Chr']}:{best_peak['Start']:,}-{best_peak['End']:,}")
        print(f"  Peak Score: {best_peak['Peak Score']:.2f}")
        print(f"  Distance to TSS: {best_peak['Distance to TSS']:,} bp")
        print(f"  Annotation: {best_peak['Annotation']}")
    else:
        print(f"  ✗ No peaks found")
    
    if len(gene_rnaseq) > 0:
        print(f"  ✓ Found in RNA-seq")
        print(f"  log2 Fold Change: {gene_rnaseq.iloc[0]['log2FoldChange']:.2f}")
        print(f"  Adjusted p-value: {gene_rnaseq.iloc[0]['padj']:.2e}")
        
        is_sig = (abs(gene_rnaseq.iloc[0]['log2FoldChange']) > 1 and 
                 gene_rnaseq.iloc[0]['padj'] < 0.01)
        print(f"  Significant DE: {'YES ✓' if is_sig else 'NO ✗'}")
        
        if len(gene_peaks) > 0 and is_sig:
            print(f"  >>> GOOD CANDIDATE FOR FIGURE 2D or 2E <<<")
    else:
        print(f"  ✗ Not found in RNA-seq")


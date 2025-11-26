import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

print("="*80)
print("PREPARING GENE LIST FOR ENRICHR - PROMOTER-BOUND GENES")
print("="*80)

# Load annotated peaks
peaks = pd.read_csv('results/homer/annotations/annotated_peaks.txt', 
                    sep='\t', comment='#')

print(f"\nTotal annotated peaks: {len(peaks):,}")

# Filter to promoter-associated peaks
promoter_peaks = peaks[
    peaks['Annotation'].str.contains('promoter', case=False, na=False)
]

print(f"Promoter-associated peaks: {len(promoter_peaks):,}")

# Extract unique gene names
promoter_genes = promoter_peaks['Gene Name'].dropna().unique()
promoter_genes = [g.strip() for g in promoter_genes if isinstance(g, str) and len(g) > 0]
promoter_genes.sort()  # Alphabetical order

print(f"\n{'='*60}")
print(f"UNIQUE GENES WITH RUNX1 BINDING AT PROMOTERS")
print(f"Total: {len(promoter_genes):,} genes")
print(f"{'='*60}")

# Save to file
output_file = 'results/genes_for_enrichr.txt'
with open(output_file, 'w') as f:
    f.write('\n'.join(promoter_genes))

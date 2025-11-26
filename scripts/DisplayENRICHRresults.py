import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Data extracted from my ENRICHR files
enrichment_data = {
    'Term': [
        # MSigDB Hallmark (top 5)
        'Unfolded Protein Response',
        'E2F Targets',
        'Myc Targets V1',
        'Estrogen Response Late',
        'DNA Repair',
        # GO Biological Process (top 5)
        'Translation',
        'DNA Damage Response',
        'DNA Repair',
        'Gene Expression',
        'Positive Regulation of\nDNA Polymerase Activity',
        # KEGG (top 5)
        'Lysine degradation',
        'Central carbon\nmetabolism in cancer',
        'Citrate cycle\n(TCA cycle)',
        'Glycolysis/\nGluconeogenesis',
        'RNA degradation'
    ],
    'Category': [
        'Hallmark', 'Hallmark', 'Hallmark', 'Hallmark', 'Hallmark',
        'GO Process', 'GO Process', 'GO Process', 'GO Process', 'GO Process',
        'KEGG', 'KEGG', 'KEGG', 'KEGG', 'KEGG'
    ],
    'Adjusted_P_value': [
        # Hallmark
        0.003782349351695863,
        0.04074324492564453,
        0.04074324492564453,
        0.0862101929570542,
        0.0862101929570542,
        # GO Biological Process
        0.02031048448121859,
        0.02031048448121859,
        0.02031048448121859,
        0.07153291231187169,
        0.11503848619384174,
        # KEGG
        0.9734537756006275,
        0.9734537756006275,
        0.9734537756006275,
        0.9734537756006275,
        0.9734537756006275
    ],
    'Combined_Score': [
        # Hallmark
        22.96407112352249,
        10.229767852006498,
        10.229767852006498,
        7.865925946931673,
        7.576985434059747,
        # GO Biological Process
        23.58318599078917,
        19.691030812243948,
        20.824905527394744,
        16.14581294288802,
        170.01801257566038,
        # KEGG
        10.293002437469585,
        9.204266035424881,
        11.424366134981568,
        8.291062649239912,
        7.599967318171916
    ]
}

# Create DataFrame
df = pd.DataFrame(enrichment_data)

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))

# Define colors for categories
category_colors = {
    'Hallmark': '#d62728',
    'GO Process': '#1f77b4',
    'KEGG': '#2ca02c'
}

colors = [category_colors[cat] for cat in df['Category']]

# Sort by combined score for better visualization
df_sorted = df.sort_values('Combined_Score', ascending=True)
colors_sorted = [category_colors[cat] for cat in df_sorted['Category']]

# Plot 1: Combined Score
y_pos = np.arange(len(df_sorted))
bars1 = ax1.barh(y_pos, df_sorted['Combined_Score'], color=colors_sorted, 
                 alpha=0.75, edgecolor='black', linewidth=1.5)

ax1.set_yticks(y_pos)
ax1.set_yticklabels(df_sorted['Term'], fontsize=10)
ax1.set_xlabel('Combined Score', fontsize=13, fontweight='bold')
ax1.set_title('Enrichment Strength', fontsize=14, fontweight='bold', pad=15)
ax1.grid(axis='x', alpha=0.3, linestyle='--')

# Add value labels (only for values < 50 to avoid clutter)
for i, (bar, val) in enumerate(zip(bars1, df_sorted['Combined_Score'])):
    if val < 50:  # Only label reasonable values
        ax1.text(val + 1, i, f'{val:.1f}', 
                 va='center', fontsize=9, fontweight='bold')

# Plot 2: -log10(Adjusted P-value)
df_sorted['-log10_P'] = -np.log10(df_sorted['Adjusted_P_value'])

bars2 = ax2.barh(y_pos, df_sorted['-log10_P'], color=colors_sorted, 
                 alpha=0.75, edgecolor='black', linewidth=1.5)

ax2.set_yticks(y_pos)
ax2.set_yticklabels(df_sorted['Term'], fontsize=10)
ax2.set_xlabel('-log₁₀(Adjusted P-value)', fontsize=13, fontweight='bold')
ax2.set_title('Statistical Significance', fontsize=14, fontweight='bold', pad=15)
ax2.grid(axis='x', alpha=0.3, linestyle='--')

# Add significance threshold line
ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', 
            linewidth=2, alpha=0.6, label='p = 0.05', zorder=0)
ax2.legend(loc='lower right', fontsize=10)

# Add value labels
for i, (bar, val) in enumerate(zip(bars2, df_sorted['-log10_P'])):
    ax2.text(val + 0.05, i, f'{val:.1f}', 
             va='center', fontsize=9, fontweight='bold')

# Add category legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=color, label=cat, alpha=0.75, edgecolor='black') 
                   for cat, color in category_colors.items()]
fig.legend(handles=legend_elements, loc='upper center', 
           ncol=3, fontsize=12, frameon=True, 
           bbox_to_anchor=(0.5, 0.98))

# Main title
fig.suptitle('Top Enriched Pathways in RUNX1 Promoter-Bound Genes', 
             fontsize=16, fontweight='bold', y=0.995)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('results/enrichment_top_pathways_figure.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Saved figure: results/enrichment_top_pathways_figure.png")

# Create summary table
print("\n" + "="*100)
print("TOP ENRICHED PATHWAYS - SUMMARY TABLE")
print("="*100)

summary_table = df.copy()
summary_table['Adjusted P-value'] = summary_table['Adjusted_P_value'].apply(lambda x: f"{x:.2e}")
summary_table['Combined Score'] = summary_table['Combined_Score'].apply(lambda x: f"{x:.1f}")
summary_table = summary_table[['Category', 'Term', 'Adjusted P-value', 'Combined Score']]

# Clean up term names for table
summary_table['Term'] = summary_table['Term'].str.replace('\n', ' ')

print(summary_table.to_string(index=False))
print("\n" + "="*100)

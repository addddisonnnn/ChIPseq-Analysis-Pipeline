import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import os
import re

print("="*80)
print("SUPPLEMENTARY FIGURE S2A - ALIGNMENT STATISTICS")
print("="*80)

# Read flagstat files
flagstat_dir = 'results/flagstat'
flagstat_files = {
    'INPUT Rep1': os.path.join(flagstat_dir, 'INPUT_rep1_flagstat.txt'),
    'INPUT Rep2': os.path.join(flagstat_dir, 'INPUT_rep2_flagstat.txt'),
    'RUNX1 IP Rep1': os.path.join(flagstat_dir, 'IP_rep1_flagstat.txt'),
    'RUNX1 IP Rep2': os.path.join(flagstat_dir, 'IP_rep2_flagstat.txt')
}

stats_data = []

for sample_name, filepath in flagstat_files.items():
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            lines = f.readlines()
            # First line: total reads
            total_line = lines[0]
            total_reads = int(total_line.split()[0])
            
            # Find mapped reads line (usually "X + Y mapped")
            mapped_reads = 0
            for line in lines:
                if 'mapped (' in line and 'primary' not in line:
                    mapped_reads = int(line.split()[0])
                    break
            
            mapping_pct = (mapped_reads / total_reads * 100) if total_reads > 0 else 0
            
            stats_data.append({
                'Sample': sample_name,
                'Total Reads': f"{total_reads:,}",
                'Mapped Reads': f"{mapped_reads:,}",
                'Mapping Rate (%)': f"{mapping_pct:.2f}"
            })
    else:
        print(f"Warning: {filepath} not found")

# Create table
stats_df = pd.DataFrame(stats_data)

print("\nYour Alignment Statistics:")
print(stats_df.to_string(index=False))

# Paper's reported statistics (you'll need to get these from Supplementary Table S2)
print("\n" + "="*80)
print("PAPER'S REPORTED STATISTICS (from Supplementary Table S2)")
print("="*80)
print("""
According to the paper's methods:
- They used single-end 100bp reads
- ChIP-seq samples had biological duplicates
- Expected high mapping rates (>70%) for good quality data

Compare your results above with the paper's supplementary table.
""")

# Create a formatted table figure
fig, ax = plt.subplots(figsize=(10, 4))
ax.axis('tight')
ax.axis('off')

table_data = []
for _, row in stats_df.iterrows():
    table_data.append([row['Sample'], row['Total Reads'], 
                      row['Mapped Reads'], row['Mapping Rate (%)']])

table = ax.table(cellText=table_data,
                colLabels=['Sample', 'Total Reads', 'Mapped Reads', 'Mapping Rate (%)'],
                cellLoc='center',
                loc='center',
                colWidths=[0.25, 0.25, 0.25, 0.25])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# Style header
for i in range(4):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Alternate row colors
for i in range(1, len(table_data) + 1):
    for j in range(4):
        if i % 2 == 0:
            table[(i, j)].set_facecolor('#D9E2F3')

plt.title('Supplementary Table S2A - Alignment Statistics', 
          fontsize=14, fontweight='bold', pad=20)
plt.savefig('results/supp_figure_S2A_table.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✓ Saved table as 'results/supp_figure_S2A_table.png'")

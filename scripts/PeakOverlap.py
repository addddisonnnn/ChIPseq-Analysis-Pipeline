print("\n" + "="*80)
print("SUPPLEMENTARY FIGURE S2C - PEAK OVERLAP BETWEEN REPLICATES")
print("="*80)

# Read peak files for each replicate
rep1_peaks_path = 'results/IP_rep1_peaks.bed'
rep2_peaks_path = 'results/IP_rep2_peaks.bed'

def read_peaks_bed(filepath):
    """Read BED file and return list of peak regions"""
    peaks = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                peaks.append((chrom, start, end))
    return peaks

def peaks_overlap(peak1, peak2):
    """Check if two peaks overlap (same chr and any bp overlap)"""
    chr1, start1, end1 = peak1
    chr2, start2, end2 = peak2
    
    if chr1 != chr2:
        return False
    
    # Check for any overlap
    return not (end1 <= start2 or end2 <= start1)

def count_overlapping_peaks(peaks1, peaks2):
    """Count how many peaks from peaks1 overlap with any peak in peaks2"""
    overlapping = 0
    for peak1 in peaks1:
        for peak2 in peaks2:
            if peaks_overlap(peak1, peak2):
                overlapping += 1
                break  # Count each peak1 only once
    return overlapping

# Read peaks
rep1_peaks = read_peaks_bed(rep1_peaks_path)
rep2_peaks = read_peaks_bed(rep2_peaks_path)

print(f"\nPeak counts:")
print(f"  IP Replicate 1: {len(rep1_peaks):,} peaks")
print(f"  IP Replicate 2: {len(rep2_peaks):,} peaks")

# Calculate overlaps (peaks that have any bp overlap)
rep1_overlap_rep2 = count_overlapping_peaks(rep1_peaks, rep2_peaks)
rep2_overlap_rep1 = count_overlapping_peaks(rep2_peaks, rep1_peaks)

# For Venn diagram, use the counts
rep1_only = len(rep1_peaks) - rep1_overlap_rep2
rep2_only = len(rep2_peaks) - rep2_overlap_rep1
both = rep1_overlap_rep2  # Use rep1's overlapping count

print(f"\nOverlap analysis:")
print(f"  Replicate 1 peaks overlapping Rep2: {rep1_overlap_rep2:,}")
print(f"  Replicate 2 peaks overlapping Rep1: {rep2_overlap_rep1:,}")
print(f"  Replicate 1 only (no overlap): {rep1_only:,}")
print(f"  Replicate 2 only (no overlap): {rep2_only:,}")

# Calculate overlap percentage
overlap_pct_rep1 = (rep1_overlap_rep2 / len(rep1_peaks) * 100) if len(rep1_peaks) > 0 else 0
overlap_pct_rep2 = (rep2_overlap_rep1 / len(rep2_peaks) * 100) if len(rep2_peaks) > 0 else 0

print(f"\nOverlap percentages:")
print(f"  {overlap_pct_rep1:.1f}% of Rep1 peaks overlap with Rep2")
print(f"  {overlap_pct_rep2:.1f}% of Rep2 peaks overlap with Rep1")

# Create Venn diagram
fig, ax = plt.subplots(figsize=(10, 8))

# Note: For Venn diagram, we need to be careful with the numbers
# The intersection should be the average of the two overlap counts
intersection_count = int((rep1_overlap_rep2 + rep2_overlap_rep1) / 2)

venn = venn2(subsets=(rep1_only, rep2_only, intersection_count),
             set_labels=('IP Replicate 1', 'IP Replicate 2'),
             set_colors=('#d62728', '#1f77b4'),
             alpha=0.7)

# Customize labels
for text in venn.set_labels:
    text.set_fontsize(14)
    text.set_fontweight('bold')

for text in venn.subset_labels:
    if text:
        text.set_fontsize(12)
        text.set_fontweight('bold')

plt.title('Peak Overlap Between Replicates (Supplementary Figure S2C)', 
          fontsize=16, fontweight='bold', pad=20)

# Add statistics text box
stats_text = f"""
Total peaks:
  Rep1: {len(rep1_peaks):,}
  Rep2: {len(rep2_peaks):,}

Overlapping:
  {overlap_pct_rep1:.1f}% of Rep1
  {overlap_pct_rep2:.1f}% of Rep2
"""

plt.text(0.02, 0.02, stats_text, transform=ax.transAxes,
         fontsize=11, verticalalignment='bottom',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('results/supp_figure_S2C_venn.png', dpi=300, bbox_inches='tight')
plt.show()

print("\n✓ Saved Venn diagram as 'results/supp_figure_S2C_venn.png'")

# Compare with reproducible peaks
reproducible_peaks_path = 'results/peaks/reproducible_peaks.bed'
if os.path.exists(reproducible_peaks_path):
    reproducible_peaks = read_peaks_bed(reproducible_peaks_path)
    print(f"\nReproducible peaks (after bedtools intersect): {len(reproducible_peaks):,}")
    
    # Percentage relative to smaller set
    smaller_set = min(len(rep1_peaks), len(rep2_peaks))
    print(f"This represents {len(reproducible_peaks)/smaller_set*100:.1f}% of the smaller replicate")
    
    # Percentage relative to overlapping peaks
    avg_overlap = (rep1_overlap_rep2 + rep2_overlap_rep1) / 2
    if avg_overlap > 0:
        print(f"Bedtools intersect retained {len(reproducible_peaks)/avg_overlap*100:.1f}% of overlapping peaks")
    
    print("\nNote: Bedtools intersect is more stringent than simple overlap counting")
    print("It requires peaks to meet specific overlap criteria (default: 1bp minimum)")

# Summary statistics table
print("\n" + "="*80)
print("SUMMARY TABLE")
print("="*80)

summary_data = {
    'Metric': [
        'Total Rep1 peaks',
        'Total Rep2 peaks',
        'Rep1 peaks overlapping Rep2',
        'Rep2 peaks overlapping Rep1',
        'Reproducible peaks (bedtools)',
        'Rep1 overlap %',
        'Rep2 overlap %'
    ],
    'Value': [
        f"{len(rep1_peaks):,}",
        f"{len(rep2_peaks):,}",
        f"{rep1_overlap_rep2:,}",
        f"{rep2_overlap_rep1:,}",
        f"{len(reproducible_peaks):,}" if os.path.exists(reproducible_peaks_path) else "N/A",
        f"{overlap_pct_rep1:.1f}%",
        f"{overlap_pct_rep2:.1f}%"
    ]
}

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

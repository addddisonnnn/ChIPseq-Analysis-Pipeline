print("\n" + "="*80)
print("SUPPLEMENTARY FIGURE S2B - CORRELATION ANALYSIS")
print("="*80)

# Your correlation plot is already generated at:
correlation_plot_path = 'results/correlation_heatmap.png'
correlation_matrix_path = 'results/correlation_matrix.tab'

if os.path.exists(correlation_matrix_path):
    # Read the correlation matrix
    corr_matrix = pd.read_csv(correlation_matrix_path, sep='\t', index_col=0)
    
    print("\nYour Spearman Correlation Matrix:")
    print(corr_matrix.round(3))
    
    # Calculate average correlations
    print("\n" + "="*80)
    print("CORRELATION SUMMARY")
    print("="*80)
    
    # Get sample names
    samples = corr_matrix.columns.tolist()
    
    # Identify IP and INPUT samples
    ip_samples = [s for s in samples if 'IP_' in s and 'INPUT' not in s]
    input_samples = [s for s in samples if 'INPUT' in s]
    
    print(f"\nIP samples: {ip_samples}")
    print(f"INPUT samples: {input_samples}")
    
    # IP-IP correlation
    if len(ip_samples) >= 2:
        ip_corr = corr_matrix.loc[ip_samples[0], ip_samples[1]]
        print(f"\nIP replicate correlation: {ip_corr:.3f}")
    
    # INPUT-INPUT correlation
    if len(input_samples) >= 2:
        input_corr = corr_matrix.loc[input_samples[0], input_samples[1]]
        print(f"INPUT replicate correlation: {input_corr:.3f}")
    
    # IP-INPUT correlation (should be lower)
    if len(ip_samples) > 0 and len(input_samples) > 0:
        ip_input_corrs = []
        for ip in ip_samples:
            for inp in input_samples:
                ip_input_corrs.append(corr_matrix.loc[ip, inp])
        avg_ip_input = np.mean(ip_input_corrs)
        print(f"Average IP-INPUT correlation: {avg_ip_input:.3f}")
    
    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)
    print("""
Expected patterns:
- High correlation between replicates (>0.9): Indicates good reproducibility
- Lower correlation between IP and INPUT (<0.8): Expected, as they represent different samples
- IP replicates should correlate better with each other than with INPUT

Paper's findings (Supplementary Figure S2B):
- The authors showed high correlation between replicates
- This validates the quality and reproducibility of the ChIP-seq experiment
    """)
    
    # Display your plot
    if os.path.exists(correlation_plot_path):
        from IPython.display import Image, display
        print("\nYour Correlation Heatmap:")
        display(Image(correlation_plot_path))
    
else:
    print(f"Warning: Correlation matrix not found at {correlation_matrix_path}")

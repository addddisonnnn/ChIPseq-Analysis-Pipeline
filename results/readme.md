## Key Results

### Quality Control Metrics
| Sample | Total Reads (M) | Mapping Rate | Duplication Rate | GC Content |
|--------|-----------------|--------------|------------------|------------|
| INPUT_rep1 | ~30 | 89.1% | 1-4% | 43-47% |
| INPUT_rep2 | ~10.7 | 74.3% | 1-4% | 43-47% |
| IP_rep1 | ~30 | 89.1% | 12-13% | 43-47% |
| IP_rep2 | ~29 | 74.3% | 12-13% | 43-47% |

**QC Assessment**: All samples showed high-quality sequencing data with Q30+ base quality scores, minimal adapter contamination, and appropriate GC content. Higher duplication rates in IP samples are expected due to enrichment.

### RUNX1 Binding Landscape
- **Total peaks called**: 89,360 (Rep1) vs 22,897 (Rep2)
- **Reproducible peaks**: 6,462 high-confidence binding sites
- **Overlap efficiency**: 7.2% of Rep1 peaks, 28.2% of Rep2 peaks
- **Promoter binding**: Strong enrichment at transcription start sites

### Motif Enrichment Analysis
- **Primary motif**: Canonical RUNX family binding motif (p < 1e-50)
- **Validation**: Confirms immunoprecipitation specificity
- **Co-factors**: Secondary motifs suggest potential RUNX1 interaction partners

### Integration with Gene Expression
| Regulation | Promoter Peaks (±5kb) | Gene Body Peaks |
|------------|----------------------|-----------------|
| Up-regulated | 14.7% | 18.3% |
| Down-regulated | 15.2% | 18.7% |

**Key Direct Targets**:
- **MALAT1**: Direct promoter binding (0bp from TSS), log2FC = -1.67
- **PIDD1**: Promoter-proximal binding (59bp from TSS), log2FC = -2.34

### Functional Enrichment Analysis

#### Hallmark Pathways (MSigDB)
![Hallmark Pathways](results/MSigDB_Hallmark_2020_bar_graph.png)

**Top Enriched Pathways**:
1. **E2F Targets** - Cell cycle regulation
2. **G2-M Checkpoint** - Mitotic control
3. **Myc Targets** - Proliferation signaling
4. **Estrogen Response** - Hormone signaling
5. **DNA Repair** - Genomic stability

#### Biological Processes (GO)
![Biological Processes](results/GO_Biological_Process_2025_bar_graph.png)

**Key Processes**:
- Cell cycle phase transition
- DNA damage response
- Double-strand break repair
- Chromatin organization
- Regulation of transcription

#### KEGG Pathways
![KEGG Pathways](results/KEGG_2021_Human_bar_graph.png)

**Signaling Networks**:
- Cell cycle
- p53 signaling pathway
- Metabolic pathways
- PI3K-Akt signaling

### Correlation Analysis
- **INPUT replicates**: High correlation (r = 0.842) - technical reproducibility
- **IP replicates**: Moderate correlation (r = 0.574) - biological consistency
- **Sample clustering**: Clear separation between INPUT and IP conditions

### Signal Profile Characteristics
![Signal Profile](results/profile_plot.png)

- **TSS enrichment**: Characteristic peak immediately upstream of transcription start sites
- **Gene body signal**: Consistent binding pattern across both replicates
- **Reproducibility**: Concordant profiles between biological replicates

### Comparative Analysis with Original Study
- **Methodological alignment**: Conservative peak calling approach
- **Biological validation**: Confirmed RUNX1's role in breast cancer pathogenesis
- **Extended insights**: Broader regulatory network beyond extracellular matrix genes

### Technical Performance
- **Pipeline reproducibility**: Successful implementation of modular Nextflow workflow
- **Resource efficiency**: 8GB RAM and 4 CPUs per alignment/peak calling job
- **Data integrity**: Comprehensive QC throughout processing steps

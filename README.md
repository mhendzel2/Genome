# Genomics Data Analysis Platform

A comprehensive Streamlit-based web application for genomics data analysis with integrated statistical methods, multiple data source support, and advanced analytical capabilities.

## Features

### Dataset Discovery & Download
- **ENCODE Database**: Search and browse functional genomics experiments (ChIP-seq, RNA-seq, ATAC-seq, etc.)
- **GEO (Gene Expression Omnibus)**: Access gene expression datasets with series matrix downloads
- **TCGA/GDC**: Cancer genomics data with clinical annotations
- **DepMap**: Cancer cell line dependency data from Figshare releases
- Direct download capabilities for all supported databases

### Data Analysis Tools

#### Statistical Analysis
- **Parametric Tests**: t-tests (paired/unpaired), ANOVA, correlation analysis
- **Non-parametric Tests**: Mann-Whitney U, Kruskal-Wallis, Wilcoxon signed-rank
- **Multiple Testing Correction**: Bonferroni, FDR (Benjamini-Hochberg), Holm methods
- **Survival Analysis**: Log-rank tests for survival data
- **Clustering**: Hierarchical clustering, K-means, PCA
- **Normality & Variance Tests**: Shapiro-Wilk, Kolmogorov-Smirnov, Levene's test

#### Genomics-Specific Analysis
- **Basic Statistics**: Region counting, score distributions, length statistics
- **Peak Calling**: ChIP-seq peak identification with customizable thresholds
- **Differential Expression**: DESeq2-compatible analysis with FDR correction
- **Enrichment Analysis**: GO/pathway enrichment via g:Profiler
- **Quality Control**: Format validation, coordinate checking, coverage metrics
- **Tissue Comparison**: Cross-tissue genomic feature comparison
- **Chromatin Interactions**: Hi-C data analysis (experimental)
- **Multi-omics Integration**: Correlation-based integration across data types

#### Specialized Data Types
- **Proteomics**: MaxQuant data loading, protein quantification, differential expression
- **Methylation**: Bisulfite-seq analysis, DMR detection, CpG island annotation
- **Gene Expression**: RNA-seq, microarray, with header detection and normalization
- **Epigenomics**: Histone modifications, accessibility data (ATAC-seq, DNase-seq)

### Data Validation & QC
- Comprehensive file format validation (BED, VCF, GTF/GFF, gene expression)
- Coordinate validity checking
- Missing value detection and reporting
- File checksum calculation
- Automatic format detection
- Quality metrics reporting

### Interactive Visualizations
- **Enrichment Plots**: Bar charts for pathway/GO enrichment
- **Heatmaps**: Correlation matrices and expression heatmaps
- **PCA Plots**: Principal component analysis with variance explained
- **Volcano Plots**: Differential expression visualization (coming soon)
- **QC Dashboards**: Quality control metrics visualization

## Installation

### System Requirements
- Python 3.11+
- R (optional, for advanced genomics analyses)
- Fortran compiler and OpenBLAS (for scipy on some systems)

### Installation Steps

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. **Install system dependencies (Linux/Debian-based):**
   ```bash
   sudo apt-get update && sudo apt-get install -y gfortran libopenblas-dev
   ```

3. **Create a virtual environment (recommended):**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

4. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

5. **Install R packages (optional, for R integration):**
   ```R
   install.packages("BiocManager")
   BiocManager::install(c("GenomicRanges", "DESeq2", "ChIPseeker"))
   ```

6. **Run the application:**
   ```bash
   streamlit run app.py
   ```

## Usage Guide

### Dataset Discovery
1. Navigate to the "Dataset Discovery" page
2. Select a data source (ENCODE, GEO, DepMap, TCGA)
3. Enter search terms (e.g., "H3K4me3", "GSE12345", "breast cancer")
4. Browse results and click "Download" to retrieve datasets

### Running Analyses
1. Go to the "Analysis" page
2. Select an analysis type:
   - **Basic Statistics**: Upload BED/bedGraph files
   - **Enrichment Analysis**: Upload gene expression data
   - **Differential Expression**: Upload control and treatment files
   - **Quality Control**: Upload files for validation
3. Upload your data files
4. Configure parameters if needed
5. Click "Run Analysis"

### Visualization
1. After running an analysis, go to the "Visualization" page
2. Results are automatically loaded
3. Select visualization types:
   - Enrichment bar charts
   - Heatmaps
   - PCA plots
4. Customize and export plots

## Supported File Formats

- **BED/bedGraph**: Genomic intervals and scores
- **Gene Expression**: Tab-delimited with optional header
- **VCF**: Variant calls
- **GTF/GFF**: Gene annotations
- **MaxQuant**: Proteomics proteinGroups.txt
- **Bismark**: Methylation coverage files
- **Generic TSV**: Tab-separated values with flexible parsing

## Advanced Features

### Statistical Analysis API
```python
from utils.statistical_analysis import StatisticalAnalyzer

analyzer = StatisticalAnalyzer()

# T-test
result = analyzer.t_test_analysis(group1, group2, paired=False)

# Multiple testing correction
corrected = analyzer.multiple_testing_correction(pvalues, method='fdr_bh')

# PCA
pca_result = analyzer.pca_analysis(data_matrix)
```

### Direct Data Download
```python
from utils.geo_client import GEOClient
from utils.tcga_client import TCGAClient

# Download from GEO
geo = GEOClient()
series_file = geo.download_series_matrix('GSE12345')

# Download from TCGA
tcga = TCGAClient()
projects = tcga.list_projects(program='TCGA')
```

### R Integration
The platform integrates with R for specialized genomics analyses. Set `r_available=True` and ensure R packages are installed.

## Configuration

### Database Setup (Optional)
For persistent storage, configure PostgreSQL:
```bash
export DATABASE_URL="postgresql://user:password@localhost/genomics_db"
```

For development, SQLite is used automatically if DATABASE_URL is not set.

### Logging
Logs are written to `app.log` and `streamlit.log`. Configure logging level in `app.py`:
```python
logging.basicConfig(level=logging.INFO)
```

## API Documentation

### Statistical Methods
- `t_test_analysis()`: Parametric two-group comparison
- `anova_analysis()`: Multi-group comparison
- `multiple_testing_correction()`: P-value adjustment
- `fold_change_analysis()`: Calculate log2 fold changes
- `correlation_analysis()`: Pearson, Spearman, Kendall
- `hierarchical_clustering()`: Dendrogram-based clustering
- `kmeans_clustering()`: Partition-based clustering
- `pca_analysis()`: Dimensionality reduction
- `survival_analysis_logrank()`: Survival curve comparison

### Genomics Analyses
- `basic_statistics()`: Summary statistics for genomic regions
- `peak_calling()`: Identify significant peaks in ChIP-seq
- `differential_expression()`: DE analysis with FDR correction
- `enrichment_analysis()`: GO/pathway enrichment
- `quality_control()`: Validate and assess data quality

## Development

### Running Tests
```bash
pytest utils/tests/
pytest analysis/tests/
```

### Code Structure
```
├── app.py                      # Main Streamlit application
├── analysis/                   # Analysis modules
│   ├── genomics_analysis.py   # Core genomics analyses
│   ├── analyzer_adapter.py    # Adapter pattern for analyses
│   └── tests/                  # Unit tests
├── utils/                      # Utility modules
│   ├── statistical_analysis.py # Statistical methods
│   ├── geo_client.py          # GEO database client
│   ├── tcga_client.py         # TCGA/GDC client
│   ├── proteomics_analysis.py # Proteomics support
│   ├── methylation_analysis.py # Methylation analysis
│   ├── data_validation.py     # File validation
│   ├── dataset_discovery.py   # ENCODE search
│   ├── depmap_client.py       # DepMap data
│   ├── gprofiler_client.py    # Enrichment analysis
│   ├── visualization.py       # Plotting functions
│   └── tests/                  # Unit tests
└── requirements.txt            # Python dependencies
```

## Troubleshooting

### Common Issues

1. **Import errors**: Ensure all dependencies are installed: `pip install -r requirements.txt`
2. **R integration fails**: Install R and required Bioconductor packages
3. **Memory errors**: Reduce dataset size or increase system memory
4. **File format errors**: Use data validation to check file format

### Getting Help
- Check logs in `app.log` and `streamlit.log`
- Review error messages in the Streamlit interface
- Ensure file formats match expected specifications

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

This project is provided for research and educational purposes.

## Citation

If you use this platform in your research, please cite:
```
[Citation information to be added]
```

## Acknowledgments

- ENCODE Consortium for genomics data
- NCBI GEO for gene expression data
- NCI GDC for cancer genomics data
- Broad Institute DepMap for cell line data
- Bioconductor project for R packages

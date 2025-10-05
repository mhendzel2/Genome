# Quick Start Guide - Genomics Data Analysis Platform

## Installation (5 minutes)

```bash
# 1. Clone and navigate
git clone <repository-url>
cd Genome

# 2. Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the application
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

## Quick Examples

### 1. Search and Download from GEO (2 minutes)

```python
from utils.geo_client import GEOClient

# Initialize client
geo = GEOClient()

# Search for datasets
results = geo.search_datasets('breast cancer')
print(f"Found {len(results)} datasets")

# Download a series matrix
file_path = geo.download_series_matrix('GSE12345')
print(f"Downloaded to: {file_path}")
```

### 2. Run Statistical Analysis (1 minute)

```python
from utils.statistical_analysis import StatisticalAnalyzer
import numpy as np

# Initialize analyzer
stats = StatisticalAnalyzer()

# Generate sample data
control = np.random.normal(10, 2, 30)
treatment = np.random.normal(12, 2, 30)

# Run t-test
result = stats.t_test_analysis(control, treatment)
print(f"p-value: {result['p_value']}")
print(f"Cohen's d: {result['cohens_d']}")

# Multiple testing correction
pvalues = [0.001, 0.01, 0.05, 0.1, 0.5]
corrected = stats.multiple_testing_correction(pvalues, method='fdr_bh')
print(f"Significant: {corrected['num_significant']}")
```

### 3. Analyze Proteomics Data (3 minutes)

```python
from utils.proteomics_analysis import ProteomicsAnalyzer

# Initialize analyzer
proteomics = ProteomicsAnalyzer()

# Load MaxQuant data
df = proteomics.load_maxquant_data('proteinGroups.txt')

# Normalize intensities
intensity_cols = ['Intensity sample 1', 'Intensity sample 2']
df_norm = proteomics.normalize_intensities(df, intensity_cols, method='median')

# Differential expression
control_cols = ['Intensity control 1', 'Intensity control 2']
treatment_cols = ['Intensity treatment 1', 'Intensity treatment 2']
results = proteomics.differential_expression_proteomics(
    df_norm, control_cols, treatment_cols
)
print(f"Significant proteins: {results['Significant'].sum()}")
```

### 4. Download from TCGA (2 minutes)

```python
from utils.tcga_client import TCGAClient

# Initialize client
tcga = TCGAClient()

# List available projects
projects = tcga.list_projects(program='TCGA')
print(f"Found {len(projects)} TCGA projects")

# Search for RNA-seq data
files = tcga.search_files(
    project_id='TCGA-BRCA',
    data_type='Gene Expression Quantification',
    experimental_strategy='RNA-Seq',
    size=10
)
print(f"Found {len(files)} files")

# Download a file
if files:
    file_path = tcga.download_file(files[0]['file_id'])
    print(f"Downloaded to: {file_path}")
```

### 5. Validate and QC Data (1 minute)

```python
from utils.data_validation import DataValidator

# Initialize validator
validator = DataValidator()

# Validate a BED file
with open('peaks.bed', 'rb') as f:
    validation = validator.validate_file(f, 'BED')
    
print(f"Valid: {validation['valid']}")
print(f"Errors: {validation['errors']}")
print(f"Metrics: {validation['metrics']}")

# Calculate checksum
checksum = validator.calculate_file_checksum(f)
print(f"MD5: {checksum}")
```

## Common Use Cases

### Differential Gene Expression Analysis

1. Upload control and treatment RNA-seq count files
2. Navigate to "Analysis" → "Differential Expression"
3. Upload files and click "Run Analysis"
4. View results with p-values, fold changes, and FDR
5. Visualize with volcano plots

### Enrichment Analysis

1. Upload a gene list (one gene per line or gene expression data)
2. Select "Enrichment Analysis"
3. Configure thresholds (default: |fold change| > 2)
4. View GO terms and pathways
5. Generate enrichment bar charts

### ChIP-seq Peak Analysis

1. Upload BED file with ChIP-seq peaks
2. Select "Basic Statistics" or "Peak Calling"
3. Configure fold change threshold
4. View peak statistics and distributions
5. Export results

### Quality Control Pipeline

1. Upload multiple files
2. Select "Quality Control"
3. Specify file type
4. View comprehensive QC metrics:
   - Coordinate validity
   - Missing values
   - Format compliance
   - Coverage statistics

## Tips and Tricks

### Performance Optimization
- Filter data before analysis to reduce memory usage
- Use appropriate significance thresholds
- Close unused browser tabs
- Restart the app if memory issues occur

### Data Preparation
- Ensure files are tab-delimited
- Remove special characters from headers
- Check coordinate formats (0-based vs 1-based)
- Compress large files with gzip

### Troubleshooting
```bash
# Check logs
tail -f app.log

# Clear cache
rm -rf .streamlit/cache

# Reinstall dependencies
pip install -r requirements.txt --force-reinstall

# Test R integration
Rscript -e "library(GenomicRanges)"
```

## Next Steps

1. **Explore the Documentation**: Read the full README.md
2. **Try the Examples**: Run the example scripts in `examples/`
3. **Run Tests**: `pytest utils/tests/ analysis/tests/`
4. **Customize**: Modify analysis parameters and thresholds
5. **Contribute**: Submit issues and pull requests

## Getting Help

- Check `docs/IMPROVEMENTS_SUMMARY.md` for detailed feature documentation
- Review `app.log` and `streamlit.log` for error messages
- Consult the README.md for API documentation
- Open an issue on GitHub for bugs or feature requests

## Key Features to Try

1. ✅ Download data from GEO, TCGA, DepMap, ENCODE
2. ✅ Run 15+ statistical tests
3. ✅ Analyze proteomics and methylation data
4. ✅ Validate data quality
5. ✅ Create interactive visualizations
6. ✅ Export results and figures

Happy analyzing! 🧬

# Genomics Data Analysis Platform

A comprehensive Streamlit-based web application for genomics data analysis, designed to provide researchers with tools for uploading, validating, analyzing, and visualizing genomics datasets.

## Features

### Core Functionality
- **Data Upload & Validation**: Support for multiple genomics file formats (BED, bedGraph, BigWig, WIG, etc.)
- **Dataset Discovery**: Search and browse public genomics databases (ENCODE, GEO, ArrayExpress, TCGA)
- **Analysis Tools**: Basic statistics, peak calling, differential expression, tissue comparison
- **Interactive Visualizations**: Genome browser views, heatmaps, volcano plots, PCA plots
- **Database Integration**: PostgreSQL backend for persistent data storage

### Supported Data Types
- **Histone Marks**: BED, bedGraph, BigWig, WIG formats
- **Hi-C**: TXT, BEDPE, HIC, COOL formats  
- **Gene Expression**: CSV, TSV, TXT, FPKM, TPM formats
- **ChIP-seq**: BED, bedGraph, narrowPeak, broadPeak, BigWig formats

### Analysis Capabilities
- Basic statistics computation
- Peak calling for ChIP-seq data
- Differential expression analysis
- Tissue comparison studies
- Enrichment analysis
- Quality control assessment

## Installation

### Prerequisites
- Python 3.11+
- PostgreSQL database
- R (optional, for advanced analysis)

### Local Installation

1. **Clone or extract the project**
```bash
cd genomics-analysis-platform
```

2. **Install Python dependencies**
```bash
pip install -r requirements.txt
```

3. **Set up PostgreSQL database**
   - Install PostgreSQL on your system
   - Create a new database for the application
   - Set the following environment variables:
   ```bash
   export DATABASE_URL="postgresql://username:password@localhost:5432/genomics_db"
   export PGUSER="your_username"
   export PGPASSWORD="your_password"
   export PGDATABASE="genomics_db"
   export PGHOST="localhost"
   export PGPORT="5432"
   ```

4. **Install R packages (optional)**
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("GenomicRanges", "DESeq2", "ChIPseeker"))
   install.packages(c("dplyr", "ggplot2"))
   ```

5. **Run the application**
```bash
streamlit run app.py --server.port 5000
```

## Usage

### Data Upload
1. Navigate to "Data Upload & Validation"
2. Select your data type (Histone Marks, Hi-C, Gene Expression, ChIP-seq)
3. Upload your genomics files
4. Review validation results and file previews

### Dataset Discovery
1. Go to "Dataset Discovery"
2. Use search terms, filters, or browse popular datasets
3. Select datasets for analysis
4. View detailed dataset information

### Data Analysis
1. Access "Data Analysis" with uploaded files or selected datasets
2. Choose analysis type:
   - Basic Statistics
   - Peak Calling (ChIP-seq)
   - Differential Expression
   - Chromatin Interaction Analysis
   - Enrichment Analysis
3. Configure parameters and run analysis

### Tissue Comparison
1. Select "Tissue Comparison"
2. Choose tissues to compare
3. Select comparison type
4. View comparative analysis results

### Visualization
1. Navigate to "Visualization"
2. Select analysis results to visualize
3. Choose visualization type
4. Generate and download interactive plots

## Database Features

The application includes a PostgreSQL database that stores:
- Public genomics datasets with metadata
- Analysis results and parameters
- Uploaded file information
- User session data

Access "Database Status" to view:
- Dataset statistics
- Analysis history
- Database management tools

## Project Structure

```
genomics-analysis-platform/
├── app.py                          # Main Streamlit application
├── requirements.txt                # Python dependencies
├── README.md                       # This file
├── .streamlit/
│   └── config.toml                # Streamlit configuration
├── utils/                         # Utility modules
│   ├── data_validation.py         # File validation logic
│   ├── r_integration.py           # R integration for analysis
│   ├── visualization.py           # Plotting and visualization
│   ├── dataset_discovery.py       # Dataset search functionality
│   └── database.py                # Database models and operations
└── analysis/                      # Analysis modules
    ├── genomics_analysis.py       # Main analysis coordinator
    └── r_scripts.py               # R script templates
```

## Configuration

### Streamlit Configuration
The application uses custom Streamlit settings in `.streamlit/config.toml`:
- Headless server mode
- Custom port binding (5000)
- Light theme

### Database Configuration
Configure your PostgreSQL connection using environment variables:
- `DATABASE_URL`: Complete connection string
- `PGUSER`, `PGPASSWORD`, `PGDATABASE`, `PGHOST`, `PGPORT`: Individual connection parameters

## Development

### Adding New Analysis Types
1. Implement analysis logic in `analysis/genomics_analysis.py`
2. Add R script templates in `analysis/r_scripts.py`
3. Update the UI in `app.py`

### Adding New Visualizations
1. Add visualization methods to `utils/visualization.py`
2. Update the visualization page in `app.py`

### Database Schema Updates
1. Modify models in `utils/database.py`
2. Handle migrations appropriately
3. Update the database manager methods

## Troubleshooting

### Database Connection Issues
- Verify PostgreSQL is running
- Check environment variables
- Ensure database exists and user has permissions

### R Integration Issues
- Install required R packages
- Check R is available in system PATH
- Review R script execution logs

### File Upload Problems
- Check file format compatibility
- Verify file size limits
- Review validation error messages

## License

This project is provided as-is for research and educational purposes.

## Support

For issues or questions:
1. Check the troubleshooting section
2. Review error messages in the application
3. Check the database status page for system health
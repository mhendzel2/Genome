# Genomics Data Analysis Platform (New Design)

This is a Streamlit-based web application for genomics data analysis, redesigned from the ground up for stability and maintainability. It provides tools for searching public genomics databases, performing basic analysis, and visualizing the results.

## Features

- **Dataset Discovery**: Search and browse public genomics databases (ENCODE).
- **Data Download**: Download datasets directly from the ENCODE database.
- **Analysis Tools**:
  - Basic statistics calculation.
  - Gene set enrichment analysis.
  - Chromatin interaction analysis (insulation score).
  - Basic multi-omics integration (correlation).
- **Interactive Visualizations**:
  - Volcano plots.
  - Heatmaps.

## Setup and Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ```

2.  **Install system dependencies (if not already present):**
    This application requires a Fortran compiler and the OpenBLAS library. On Debian-based systems, you can install them with:
    ```bash
    sudo apt-get update && sudo apt-get install -y gfortran libopenblas-dev
    ```

3.  **Install Python dependencies:**
    It is recommended to use a virtual environment.
    (e.g., `python -m venv venv` and then `source venv/bin/activate`)
    ```bash
    pip install -r requirements.txt
    ```

4.  **Run the application:**
    ```bash
    streamlit run app.py
    ```

## Usage

-   **Dataset Discovery**: Use the search bar on the "Dataset Discovery" page to find datasets from ENCODE.
-   **Data Download**: Click the "Download" button on a search result to download the dataset.
-   **Analysis**: Go to the "Analysis" page, upload a file, and select an analysis to run.
-   **Visualization**: After running an analysis, go to the "Visualization" page to create plots from the results.

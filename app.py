import streamlit as st
from utils.dataset_discovery import DatasetDiscovery
from utils.depmap_client import DepMapClient
from utils.visualization import GenomicsVisualizer
from analysis.analyzer import GenomicsAnalyzer
import pandas as pd
import numpy as np
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    filename='app.log',
                    filemode='a')

def dataset_discovery_page(dataset_discovery, depmap_client):
    st.header("🔍 Dataset Discovery")

    source = st.selectbox("Select Data Source", ["ENCODE", "DepMap"])

    if source == "ENCODE":
        try:
            search_term = st.text_input("Search ENCODE", placeholder="e.g., H3K4me3, GSE12345")

            if st.button("Search"):
                with st.spinner("Searching..."):
                    results = dataset_discovery.search_datasets(search_term=search_term)
                    if results:
                        st.success(f"Found {len(results)} datasets.")
                        st.session_state['search_results'] = results
                        for result in results:
                            with st.expander(f"{result['accession']}: {result['title']}"):
                                st.write(f"**Description:** {result['description']}")
                                st.write(f"**Data Type:** {result['data_type']}")
                                st.write(f"**Tissue:** {result['tissue']}")
                                st.write(f"**Organism:** {result['organism']}")
                                if st.button("Download", key=f"download_{result['accession']}"):
                                    with st.spinner(f"Downloading {result['accession']}..."):
                                        try:
                                            download_path = dataset_discovery.download_dataset(result['accession'])
                                            st.success(f"Dataset downloaded to `{download_path}`")
                                        except Exception as e:
                                            st.error(f"Download failed: {e}")
                                            logging.error(f"Download failed for {result['accession']}: {e}", exc_info=True)
                    else:
                        st.warning("No datasets found.")
        except Exception as e:
            st.error("An error occurred on the Dataset Discovery page.")
            logging.error(f"Error on Dataset Discovery page: {e}", exc_info=True)

    elif source == "DepMap":
        try:
            st.subheader("DepMap Data Files (Release 24Q4)")
            files = depmap_client.list_files_in_release()
            if files:
                for f in files:
                    with st.expander(f"{f['name']}"):
                        st.write(f"**File Size:** {f['size'] / 1e6:.2f} MB")
                        if st.button("Download", key=f"download_{f['id']}"):
                            with st.spinner(f"Downloading {f['name']}..."):
                                download_dir = "downloads/depmap"
                                os.makedirs(download_dir, exist_ok=True)
                                local_filename = os.path.join(download_dir, f['name'])
                                success = depmap_client.download_file(f['download_url'], local_filename)
                                if success:
                                    st.success(f"File downloaded to `{local_filename}`")
                                else:
                                    st.error("Download failed.")
            else:
                st.warning("Could not retrieve file list from DepMap/Figshare.")
        except Exception as e:
            st.error("An error occurred while fetching DepMap data.")
            logging.error(f"Error on DepMap section: {e}", exc_info=True)

def analysis_page(analyzer):
    st.header("🔬 Analysis")
    try:
        uploaded_file = st.file_uploader("Upload a file for analysis")

        if uploaded_file:
            analysis_type = st.selectbox("Choose Analysis", ["Basic Statistics", "Enrichment Analysis"])
            if st.button("Run Analysis"):
                with st.spinner(f"Running {analysis_type}..."):
                    if analysis_type == "Basic Statistics":
                        st.session_state['analysis_results'] = analyzer.basic_statistics(uploaded_file, "BED")
                    elif analysis_type == "Enrichment Analysis":
                        # Wrap the uploaded file in the structure expected by the analyzer
                        wrapped_file = {
                            uploaded_file.name: {
                                'file': uploaded_file,
                                'type': 'Gene Expression'
                            }
                        }
                        st.session_state['analysis_results'] = analyzer.enrichment_analysis(wrapped_file)

                st.success("Analysis complete!")
                st.write(st.session_state['analysis_results'])
    except Exception as e:
        st.error("An error occurred on the Analysis page.")
        logging.error(f"Error on Analysis page: {e}", exc_info=True)

def visualization_page(visualizer):
    st.header("📈 Visualization")
    try:
        if 'analysis_results' in st.session_state:
            results_df = st.session_state['analysis_results']

            if isinstance(results_df, pd.DataFrame) and not results_df.empty:
                # Check if the results are from enrichment analysis
                if 'p_value' in results_df.columns and 'name' in results_df.columns:
                    st.subheader("Enrichment Analysis Results")
                    st.dataframe(results_df) # Display the full results table
                    if st.button("Generate Enrichment Plot"):
                        with st.spinner("Generating plot..."):
                            fig = visualizer.create_enrichment_bar_chart(results_df)
                            st.plotly_chart(fig)
                else:
                    st.subheader("General Data Visualization")
                    st.dataframe(results_df)
                    if st.button("Generate Heatmap"):
                        with st.spinner("Generating heatmap..."):
                            # Create a heatmap of numeric columns
                            numeric_df = results_df.select_dtypes(include='number')
                            if not numeric_df.empty:
                                fig = visualizer.create_heatmap(numeric_df)
                                st.plotly_chart(fig)
                            else:
                                st.warning("No numeric data available for a heatmap.")
            else:
                st.info("Run an analysis to generate results to visualize.")
        else:
            st.warning("No analysis results to visualize. Please run an analysis first.")
    except Exception as e:
        st.error("An error occurred on the Visualization page.")
        logging.error(f"Error on Visualization page: {e}", exc_info=True)

def main():
    st.set_page_config(
        page_title="Genomics Data Analysis Platform",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 Genomics Data Analysis Platform (New Design)")
    st.sidebar.title("Navigation")
    
    dataset_discovery = DatasetDiscovery()
    depmap_client = DepMapClient()
    visualizer = GenomicsVisualizer()
    analyzer = GenomicsAnalyzer()
    
    page = st.sidebar.selectbox(
        "Choose a section",
        ["Home", "Dataset Discovery", "Analysis", "Visualization"]
    )
    
    if page == "Home":
        st.write("Welcome to the new and improved Genomics Data Analysis Platform!")
    elif page == "Dataset Discovery":
        dataset_discovery_page(dataset_discovery, depmap_client)
    elif page == "Analysis":
        analysis_page(analyzer)
    elif page == "Visualization":
        visualization_page(visualizer)

if __name__ == "__main__":
    main()

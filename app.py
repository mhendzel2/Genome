import streamlit as st
from utils.dataset_discovery import DatasetDiscovery
from utils.depmap_client import DepMapClient
from utils.geo_client import GEOClient
from utils.visualization import GenomicsVisualizer
from analysis.analyzer_adapter import AnalyzerAdapter
from utils.r_integration import RIntegration
import pandas as pd
import numpy as np
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    filename='app.log',
                    filemode='a')

def dataset_discovery_page(dataset_discovery, depmap_client, geo_client):
    st.header("🔍 Dataset Discovery")

    source = st.selectbox("Select Data Source", ["ENCODE", "GEO", "DepMap"])

    if source == "ENCODE":
        try:
            search_term = st.text_input("Search ENCODE", placeholder="e.g., H4K5ac, ENCSR675FLJ, ChIP-seq")
            
            col1, col2 = st.columns(2)
            with col1:
                organism = st.selectbox("Organism", ["Homo sapiens", "Mus musculus", ""], index=0)
            with col2:
                use_semantic = st.checkbox("Use semantic expansion", value=False)

            if st.button("Search ENCODE"):
                with st.spinner("Searching ENCODE..."):
                    results = dataset_discovery.search_datasets(
                        search_term=search_term,
                        organism=organism if organism else None,
                        use_semantic_search=use_semantic
                    )
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

    elif source == "GEO":
        try:
            search_term = st.text_input("Search GEO", placeholder="e.g., GSE135771, KANSL1, MSL1 ChIP-seq")
            
            col1, col2 = st.columns(2)
            with col1:
                max_results = st.slider("Max results", 5, 50, 20)
            with col2:
                dataset_type = st.selectbox("Dataset type", ["gse", "gds", "gpl"], index=0)

            if st.button("Search GEO"):
                with st.spinner("Searching GEO..."):
                    results = geo_client.search_datasets(
                        search_term=search_term,
                        dataset_type=dataset_type,
                        max_results=max_results
                    )
                    if results:
                        st.success(f"Found {len(results)} datasets.")
                        st.session_state['search_results'] = results
                        st.session_state['geo_results'] = results
                        for result in results:
                            with st.expander(f"{result['accession']}: {result['title']}"):
                                st.write(f"**Description:** {result.get('description', 'N/A')[:500]}")
                                st.write(f"**Organism:** {result.get('organism', 'N/A')}")
                                st.write(f"**Samples:** {result.get('n_samples', 'N/A')}")
                                st.write(f"**Platform:** {result.get('platform', 'N/A')}")
                                st.write(f"**Submission Date:** {result.get('submission_date', 'N/A')}")
                                if result.get('pubmed_id'):
                                    st.write(f"**PubMed ID:** {result['pubmed_id']}")
                                
                                # Download buttons
                                dl_col1, dl_col2 = st.columns(2)
                                with dl_col1:
                                    if st.button("Download Matrix", key=f"matrix_{result['accession']}"):
                                        with st.spinner(f"Downloading series matrix for {result['accession']}..."):
                                            try:
                                                files = geo_client.download_series_matrix(result['accession'])
                                                st.success(f"Downloaded {len(files)} file(s) to downloads/geo/")
                                                for f in files:
                                                    st.write(f"  - `{f}`")
                                            except Exception as e:
                                                st.error(f"Download failed: {e}")
                                                logging.error(f"GEO matrix download failed: {e}", exc_info=True)
                                with dl_col2:
                                    if st.button("Download Supplementary", key=f"suppl_{result['accession']}"):
                                        with st.spinner(f"Downloading supplementary files for {result['accession']}..."):
                                            try:
                                                files = geo_client.download_supplementary_files(result['accession'])
                                                if files:
                                                    st.success(f"Downloaded {len(files)} file(s) to downloads/geo/")
                                                    for f in files:
                                                        st.write(f"  - `{f}`")
                                                else:
                                                    st.warning("No supplementary files found.")
                                            except Exception as e:
                                                st.error(f"Download failed: {e}")
                                                logging.error(f"GEO supplementary download failed: {e}", exc_info=True)
                    else:
                        st.warning("No datasets found.")
        except Exception as e:
            st.error("An error occurred while searching GEO.")
            logging.error(f"Error on GEO search: {e}", exc_info=True)

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
        analysis_type = st.selectbox("Choose Analysis", ["Basic Statistics", "Enrichment Analysis", "Differential Expression", "Quality Control"])

        if analysis_type == "Differential Expression":
            st.subheader("Differential Expression Analysis")
            col1, col2 = st.columns(2)
            with col1:
                control_files = st.file_uploader("Upload Control Group Files", accept_multiple_files=True)
            with col2:
                treatment_files = st.file_uploader("Upload Treatment Group Files", accept_multiple_files=True)

            if st.button("Run Analysis"):
                if control_files and treatment_files:
                    with st.spinner("Running Differential Expression Analysis..."):
                        st.session_state['uploaded_data'] = {'control': control_files, 'treatment': treatment_files}
                        result = analyzer.differential_expression(
                            st.session_state['uploaded_data']
                        )
                        
                        # Extract actual result from adapter wrapper
                        if isinstance(result, dict) and 'result' in result:
                            st.session_state['analysis_results'] = result['result']
                        else:
                            st.session_state['analysis_results'] = result
                        
                        st.success("Analysis complete!")
                        st.write(st.session_state['analysis_results'])
                else:
                    st.warning("Please upload files for both control and treatment groups.")

        elif analysis_type == "Quality Control":
            st.subheader("Quality Control")
            uploaded_files = st.file_uploader("Upload files for QC", accept_multiple_files=True)
            file_type = st.selectbox("Select File Type", ["Gene Expression", "ChIP-seq", "Other"])

            if st.button("Run QC"):
                if uploaded_files:
                    with st.spinner("Running Quality Control..."):
                        st.session_state['uploaded_data'] = uploaded_files
                        wrapped_files = {
                            f.name: {'file': f, 'type': file_type} for f in uploaded_files
                        }
                        result = analyzer.quality_control(wrapped_files)
                        
                        # Extract actual result from adapter wrapper
                        if isinstance(result, dict) and 'result' in result:
                            qc_results = result['result']
                        else:
                            qc_results = result
                        
                        if isinstance(qc_results, dict):
                            st.session_state['analysis_results'] = pd.DataFrame.from_dict(qc_results, orient='index')
                        else:
                            st.session_state['analysis_results'] = qc_results
                        st.success("QC complete!")
                        st.write(st.session_state['analysis_results'])
                else:
                    st.warning("Please upload at least one file.")

        else:
            uploaded_file = st.file_uploader("Upload a file for analysis")
            if uploaded_file:
                if st.button("Run Analysis"):
                    with st.spinner(f"Running {analysis_type}..."):
                        st.session_state['uploaded_data'] = [uploaded_file]
                        if analysis_type == "Basic Statistics":
                            result = analyzer.basic_statistics(
                                {uploaded_file.name: {'file': uploaded_file, 'type': 'BED'}}
                            )
                        elif analysis_type == "Enrichment Analysis":
                            wrapped_file = {
                                uploaded_file.name: {
                                    'file': uploaded_file,
                                    'type': 'Gene Expression'
                                }
                            }
                            result = analyzer.enrichment_analysis(wrapped_file)
                        
                        # Extract actual result from adapter wrapper
                        if isinstance(result, dict) and 'result' in result:
                            st.session_state['analysis_results'] = result['result']
                        else:
                            st.session_state['analysis_results'] = result
                    
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
                            numeric_df = results_df.select_dtypes(include='number')
                            if not numeric_df.empty:
                                fig = visualizer.create_heatmap(numeric_df)
                                st.plotly_chart(fig)
                            else:
                                st.warning("No numeric data available for a heatmap.")

                if st.button("Generate PCA Plot"):
                    with st.spinner("Generating PCA plot..."):
                        if 'uploaded_data' in st.session_state:
                            uploaded_data = st.session_state['uploaded_data']
                            if isinstance(uploaded_data, dict):
                                files_to_plot = [f for group in uploaded_data.values() for f in group]
                            else:
                                files_to_plot = uploaded_data

                            fig = visualizer.create_pca_plot(files_to_plot)
                            st.plotly_chart(fig)
                        else:
                            st.warning("Uploaded data not found. Please run an analysis first.")
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
    geo_client = GEOClient()
    visualizer = GenomicsVisualizer()
    r_integration = RIntegration()
    analyzer = AnalyzerAdapter(r_integration)
    
    page = st.sidebar.selectbox(
        "Choose a section",
        ["Home", "Dataset Discovery", "Analysis", "Visualization"]
    )
    
    if page == "Home":
        st.write("Welcome to the new and improved Genomics Data Analysis Platform!")
    elif page == "Dataset Discovery":
        dataset_discovery_page(dataset_discovery, depmap_client, geo_client)
    elif page == "Analysis":
        analysis_page(analyzer)
    elif page == "Visualization":
        visualization_page(visualizer)

if __name__ == "__main__":
    main()

import streamlit as st
import streamlit.components.v1 as components
from utils.dataset_discovery import DatasetDiscovery
from utils.depmap_client import DepMapClient
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
        analysis_options = ["Basic Statistics", "Enrichment Analysis", "Differential Expression", "Quality Control", "Multi-omics Integration", "Chromatin Interaction"]
        analysis_type = st.selectbox("Choose Analysis", analysis_options)

        if analysis_type == "Chromatin Interaction":
            st.subheader("Chromatin Interaction Analysis (BEDPE)")
            uploaded_file = st.file_uploader("Upload Hi-C Data File (.bedpe)", accept_multiple_files=False)
            if st.button("Run Analysis"):
                if uploaded_file:
                    with st.spinner("Parsing Hi-C data..."):
                        wrapped_file = {uploaded_file.name: {'file': uploaded_file, 'type': 'HiC'}}
                        st.session_state['uploaded_data'] = wrapped_file
                        results = analyzer.chromatin_interaction_analysis(wrapped_file)
                        st.session_state['analysis_results'] = results
                        st.session_state['analysis_type'] = "Hi-C"
                        st.success("Analysis complete!")
                        st.write(f"Parsed {len(results)} interactions.")
                        st.dataframe(results.head())
                else:
                    st.warning("Please upload a file for analysis.")

        elif analysis_type == "Multi-omics Integration":
            st.subheader("Multi-omics Integration Analysis")
            uploaded_files = st.file_uploader("Upload Omics Files (at least 2)", accept_multiple_files=True)
            corr_threshold = st.slider("Correlation Threshold", 0.5, 1.0, 0.7, 0.05)

            if st.button("Run Analysis"):
                if uploaded_files and len(uploaded_files) >= 2:
                    with st.spinner("Running Multi-omics Integration..."):
                        wrapped_files = {
                            f.name: {'file': f, 'type': 'Omics'} for f in uploaded_files
                        }
                        st.session_state['uploaded_data'] = wrapped_files
                        results = analyzer.multi_omics_integration(wrapped_files, correlation_threshold=corr_threshold)
                        st.session_state['analysis_results'] = results
                        st.session_state['analysis_type'] = "Multi-omics"
                        st.success("Analysis complete!")
                        st.write(f"Found {len(results)} significant correlations.")
                        st.dataframe(results)
                else:
                    st.warning("Please upload at least two files for integration.")

        elif analysis_type == "Differential Expression":
            st.subheader("Differential Expression Analysis")
            col1, col2 = st.columns(2)
            with col1:
                control_files = st.file_uploader("Upload Control Group Files", accept_multiple_files=True)
            with col2:
                treatment_files = st.file_uploader("Upload Treatment Group Files", accept_multiple_files=True)

            if st.button("Run Analysis"):
                if control_files and treatment_files:
                    with st.spinner("Running Differential Expression Analysis..."):
                        # For simplicity, we'll just use the first control file for this example
                        # A real implementation would handle both groups properly in the R script
                        st.session_state['uploaded_data'] = {'control': control_files, 'treatment': treatment_files}
                        wrapped_files = {f.name: {'file': f, 'type': 'Gene Expression'} for f in control_files}

                        results = analyzer.differential_expression(wrapped_files)
                        st.session_state['analysis_results'] = results
                        st.session_state['analysis_type'] = "DE"
                        st.success("Analysis complete!")
                        st.dataframe(results.head())
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
                        results = analyzer.quality_control(wrapped_files)
                        st.session_state['analysis_results'] = pd.DataFrame.from_dict(results, orient='index')
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
                            st.session_state['analysis_results'] = analyzer.basic_statistics(
                                {uploaded_file.name: {'file': uploaded_file, 'type': 'BED'}}
                            )
                        elif analysis_type == "Enrichment Analysis":
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
                analysis_type = st.session_state.get('analysis_type')

                if analysis_type == "Hi-C":
                    st.subheader("Chromatin Interaction (Circos Plot)")
                    st.dataframe(results_df.head())
                    if st.button("Generate Circos Plot"):
                        with st.spinner("Generating Circos plot... This may take a moment."):
                            circos_path = visualizer.create_circos_plot(results_df)
                            if circos_path:
                                st.image(circos_path)
                                os.remove(circos_path) # Clean up temp file
                            else:
                                st.warning("Could not generate Circos plot.")

                elif analysis_type == "Multi-omics":
                    st.subheader("Multi-omics Correlation Network")
                    st.dataframe(results_df)
                    if st.button("Generate Network Graph"):
                        with st.spinner("Generating network graph..."):
                            network_html_path = visualizer.create_network_graph(results_df)
                            if network_html_path:
                                with open(network_html_path, 'r', encoding='utf-8') as f:
                                    html_code = f.read()
                                components.html(html_code, height=800)
                                os.remove(network_html_path) # Clean up the temp file
                            else:
                                st.warning("Could not generate network graph.")

                # Check if the results are from differential expression
                elif analysis_type == "DE" and 'log2FoldChange' in results_df.columns and 'pvalue' in results_df.columns:
                    st.subheader("Differential Expression Results")
                    st.dataframe(results_df)
                    lfc_thresh = st.slider("Log Fold Change Threshold", 0.0, 5.0, 1.0, 0.1)
                    pval_thresh = st.slider("P-value Threshold", 0.0, 0.1, 0.05, 0.001)
                    if st.button("Generate Volcano Plot"):
                        with st.spinner("Generating volcano plot..."):
                            fig = visualizer.create_volcano_plot(results_df, lfc_threshold=lfc_thresh, p_value_threshold=pval_thresh)
                            st.plotly_chart(fig)

                # Check if the results are from enrichment analysis
                elif 'pvalue' in results_df.columns and 'name' in results_df.columns:
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

                st.subheader("Principal Component Analysis (PCA)")
                pca_dims = st.radio("Select PCA Dimensions", (2, 3), horizontal=True)
                if st.button("Generate PCA Plot"):
                    with st.spinner(f"Generating {pca_dims}D PCA plot..."):
                        if 'uploaded_data' in st.session_state:
                            uploaded_data = st.session_state['uploaded_data']
                            if isinstance(uploaded_data, dict): # Handle dict of lists of files
                                files_to_plot = [f for group in uploaded_data.values() for f in group]
                            else: # Handle list of files
                                files_to_plot = uploaded_data

                            fig = visualizer.create_pca_plot(files_to_plot, dimensions=pca_dims)
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
        dataset_discovery_page(dataset_discovery, depmap_client)
    elif page == "Analysis":
        analysis_page(analyzer)
    elif page == "Visualization":
        visualization_page(visualizer)

if __name__ == "__main__":
    main()

import streamlit as st
from utils.dataset_discovery import DatasetDiscovery
from utils.encode_client import ENCODEClient
from utils.visualization import GenomicsVisualizer
from analysis.analyzer import GenomicsAnalyzer
import pandas as pd
import logging

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    filename='app.log',
                    filemode='a')

def dataset_discovery_page(dataset_discovery):
    st.header("🔍 Dataset Discovery")
    try:
        db_options = list(dataset_discovery.clients.keys())
        selected_dbs = st.multiselect("Select databases to search", db_options, default=db_options)

        search_term = st.text_input("Search", placeholder="e.g., H3K4me3, cancer")

        if st.button("Search"):
            if not selected_dbs:
                st.warning("Please select at least one database to search.")
                return

            with st.spinner("Searching..."):
                results = dataset_discovery.search(databases=selected_dbs, search_term=search_term)
                if results:
                    st.success(f"Found {len(results)} datasets across {len(selected_dbs)} database(s).")
                    st.session_state['search_results'] = results
                    for result in results:
                        with st.expander(f"[{result['database']}] {result['accession']}: {result['title']}"):
                            st.write(f"**Description:** {result['description']}")
                            st.write(f"**Data Type:** {result['data_type']}")
                            st.write(f"**Tissue:** {result['tissue']}")
                            st.write(f"**Organism:** {result['organism']}")
                            if st.button("Download", key=f"download_{result['database']}_{result['accession']}"):
                                with st.spinner(f"Downloading {result['accession']} from {result['database']}..."):
                                    try:
                                        download_path = dataset_discovery.download(result['database'], result['accession'])
                                        st.success(f"Dataset downloaded to `{download_path}`")
                                    except Exception as e:
                                        st.error(f"Download failed: {e}")
                                        logging.error(f"Download failed for {result['accession']}: {e}", exc_info=True)
                else:
                    st.warning("No datasets found.")
    except Exception as e:
        st.error("An error occurred on the Dataset Discovery page.")
        logging.error(f"Error on Dataset Discovery page: {e}", exc_info=True)

def analysis_page(analyzer):
    st.header("🔬 Analysis")
    try:
        analysis_type = st.selectbox("Choose Analysis Type", ["Single File Analysis", "Comparative Analysis"])

        if analysis_type == "Single File Analysis":
            uploaded_file = st.file_uploader("Upload a file for analysis")
            if uploaded_file:
                file_type = st.selectbox("Select File Type", ["Gene Expression", "BED", "Proteomics"])
                analysis = st.selectbox("Choose Analysis", ["Basic Statistics", "Enrichment Analysis"])
                if st.button("Run Analysis"):
                    with st.spinner(f"Running {analysis}..."):
                        if analysis == "Basic Statistics":
                            st.session_state['analysis_results'] = analyzer.basic_statistics(uploaded_file, file_type)
                        elif analysis == "Enrichment Analysis":
                            st.session_state['analysis_results'] = analyzer.enrichment_analysis(uploaded_file)
                    st.success("Analysis complete!")
                    st.write(st.session_state['analysis_results'])

        elif analysis_type == "Comparative Analysis":
            st.subheader("Gene Expression vs. Proteomics")
            gene_expr_file = st.file_uploader("Upload Gene Expression File")
            proteomics_file = st.file_uploader("Upload Proteomics File")
            if gene_expr_file and proteomics_file:
                if st.button("Run Comparative Analysis"):
                    with st.spinner("Running comparative analysis..."):
                        st.session_state['analysis_results'] = analyzer.comparative_analysis(gene_expr_file, proteomics_file)
                    st.success("Analysis complete!")
                    st.write(st.session_state['analysis_results'])

    except Exception as e:
        st.error("An error occurred on the Analysis page.")
        logging.error(f"Error on Analysis page: {e}", exc_info=True)

def visualization_page(visualizer):
    st.header("📈 Visualization")
    try:
        if 'analysis_results' in st.session_state:
            results = st.session_state['analysis_results']
            if isinstance(results, list) and results:
                df = pd.DataFrame(results)

                # Check if the results are from enrichment analysis
                if 'Term' in df.columns and 'Combined Score' in df.columns:
                    st.subheader("Enrichment Analysis Results")
                    if st.button("Generate Enrichment Plot"):
                        with st.spinner("Generating plot..."):
                            fig = visualizer.create_enrichment_bar_chart(df)
                            st.plotly_chart(fig)

                # Check if the results are from comparative analysis
                elif 'correlation_coefficient' in results:
                    st.subheader("Comparative Analysis Results")
                    st.write(f"Correlation Coefficient: {results['correlation_coefficient']:.4f}")
                    st.write(f"P-value: {results['p_value']:.4f}")
                    if st.button("Generate Correlation Plot"):
                        with st.spinner("Generating plot..."):
                            df = results['dataframe']
                            fig = visualizer.create_correlation_scatter_plot(df, df.columns[1], df.columns[2])
                            st.plotly_chart(fig)
                else:
                    st.subheader("General Data Visualization")
                    if st.button("Generate Heatmap"):
                        with st.spinner("Generating heatmap..."):
                            # Create a heatmap of numeric columns
                            numeric_df = df.select_dtypes(include=np.number)
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
    visualizer = GenomicsVisualizer()
    analyzer = GenomicsAnalyzer()
    
    page = st.sidebar.selectbox(
        "Choose a section",
        ["Home", "Dataset Discovery", "Analysis", "Visualization"]
    )
    
    if page == "Home":
        st.write("Welcome to the new and improved Genomics Data Analysis Platform!")
    elif page == "Dataset Discovery":
        dataset_discovery_page(dataset_discovery)
    elif page == "Analysis":
        analysis_page(analyzer)
    elif page == "Visualization":
        visualization_page(visualizer)

if __name__ == "__main__":
    main()

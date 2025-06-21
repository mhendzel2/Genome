import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO
import os
import sys

# Add utils to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from utils.data_validation import DataValidator
from utils.r_integration import RIntegration
from utils.visualization import GenomicsVisualizer
from utils.dataset_discovery import DatasetDiscovery
from utils.database import DatabaseManager
from analysis.genomics_analysis import GenomicsAnalyzer

# Initialize session state
if 'uploaded_files' not in st.session_state:
    st.session_state.uploaded_files = {}
if 'analysis_results' not in st.session_state:
    st.session_state.analysis_results = {}
if 'selected_datasets' not in st.session_state:
    st.session_state.selected_datasets = []
if 'session_id' not in st.session_state:
    import uuid
    st.session_state.session_id = str(uuid.uuid4())

@st.cache_resource
def get_db_manager():
    """Initialize database manager with caching to prevent repeated initialization"""
    try:
        db_manager = DatabaseManager()
        # Populate only if the database is truly empty
        if db_manager.get_database_stats()['total_datasets'] == 0:
            count = db_manager.populate_sample_datasets()
            st.toast(f"Database initialized with {count} sample datasets.")
        return db_manager
    except Exception as e:
        st.error(f"Database connection failed: {e}")
        return None

def main():
    st.set_page_config(
        page_title="Genomics Data Analysis Platform",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 Genomics Data Analysis Platform")
    st.sidebar.title("Navigation")
    
    # Initialize components
    validator = DataValidator()
    r_integration = RIntegration()
    visualizer = GenomicsVisualizer()
    dataset_discovery = DatasetDiscovery()
    analyzer = GenomicsAnalyzer(r_integration)
    
    # Initialize database with caching
    db_manager = get_db_manager()
    if not db_manager:
        st.stop()  # Stop execution if DB is not available
    
    # Sidebar navigation
    page = st.sidebar.selectbox(
        "Choose Analysis Module",
        ["Data Upload & Validation", "Dataset Discovery", "Data Analysis", "Tissue Comparison", "Visualization", "Database Status", "Export Results"]
    )
    
    if page == "Data Upload & Validation":
        data_upload_page(validator, db_manager)
    elif page == "Dataset Discovery":
        dataset_discovery_page(dataset_discovery, db_manager)
    elif page == "Data Analysis":
        data_analysis_page(analyzer, db_manager)
    elif page == "Tissue Comparison":
        tissue_comparison_page(analyzer, visualizer)
    elif page == "Visualization":
        visualization_page(visualizer)
    elif page == "Database Status":
        database_status_page(db_manager)
    elif page == "Export Results":
        export_results_page()

def database_status_page(db_manager):
    st.header("🗄️ Database Status")
    
    if not db_manager:
        st.error("Database is not available")
        return
    
    try:
        # Get database statistics
        stats = db_manager.get_database_stats()
        
        # Display overview
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Datasets", stats['total_datasets'])
        with col2:
            st.metric("Total Analyses", stats['total_analyses'])
        with col3:
            st.metric("Uploaded Files", stats['total_uploaded_files'])
        
        # Display breakdowns
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Datasets by Organism")
            if stats['datasets_by_organism']:
                organism_df = pd.DataFrame(list(stats['datasets_by_organism'].items()), 
                                         columns=['Organism', 'Count'])
                st.dataframe(organism_df)
            
            st.subheader("Datasets by Data Type")
            if stats['datasets_by_data_type']:
                datatype_df = pd.DataFrame(list(stats['datasets_by_data_type'].items()), 
                                         columns=['Data Type', 'Count'])
                st.dataframe(datatype_df)
        
        with col2:
            st.subheader("Datasets by Database")
            if stats['datasets_by_database']:
                database_df = pd.DataFrame(list(stats['datasets_by_database'].items()), 
                                         columns=['Database', 'Count'])
                st.dataframe(database_df)
        
        # Analysis history
        if st.checkbox("Show Analysis History"):
            history = db_manager.get_analysis_history(st.session_state.session_id)
            if history:
                st.subheader("Your Analysis History")
                for analysis in history[:10]:  # Show last 10
                    with st.expander(f"{analysis['analysis_type']} - {analysis['created_at']}"):
                        st.json(analysis['results'])
            else:
                st.info("No analysis history found")
        
        # Database management
        st.subheader("Database Management")
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Refresh Sample Data"):
                count = db_manager.populate_sample_datasets()
                st.success(f"Added {count} new sample datasets")
        with col2:
            if st.button("Export Database"):
                st.info("Database export functionality would be implemented here")
                
    except Exception as e:
        st.error(f"Error accessing database: {e}")

def data_upload_page(validator, db_manager):
    st.header("📁 Data Upload & Validation")
    
    # Data type selection
    data_types = ["Histone Marks", "HiC", "Gene Expression", "ChIP-seq"]
    selected_type = st.selectbox("Select Data Type", data_types)
    
    # File upload
    uploaded_file = st.file_uploader(
        f"Upload {selected_type} data",
        type=['bed', 'bedgraph', 'txt', 'csv', 'tsv', 'gz', 'bigwig', 'wig'],
        help=f"Supported formats for {selected_type}: BED, bedGraph, TXT, CSV, TSV, GZ, BigWig, WIG"
    )
    
    if uploaded_file is not None:
        # Validate file
        with st.spinner("Validating file..."):
            validation_result = validator.validate_file(uploaded_file, selected_type)
        
        if validation_result['valid']:
            st.success(f"✅ File validated successfully as {selected_type} data")
            st.info(f"**File Info:** {validation_result['info']}")
            
            # Store in session state
            st.session_state.uploaded_files[uploaded_file.name] = {
                'file': uploaded_file,
                'type': selected_type,
                'validation': validation_result
            }
            
            # Display preview
            if st.checkbox("Show data preview"):
                try:
                    df = validator.load_as_dataframe(uploaded_file, selected_type)
                    st.dataframe(df.head(10))
                    st.info(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
                except Exception as e:
                    st.error(f"Error loading preview: {str(e)}")
        else:
            st.error(f"❌ Validation failed: {validation_result['error']}")
    
    # Display uploaded files
    if st.session_state.uploaded_files:
        st.subheader("Uploaded Files")
        for filename, file_info in st.session_state.uploaded_files.items():
            col1, col2, col3 = st.columns([3, 2, 1])
            with col1:
                st.text(filename)
            with col2:
                st.text(file_info['type'])
            with col3:
                if st.button("Remove", key=f"remove_{filename}"):
                    del st.session_state.uploaded_files[filename]
                    st.rerun()

def dataset_discovery_page(dataset_discovery, db_manager):
    st.header("🔍 Dataset Discovery")
    
    # Search parameters
    col1, col2, col3 = st.columns(3)
    with col1:
        search_term = st.text_input("Search Term", placeholder="e.g., H3K4me3, GSE12345, ENCODE")
        data_types_filter = st.multiselect(
            "Data Types", 
            ["Histone Marks", "HiC", "Gene Expression", "ChIP-seq"],
            default=[]
        )
    
    with col2:
        tissues_filter = st.multiselect(
            "Tissues",
            ["Liver", "Heart", "Brain", "Lung", "Kidney", "Muscle", "Skin", "Blood", "Spleen", "Pancreas", "Intestine", "Stomach", "Thymus", "Bone Marrow"],
            default=[]
        )
        organism = st.selectbox("Organism", ["Human", "Mouse"], index=0)
    
    with col3:
        cell_types_filter = st.multiselect(
            "Cell Types",
            ["Any", "Primary cells", "Cell lines", "Stem cells", "Immune cells"],
            default=["Any"]
        )
        # Add action buttons
        col_btn1, col_btn2, col_btn3 = st.columns(3)
        with col_btn1:
            search_button = st.button("Search Datasets", type="primary")
        with col_btn2:
            browse_button = st.button("Browse Popular")
        with col_btn3:
            show_all_button = st.button("Show All Datasets")
    
    # Handle clear filters
    if 'clear_filters' in st.session_state and st.session_state.clear_filters:
        st.session_state.clear_filters = False
        st.rerun()
    
    if search_button or browse_button or show_all_button:
        with st.spinner("Searching public databases..."):
            if browse_button:
                # Get popular datasets from database
                if db_manager:
                    results = db_manager.search_datasets("", None, None, None, organism, 20)
                else:
                    results = dataset_discovery.get_popular_datasets(20)
            elif show_all_button:
                # Show all datasets with current filters from database
                if db_manager:
                    results = db_manager.search_datasets("", data_types_filter, tissues_filter, cell_types_filter, organism)
                else:
                    results = dataset_discovery.search_datasets("", data_types_filter, tissues_filter, cell_types_filter, organism)
            else:
                # Regular search with filters from database
                if db_manager:
                    results = db_manager.search_datasets(search_term, data_types_filter, tissues_filter, cell_types_filter, organism)
                else:
                    results = dataset_discovery.search_datasets(search_term, data_types_filter, tissues_filter, cell_types_filter, organism)
        
        if results:
            st.success(f"Found {len(results)} datasets")
            
            # Display results
            for i, dataset in enumerate(results):
                relevance_score = dataset.get('relevance_score', 0)
                score_indicator = "🔥" if relevance_score > 20 else "⭐" if relevance_score > 10 else ""
                
                with st.expander(f"{score_indicator} {dataset['title']} - {dataset['database']}"):
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        st.write(f"**Description:** {dataset['description']}")
                        st.write(f"**Data Type:** {dataset['data_type']}")
                        st.write(f"**Tissue:** {dataset['tissue']}")
                        st.write(f"**Organism:** {dataset['organism']}")
                        st.write(f"**Accession:** {dataset['accession']}")
                        
                        # Show cell type information
                        if dataset.get('cell_type'):
                            st.write(f"**Cell Type:** {dataset['cell_type']}")
                        if dataset.get('cell_line'):
                            st.write(f"**Cell Line:** {dataset['cell_line']}")
                        
                        # Show file information
                        st.write(f"**Files:** {dataset['file_count']} files ({dataset['file_size_gb']} GB)")
                        st.write(f"**Published:** {dataset['publication_date']}")
                        
                        if relevance_score > 0:
                            st.caption(f"Relevance score: {relevance_score}")
                    
                    with col2:
                        if st.button("Select", key=f"select_{i}"):
                            if dataset not in st.session_state.selected_datasets:
                                st.session_state.selected_datasets.append(dataset)
                                st.success("Dataset selected!")
                            else:
                                st.info("Already selected")
        else:
            st.warning("No datasets found matching your criteria")
    
    # Selected datasets
    if st.session_state.selected_datasets:
        st.subheader("Selected Datasets")
        for i, dataset in enumerate(st.session_state.selected_datasets):
            col1, col2, col3 = st.columns([4, 2, 1])
            with col1:
                st.text(dataset['title'])
            with col2:
                st.text(f"{dataset['data_type']} - {dataset['tissue']}")
            with col3:
                if st.button("Remove", key=f"remove_dataset_{i}"):
                    st.session_state.selected_datasets.pop(i)
                    st.rerun()

def data_analysis_page(analyzer, db_manager):
    st.header("📊 Data Analysis")
    
    if not st.session_state.uploaded_files and not st.session_state.selected_datasets:
        st.warning("Please upload files or select datasets first")
        return
    
    # Analysis type selection
    analysis_types = [
        "Basic Statistics",
        "Peak Calling (ChIP-seq)",
        "Differential Expression",
        "Chromatin Interaction Analysis",
        "Enrichment Analysis"
    ]
    
    selected_analysis = st.selectbox("Select Analysis Type", analysis_types)
    
    # Analysis parameters
    st.subheader("Analysis Parameters")
    
    if selected_analysis == "Basic Statistics":
        if st.button("Run Basic Statistics"):
            with st.spinner("Computing statistics..."):
                results = analyzer.basic_statistics(st.session_state.uploaded_files)
            
            if results:
                st.session_state.analysis_results['basic_stats'] = results
                st.success("Analysis completed!")
                
                for filename, stats in results.items():
                    st.subheader(f"Statistics for {filename}")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Total Regions", stats.get('total_regions', 'N/A'))
                        st.metric("Mean Score", f"{stats.get('mean_score', 0):.2f}")
                    with col2:
                        st.metric("Median Score", f"{stats.get('median_score', 0):.2f}")
                        st.metric("Max Score", f"{stats.get('max_score', 0):.2f}")
    
    elif selected_analysis == "Peak Calling (ChIP-seq)":
        col1, col2 = st.columns(2)
        with col1:
            p_value = st.number_input("P-value threshold", value=0.05, format="%.3f")
            fold_change = st.number_input("Fold change threshold", value=2.0)
        with col2:
            window_size = st.number_input("Window size", value=200, step=50)
            min_length = st.number_input("Minimum peak length", value=50, step=10)
        
        if st.button("Run Peak Calling"):
            with st.spinner("Calling peaks..."):
                results = analyzer.peak_calling(
                    st.session_state.uploaded_files,
                    p_value=p_value,
                    fold_change=fold_change,
                    window_size=window_size,
                    min_length=min_length
                )
            
            if results:
                st.session_state.analysis_results['peak_calling'] = results
                st.success(f"Peak calling completed! Found {len(results)} significant peaks")
                
                # Display peak summary
                if len(results) > 0:
                    df_peaks = pd.DataFrame(results)
                    st.dataframe(df_peaks.head(20))
    
    elif selected_analysis == "Differential Expression":
        # Parameters for differential expression
        col1, col2 = st.columns(2)
        with col1:
            padj_threshold = st.number_input("Adjusted p-value", value=0.05, format="%.3f")
            lfc_threshold = st.number_input("Log2 fold change", value=1.0)
        
        if st.button("Run Differential Expression"):
            with st.spinner("Running DESeq2 analysis..."):
                results = analyzer.differential_expression(
                    st.session_state.uploaded_files,
                    padj_threshold=padj_threshold,
                    lfc_threshold=lfc_threshold
                )
            
            if results:
                st.session_state.analysis_results['diff_expression'] = results
                st.success("Differential expression analysis completed!")
                
                # Display results summary
                st.write(f"**Upregulated genes:** {results.get('upregulated', 0)}")
                st.write(f"**Downregulated genes:** {results.get('downregulated', 0)}")

def tissue_comparison_page(analyzer, visualizer):
    st.header("🔬 Tissue Comparison Analysis")
    
    if not st.session_state.uploaded_files and not st.session_state.selected_datasets:
        st.warning("Please upload files or select datasets first")
        return
    
    # Get available tissues from uploaded data
    available_tissues = []
    for filename, file_info in st.session_state.uploaded_files.items():
        # Extract tissue info from filename or metadata
        if 'liver' in filename.lower():
            available_tissues.append('Liver')
        elif 'heart' in filename.lower():
            available_tissues.append('Heart')
        elif 'brain' in filename.lower():
            available_tissues.append('Brain')
    
    # Add tissues from selected datasets
    for dataset in st.session_state.selected_datasets:
        if dataset['tissue'] not in available_tissues:
            available_tissues.append(dataset['tissue'])
    
    if len(available_tissues) < 2:
        st.warning("Need at least 2 different tissues for comparison")
        available_tissues = ['Liver', 'Heart', 'Brain']  # Default for demo
    
    # Tissue selection for comparison
    tissues_to_compare = st.multiselect(
        "Select tissues to compare",
        available_tissues,
        default=available_tissues[:2] if len(available_tissues) >= 2 else available_tissues
    )
    
    comparison_type = st.selectbox(
        "Comparison Type",
        ["Expression Levels", "Peak Overlap", "Chromatin Accessibility", "Histone Modifications"]
    )
    
    if len(tissues_to_compare) >= 2 and st.button("Run Tissue Comparison"):
        with st.spinner("Performing tissue comparison..."):
            results = analyzer.tissue_comparison(
                tissues_to_compare,
                comparison_type,
                st.session_state.uploaded_files
            )
        
        if results:
            st.session_state.analysis_results['tissue_comparison'] = results
            st.success("Tissue comparison completed!")
            
            # Display comparison results
            st.subheader("Comparison Results")
            
            # Statistical summary
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Common Features", results.get('common_features', 0))
            with col2:
                st.metric("Unique to " + tissues_to_compare[0], results.get('unique_tissue1', 0))
            with col3:
                st.metric("Unique to " + tissues_to_compare[1], results.get('unique_tissue2', 0))
            
            # Visualization
            if st.checkbox("Show comparison plot"):
                fig = visualizer.create_tissue_comparison_plot(results, tissues_to_compare)
                st.plotly_chart(fig, use_container_width=True)

def visualization_page(visualizer):
    st.header("📈 Data Visualization")
    
    if not st.session_state.analysis_results:
        st.warning("No analysis results available. Please run some analyses first.")
        return
    
    # Visualization type selection
    viz_types = [
        "Genome Browser View",
        "Heatmap",
        "Volcano Plot",
        "MA Plot",
        "PCA Plot",
        "Peak Distribution"
    ]
    
    selected_viz = st.selectbox("Select Visualization Type", viz_types)
    
    # Get available data for visualization
    available_results = list(st.session_state.analysis_results.keys())
    selected_result = st.selectbox("Select Analysis Result", available_results)
    
    if st.button("Generate Visualization"):
        data = st.session_state.analysis_results[selected_result]
        
        with st.spinner("Generating visualization..."):
            fig = None
            if selected_viz == "Genome Browser View":
                fig = visualizer.create_genome_browser_view(data)
            elif selected_viz == "Heatmap":
                fig = visualizer.create_heatmap(data)
            elif selected_viz == "Volcano Plot":
                fig = visualizer.create_volcano_plot(data)
            elif selected_viz == "MA Plot":
                fig = visualizer.create_ma_plot(data)
            elif selected_viz == "PCA Plot":
                fig = visualizer.create_pca_plot(data)
            elif selected_viz == "Peak Distribution":
                fig = visualizer.create_peak_distribution(data)
        
        if fig:
            st.plotly_chart(fig, use_container_width=True)
            
            # Download option
            if st.button("Download Plot"):
                # Create download link for the plot
                fig.write_html("plot.html")
                with open("plot.html", "rb") as file:
                    st.download_button(
                        label="Download HTML",
                        data=file,
                        file_name=f"{selected_viz.lower().replace(' ', '_')}.html",
                        mime="text/html"
                    )

def export_results_page():
    st.header("💾 Export Results")
    
    if not st.session_state.analysis_results:
        st.warning("No analysis results available for export.")
        return
    
    # Export format selection
    export_formats = ["CSV", "TSV", "Excel", "JSON"]
    selected_format = st.selectbox("Select Export Format", export_formats)
    
    # Results selection
    available_results = list(st.session_state.analysis_results.keys())
    results_to_export = st.multiselect(
        "Select Results to Export",
        available_results,
        default=available_results
    )
    
    if st.button("Export Results"):
        try:
            # Create export data
            export_data = {}
            for result_key in results_to_export:
                export_data[result_key] = st.session_state.analysis_results[result_key]
            
            if selected_format == "JSON":
                import json
                export_string = json.dumps(export_data, indent=2, default=str)
                st.download_button(
                    label="Download JSON",
                    data=export_string,
                    file_name="genomics_analysis_results.json",
                    mime="application/json"
                )
            
            elif selected_format in ["CSV", "TSV"]:
                export_files = {}
                
                for result_key, result_data in export_data.items():
                    if isinstance(result_data, dict):
                        if 'peaks' in str(result_data) and isinstance(result_data, list):
                            df = pd.DataFrame(result_data)
                            export_files[f"{result_key}_peaks"] = df
                        elif any(key in result_data for key in ['upregulated', 'downregulated', 'total_genes']):
                            # Differential expression results
                            df = pd.DataFrame([result_data])
                            export_files[f"{result_key}_summary"] = df
                        elif any(key in result_data for key in ['common_features', 'unique_tissue1']):
                            # Tissue comparison results
                            df = pd.DataFrame([result_data])
                            export_files[f"{result_key}_comparison"] = df
                        else:
                            df = pd.DataFrame([result_data])
                            export_files[f"{result_key}"] = df
                
                import zipfile
                import io
                
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                    for file_key, df in export_files.items():
                        separator = ',' if selected_format == "CSV" else '\t'
                        csv_string = df.to_csv(sep=separator, index=False)
                        zip_file.writestr(f"{file_key}.{selected_format.lower()}", csv_string)
                
                st.download_button(
                    label=f"Download {selected_format} Archive",
                    data=zip_buffer.getvalue(),
                    file_name=f"genomics_analysis_results.zip",
                    mime="application/zip"
                )
            
            st.success("Export file generated successfully!")
            
        except Exception as e:
            st.error(f"Export failed: {str(e)}")
    
    # Summary of available results
    st.subheader("Available Results Summary")
    for result_key, result_data in st.session_state.analysis_results.items():
        with st.expander(f"{result_key}"):
            if isinstance(result_data, dict):
                for key, value in result_data.items():
                    st.write(f"**{key}:** {value}")
            else:
                st.write(str(result_data))

if __name__ == "__main__":
    main()

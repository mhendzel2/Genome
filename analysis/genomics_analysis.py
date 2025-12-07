import pandas as pd
import numpy as np
from typing import Dict, Any, List, Optional
import tempfile
import os
from utils.gprofiler_client import GProfilerClient
from statsmodels.stats.multitest import fdrcorrection

# Try to import cooler for Hi-C analysis
try:
    import cooler
    import cooltools
    import cooltools.lib.numutils
except ImportError:
    cooler = None
    cooltools = None

class GenomicsAnalyzer:
    """Main genomics analysis coordinator"""
    
    def __init__(self, r_integration):
        self.r_integration = r_integration
        self.supported_analyses = [
            'basic_statistics',
            'peak_calling', 
            'differential_expression',
            'tissue_comparison',
            'enrichment_analysis',
            'chromatin_interaction_analysis',
            'ai_interpretation'
        ]
    
    def basic_statistics(self, uploaded_files: Dict[str, Any]) -> Dict[str, Any]:
        """Compute basic statistics for uploaded files"""
        if not uploaded_files:
            return {}
        
        # Use R integration for statistics
        return self.r_integration.basic_statistics(uploaded_files)
    
    def peak_calling(self, uploaded_files: Dict[str, Any], **params) -> List[Dict[str, Any]]:
        """Perform peak calling on ChIP-seq data"""
        if not uploaded_files:
            return []
        
        # Filter for ChIP-seq files
        chip_files = {
            filename: file_info 
            for filename, file_info in uploaded_files.items()
            if file_info['type'] in ['ChIP-seq', 'Histone Marks']
        }
        
        if not chip_files:
            return []
        
        return self.r_integration.peak_calling(chip_files, **params)
    
    def differential_expression(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """Perform differential expression analysis"""
        if not uploaded_files:
            return {}
        
        # Filter for expression files
        expr_files = {
            filename: file_info 
            for filename, file_info in uploaded_files.items()
            if file_info['type'] == 'Gene Expression'
        }
        
        if not expr_files:
            return {'error': 'No gene expression data found'}
        
        results = self.r_integration.differential_expression(expr_files, **params)

        if isinstance(results, dict) and 'pvalue' in results:
            pvals = results['pvalue']
            if isinstance(pvals, list) and len(pvals) > 0:
                _, padj = fdrcorrection(pvals)
                results['padj'] = padj.tolist()

        return results
    
    def tissue_comparison(self, tissues: List[str], comparison_type: str, 
                         uploaded_files: Dict[str, Any]) -> Dict[str, Any]:
        """Compare genomics data across tissues"""
        if len(tissues) < 2:
            return {'error': 'At least 2 tissues required for comparison'}
        
        return self.r_integration.tissue_comparison(tissues, comparison_type, uploaded_files)
    
    def enrichment_analysis(self, uploaded_files: Dict[str, Any], 
                          gene_sets: Optional[List[str]] = None, **params) -> Dict[str, Any]:
        """Perform pathway/GO enrichment analysis"""
        if not uploaded_files:
            return {}
        
        # Extract genes from uploaded expression files
        all_genes = set()
        significant_genes = set()
        
        for filename, file_info in uploaded_files.items():
            if file_info['type'] == 'Gene Expression':
                try:
                    file_info['file'].seek(0)
                    content = file_info['file'].read().decode('utf-8')
                    lines = [line for line in content.split('\n') 
                            if line.strip() and not line.startswith('#')]
                    
                    for line in lines:
                        fields = line.split('\t')
                        if len(fields) >= 2:
                            gene_id = fields[0]
                            all_genes.add(gene_id)
                            
                            try:
                                expr_value = float(fields[1])
                                if abs(expr_value) > params.get('expression_threshold', 2.0):
                                    significant_genes.add(gene_id)
                            except ValueError:
                                continue
                                
                except Exception as e:
                    print(f"Error processing {filename}: {e}")
        
        if significant_genes:
            gprofiler_client = GProfilerClient()
            enrichment_results = gprofiler_client.perform_enrichment(list(significant_genes))
            return enrichment_results
        else:
            return pd.DataFrame()
    
    def chromatin_interaction_analysis(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """Analyze chromatin interaction data (Hi-C) using cooler"""
        if not uploaded_files:
            return {}
        
        # Filter for Hi-C files (.cool or .mcool)
        hic_files = {
            filename: file_info 
            for filename, file_info in uploaded_files.items()
            if file_info['type'] == 'HiC' or filename.endswith('.cool') or filename.endswith('.mcool')
        }
        
        if not hic_files:
            return {'error': 'No Hi-C data found (expecting .cool files)'}
        
        if not cooler:
            return {'error': 'cooler library not installed'}

        results = {}
        
        for filename, file_info in hic_files.items():
            try:
                # Cooler requires a file path, write to temp
                with tempfile.NamedTemporaryFile(suffix='.cool', delete=False) as tmp:
                    file_info['file'].seek(0)
                    tmp.write(file_info['file'].read())
                    tmp_path = tmp.name
                
                try:
                    clr = cooler.Cooler(tmp_path)

                    # Basic stats
                    n_interactions = clr.info['sum']
                    n_bins = clr.info['nbins']
                    chromosomes = clr.chromnames

                    # Calculate simple stats on pixels (sparse matrix)
                    pixels = clr.pixels()[:]
                    mean_counts = pixels['count'].mean()
                    max_counts = pixels['count'].max()

                    # Insulation score (TAD boundaries approximation) using cooltools if available
                    tad_info = "Insulation score calculation skipped (cooltools not fully configured)"
                    if cooltools:
                         # Simplified: just mention availability
                         tad_info = "Cooltools available for TAD analysis"

                    results[filename] = {
                        'total_interactions': int(n_interactions),
                        'number_of_bins': int(n_bins),
                        'chromosomes': chromosomes,
                        'mean_interactions_per_bin': float(mean_counts),
                        'max_interactions': int(max_counts),
                        'tad_analysis': tad_info
                    }

                finally:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
                
            except Exception as e:
                results[f'error_{filename}'] = str(e)
        
        return results
    
    def multi_omics_integration(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """
        Integrate multiple omics data types using established integration methods.
        
        Supports integration of gene expression, proteomics, ChIP-seq, and other
        omics data types using dimensionality reduction and correlation analysis.
        
        Args:
            uploaded_files: Dictionary of filename -> file_info with 'file' and 'type' keys
            **params: Additional parameters:
                - method: Integration method ('pca', 'cca', 'correlation', 'mofa')
                - n_components: Number of components for dimensionality reduction
                - scale: Whether to scale data before integration
        
        Returns:
            Dictionary with:
            - data_types_available: List of omics types found
            - integration_possible: Whether integration was successful
            - method_used: Integration method applied
            - joint_embedding: DataFrame with integrated low-dimensional representation
            - correlation_matrix: Cross-omics correlation matrix
            - feature_contributions: Top contributing features per component
            - common_features: Shared features across data types
            - quality_metrics: Integration quality scores
        """
        if not uploaded_files:
            return {'error': 'No files provided for integration'}
        
        # Import required libraries with fallbacks
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            from sklearn.cross_decomposition import CCA
        except ImportError:
            return {'error': 'scikit-learn required for multi-omics integration'}
        
        method = params.get('method', 'pca')
        n_components = params.get('n_components', 10)
        scale = params.get('scale', True)
        
        # Categorize files by data type
        data_matrices = {}
        data_types = {}
        
        for filename, file_info in uploaded_files.items():
            data_type = file_info.get('type', 'Unknown')
            if data_type not in data_types:
                data_types[data_type] = []
            data_types[data_type].append(filename)
            
            try:
                # Load data into DataFrame
                file_obj = file_info.get('file')
                if file_obj is not None:
                    file_obj.seek(0)
                    content = file_obj.read()
                    if isinstance(content, bytes):
                        content = content.decode('utf-8')
                    
                    # Parse as TSV
                    from io import StringIO
                    df = pd.read_csv(StringIO(content), sep='\t', index_col=0)
                    data_matrices[data_type] = df
                    
            except Exception as e:
                logger.warning(f"Could not load {filename}: {e}")
        
        results = {
            'data_types_available': list(data_types.keys()),
            'n_data_types': len(data_types),
            'integration_possible': len(data_types) >= 2,
            'method_used': method,
            'joint_embedding': None,
            'correlation_matrix': {},
            'feature_contributions': {},
            'common_features': {},
            'quality_metrics': {}
        }
        
        if len(data_matrices) < 2:
            results['integration_possible'] = False
            results['error'] = 'Need at least 2 data types for integration'
            return results
        
        # Find common samples/features
        matrix_list = list(data_matrices.items())
        
        # Try sample-based integration (samples as rows)
        common_samples = None
        for data_type, df in matrix_list:
            if common_samples is None:
                common_samples = set(df.index)
            else:
                common_samples = common_samples & set(df.index)
        
        # Also try feature-based integration (features as rows)
        common_features = None
        for data_type, df in matrix_list:
            if common_features is None:
                common_features = set(df.columns)
            else:
                common_features = common_features & set(df.columns)
        
        results['common_samples'] = len(common_samples) if common_samples else 0
        results['common_features_count'] = len(common_features) if common_features else 0
        
        # Determine integration axis
        if common_samples and len(common_samples) >= 5:
            # Sample-based integration
            integration_axis = 'samples'
            aligned_matrices = {}
            
            for data_type, df in matrix_list:
                aligned_df = df.loc[list(common_samples)].copy()
                # Remove constant features
                aligned_df = aligned_df.loc[:, aligned_df.std() > 0]
                aligned_matrices[data_type] = aligned_df
                results['common_features'][data_type] = len(aligned_df.columns)
            
        elif common_features and len(common_features) >= 5:
            # Feature-based integration (transpose)
            integration_axis = 'features'
            aligned_matrices = {}
            
            for data_type, df in matrix_list:
                aligned_df = df[list(common_features)].T.copy()
                aligned_df = aligned_df.loc[:, aligned_df.std() > 0]
                aligned_matrices[data_type] = aligned_df
                results['common_features'][data_type] = len(aligned_df.columns)
        else:
            results['integration_possible'] = False
            results['error'] = 'Insufficient common samples or features for integration'
            return results
        
        results['integration_axis'] = integration_axis
        
        # Perform integration based on method
        scaler = StandardScaler()
        
        if method == 'pca':
            # Concatenate all matrices and perform joint PCA
            try:
                combined_data = []
                feature_names = []
                type_labels = []
                
                for data_type, df in aligned_matrices.items():
                    if scale:
                        scaled_data = scaler.fit_transform(df.values)
                    else:
                        scaled_data = df.values
                    combined_data.append(scaled_data)
                    feature_names.extend([f"{data_type}_{c}" for c in df.columns])
                    type_labels.extend([data_type] * df.shape[1])
                
                combined_matrix = np.hstack(combined_data)
                
                # Handle case where n_components exceeds dimensions
                n_comp = min(n_components, combined_matrix.shape[0], combined_matrix.shape[1])
                
                pca = PCA(n_components=n_comp)
                embedding = pca.fit_transform(combined_matrix)
                
                # Create embedding DataFrame
                sample_names = list(list(aligned_matrices.values())[0].index)
                embedding_df = pd.DataFrame(
                    embedding,
                    index=sample_names,
                    columns=[f'PC{i+1}' for i in range(n_comp)]
                )
                results['joint_embedding'] = embedding_df.to_dict()
                
                # Feature contributions (loadings)
                loadings = pd.DataFrame(
                    pca.components_.T,
                    index=feature_names,
                    columns=[f'PC{i+1}' for i in range(n_comp)]
                )
                
                # Top contributing features per component
                for i in range(min(3, n_comp)):
                    pc_name = f'PC{i+1}'
                    top_features = loadings[pc_name].abs().nlargest(10).index.tolist()
                    results['feature_contributions'][pc_name] = top_features
                
                # Quality metrics
                results['quality_metrics']['variance_explained'] = pca.explained_variance_ratio_.tolist()
                results['quality_metrics']['total_variance_explained'] = float(sum(pca.explained_variance_ratio_))
                
            except Exception as e:
                logger.exception("PCA integration failed: %s", e)
                results['error'] = f'PCA integration failed: {e}'
        
        elif method == 'cca':
            # Canonical Correlation Analysis (for 2 data types)
            if len(aligned_matrices) != 2:
                results['error'] = 'CCA requires exactly 2 data types'
                return results
            
            try:
                type_names = list(aligned_matrices.keys())
                X = aligned_matrices[type_names[0]].values
                Y = aligned_matrices[type_names[1]].values
                
                if scale:
                    X = scaler.fit_transform(X)
                    Y = scaler.fit_transform(Y)
                
                n_comp = min(n_components, X.shape[0], X.shape[1], Y.shape[1])
                
                cca = CCA(n_components=n_comp)
                X_c, Y_c = cca.fit_transform(X, Y)
                
                # Combined embedding
                sample_names = list(aligned_matrices[type_names[0]].index)
                embedding_df = pd.DataFrame(
                    np.hstack([X_c, Y_c]),
                    index=sample_names,
                    columns=[f'{type_names[0]}_CC{i+1}' for i in range(n_comp)] + 
                           [f'{type_names[1]}_CC{i+1}' for i in range(n_comp)]
                )
                results['joint_embedding'] = embedding_df.to_dict()
                
                # Canonical correlations
                correlations = []
                for i in range(n_comp):
                    corr = np.corrcoef(X_c[:, i], Y_c[:, i])[0, 1]
                    correlations.append(float(corr))
                
                results['quality_metrics']['canonical_correlations'] = correlations
                
            except Exception as e:
                logger.exception("CCA integration failed: %s", e)
                results['error'] = f'CCA integration failed: {e}'
        
        elif method == 'correlation':
            # Simple correlation-based integration
            try:
                type_names = list(aligned_matrices.keys())
                
                # Cross-correlation matrix between data types
                for i, type1 in enumerate(type_names):
                    for type2 in type_names[i+1:]:
                        df1 = aligned_matrices[type1]
                        df2 = aligned_matrices[type2]
                        
                        # Sample-wise correlations
                        correlations = []
                        for sample in df1.index:
                            v1 = df1.loc[sample].values
                            v2 = df2.loc[sample].values
                            if len(v1) > 0 and len(v2) > 0:
                                # Handle different feature dimensions
                                min_len = min(len(v1), len(v2))
                                corr = np.corrcoef(v1[:min_len], v2[:min_len])[0, 1]
                                if not np.isnan(corr):
                                    correlations.append(corr)
                        
                        key = f'{type1}_vs_{type2}'
                        results['correlation_matrix'][key] = {
                            'mean_correlation': float(np.mean(correlations)) if correlations else 0,
                            'std_correlation': float(np.std(correlations)) if correlations else 0,
                            'n_samples': len(correlations)
                        }
                
                results['quality_metrics']['integration_type'] = 'correlation-based'
                
            except Exception as e:
                logger.exception("Correlation integration failed: %s", e)
                results['error'] = f'Correlation integration failed: {e}'
        
        else:
            results['error'] = f'Unknown integration method: {method}'
        
        # Calculate integration quality score
        if 'error' not in results:
            quality_score = 0.0
            
            if results.get('quality_metrics', {}).get('variance_explained'):
                # For PCA: use cumulative variance explained
                quality_score = results['quality_metrics']['total_variance_explained']
            elif results.get('quality_metrics', {}).get('canonical_correlations'):
                # For CCA: use mean canonical correlation
                quality_score = np.mean(results['quality_metrics']['canonical_correlations'])
            elif results.get('correlation_matrix'):
                # For correlation: use mean cross-correlation
                cors = [v.get('mean_correlation', 0) for v in results['correlation_matrix'].values()]
                quality_score = np.mean(cors) if cors else 0
            
            results['quality_metrics']['integration_score'] = float(quality_score)
        
        logger.info(f"Multi-omics integration completed: {len(data_types)} data types, "
                   f"method={method}, score={results.get('quality_metrics', {}).get('integration_score', 'N/A')}")
        
        return results
    
    def quality_control(self, uploaded_files: Dict[str, Any]) -> Dict[str, Any]:
        """Perform quality control analysis on uploaded files"""
        if not uploaded_files:
            return {}
        
        qc_results = {}
        
        for filename, file_info in uploaded_files.items():
            try:
                # Reset pointer
                try:
                    file_info['file'].seek(0)
                except Exception:
                    pass

                # Read first chunk for QC to avoid memory issues
                content_chunk = file_info['file'].read(10000)
                if isinstance(content_chunk, bytes):
                    content_chunk = content_chunk.decode('utf-8', errors='ignore')

                lines = content_chunk.split('\n')
                
                # Basic QC metrics (estimated from chunk)
                total_lines = len(lines)
                data_lines = len([l for l in lines if l.strip() and not l.startswith('#')])
                
                qc_metrics = {
                    'preview_lines': total_lines,
                    'preview_data_lines': data_lines,
                    'quality_score': 'PASS' if data_lines > 0 else 'FAIL'
                }
                
                # Add data type-specific QC
                data_type = file_info['type']
                
                if data_type in ['ChIP-seq', 'Histone Marks', 'BED']:
                    # Check for proper BED format
                    bed_format_ok = True
                    invalid_coords = 0
                    
                    for line in lines[:100]:  # Check first 100 lines
                        if line.strip() and not line.startswith('#'):
                            fields = line.split('\t')
                            if len(fields) >= 3:
                                try:
                                    start = int(fields[1])
                                    end = int(fields[2])
                                    if start >= end or start < 0:
                                        invalid_coords += 1
                                except ValueError:
                                    invalid_coords += 1
                            else:
                                bed_format_ok = False
                    
                    qc_metrics.update({
                        'bed_format_valid': bed_format_ok,
                        'invalid_coordinates': invalid_coords,
                    })
                
                elif data_type == 'Gene Expression':
                    # Check for consistent column counts
                    column_counts = []
                    for line in lines[:100]:
                        if line.strip() and not line.startswith('#'):
                            column_counts.append(len(line.split('\t')))
                    
                    consistent_columns = len(set(column_counts)) <= 1
                    
                    qc_metrics.update({
                        'consistent_columns': consistent_columns,
                        'column_count': max(column_counts) if column_counts else 0,
                    })
                
                qc_results[filename] = qc_metrics
                
            except Exception as e:
                qc_results[filename] = {'quality_score': 'ERROR', 'error': str(e)}
        
        return qc_results
    
    def generate_ai_interpretation(self, enrichment_results: pd.DataFrame, tissue: str = "Unknown") -> str:
        """
        Generate an AI-driven interpretation of enrichment results using LLM (RAG).

        Args:
            enrichment_results: DataFrame from enrichment analysis.
            tissue: Tissue type for context.

        Returns:
            String containing the AI interpretation.
        """
        if enrichment_results.empty:
            return "No significant pathways found to interpret."

        # Extract top pathways
        if 'name' in enrichment_results.columns:
            pathways = enrichment_results['name'].head(10).tolist()
        elif 'term_name' in enrichment_results.columns:
            pathways = enrichment_results['term_name'].head(10).tolist()
        else:
            return "Could not identify pathway names in results."

        pathway_list_str = ", ".join(pathways)

        # Construct Prompt
        prompt = f"""
        You are a helpful assistant for genomic data analysis.
        Given these significant pathways: [{pathway_list_str}] found in a {tissue} sample,
        summarize potentially biological implications citing recent literature (2024-2025 context).
        Focus on clinical relevance and mechanism.
        """

        # Mock LLM Call (Placeholder for GPT-4o / BioMistral API)
        # In a real implementation, we would call:
        # response = requests.post(LLM_API_URL, json={"prompt": prompt})

        interpretation = f"""
        ### AI Interpretation (GeneChat)

        **Context:** {tissue} Tissue Analysis
        **Top Pathways:** {pathway_list_str}

        **Biological Implications:**
        Based on the identified pathways, the sample shows strong enrichment in processes related to cellular signaling and metabolism.
        Recent literature suggests these pathways are critical in {tissue} homeostasis and disease progression.

        *Note: This is a generated interpretation based on pathway names. For clinical decisions, verify with targeted assays.*

        (This feature is currently running in mock mode. Integrate an LLM API key to enable full RAG capabilities.)
        """

        return interpretation

    def generate_analysis_report(self, analysis_results: Dict[str, Any], ai_interpretation: Optional[str] = None) -> str:
        """Generate a comprehensive analysis report"""
        report_lines = []
        report_lines.append("# Genomics Analysis Report")
        report_lines.append(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append("")
        
        if ai_interpretation:
            report_lines.append(ai_interpretation)
            report_lines.append("")
            report_lines.append("---")
            report_lines.append("")

        for analysis_type, results in analysis_results.items():
            report_lines.append(f"## {analysis_type.replace('_', ' ').title()}")
            
            if isinstance(results, dict):
                for key, value in results.items():
                    if isinstance(value, (int, float)):
                        report_lines.append(f"- {key}: {value:,.2f}" if isinstance(value, float) else f"- {key}: {value:,}")
                    else:
                        report_lines.append(f"- {key}: {value}")
            elif isinstance(results, list):
                report_lines.append(f"- Total items: {len(results)}")
                if results and isinstance(results[0], dict):
                    report_lines.append("- Sample entries:")
                    for item in results[:3]:  # Show first 3 items
                        report_lines.append(f"  - {item}")
            
            report_lines.append("")
        
        return "\n".join(report_lines)
    
    def generate_ai_interpretation(self, analysis_results: Dict[str, Any],
                                  tissue_context: Optional[str] = None,
                                  api_provider: str = 'openai',
                                  include_literature: bool = True,
                                  **params) -> Dict[str, Any]:
        """
        Generate AI-powered interpretation of genomics analysis results.
        
        Uses language models with optional retrieval-augmented generation (RAG)
        to provide clinically relevant interpretations of enrichment results,
        differential expression, and other genomics findings.
        
        Args:
            analysis_results: Dictionary containing analysis outputs (e.g., enriched pathways)
            tissue_context: Tissue/cell type context for interpretation
            api_provider: LLM provider ('openai', 'anthropic', or 'local')
            include_literature: Whether to include PubMed literature search
            **params: Additional parameters:
                - api_key: API key (can also use environment variable)
                - model: Specific model to use
                - max_tokens: Maximum response length
        
        Returns:
            Dictionary with:
            - interpretation: AI-generated interpretation text
            - key_findings: List of main findings
            - clinical_relevance: Clinical implications (if applicable)
            - literature_references: Relevant recent publications
            - confidence_notes: Caveats about interpretation
            - status: Success/error status
        """
        import os
        import json
        
        result = {
            'interpretation': '',
            'key_findings': [],
            'clinical_relevance': '',
            'literature_references': [],
            'confidence_notes': [],
            'status': 'pending'
        }
        
        # Extract key information from analysis results
        enriched_pathways = []
        significant_genes = []
        fold_changes = {}
        
        if isinstance(analysis_results, dict):
            # Handle enrichment results
            if 'enrichment' in analysis_results:
                enr = analysis_results['enrichment']
                if isinstance(enr, pd.DataFrame) and not enr.empty:
                    enriched_pathways = enr.head(20).to_dict('records')
                elif isinstance(enr, list):
                    enriched_pathways = enr[:20]
            
            # Handle differential expression
            if 'significant_genes' in analysis_results:
                significant_genes = analysis_results['significant_genes'][:100]
            if 'fold_changes' in analysis_results:
                fold_changes = analysis_results['fold_changes']
            
            # Direct pathway data
            if 'name' in analysis_results or 'term_name' in analysis_results:
                enriched_pathways = [analysis_results]
        
        elif isinstance(analysis_results, pd.DataFrame) and not analysis_results.empty:
            # DataFrame with enrichment results
            enriched_pathways = analysis_results.head(20).to_dict('records')
        
        elif isinstance(analysis_results, list):
            enriched_pathways = analysis_results[:20]
        
        # Build context for LLM
        context_parts = []
        
        if enriched_pathways:
            pathway_text = "Enriched Pathways/GO Terms:\n"
            for i, p in enumerate(enriched_pathways[:15], 1):
                name = p.get('name', p.get('term_name', p.get('description', 'Unknown')))
                pval = p.get('p_value', p.get('pvalue', p.get('adjusted_p_value', 'N/A')))
                pathway_text += f"{i}. {name} (p={pval})\n"
            context_parts.append(pathway_text)
        
        if significant_genes:
            gene_text = f"Significant genes ({len(significant_genes)} total): {', '.join(significant_genes[:30])}"
            if len(significant_genes) > 30:
                gene_text += f"... and {len(significant_genes) - 30} more"
            context_parts.append(gene_text)
        
        if tissue_context:
            context_parts.append(f"Tissue/Cell Type Context: {tissue_context}")
        
        if not context_parts:
            result['status'] = 'error'
            result['interpretation'] = 'Insufficient data for interpretation. Please provide enrichment results or significant gene lists.'
            return result
        
        full_context = "\n\n".join(context_parts)
        
        # Fetch relevant literature if requested
        if include_literature:
            literature_refs = self._search_pubmed_literature(
                enriched_pathways[:5] if enriched_pathways else [],
                significant_genes[:10] if significant_genes else [],
                tissue_context
            )
            result['literature_references'] = literature_refs
        
        # Build prompt for LLM
        prompt = f"""You are an expert genomics researcher and bioinformatician. Analyze the following genomics analysis results and provide a clear, scientifically accurate interpretation.

## Analysis Results:
{full_context}

## Task:
1. Summarize the main biological themes and pathways identified
2. Explain the potential biological significance of these findings
3. Identify any notable gene interactions or pathway crosstalk
4. If applicable, discuss clinical or therapeutic implications
5. Note any limitations or caveats in interpreting these results

Please structure your response with clear sections and be specific about the biological mechanisms involved."""

        if result['literature_references']:
            lit_text = "\n## Relevant Literature:\n"
            for ref in result['literature_references'][:5]:
                lit_text += f"- {ref.get('title', 'N/A')} ({ref.get('year', 'N/A')})\n"
            prompt += lit_text
        
        # Try to call LLM API
        api_key = params.get('api_key') or os.environ.get('OPENAI_API_KEY') or os.environ.get('ANTHROPIC_API_KEY')
        
        if api_key:
            try:
                interpretation = self._call_llm_api(
                    prompt=prompt,
                    api_provider=api_provider,
                    api_key=api_key,
                    model=params.get('model'),
                    max_tokens=params.get('max_tokens', 2000)
                )
                result['interpretation'] = interpretation
                result['status'] = 'success'
                
                # Extract key findings from interpretation
                result['key_findings'] = self._extract_key_findings(interpretation)
                
            except Exception as e:
                logger.exception("LLM API call failed: %s", e)
                result['status'] = 'api_error'
                result['interpretation'] = self._generate_fallback_interpretation(
                    enriched_pathways, significant_genes, tissue_context
                )
                result['confidence_notes'].append(f"API call failed, using rule-based interpretation: {e}")
        else:
            # Fallback to rule-based interpretation
            result['status'] = 'no_api_key'
            result['interpretation'] = self._generate_fallback_interpretation(
                enriched_pathways, significant_genes, tissue_context
            )
            result['confidence_notes'].append(
                "No API key provided. Using rule-based interpretation. "
                "For AI-powered insights, set OPENAI_API_KEY or ANTHROPIC_API_KEY environment variable."
            )
        
        return result
    
    def _call_llm_api(self, prompt: str, api_provider: str, api_key: str,
                     model: Optional[str] = None, max_tokens: int = 2000) -> str:
        """Call LLM API for interpretation generation."""
        import requests
        
        if api_provider == 'openai':
            model = model or 'gpt-4-turbo-preview'
            headers = {
                'Authorization': f'Bearer {api_key}',
                'Content-Type': 'application/json'
            }
            payload = {
                'model': model,
                'messages': [
                    {'role': 'system', 'content': 'You are an expert genomics researcher providing scientific interpretations of analysis results.'},
                    {'role': 'user', 'content': prompt}
                ],
                'max_tokens': max_tokens,
                'temperature': 0.3
            }
            
            response = requests.post(
                'https://api.openai.com/v1/chat/completions',
                headers=headers,
                json=payload,
                timeout=120
            )
            response.raise_for_status()
            return response.json()['choices'][0]['message']['content']
        
        elif api_provider == 'anthropic':
            model = model or 'claude-3-sonnet-20240229'
            headers = {
                'x-api-key': api_key,
                'Content-Type': 'application/json',
                'anthropic-version': '2023-06-01'
            }
            payload = {
                'model': model,
                'max_tokens': max_tokens,
                'messages': [
                    {'role': 'user', 'content': prompt}
                ]
            }
            
            response = requests.post(
                'https://api.anthropic.com/v1/messages',
                headers=headers,
                json=payload,
                timeout=120
            )
            response.raise_for_status()
            return response.json()['content'][0]['text']
        
        else:
            raise ValueError(f"Unsupported API provider: {api_provider}")
    
    def _search_pubmed_literature(self, pathways: List[Dict], genes: List[str],
                                 tissue: Optional[str] = None) -> List[Dict]:
        """Search PubMed for relevant recent literature."""
        import requests
        
        references = []
        
        # Build search queries from top pathways and genes
        search_terms = []
        
        for p in pathways[:3]:
            name = p.get('name', p.get('term_name', ''))
            if name:
                search_terms.append(name)
        
        if genes:
            search_terms.append(' OR '.join(genes[:5]))
        
        if tissue:
            search_terms.append(tissue)
        
        if not search_terms:
            return references
        
        # Combine search terms
        query = ' AND '.join([f'({term})' for term in search_terms[:3]])
        query += ' AND (2024[pdat] OR 2025[pdat])'  # Recent papers
        
        try:
            # Search PubMed via E-utilities
            search_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
            params = {
                'db': 'pubmed',
                'term': query,
                'retmax': 10,
                'retmode': 'json',
                'sort': 'relevance'
            }
            
            response = requests.get(search_url, params=params, timeout=30)
            response.raise_for_status()
            search_result = response.json()
            
            pmids = search_result.get('esearchresult', {}).get('idlist', [])
            
            if pmids:
                # Fetch article details
                fetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
                fetch_params = {
                    'db': 'pubmed',
                    'id': ','.join(pmids),
                    'retmode': 'json'
                }
                
                fetch_response = requests.get(fetch_url, params=fetch_params, timeout=30)
                fetch_response.raise_for_status()
                articles = fetch_response.json().get('result', {})
                
                for pmid in pmids:
                    if pmid in articles:
                        article = articles[pmid]
                        references.append({
                            'pmid': pmid,
                            'title': article.get('title', 'N/A'),
                            'authors': article.get('authors', [])[:3],
                            'journal': article.get('fulljournalname', 'N/A'),
                            'year': article.get('pubdate', 'N/A')[:4],
                            'url': f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
                        })
        
        except Exception as e:
            logger.warning(f"PubMed search failed: {e}")
        
        return references
    
    def _generate_fallback_interpretation(self, pathways: List[Dict], genes: List[str],
                                         tissue: Optional[str] = None) -> str:
        """Generate rule-based interpretation when LLM is unavailable."""
        interpretation_parts = []
        
        interpretation_parts.append("## Analysis Interpretation\n")
        
        if pathways:
            interpretation_parts.append("### Enriched Biological Processes\n")
            
            # Categorize pathways
            immune_pathways = []
            metabolic_pathways = []
            signaling_pathways = []
            other_pathways = []
            
            for p in pathways:
                name = p.get('name', p.get('term_name', '')).lower()
                if any(term in name for term in ['immune', 'inflammatory', 'cytokine', 'interferon', 'interleukin']):
                    immune_pathways.append(p)
                elif any(term in name for term in ['metabol', 'glycoly', 'lipid', 'amino acid', 'oxidative']):
                    metabolic_pathways.append(p)
                elif any(term in name for term in ['signaling', 'pathway', 'receptor', 'kinase']):
                    signaling_pathways.append(p)
                else:
                    other_pathways.append(p)
            
            if immune_pathways:
                interpretation_parts.append(f"**Immune/Inflammatory Response**: {len(immune_pathways)} pathways enriched, "
                                          f"suggesting involvement in immune processes. ")
            if metabolic_pathways:
                interpretation_parts.append(f"**Metabolic Changes**: {len(metabolic_pathways)} metabolic pathways identified, "
                                          f"indicating altered cellular metabolism. ")
            if signaling_pathways:
                interpretation_parts.append(f"**Cell Signaling**: {len(signaling_pathways)} signaling pathways enriched, "
                                          f"suggesting active signal transduction. ")
            
            interpretation_parts.append("\n")
        
        if genes:
            interpretation_parts.append(f"### Significant Genes\n")
            interpretation_parts.append(f"A total of {len(genes)} genes showed significant changes. ")
            interpretation_parts.append(f"Notable genes include: {', '.join(genes[:10])}.\n")
        
        if tissue:
            interpretation_parts.append(f"\n### Tissue Context\n")
            interpretation_parts.append(f"Analysis performed in context of {tissue}. "
                                       f"Interpret results considering tissue-specific expression patterns.\n")
        
        interpretation_parts.append("\n### Limitations\n")
        interpretation_parts.append("This is an automated rule-based interpretation. "
                                   "For more detailed analysis, consider using AI-powered interpretation "
                                   "with an API key or consulting domain experts.\n")
        
        return ''.join(interpretation_parts)
    
    def _extract_key_findings(self, interpretation: str) -> List[str]:
        """Extract key findings from LLM interpretation."""
        findings = []
        
        # Look for numbered lists or bullet points
        lines = interpretation.split('\n')
        for line in lines:
            line = line.strip()
            # Match numbered items or bullet points
            if line and (line[0].isdigit() or line.startswith('- ') or line.startswith('• ')):
                # Clean up the line
                clean_line = line.lstrip('0123456789.-•) ').strip()
                if len(clean_line) > 20 and len(clean_line) < 300:
                    findings.append(clean_line)
        
        return findings[:10]  # Return top 10 findings

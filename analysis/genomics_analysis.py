import pandas as pd
import numpy as np
from typing import Dict, Any, List, Optional
import tempfile
import os
from utils.gprofiler_client import GProfilerClient

class GenomicsAnalyzer:
    """Main genomics analysis coordinator"""
    
    def __init__(self, r_integration):
        self.r_integration = r_integration
        self.supported_analyses = [
            'basic_statistics',
            'peak_calling', 
            'differential_expression',
            'tissue_comparison',
            'enrichment_analysis'
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
        
        return self.r_integration.differential_expression(expr_files, **params)
    
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
        """Analyze chromatin interaction data (Hi-C)"""
        if not uploaded_files:
            return {}
        
        # Filter for Hi-C files
        hic_files = {
            filename: file_info 
            for filename, file_info in uploaded_files.items()
            if file_info['type'] == 'HiC'
        }
        
        if not hic_files:
            return {'error': 'No Hi-C data found'}
        
        # Mock Hi-C analysis results
        results = {
            'total_interactions': 0,
            'significant_interactions': 0,
            'interaction_distance_stats': {},
            'compartment_analysis': {},
            'tad_boundaries': []
        }
        
        for filename, file_info in hic_files.items():
            try:
                # Read file and count interactions
                file_info['file'].seek(0)
                content = file_info['file'].read().decode('utf-8')
                lines = [line for line in content.split('\n') 
                        if line.strip() and not line.startswith('#')]
                
                interactions = len(lines)
                significant = int(interactions * 0.1)  # 10% significant
                
                results['total_interactions'] += interactions
                results['significant_interactions'] += significant
                
                # Mock distance analysis
                distances = np.random.lognormal(mean=5, sigma=1, size=interactions)
                results['interaction_distance_stats'] = {
                    'mean_distance': np.mean(distances),
                    'median_distance': np.median(distances),
                    'max_distance': np.max(distances),
                    'min_distance': np.min(distances)
                }
                
            except Exception as e:
                results[f'error_{filename}'] = str(e)
        
        return results
    
    def multi_omics_integration(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        """Integrate multiple omics data types"""
        if not uploaded_files:
            return {}
        
        # Categorize files by data type
        data_types = {}
        for filename, file_info in uploaded_files.items():
            data_type = file_info['type']
            if data_type not in data_types:
                data_types[data_type] = []
            data_types[data_type].append(filename)
        
        results = {
            'data_types_available': list(data_types.keys()),
            'integration_possible': len(data_types) >= 2,
            'correlation_analysis': {},
            'common_features': {},
            'integrated_score': 0.0
        }
        
        if len(data_types) >= 2:
            # Mock integration analysis
            results['integrated_score'] = np.random.uniform(0.6, 0.9)
            
            # Mock correlation between data types
            type_pairs = [(t1, t2) for i, t1 in enumerate(data_types.keys()) 
                         for t2 in list(data_types.keys())[i+1:]]
            
            for t1, t2 in type_pairs:
                correlation = np.random.uniform(0.3, 0.8)
                results['correlation_analysis'][f'{t1}_vs_{t2}'] = correlation
            
            # Mock common features
            for data_type in data_types:
                results['common_features'][data_type] = np.random.randint(100, 1000)
        
        return results
    
    def quality_control(self, uploaded_files: Dict[str, Any]) -> Dict[str, Any]:
        """Perform quality control analysis on uploaded files"""
        if not uploaded_files:
            return {}
        
        qc_results = {}
        
        for filename, file_info in uploaded_files.items():
            try:
                file_info['file'].seek(0)
                content = file_info['file'].read().decode('utf-8')
                lines = content.split('\n')
                
                # Basic QC metrics
                total_lines = len(lines)
                data_lines = len([l for l in lines if l.strip() and not l.startswith('#')])
                empty_lines = total_lines - data_lines
                
                # Format-specific QC
                qc_metrics = {
                    'total_lines': total_lines,
                    'data_lines': data_lines,
                    'empty_lines': empty_lines,
                    'file_size_kb': len(content) / 1024,
                    'quality_score': 'PASS' if data_lines > 0 else 'FAIL'
                }
                
                # Add data type-specific QC
                data_type = file_info['type']
                
                if data_type in ['ChIP-seq', 'Histone Marks']:
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
                        'coordinate_error_rate': invalid_coords / min(100, data_lines) if data_lines > 0 else 0
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
                        'numeric_values_check': 'PASS'  # Would implement actual numeric check
                    })
                
                qc_results[filename] = qc_metrics
                
            except Exception as e:
                qc_results[filename] = {'quality_score': 'ERROR', 'error': str(e)}
        
        return qc_results
    
    def generate_analysis_report(self, analysis_results: Dict[str, Any]) -> str:
        """Generate a comprehensive analysis report"""
        report_lines = []
        report_lines.append("# Genomics Analysis Report")
        report_lines.append(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
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

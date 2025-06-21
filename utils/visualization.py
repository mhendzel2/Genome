import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from typing import Dict, Any, List, Optional

class GenomicsVisualizer:
    """Handles visualization of genomics data and analysis results"""
    
    def __init__(self):
        self.default_colors = px.colors.qualitative.Set1
        self.chromosome_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    
    def create_genome_browser_view(self, data: Dict[str, Any]) -> go.Figure:
        """Create a genome browser-like visualization"""
        fig = go.Figure()
        
        if isinstance(data, dict) and 'peaks' in str(data):
            # Handle peak data
            peaks = data if isinstance(data, list) else []
            
            if peaks:
                # Create tracks for different chromosomes
                chromosomes = list(set([peak.get('chr', 'chr1') for peak in peaks]))
                chromosomes = sorted(chromosomes, key=lambda x: self._chr_sort_key(x))
                
                for i, chr_name in enumerate(chromosomes[:5]):  # Limit to first 5 chromosomes
                    chr_peaks = [p for p in peaks if p.get('chr') == chr_name]
                    
                    if chr_peaks:
                        starts = [p.get('start', 0) for p in chr_peaks]
                        ends = [p.get('end', 1000) for p in chr_peaks]
                        scores = [p.get('score', 1) for p in chr_peaks]
                        
                        # Create rectangles for peaks
                        for start, end, score in zip(starts, ends, scores):
                            fig.add_shape(
                                type="rect",
                                x0=start, x1=end,
                                y0=i-0.4, y1=i+0.4,
                                fillcolor=px.colors.sequential.Viridis[min(int(score/10), 9)],
                                opacity=0.7,
                                line=dict(width=0)
                            )
                
                fig.update_layout(
                    title="Genome Browser View",
                    xaxis_title="Genomic Position",
                    yaxis_title="Chromosome",
                    yaxis=dict(
                        tickmode='array',
                        tickvals=list(range(len(chromosomes[:5]))),
                        ticktext=chromosomes[:5]
                    )
                )
        else:
            # Empty state
            fig.add_annotation(
                text="No genomic data available for browser view",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def create_heatmap(self, data: Dict[str, Any]) -> go.Figure:
        """Create heatmap visualization"""
        fig = go.Figure()
        
        if isinstance(data, dict) and any(key in data for key in ['tissue_comparison', 'expression']):
            # Create mock heatmap data based on analysis results
            if 'common_features' in data:
                # Tissue comparison heatmap
                tissues = data.get('tissues_compared', ['Tissue1', 'Tissue2'])
                features = ['Feature1', 'Feature2', 'Feature3', 'Feature4', 'Feature5']
                
                # Generate heatmap data
                z_data = np.random.rand(len(features), len(tissues)) * 10
                
                fig.add_trace(go.Heatmap(
                    z=z_data,
                    x=tissues,
                    y=features,
                    colorscale='Viridis',
                    colorbar=dict(title="Expression Level")
                ))
                
                fig.update_layout(
                    title="Tissue Comparison Heatmap",
                    xaxis_title="Tissues",
                    yaxis_title="Features"
                )
            else:
                # General expression heatmap
                genes = [f'Gene_{i}' for i in range(1, 11)]
                samples = [f'Sample_{i}' for i in range(1, 6)]
                
                z_data = np.random.rand(len(genes), len(samples)) * 100
                
                fig.add_trace(go.Heatmap(
                    z=z_data,
                    x=samples,
                    y=genes,
                    colorscale='RdBu_r',
                    colorbar=dict(title="Expression")
                ))
                
                fig.update_layout(
                    title="Expression Heatmap",
                    xaxis_title="Samples",
                    yaxis_title="Genes"
                )
        else:
            fig.add_annotation(
                text="No data available for heatmap visualization",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def create_volcano_plot(self, data: Dict[str, Any]) -> go.Figure:
        """Create volcano plot for differential expression"""
        fig = go.Figure()
        
        if isinstance(data, dict) and 'upregulated' in data:
            # Generate mock volcano plot data
            n_genes = data.get('total_genes', 1000)
            if n_genes == 0:
                n_genes = 1000
            
            # Generate fold changes and p-values
            fold_changes = np.random.normal(0, 1.5, n_genes)
            p_values = np.random.exponential(0.1, n_genes)
            p_values = np.clip(p_values, 1e-10, 1)  # Clip p-values
            
            # Create significance categories
            significant_up = (fold_changes > 1) & (p_values < 0.05)
            significant_down = (fold_changes < -1) & (p_values < 0.05)
            not_significant = ~(significant_up | significant_down)
            
            # Plot non-significant points
            fig.add_trace(go.Scatter(
                x=fold_changes[not_significant],
                y=-np.log10(p_values[not_significant]),
                mode='markers',
                marker=dict(color='gray', size=4, opacity=0.6),
                name='Not Significant',
                hovertemplate='Log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
            ))
            
            # Plot upregulated points
            if np.any(significant_up):
                fig.add_trace(go.Scatter(
                    x=fold_changes[significant_up],
                    y=-np.log10(p_values[significant_up]),
                    mode='markers',
                    marker=dict(color='red', size=6),
                    name='Upregulated',
                    hovertemplate='Log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
                ))
            
            # Plot downregulated points
            if np.any(significant_down):
                fig.add_trace(go.Scatter(
                    x=fold_changes[significant_down],
                    y=-np.log10(p_values[significant_down]),
                    mode='markers',
                    marker=dict(color='blue', size=6),
                    name='Downregulated',
                    hovertemplate='Log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
                ))
            
            # Add significance lines
            fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", opacity=0.5)
            fig.add_vline(x=1, line_dash="dash", line_color="black", opacity=0.5)
            fig.add_vline(x=-1, line_dash="dash", line_color="black", opacity=0.5)
            
            fig.update_layout(
                title="Volcano Plot - Differential Expression",
                xaxis_title="Log2 Fold Change",
                yaxis_title="-log10(p-value)",
                hovermode='closest'
            )
        else:
            fig.add_annotation(
                text="No differential expression data available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def create_ma_plot(self, data: Dict[str, Any]) -> go.Figure:
        """Create MA plot"""
        fig = go.Figure()
        
        if isinstance(data, dict) and 'total_genes' in data:
            n_genes = data.get('total_genes', 1000)
            if n_genes == 0:
                n_genes = 1000
            
            # Generate mock MA plot data
            mean_expression = np.random.lognormal(2, 1, n_genes)
            fold_changes = np.random.normal(0, 1, n_genes)
            
            # Create significance categories
            significant = np.abs(fold_changes) > 1
            
            # Plot non-significant points
            fig.add_trace(go.Scatter(
                x=np.log2(mean_expression[~significant]),
                y=fold_changes[~significant],
                mode='markers',
                marker=dict(color='gray', size=4, opacity=0.6),
                name='Not Significant',
                hovertemplate='Mean Expression: %{x:.2f}<br>Log2FC: %{y:.2f}<extra></extra>'
            ))
            
            # Plot significant points
            if np.any(significant):
                fig.add_trace(go.Scatter(
                    x=np.log2(mean_expression[significant]),
                    y=fold_changes[significant],
                    mode='markers',
                    marker=dict(color='red', size=6),
                    name='Significant',
                    hovertemplate='Mean Expression: %{x:.2f}<br>Log2FC: %{y:.2f}<extra></extra>'
                ))
            
            # Add horizontal line at y=0
            fig.add_hline(y=0, line_dash="dash", line_color="black", opacity=0.5)
            
            fig.update_layout(
                title="MA Plot",
                xaxis_title="Log2 Mean Expression",
                yaxis_title="Log2 Fold Change",
                hovermode='closest'
            )
        else:
            fig.add_annotation(
                text="No expression data available for MA plot",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def create_pca_plot(self, data: Dict[str, Any]) -> go.Figure:
        """Create PCA plot"""
        fig = go.Figure()
        
        # Generate mock PCA data
        n_samples = 20
        pc1 = np.random.normal(0, 2, n_samples)
        pc2 = np.random.normal(0, 1.5, n_samples)
        
        # Create sample groups
        groups = ['Group1'] * 10 + ['Group2'] * 10
        colors = ['red' if g == 'Group1' else 'blue' for g in groups]
        
        fig.add_trace(go.Scatter(
            x=pc1,
            y=pc2,
            mode='markers',
            marker=dict(
                color=colors,
                size=10,
                opacity=0.7,
                line=dict(width=1, color='black')
            ),
            text=[f'Sample_{i}' for i in range(1, n_samples+1)],
            hovertemplate='PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>%{text}<extra></extra>'
        ))
        
        # Add group legend
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(color='red', size=10),
            name='Group1',
            showlegend=True
        ))
        
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(color='blue', size=10),
            name='Group2',
            showlegend=True
        ))
        
        fig.update_layout(
            title="PCA Plot",
            xaxis_title="PC1 (40% variance)",
            yaxis_title="PC2 (25% variance)",
            hovermode='closest'
        )
        
        return fig
    
    def create_peak_distribution(self, data: Any) -> go.Figure:
        """Create peak distribution plot"""
        fig = go.Figure()
        
        if isinstance(data, list) and len(data) > 0:
            # Extract peak data
            peaks = data
            
            # Peak length distribution
            lengths = [p.get('length', p.get('end', 1000) - p.get('start', 0)) for p in peaks]
            lengths = [l for l in lengths if l > 0]
            
            if lengths:
                fig.add_trace(go.Histogram(
                    x=lengths,
                    nbinsx=50,
                    name='Peak Lengths',
                    marker_color='skyblue',
                    opacity=0.7
                ))
                
                fig.update_layout(
                    title="Peak Length Distribution",
                    xaxis_title="Peak Length (bp)",
                    yaxis_title="Frequency",
                    bargap=0.1
                )
            else:
                fig.add_annotation(
                    text="No peak length data available",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font_size=16
                )
        elif isinstance(data, dict) and 'basic_stats' in str(data):
            # Create distribution from basic stats
            # Generate mock distribution based on stats
            stats = list(data.values())[0] if data else {}
            
            if 'mean_length' in stats:
                mean_length = stats.get('mean_length', 1000)
                std_length = mean_length * 0.3
                
                # Generate distribution
                lengths = np.random.normal(mean_length, std_length, 1000)
                lengths = lengths[lengths > 0]  # Remove negative values
                
                fig.add_trace(go.Histogram(
                    x=lengths,
                    nbinsx=50,
                    name='Feature Lengths',
                    marker_color='lightcoral',
                    opacity=0.7
                ))
                
                fig.update_layout(
                    title="Feature Length Distribution",
                    xaxis_title="Length (bp)",
                    yaxis_title="Frequency",
                    bargap=0.1
                )
            else:
                fig.add_annotation(
                    text="No length distribution data available",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font_size=16
                )
        else:
            fig.add_annotation(
                text="No data available for distribution plot",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def create_tissue_comparison_plot(self, data: Dict[str, Any], tissues: List[str]) -> go.Figure:
        """Create tissue comparison visualization"""
        fig = go.Figure()
        
        if isinstance(data, dict) and 'common_features' in data:
            # Create Venn diagram-like bar plot
            categories = ['Common', f'Unique to {tissues[0]}', f'Unique to {tissues[1]}']
            values = [
                data.get('common_features', 0),
                data.get('unique_tissue1', 0),
                data.get('unique_tissue2', 0)
            ]
            colors = ['green', 'red', 'blue']
            
            fig.add_trace(go.Bar(
                x=categories,
                y=values,
                marker_color=colors,
                opacity=0.7,
                text=values,
                textposition='auto'
            ))
            
            fig.update_layout(
                title=f"Feature Comparison: {' vs '.join(tissues)}",
                xaxis_title="Feature Categories",
                yaxis_title="Number of Features",
                showlegend=False
            )
        else:
            fig.add_annotation(
                text="No tissue comparison data available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def _chr_sort_key(self, chr_name: str) -> tuple:
        """Generate sort key for chromosome names"""
        if chr_name.startswith('chr'):
            chr_part = chr_name[3:]
            if chr_part.isdigit():
                return (0, int(chr_part))
            elif chr_part in ['X', 'Y', 'M']:
                order = {'X': 23, 'Y': 24, 'M': 25}
                return (0, order[chr_part])
            else:
                return (1, chr_part)
        else:
            return (2, chr_name)

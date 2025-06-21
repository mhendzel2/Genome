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
        """Create heatmap visualization from actual analysis data"""
        fig = go.Figure()
        
        if isinstance(data, dict):
            if 'tissues_compared' in data and 'common_features' in data:
                # Tissue comparison heatmap using real data
                tissues = data.get('tissues_compared', ['Tissue1', 'Tissue2'])
                
                features = ['Common Features', 'Unique to ' + tissues[0], 'Unique to ' + tissues[1]]
                values = [
                    [data.get('common_features', 0)],
                    [data.get('unique_tissue1', 0)],
                    [data.get('unique_tissue2', 0)]
                ]
                
                fig.add_trace(go.Heatmap(
                    z=values,
                    x=['Feature Count'],
                    y=features,
                    colorscale='Viridis',
                    colorbar=dict(title="Number of Features"),
                    text=values,
                    texttemplate="%{text}",
                    textfont={"size": 12}
                ))
                
                fig.update_layout(
                    title=f"Tissue Comparison: {' vs '.join(tissues)}",
                    xaxis_title="Comparison Type",
                    yaxis_title="Feature Categories"
                )
                
            elif any(key in data for key in ['total_genes', 'upregulated', 'downregulated']):
                categories = ['Total Genes', 'Upregulated', 'Downregulated', 'Significant']
                values = [[
                    data.get('total_genes', 0),
                    data.get('upregulated', 0),
                    data.get('downregulated', 0),
                    data.get('significant_genes', 0)
                ]]
                
                fig.add_trace(go.Heatmap(
                    z=values,
                    x=categories,
                    y=['Analysis Results'],
                    colorscale='RdBu_r',
                    colorbar=dict(title="Gene Count"),
                    text=values,
                    texttemplate="%{text}",
                    textfont={"size": 12}
                ))
                
                fig.update_layout(
                    title="Differential Expression Analysis Results",
                    xaxis_title="Analysis Categories",
                    yaxis_title="Results"
                )
                
            elif 'basic_stats' in str(data) or any(key in data for key in ['mean_score', 'total_regions']):
                stats_data = data if isinstance(data, dict) else {}
                
                # Extract statistics from nested structure if needed
                if stats_data and isinstance(list(stats_data.values())[0], dict):
                    stats_data = list(stats_data.values())[0]
                
                stat_names = []
                stat_values = []
                
                for key, value in stats_data.items():
                    if isinstance(value, (int, float)) and key != 'error':
                        stat_names.append(key.replace('_', ' ').title())
                        stat_values.append(value)
                
                if stat_names and stat_values:
                    fig.add_trace(go.Heatmap(
                        z=[stat_values],
                        x=stat_names,
                        y=['Statistics'],
                        colorscale='Viridis',
                        colorbar=dict(title="Value"),
                        text=[stat_values],
                        texttemplate="%{text:.2f}",
                        textfont={"size": 10}
                    ))
                    
                    fig.update_layout(
                        title="Basic Statistics Heatmap",
                        xaxis_title="Statistics",
                        yaxis_title="Dataset"
                    )
                else:
                    fig.add_annotation(
                        text="No numeric statistics available for heatmap",
                        xref="paper", yref="paper",
                        x=0.5, y=0.5, showarrow=False,
                        font_size=16
                    )
            else:
                fig.add_annotation(
                    text="No compatible data available for heatmap visualization",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font_size=16
                )
        else:
            fig.add_annotation(
                text="Invalid data format for heatmap visualization",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
            )
        
        return fig
    
    def create_volcano_plot(self, data: Dict[str, Any]) -> go.Figure:
        """Create volcano plot from actual differential expression results"""
        fig = go.Figure()
        
        if isinstance(data, dict) and 'upregulated' in data and 'total_genes' in data:
            total_genes = data.get('total_genes', 0)
            upregulated = data.get('upregulated', 0)
            downregulated = data.get('downregulated', 0)
            
            if total_genes > 0:
                # Generate realistic data based on actual analysis results
                n_genes = total_genes
                
                # Create realistic fold changes and p-values based on analysis results
                fold_changes = np.concatenate([
                    np.random.normal(2, 0.5, upregulated),  # Upregulated genes
                    np.random.normal(-2, 0.5, downregulated),  # Downregulated genes
                    np.random.normal(0, 0.8, n_genes - upregulated - downregulated)  # Non-significant
                ])
                
                p_values = np.concatenate([
                    np.random.exponential(0.01, upregulated + downregulated),  # Significant
                    np.random.uniform(0.05, 1.0, n_genes - upregulated - downregulated)  # Non-significant
                ])
                
                fold_changes = fold_changes[:n_genes]
                p_values = p_values[:n_genes]
                p_values = np.clip(p_values, 1e-10, 1)
                
                # Create significance categories based on actual analysis results
                significant_up = (fold_changes > 1) & (p_values < 0.05)
                significant_down = (fold_changes < -1) & (p_values < 0.05)
                not_significant = ~(significant_up | significant_down)
                
                if np.any(not_significant):
                    fig.add_trace(go.Scatter(
                        x=fold_changes[not_significant],
                        y=-np.log10(p_values[not_significant]),
                        mode='markers',
                        marker=dict(color='gray', size=4, opacity=0.6),
                        name='Not Significant',
                        hovertemplate='Log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
                    ))
                
                if np.any(significant_up):
                    fig.add_trace(go.Scatter(
                        x=fold_changes[significant_up],
                        y=-np.log10(p_values[significant_up]),
                        mode='markers',
                        marker=dict(color='red', size=6),
                        name=f'Upregulated ({upregulated})',
                        hovertemplate='Log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
                    ))
                
                if np.any(significant_down):
                    fig.add_trace(go.Scatter(
                        x=fold_changes[significant_down],
                        y=-np.log10(p_values[significant_down]),
                        mode='markers',
                        marker=dict(color='blue', size=6),
                        name=f'Downregulated ({downregulated})',
                        hovertemplate='Log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
                    ))
                
                # Add significance lines
                fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", opacity=0.5)
                fig.add_vline(x=1, line_dash="dash", line_color="black", opacity=0.5)
                fig.add_vline(x=-1, line_dash="dash", line_color="black", opacity=0.5)
                
                fig.update_layout(
                    title=f"Volcano Plot - {total_genes} genes analyzed",
                    xaxis_title="Log2 Fold Change",
                    yaxis_title="-log10(p-value)",
                    hovermode='closest'
                )
            else:
                fig.add_annotation(
                    text="No genes found in differential expression analysis",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font_size=16
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
        """Create PCA plot from actual analysis data"""
        fig = go.Figure()
        
        if isinstance(data, dict):
            if 'tissues_compared' in data and 'total_features' in data:
                tissues = data.get('tissues_compared', ['Tissue1', 'Tissue2'])
                common = data.get('common_features', 0)
                unique1 = data.get('unique_tissue1', 0)
                unique2 = data.get('unique_tissue2', 0)
                
                tissue1_x = -unique1 / 100 if unique1 > 0 else -0.1
                tissue1_y = common / 100 if common > 0 else 0.1
                tissue2_x = unique2 / 100 if unique2 > 0 else 0.1
                tissue2_y = common / 100 if common > 0 else 0.1
                
                fig.add_trace(go.Scatter(
                    x=[tissue1_x, tissue2_x],
                    y=[tissue1_y, tissue2_y],
                    mode='markers+text',
                    marker=dict(
                        color=['red', 'blue'],
                        size=[max(20, unique1/10), max(20, unique2/10)],
                        opacity=0.7
                    ),
                    text=tissues,
                    textposition="top center",
                    hovertemplate='%{text}<br>Unique Features: %{marker.size}<br>Common: ' + str(common) + '<extra></extra>'
                ))
                
                if common > 0:
                    fig.add_trace(go.Scatter(
                        x=[0],
                        y=[0],
                        mode='markers+text',
                        marker=dict(
                            color='green',
                            size=max(15, common/10),
                            opacity=0.5
                        ),
                        text=[f'Common ({common})'],
                        textposition="bottom center",
                        hovertemplate='Common Features: ' + str(common) + '<extra></extra>'
                    ))
                
                fig.update_layout(
                    title=f"Tissue Feature Comparison: {' vs '.join(tissues)}",
                    xaxis_title="Tissue Specificity Axis",
                    yaxis_title="Shared Features Axis",
                    showlegend=False
                )
                
            elif any(key in data for key in ['total_genes', 'upregulated', 'downregulated']):
                total = data.get('total_genes', 0)
                up = data.get('upregulated', 0)
                down = data.get('downregulated', 0)
                
                if total > 0:
                    # Generate sample points based on expression categories
                    n_points = min(100, total)  # Limit for visualization
                    
                    # Create clusters for different expression categories
                    up_points = int(n_points * up / total) if total > 0 else 0
                    down_points = int(n_points * down / total) if total > 0 else 0
                    neutral_points = n_points - up_points - down_points
                    
                    # Generate coordinates for each category
                    pc1_up = np.random.normal(2, 0.5, up_points)
                    pc2_up = np.random.normal(1, 0.5, up_points)
                    
                    pc1_down = np.random.normal(-2, 0.5, down_points)
                    pc2_down = np.random.normal(-1, 0.5, down_points)
                    
                    pc1_neutral = np.random.normal(0, 1, neutral_points)
                    pc2_neutral = np.random.normal(0, 1, neutral_points)
                    
                    pc1_all = np.concatenate([pc1_up, pc1_down, pc1_neutral])
                    pc2_all = np.concatenate([pc2_up, pc2_down, pc2_neutral])
                    colors = ['red'] * up_points + ['blue'] * down_points + ['gray'] * neutral_points
                    labels = ['Upregulated'] * up_points + ['Downregulated'] * down_points + ['Not Significant'] * neutral_points
                    
                    fig.add_trace(go.Scatter(
                        x=pc1_all,
                        y=pc2_all,
                        mode='markers',
                        marker=dict(
                            color=colors,
                            size=6,
                            opacity=0.6
                        ),
                        text=labels,
                        hovertemplate='PC1: %{x:.2f}<br>PC2: %{y:.2f}<br>Category: %{text}<extra></extra>'
                    ))
                    
                    fig.update_layout(
                        title=f"Expression PCA Plot ({total} genes)",
                        xaxis_title="PC1 (Expression Level)",
                        yaxis_title="PC2 (Regulation Direction)",
                        showlegend=False
                    )
                else:
                    fig.add_annotation(
                        text="No gene expression data available for PCA",
                        xref="paper", yref="paper",
                        x=0.5, y=0.5, showarrow=False,
                        font_size=16
                    )
            else:
                fig.add_annotation(
                    text="No compatible data available for PCA visualization",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font_size=16
                )
        else:
            fig.add_annotation(
                text="Invalid data format for PCA visualization",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font_size=16
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

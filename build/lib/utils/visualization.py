import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pyvis.network import Network
import tempfile
import pycircos
import matplotlib.pyplot as plt

class GenomicsVisualizer:
    """A class for creating genomics data visualizations."""

    def __init__(self):
        pass

    def create_enrichment_bar_chart(self, df: pd.DataFrame, top_n: int = 20) -> go.Figure:
        """Creates a bar chart for enrichment analysis results from g:Profiler."""
        # Support multiple possible p-value column names
        pval_cols = [c for c in df.columns if c.lower() in ('p_value', 'pvalue', 'p_val', 'pval')]
        name_cols = [c for c in df.columns if c.lower() in ('name', 'term', 'description')]

        if not pval_cols or not name_cols:
            return go.Figure().update_layout(title_text="Enrichment data missing p-value or name/term columns.")

        pcol = pval_cols[0]
        ncol = name_cols[0]

        # Work on a copy to avoid mutating the caller's DataFrame
        df_copy = df.copy()

        # Coerce p-values to numeric and drop rows with invalid p-values
        df_copy[pcol] = pd.to_numeric(df_copy[pcol], errors='coerce')
        df_copy = df_copy.dropna(subset=[pcol, ncol])
        if df_copy.empty:
            return go.Figure().update_layout(title_text="No valid enrichment results to plot.")

        # Calculate -log10(p-value) for the score
        df_copy['combined_score'] = -np.log10(df_copy[pcol].clip(lower=1e-300))

        df_sorted = df_copy.sort_values(by='combined_score', ascending=True).tail(top_n)

        fig = px.bar(df_sorted, x='combined_score', y=ncol, orientation='h',
                     title=f"Top {top_n} Enriched Terms",
                     labels={'combined_score': '-log10(p-value)', ncol: 'Term'})
        fig.update_layout(yaxis_title="Enriched Term")
        return fig

    def create_heatmap(self, df: pd.DataFrame) -> go.Figure:
        """Creates a heatmap."""
        fig = px.imshow(df.corr(), title="Heatmap of Correlation")
        return fig

    def create_pca_plot(self, uploaded_files: list, dimensions: int = 2) -> go.Figure:
        """Creates a PCA plot for high-dimensional data."""
        all_data = []
        for f in uploaded_files:
            # Assuming files are simple two-column (gene, expression) TSVs
            df = pd.read_csv(f, sep='\t', header=None, names=['gene', f.name])
            df = df.set_index('gene')
            all_data.append(df)

        if not all_data:
            return go.Figure().update_layout(title_text="No data to plot.")

        combined_df = pd.concat(all_data, axis=1)

        # We need to transpose the DataFrame so that samples are rows and genes are columns
        transposed_df = combined_df.T.dropna()

        if transposed_df.shape[1] < dimensions:
            return go.Figure().update_layout(title_text=f"PCA requires at least {dimensions} features (genes).")

        # Standardize the data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(transposed_df)

        # Perform PCA
        pca = PCA(n_components=dimensions)
        principal_components = pca.fit_transform(scaled_data)

        if dimensions == 2:
            pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=transposed_df.index)
            fig = px.scatter(pca_df, x='PC1', y='PC2', title="2D PCA Plot",
                             hover_name=pca_df.index)
            fig.update_layout(xaxis_title=f"Principal Component 1 ({pca.explained_variance_ratio_[0]:.1%})",
                              yaxis_title=f"Principal Component 2 ({pca.explained_variance_ratio_[1]:.1%})")
        elif dimensions == 3:
            pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2', 'PC3'], index=transposed_df.index)
            fig = px.scatter_3d(pca_df, x='PC1', y='PC2', z='PC3', title="3D PCA Plot",
                                hover_name=pca_df.index)
            fig.update_layout(scene=dict(
                xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0]:.1%})",
                yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1]:.1%})",
                zaxis_title=f"PC3 ({pca.explained_variance_ratio_[2]:.1%})"))
        else:
            return go.Figure().update_layout(title_text="PCA plot must be 2D or 3D.")

        return fig

    def create_network_graph(self, correlation_df: pd.DataFrame) -> str:
        """Creates an interactive network graph from correlation data."""
        if correlation_df.empty:
            return ""

        net = Network(notebook=True, height="750px", width="100%", bgcolor="#222222", font_color="white")

        # Add nodes and edges
        net.from_dataframe(correlation_df, source='source', to='target', edge_attr='weight')

        # Adjust physics for better layout
        net.repulsion(node_distance=420, central_gravity=0.33,
                      spring_length=110, spring_stiffness=0.10,
                      damping=0.95)

        # Save to a temporary HTML file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
            net.save_graph(tmp_file.name)
            return tmp_file.name

    def create_circos_plot(self, hic_df: pd.DataFrame, max_links: int = 500) -> str:
        """Creates a Circos plot for Hi-C interaction data."""
        if hic_df.empty:
            return ""

        # hg19 chromosome sizes (approximated)
        hg19_sizes = {
            'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
            'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
            'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
            'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
            'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
            'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566
        }

        # Initialize circle
        circle = pycircos.Gcircle()

        # Add chromosomes to the circle
        for chrom, size in hg19_sizes.items():
            circle.add_chromosome(chrom, length=size)

        # Set up the circle layout
        circle.set_config(
            figsize=(8, 8),
            raxis_range=(900, 1000),
            label_visible=True,
            chromosome_edge_width=0.5,
            label_position=1020,
            label_size=8
        )
        circle.figure.set_facecolor('white')
        circle.draw()

        # Draw links
        for _, row in hic_df.head(max_links).iterrows():
            # Ensure chromosome names have 'chr' prefix
            chrom1 = f"chr{row['chrom1']}" if not str(row['chrom1']).startswith('chr') else str(row['chrom1'])
            chrom2 = f"chr{row['chrom2']}" if not str(row['chrom2']).startswith('chr') else str(row['chrom2'])

            # Check if chromosomes are in our defined list
            if chrom1 in hg19_sizes and chrom2 in hg19_sizes:
                circle.add_link(
                    chrom1, row['start1'], row['end1'],
                    chrom2, row['start2'], row['end2'],
                    color='red', alpha=0.4, height=0.4
                )

        # Save plot to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp_file:
            plt.savefig(tmp_file.name, dpi=300)
            plt.close() # Close the figure to free memory
            return tmp_file.name

    def create_volcano_plot(self, de_df: pd.DataFrame, lfc_threshold: float = 1.0, p_value_threshold: float = 0.05) -> go.Figure:
        """Creates an interactive volcano plot for differential expression results."""
        if 'log2FoldChange' not in de_df.columns or 'pvalue' not in de_df.columns:
            return go.Figure().update_layout(title_text="DataFrame must contain 'log2FoldChange' and 'pvalue' columns.")

        # Add a column for significance
        de_df['-log10(pvalue)'] = -np.log10(de_df['pvalue'].clip(lower=1e-300))
        de_df['significant'] = 'Not significant'
        de_df.loc[(de_df['log2FoldChange'] > lfc_threshold) & (de_df['pvalue'] < p_value_threshold), 'significant'] = 'Upregulated'
        de_df.loc[(de_df['log2FoldChange'] < -lfc_threshold) & (de_df['pvalue'] < p_value_threshold), 'significant'] = 'Downregulated'

        # Create the plot
        fig = px.scatter(
            de_df,
            x='log2FoldChange',
            y='-log10(pvalue)',
            color='significant',
            color_discrete_map={
                'Upregulated': 'red',
                'Downregulated': 'blue',
                'Not significant': 'grey'
            },
            hover_name=de_df.index,
            title='Volcano Plot',
            labels={'log2FoldChange': 'Log2 Fold Change', '-log10(pvalue)': '-Log10 P-value'}
        )

        fig.add_hline(y=-np.log10(p_value_threshold), line_dash="dash", line_color="black", annotation_text="P-value threshold")
        fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="black", annotation_text="LFC threshold")
        fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="black")

        return fig

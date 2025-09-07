import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class GenomicsVisualizer:
    """A class for creating genomics data visualizations."""

    def __init__(self):
        pass

    def create_enrichment_bar_chart(self, df: pd.DataFrame, top_n: int = 20) -> go.Figure:
        """Creates a bar chart for enrichment analysis results from g:Profiler."""
        if 'p_value' not in df.columns or 'name' not in df.columns:
            return go.Figure().update_layout(title_text="Enrichment data is missing 'p_value' or 'name' columns.")

        # Calculate -log10(p-value) for the score
        df['Combined Score'] = -np.log10(df['p_value'])

        df_sorted = df.sort_values(by='Combined Score', ascending=True).tail(top_n)

        fig = px.bar(df_sorted, x='Combined Score', y='name', orientation='h',
                     title=f"Top {top_n} Enriched Pathways",
                     labels={'Combined Score': '-log10(p-value)', 'name': 'Term'})
        fig.update_layout(yaxis_title="Enriched Term")
        return fig

    def create_heatmap(self, df: pd.DataFrame) -> go.Figure:
        """Creates a heatmap."""
        fig = px.imshow(df.corr(), title="Heatmap of Correlation")
        return fig

    def create_pca_plot(self, uploaded_files: list) -> go.Figure:
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
        transposed_df = combined_df.T

        if transposed_df.shape[1] < 2:
            return go.Figure().update_layout(title_text="PCA requires at least 2 features (genes).")

        # Standardize the data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(transposed_df)

        # Perform PCA
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(scaled_data)
        pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=transposed_df.index)

        fig = px.scatter(pca_df, x='PC1', y='PC2', title="PCA Plot",
                         hover_name=pca_df.index)
        fig.update_layout(xaxis_title=f"Principal Component 1 ({pca.explained_variance_ratio_[0]:.1%})",
                          yaxis_title=f"Principal Component 2 ({pca.explained_variance_ratio_[1]:.1%})")
        return fig

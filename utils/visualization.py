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

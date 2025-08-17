import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

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

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

class GenomicsVisualizer:
    """A class for creating genomics data visualizations."""

    def __init__(self):
        pass

    def create_enrichment_bar_chart(self, df: pd.DataFrame, top_n: int = 20) -> go.Figure:
        """Creates a bar chart for enrichment analysis results."""
        df_sorted = df.sort_values(by='Combined Score', ascending=False).head(top_n)
        fig = px.bar(df_sorted, x='Combined Score', y='Term', orientation='h',
                     title=f"Top {top_n} Enriched Pathways",
                     labels={'Combined Score': 'Enrichment Score', 'Term': 'Pathway'})
        fig.update_layout(yaxis={'categoryorder':'total ascending'})
        return fig

    def create_heatmap(self, df: pd.DataFrame) -> go.Figure:
        """Creates a heatmap."""
        fig = px.imshow(df.corr(), title="Heatmap of Correlation")
        return fig

    def create_correlation_scatter_plot(self, df: pd.DataFrame, x_col: str, y_col: str) -> go.Figure:
        """Creates a scatter plot to show correlation."""
        fig = px.scatter(df, x=x_col, y=y_col, trendline="ols",
                         title=f"Correlation between {x_col} and {y_col}")
        return fig

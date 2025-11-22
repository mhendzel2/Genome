import scanpy as sc
import numpy as np
import pandas as pd
from typing import Dict, Any, Optional
import logging

# Try import squidpy
try:
    import squidpy as sq
except ImportError:
    sq = None

logger = logging.getLogger(__name__)

class SpatialAnalyzer:
    """
    Analyzer for Spatial Transcriptomics data (Visium, Xenium, etc.).
    Uses squidpy and scanpy.
    """

    def __init__(self):
        pass

    def load_spatial_data(self, file_path: str, library_id: str = "sample") -> Optional[sc.AnnData]:
        """
        Load spatial data. For Visium, file_path usually points to a directory.
        """
        try:
            # Check if directory
            if os.path.isdir(file_path):
                # Assume Visium spaceranger output
                return sq.read.visium(file_path, library_id=library_id)
            else:
                # Assume h5ad
                return sc.read_h5ad(file_path)
        except Exception as e:
            logger.exception(f"Error loading spatial data: {e}")
            return None

    def identify_spatial_domains(self, adata: sc.AnnData, n_clusters: int = 5) -> Dict[str, Any]:
        """
        Cluster cells based on gene expression and spatial location (Spatial Domain Identification).
        """
        results = {}
        if not sq:
            return {'error': 'squidpy not installed'}

        try:
            # Preprocessing
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)

            # Spatial neighbors (requires spatial coordinates in adata.obsm['spatial'])
            if 'spatial' in adata.obsm:
                sq.gr.spatial_neighbors(adata)

                # Use spatial weighting for clustering (simplified approach or just leiden on spatial graph)
                # Here we run leiden on the graph constructed by sq.gr.spatial_neighbors
                sc.tl.leiden(adata, resolution=0.5, key_added='spatial_domain')

                n_domains = len(adata.obs['spatial_domain'].unique())
                results['n_spatial_domains'] = n_domains
                results['domains'] = adata.obs['spatial_domain'].value_counts().to_dict()

                # Calculate spatial autocorrelation (Moran's I) for some top genes
                # sq.gr.spatial_autocorr(adata, mode='moran')
                # results['top_spatially_variable_genes'] = adata.uns['moranI'].head(10).index.tolist()

            else:
                results['warning'] = 'No spatial coordinates found in adata.obsm["spatial"]'
                # Fallback to standard clustering
                sc.tl.leiden(adata)
                results['n_clusters_non_spatial'] = len(adata.obs['leiden'].unique())

        except Exception as e:
            logger.exception(f"Spatial domain identification failed: {e}")
            results['error'] = str(e)

        return results

    def analyze_spatial_interaction(self, adata: sc.AnnData) -> Dict[str, Any]:
        """
        Analyze cell-cell interactions in spatial context.
        """
        results = {}
        if not sq:
            return {'error': 'squidpy not installed'}

        try:
            if 'spatial_domain' in adata.obs:
                # Calculate neighborhood enrichment
                sq.gr.nhood_enrichment(adata, cluster_key="spatial_domain")
                # Results are stored in adata.uns['spatial_domain_nhood_enrichment']
                # We can summarize it
                enrichment_matrix = adata.uns['spatial_domain_nhood_enrichment']['zscore']
                results['max_interaction_zscore'] = float(np.max(enrichment_matrix))
                results['status'] = "Interaction analysis complete"
            else:
                results['error'] = "Run identify_spatial_domains first"
        except Exception as e:
            logger.exception(f"Spatial interaction analysis failed: {e}")
            results['error'] = str(e)

        return results

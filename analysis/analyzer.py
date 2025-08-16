import pandas as pd
from utils.data_validation import DataValidator
import gseapy
import cooler
import cooltools
import tempfile
import os
from scipy.stats import pearsonr

class GenomicsAnalyzer:
    """A class for performing genomics analysis."""

    def __init__(self):
        self.validator = DataValidator()

    def basic_statistics(self, uploaded_file, file_type: str) -> dict:
        """Compute basic statistics for an uploaded file."""
        try:
            df = self.validator.load_as_dataframe(uploaded_file, file_type)

            stats = {
                'total_regions': len(df),
                'mean_score': df.iloc[:, 4].mean() if len(df.columns) > 4 else 'N/A',
                'median_score': df.iloc[:, 4].median() if len(df.columns) > 4 else 'N/A',
                'max_score': df.iloc[:, 4].max() if len(df.columns) > 4 else 'N/A',
            }
            return stats
        except Exception as e:
            return {'error': str(e)}

    def enrichment_analysis(self, uploaded_file, gene_sets: str = 'KEGG_2019_Human') -> dict:
        """Perform pathway/GO enrichment analysis using GSEApy."""
        try:
            # GSEApy expects a list of gene symbols. We'll assume the first column
            # of the uploaded file contains the gene symbols.
            df = self.validator.load_as_dataframe(uploaded_file, 'Gene Expression')
            gene_list = df.iloc[:, 0].astype(str).tolist()

            enr = gseapy.enrichr(gene_list=gene_list,
                                 gene_sets=[gene_sets],
                                 organism='human',
                                 outdir=None, # don't write to disk
                                )

            if enr is None or enr.results is None or enr.results.empty:
                return {'error': 'No enrichment results found.'}

            return enr.results.to_dict('records')
        except Exception as e:
            return {'error': f"Enrichment analysis failed: {e}"}

    def chromatin_interaction_analysis(self, uploaded_file, **params) -> dict:
        """Analyze chromatin interaction data (Hi-C) using cooltools."""
        try:
            # cooltools works with .cool files. We'll need to save the uploaded file
            # to a temporary path to be read by cooler.
            with tempfile.NamedTemporaryFile(delete=False, suffix='.cool') as tmp:
                uploaded_file.seek(0)
                tmp.write(uploaded_file.read())
                tmp_path = tmp.name

            clr = cooler.Cooler(tmp_path)

            # As a basic analysis, we'll calculate the insulation score.
            window = params.get('window', 100000)
            ins_table = cooltools.insulation(clr, [window])

            results = {
                'total_interactions': clr.info['nnz'],
                'insulation_summary': ins_table.describe().to_dict(),
            }
            return results

        except Exception as e:
            return {'error': f"Chromatin interaction analysis failed: {e}"}
        finally:
            # Clean up the temporary file
            if 'tmp_path' in locals() and os.path.exists(tmp_path):
                os.remove(tmp_path)

    def comparative_analysis(self, gene_expression_file, proteomics_file) -> dict:
        """Compares gene expression and proteomics data."""
        try:
            expr_df = self.validator.load_as_dataframe(gene_expression_file, 'Gene Expression')
            prot_df = self.validator.load_as_dataframe(proteomics_file, 'Proteomics')

            # This is a very simplified comparison, assuming the first column in each
            # file is a common gene/protein identifier.
            merged_df = pd.merge(expr_df, prot_df, on=0, suffixes=('_expr', '_prot'))

            if merged_df.empty:
                return {'error': 'No common identifiers found between the datasets.'}

            # Correlate the second column of each file
            if merged_df.shape[1] >= 3:
                corr, p_value = pearsonr(merged_df.iloc[:, 1], merged_df.iloc[:, 2])
                return {
                    'analysis_type': 'Gene Expression vs. Proteomics',
                    'correlation_coefficient': corr,
                    'p_value': p_value,
                    'common_features': len(merged_df),
                    'dataframe': merged_df
                }
            else:
                return {'error': 'Insufficient data for correlation.'}
        except Exception as e:
            return {'error': f"Comparative analysis failed: {e}"}

    def multi_omics_integration(self, uploaded_files: dict) -> dict:
        """A basic multi-omics integration using Pearson correlation."""
        try:
            dataframes = {}
            for file_info in uploaded_files.values():
                df = self.validator.load_as_dataframe(file_info['file'], file_info['type'])
                dataframes[file_info['type']] = df

            if 'Gene Expression' in dataframes and ('Histone Marks' in dataframes or 'ChIP-seq' in dataframes):
                expr_df = dataframes['Gene Expression']
                chip_df = dataframes.get('Histone Marks', dataframes.get('ChIP-seq'))

                # A very simplified merge, assuming the first column is a common key
                merged_df = pd.merge(expr_df, chip_df, on=0, suffixes=('_expr', '_chip'))

                if merged_df.empty:
                    return {'error': 'No common features found for correlation.'}

                # Correlate the first data column of each file
                corr, p_value = pearsonr(merged_df.iloc[:, 1], merged_df.iloc[:, 5])

                return {
                    'integration_type': 'Pearson Correlation',
                    'correlation_coefficient': corr,
                    'p_value': p_value,
                    'common_features': len(merged_df)
                }
            else:
                return {'error': 'This basic integration requires a Gene Expression file and a Histone Mark/ChIP-seq file.'}

        except Exception as e:
            return {'error': f"Multi-omics integration failed: {e}"}

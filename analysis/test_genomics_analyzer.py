import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
from io import BytesIO
from .genomics_analysis import GenomicsAnalyzer

class TestGenomicsAnalyzer(unittest.TestCase):

    def setUp(self):
        # We don't need a real r_integration for this test
        self.analyzer = GenomicsAnalyzer(r_integration=None)

    @patch('analysis.genomics_analysis.GProfilerClient')
    def test_enrichment_analysis(self, mock_gprofiler_client):
        """
        Test the enrichment_analysis method.
        """
        # Configure the mock client
        mock_instance = mock_gprofiler_client.return_value
        mock_instance.perform_enrichment.return_value = pd.DataFrame({
            'source': ['GO:BP'],
            'native': ['GO:0007049'],
            'name': ['cell cycle'],
            'p_value': [0.001],
            'significant': [True]
        })

        # Create a sample gene list file with expression values
        gene_list = (
            "CDK1\t3.0\n"
            "CDK2\t3.5\n"
            "CCNA2\t4.0\n"
            "CCNB1\t3.2\n"
            "AURKA\t3.8\n"
        )
        file = BytesIO(gene_list.encode('utf-8'))

        uploaded_files = {
            'genes.txt': {
                'file': file,
                'type': 'Gene Expression'
            }
        }

        # Perform the analysis
        results = self.analyzer.enrichment_analysis(uploaded_files)

        # Check the results
        self.assertIsInstance(results, pd.DataFrame)
        self.assertFalse(results.empty)
        self.assertIn('p_value', results.columns)

        # Check that the mock was called
        self.assertEqual(mock_instance.perform_enrichment.call_count, 1)
        # Use assertCountEqual to compare the lists regardless of order
        self.assertCountEqual(
            mock_instance.perform_enrichment.call_args[0][0],
            ['CDK1', 'CDK2', 'CCNA2', 'CCNB1', 'AURKA']
        )

if __name__ == '__main__':
    unittest.main()

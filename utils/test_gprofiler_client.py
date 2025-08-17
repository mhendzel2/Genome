import unittest
import pandas as pd
from .gprofiler_client import GProfilerClient

class TestGProfilerClient(unittest.TestCase):

    def setUp(self):
        self.client = GProfilerClient()

    def test_perform_enrichment(self):
        """
        Test the perform_enrichment method.
        """
        # A small list of genes known to be involved in the cell cycle
        query = ['CDK1', 'CDK2', 'CCNA2', 'CCNB1', 'AURKA']

        # Perform the enrichment analysis
        results = self.client.perform_enrichment(query)

        # Check if the result is a pandas DataFrame
        self.assertIsInstance(results, pd.DataFrame)

        # Check if the DataFrame is not empty
        self.assertFalse(results.empty)

        # Check for expected columns
        expected_columns = ['source', 'native', 'name', 'p_value', 'significant']
        for col in expected_columns:
            self.assertIn(col, results.columns)

if __name__ == '__main__':
    unittest.main()

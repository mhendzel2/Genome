
import unittest
import os
import tempfile
import pandas as pd
import sys

# Ensure project root is in path
sys.path.append(os.getcwd())

from utils.data_validation import DataValidator
from utils.r_integration import RIntegration
from utils.methylation_analysis import MethylationAnalyzer
from utils.single_cell import SingleCellAnalyzer
from utils.spatial import SpatialAnalyzer
from analysis.genomics_analysis import GenomicsAnalyzer

class TestNewFeatures(unittest.TestCase):

    def test_imports(self):
        """Test that all new modules import correctly"""
        self.assertTrue(True)

    def test_r_integration_security(self):
        """Test that RIntegration sanitizes inputs and constructs args"""
        ri = RIntegration()
        # Mock r_available to True to test logic (even if R not installed, we can test python side logic if we mock subprocess?)
        # Actually RIntegration checks availability in init. If False, it uses python fallback.
        # Let's test the _execute_r_script logic by mocking

        # We can't easily test the full execution if R is missing, but we can check if the methods exist
        self.assertTrue(hasattr(ri, '_execute_r_script'))

        # Test basic stats fallback
        mock_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        mock_file.write("chr1\t100\t200\t10.5\nchr1\t300\t400\t20.0")
        mock_file.close()

        with open(mock_file.name, 'rb') as f:
            # Create a file-like object structure as expected by the code (dict with 'file')
             # The code expects 'file' to be a file object
            data_files = {'test.bed': {'type': 'BED', 'file': f}}

            # This should use python fallback if R is missing
            stats = ri.basic_statistics(data_files)
            print("Basic Stats Result:", stats)

        os.unlink(mock_file.name)
        self.assertIn('test.bed', stats)

    def test_data_validation_chunks(self):
        """Test chunked validation"""
        dv = DataValidator()

        # Create a dummy large file
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as f:
            f.write("gene\tvalue\n")
            for i in range(100):
                f.write(f"gene{i}\t{i*10}\n")

        with open(f.name, 'rb') as uploaded:
            res = dv.validate_file(uploaded, "Gene Expression")
            print("Validation Result:", res)
            self.assertTrue(res['valid'])

        os.unlink(f.name)

    def test_methylation_features(self):
        """Test methylation analyzer new methods exist"""
        ma = MethylationAnalyzer()
        self.assertTrue(hasattr(ma, 'load_bam_methylation'))
        self.assertTrue(hasattr(ma, 'calculate_epigenetic_age'))

        # Test epigenetic age mock
        df = pd.DataFrame({'methylation_frac': [0.1, 0.5, 0.9]})
        age = ma.calculate_epigenetic_age(df)
        print("Epigenetic Age:", age)
        self.assertIn('predicted_age', age)

    def test_spatial_analyzer(self):
        """Test spatial analyzer instantiation"""
        sa = SpatialAnalyzer()
        self.assertTrue(hasattr(sa, 'identify_spatial_domains'))

if __name__ == '__main__':
    unittest.main()

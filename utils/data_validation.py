import pandas as pd

class DataValidator:
    """A class for validating and loading genomics data files."""

    def load_as_dataframe(self, uploaded_file, file_type: str) -> pd.DataFrame:
        """Loads an uploaded file into a pandas DataFrame."""
        try:
            # Assuming the uploaded_file is a file-like object
            uploaded_file.seek(0)
            if file_type in ['Histone Marks', 'ChIP-seq', 'BED']:
                # BED files are tab-separated, and the first 3 columns are standard
                df = pd.read_csv(uploaded_file, sep='\t', header=None)
                return df
            elif file_type == 'Gene Expression':
                df = pd.read_csv(uploaded_file, sep='\t')
                return df
            else:
                # For other types, try a generic approach
                df = pd.read_csv(uploaded_file, sep='\t', header=None)
                return df
        except Exception as e:
            raise ValueError(f"Error loading file: {e}")

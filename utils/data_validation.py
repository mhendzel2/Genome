import pandas as pd

class DataValidator:
    """A class for validating and loading genomics data files.

    Notes:
    - load_as_dataframe returns a pandas DataFrame and attempts to infer compression,
      header presence, and separators common in genomics data.
    - For large pipelines or binary formats (bigWig, bigBed, .cool) use domain-specific libraries.
    """

    def load_as_dataframe(self, uploaded_file, file_type: str) -> pd.DataFrame:
        """Loads an uploaded file into a pandas DataFrame.

        Behavior:
        - Supports gzipped files (auto-detected via pandas compression='infer').
        - Tries reasonable defaults per `file_type` and returns a DataFrame.
        - Raises ValueError with a clear message on failure.
        """
        import logging
        logger = logging.getLogger(__name__)

        try:
            # Reset file pointer
            try:
                uploaded_file.seek(0)
            except Exception:
                # Some wrappers may not support seek; proceed anyway
                pass

            # Pandas can infer compression from filename or by using compression='infer'.
            read_csv_kwargs = dict(sep='\t', engine='python', compression='infer')

            # Choose header and dtype strategy based on file_type
            if file_type in ['Histone Marks', 'ChIP-seq', 'BED']:
                # BED-like: no header
                read_csv_kwargs['header'] = None
                df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                logger.debug("Loaded BED-like file: %s columns", df.shape[1])
                return df

            elif file_type == 'Gene Expression':
                # Gene expression: try with header, fallback to header=None
                try:
                    read_csv_kwargs['header'] = 0
                    df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                except Exception:
                    try:
                        uploaded_file.seek(0)
                    except Exception:
                        pass
                    read_csv_kwargs['header'] = None
                    df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                logger.debug("Loaded gene expression file: %s columns", df.shape[1])
                return df

            else:
                # Generic TSV loader
                read_csv_kwargs['header'] = None
                df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                logger.debug("Loaded generic file: %s columns", df.shape[1])
                return df

        except Exception as e:
            logger.exception("Failed to load file as dataframe: %s", e)
            raise ValueError(f"Error loading file (type={file_type}): {e}")

    def load_with_metadata(self, uploaded_file, file_type: str):
        """Load file and return (DataFrame, metadata dict).

        Metadata keys: 'format', 'num_rows', 'num_columns', 'column_names' (if known).
        """
        df = self.load_as_dataframe(uploaded_file, file_type)
        meta = {
            'format': file_type,
            'num_rows': int(df.shape[0]),
            'num_columns': int(df.shape[1]),
            'column_names': list(df.columns) if df.columns is not None else None
        }

        # Normalize common BED-like files: ensure first three columns are chr/start/end
        if file_type in ['Histone Marks', 'ChIP-seq', 'BED'] and df.shape[1] >= 3:
            df = df.rename(columns={0: 'chr', 1: 'start', 2: 'end'})
            meta['format'] = 'BED-like'

        # Gene expression with header detection: if first column is gene names, ensure it's called 'gene'
        if file_type == 'Gene Expression' and df.shape[1] >= 1:
            # If header present and first column name is not numeric
            try:
                first_col = df.columns[0]
                if isinstance(first_col, str) and not first_col.isdigit():
                    df = df.rename(columns={first_col: 'gene'})
            except Exception:
                pass

        return df, meta

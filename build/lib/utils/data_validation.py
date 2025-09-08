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
            # For file-like objects without a filename, try to detect gzip magic bytes.
            read_csv_kwargs = dict(sep='\t', engine='python', compression='infer')

            # If uploaded_file is a bytes buffer without a filename, peek to detect gzip
            try:
                has_name = hasattr(uploaded_file, 'name') and isinstance(getattr(uploaded_file, 'name'), str)
            except Exception:
                has_name = False

            if not has_name:
                try:
                    # Read first two bytes to check for gzip magic number
                    pos = None
                    if hasattr(uploaded_file, 'tell') and hasattr(uploaded_file, 'seek'):
                        try:
                            pos = uploaded_file.tell()
                        except Exception:
                            pos = None

                    first2 = uploaded_file.read(2)
                    # Ensure we restore pointer
                    try:
                        if pos is not None:
                            uploaded_file.seek(pos)
                        else:
                            uploaded_file.seek(0)
                    except Exception:
                        pass

                    if isinstance(first2, (bytes, bytearray)) and len(first2) >= 2 and first2[0] == 0x1f and first2[1] == 0x8b:
                        read_csv_kwargs['compression'] = 'gzip'
                except Exception:
                    # If any peek fails, continue and let pandas try to infer
                    pass

            # Choose header and dtype strategy based on file_type
            if file_type in ['Histone Marks', 'ChIP-seq', 'BED']:
                # BED-like: no header
                read_csv_kwargs['header'] = None
                df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                logger.debug("Loaded BED-like file: %s columns", df.shape[1])
                return df

            elif file_type == 'Gene Expression':
                # Gene expression: attempt a lightweight header detection by peeking first line.
                header_detected = False
                try:
                    # Preserve position
                    pos = None
                    if hasattr(uploaded_file, 'tell') and hasattr(uploaded_file, 'seek'):
                        try:
                            pos = uploaded_file.tell()
                        except Exception:
                            pos = None

                    # Peek a chunk safely (works for StringIO and BytesIO)
                    peek = uploaded_file.read(2048)
                    if isinstance(peek, (bytes, bytearray)):
                        try:
                            peek = peek.decode('utf-8')
                        except Exception:
                            peek = ''

                    # Heuristic: if first non-empty line contains non-numeric tokens, treat as header
                    first_line = ''
                    if isinstance(peek, str) and peek:
                        for line in peek.splitlines():
                            if line.strip():
                                first_line = line
                                break

                    if first_line and '\t' in first_line:
                        tokens = [t for t in first_line.strip().split('\t') if t]
                        if tokens and any(not t.replace('.', '', 1).isdigit() for t in tokens):
                            header_detected = True

                    # restore pointer
                    try:
                        if pos is not None:
                            uploaded_file.seek(pos)
                        else:
                            uploaded_file.seek(0)
                    except Exception:
                        pass
                except Exception:
                    header_detected = False

                try:
                    read_csv_kwargs['header'] = 0 if header_detected else None
                    df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                except Exception:
                    # If parsing fails, attempt fallback without header
                    try:
                        uploaded_file.seek(0)
                    except Exception:
                        pass
                    read_csv_kwargs['header'] = None
                    df = pd.read_csv(uploaded_file, **read_csv_kwargs)

                logger.debug("Loaded gene expression file: %s columns (header_detected=%s)", df.shape[1], header_detected)
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

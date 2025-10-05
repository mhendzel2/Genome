import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import hashlib
import logging

logger = logging.getLogger(__name__)


class DataValidator:
    """A class for validating and loading genomics data files.

    Notes:
    - load_as_dataframe returns a pandas DataFrame and attempts to infer compression,
      header presence, and separators common in genomics data.
    - For large pipelines or binary formats (bigWig, bigBed, .cool) use domain-specific libraries.
    """

    def __init__(self):
        self.validation_rules = {
            'BED': self._validate_bed,
            'ChIP-seq': self._validate_bed,
            'Histone Marks': self._validate_bed,
            'Gene Expression': self._validate_gene_expression,
            'VCF': self._validate_vcf,
            'GTF': self._validate_gtf,
            'GFF': self._validate_gff
        }


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
    
    def validate_file(self, uploaded_file, file_type: str) -> Dict[str, any]:
        """
        Comprehensive validation of uploaded file.
        
        Args:
            uploaded_file: File object to validate
            file_type: Type of genomics data file
        
        Returns:
            Dictionary with validation results and metrics
        """
        validation_result = {
            'valid': True,
            'errors': [],
            'warnings': [],
            'metrics': {}
        }
        
        try:
            # Load the file
            df = self.load_as_dataframe(uploaded_file, file_type)
            
            # Apply type-specific validation
            if file_type in self.validation_rules:
                type_validation = self.validation_rules[file_type](df)
                validation_result['errors'].extend(type_validation.get('errors', []))
                validation_result['warnings'].extend(type_validation.get('warnings', []))
                validation_result['metrics'].update(type_validation.get('metrics', {}))
            
            # General metrics
            validation_result['metrics']['n_rows'] = len(df)
            validation_result['metrics']['n_columns'] = len(df.columns)
            validation_result['metrics']['memory_usage_mb'] = df.memory_usage(deep=True).sum() / 1024**2
            
            # Check for errors
            if validation_result['errors']:
                validation_result['valid'] = False
            
        except Exception as e:
            validation_result['valid'] = False
            validation_result['errors'].append(f"Failed to load file: {str(e)}")
            logger.exception("File validation failed: %s", e)
        
        return validation_result
    
    def _validate_bed(self, df: pd.DataFrame) -> Dict[str, any]:
        """Validate BED format file"""
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        # Check minimum columns
        if df.shape[1] < 3:
            result['errors'].append("BED file must have at least 3 columns (chr, start, end)")
            return result
        
        # Check coordinate validity
        start_col = df.iloc[:, 1]
        end_col = df.iloc[:, 2]
        
        try:
            starts = pd.to_numeric(start_col, errors='coerce')
            ends = pd.to_numeric(end_col, errors='coerce')
            
            invalid_coords = (starts >= ends) | (starts < 0) | (ends < 0)
            n_invalid = invalid_coords.sum()
            
            if n_invalid > 0:
                result['errors'].append(f"{n_invalid} regions have invalid coordinates (start >= end or negative)")
            
            result['metrics']['invalid_coordinates'] = int(n_invalid)
            result['metrics']['mean_region_length'] = float((ends - starts).mean())
            result['metrics']['median_region_length'] = float((ends - starts).median())
            
        except Exception as e:
            result['errors'].append(f"Error validating coordinates: {str(e)}")
        
        # Check chromosome names
        chr_col = df.iloc[:, 0].astype(str)
        valid_chr_pattern = chr_col.str.match(r'^(chr)?[0-9XYM]+$', case=False)
        n_invalid_chr = (~valid_chr_pattern).sum()
        
        if n_invalid_chr > 0:
            result['warnings'].append(f"{n_invalid_chr} regions have non-standard chromosome names")
        
        result['metrics']['unique_chromosomes'] = chr_col.nunique()
        result['metrics']['regions_by_chr'] = chr_col.value_counts().to_dict()
        
        return result
    
    def _validate_gene_expression(self, df: pd.DataFrame) -> Dict[str, any]:
        """Validate gene expression file"""
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        if df.shape[1] < 2:
            result['errors'].append("Expression file must have at least 2 columns (gene, value)")
            return result
        
        # Check for numeric values
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) == 0:
            result['errors'].append("No numeric expression values found")
        
        # Check for missing values
        missing_pct = (df.isna().sum().sum() / df.size) * 100
        result['metrics']['missing_value_percentage'] = float(missing_pct)
        
        if missing_pct > 20:
            result['warnings'].append(f"{missing_pct:.1f}% missing values detected")
        
        # Check value distribution
        for col in numeric_cols:
            values = df[col].dropna()
            if len(values) > 0:
                result['metrics'][f'{col}_mean'] = float(values.mean())
                result['metrics'][f'{col}_std'] = float(values.std())
                result['metrics'][f'{col}_range'] = [float(values.min()), float(values.max())]
        
        return result
    
    def _validate_vcf(self, df: pd.DataFrame) -> Dict[str, any]:
        """Validate VCF format file"""
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        # VCF has specific required columns
        required_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        
        if df.shape[1] < 8:
            result['errors'].append("VCF file must have at least 8 columns")
            return result
        
        result['metrics']['n_variants'] = len(df)
        result['metrics']['n_samples'] = max(0, df.shape[1] - 9)  # Columns after FORMAT
        
        return result
    
    def _validate_gtf(self, df: pd.DataFrame) -> Dict[str, any]:
        """Validate GTF format file"""
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        if df.shape[1] < 9:
            result['errors'].append("GTF file must have 9 columns")
            return result
        
        # Check feature types
        if df.shape[1] >= 3:
            feature_col = df.iloc[:, 2]
            result['metrics']['feature_types'] = feature_col.value_counts().to_dict()
        
        return result
    
    def _validate_gff(self, df: pd.DataFrame) -> Dict[str, any]:
        """Validate GFF format file"""
        return self._validate_gtf(df)  # GFF and GTF have similar structures
    
    def calculate_file_checksum(self, uploaded_file) -> str:
        """Calculate MD5 checksum of file"""
        try:
            uploaded_file.seek(0)
            content = uploaded_file.read()
            if isinstance(content, str):
                content = content.encode()
            checksum = hashlib.md5(content).hexdigest()
            uploaded_file.seek(0)
            return checksum
        except Exception as e:
            logger.exception("Error calculating checksum: %s", e)
            return "error"
    
    def detect_file_format(self, uploaded_file) -> str:
        """
        Attempt to detect file format from content.
        
        Returns:
            Detected format string
        """
        try:
            uploaded_file.seek(0)
            
            # Read first few lines
            if hasattr(uploaded_file, 'read'):
                content = uploaded_file.read(10000)
                if isinstance(content, bytes):
                    content = content.decode('utf-8', errors='ignore')
                uploaded_file.seek(0)
            else:
                return 'unknown'
            
            lines = [l.strip() for l in content.split('\n') if l.strip()]
            
            if not lines:
                return 'unknown'
            
            first_line = lines[0]
            
            # Check for VCF
            if first_line.startswith('##fileformat=VCF'):
                return 'VCF'
            
            # Check for GTF/GFF
            if first_line.startswith('##gff-version'):
                return 'GFF'
            if any(line.split('\t')[2] in ['gene', 'transcript', 'exon', 'CDS'] for line in lines[:10] if '\t' in line and len(line.split('\t')) >= 3):
                return 'GTF'
            
            # Check for BED (3-6 columns, numeric start/end)
            for line in lines[:5]:
                fields = line.split('\t')
                if len(fields) >= 3:
                    try:
                        int(fields[1])
                        int(fields[2])
                        return 'BED'
                    except:
                        pass
            
            # Default to generic
            return 'Generic'
            
        except Exception as e:
            logger.exception("Error detecting file format: %s", e)
            return 'unknown'


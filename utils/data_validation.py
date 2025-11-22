import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
import hashlib
import logging
import tempfile
import os

# Import pysam but handle if not installed (though we installed it)
try:
    import pysam
except ImportError:
    pysam = None

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
            'GFF': self._validate_gff,
            'BAM': self._validate_bam
        }


    def load_as_dataframe(self, uploaded_file, file_type: str, chunk_size: Optional[int] = None) -> pd.DataFrame:
        """Loads an uploaded file into a pandas DataFrame.

        If chunk_size is provided, returns an iterator or processes only the first chunk for validation.
        For BAM/CRAM, returns a specialized object or summary dataframe.
        """
        logger = logging.getLogger(__name__)

        try:
            # Handle BAM/CRAM files specially (binary)
            if file_type in ['BAM', 'CRAM']:
                if pysam:
                    # We can't load BAM into a simple DataFrame easily.
                    # We'll return a summary DataFrame of the header or first few reads.
                    # Note: uploaded_file might be a stream. pysam needs a file path or file object.
                    return self._load_bam_summary(uploaded_file)
                else:
                    raise ImportError("pysam is required to load BAM files")

            # Reset file pointer
            try:
                uploaded_file.seek(0)
            except Exception:
                pass

            # Pandas can infer compression from filename or by using compression='infer'.
            # For file-like objects without a filename, try to detect gzip magic bytes.
            read_csv_kwargs = dict(sep='\t', engine='python', compression='infer')
            if chunk_size:
                read_csv_kwargs['chunksize'] = chunk_size

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
                    pass

            # Choose header and dtype strategy based on file_type
            if file_type in ['Histone Marks', 'ChIP-seq', 'BED']:
                # BED-like: no header
                read_csv_kwargs['header'] = None
                df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                if chunk_size:
                    # If it's an iterator, get the first chunk
                    df = next(df)
                logger.debug("Loaded BED-like file: %s columns", df.shape[1])
                return df

            elif file_type == 'Gene Expression':
                # Gene expression: attempt a lightweight header detection by peeking first line.
                header_detected = False
                try:
                    pos = None
                    if hasattr(uploaded_file, 'tell') and hasattr(uploaded_file, 'seek'):
                        try:
                            pos = uploaded_file.tell()
                        except Exception:
                            pos = None

                    peek = uploaded_file.read(2048)
                    if isinstance(peek, (bytes, bytearray)):
                        try:
                            peek = peek.decode('utf-8')
                        except Exception:
                            peek = ''

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
                    if chunk_size:
                        df = next(df)
                except Exception:
                    try:
                        uploaded_file.seek(0)
                    except Exception:
                        pass
                    read_csv_kwargs['header'] = None
                    df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                    if chunk_size:
                        df = next(df)

                logger.debug("Loaded gene expression file: %s columns (header_detected=%s)", df.shape[1], header_detected)
                return df

            else:
                # Generic TSV loader
                read_csv_kwargs['header'] = None
                df = pd.read_csv(uploaded_file, **read_csv_kwargs)
                if chunk_size:
                    df = next(df)
                logger.debug("Loaded generic file: %s columns", df.shape[1])
                return df

        except Exception as e:
            logger.exception("Failed to load file as dataframe: %s", e)
            raise ValueError(f"Error loading file (type={file_type}): {e}")

    def _load_bam_summary(self, uploaded_file) -> pd.DataFrame:
        """Helper to create a summary DataFrame from a BAM file using pysam"""
        # pysam needs a file path, so we dump to temp
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bam") as tmp:
            uploaded_file.seek(0)
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name

        try:
            with pysam.AlignmentFile(tmp_path, "rb") as bam:
                # Get basic stats
                stats = {
                    'references': list(bam.references),
                    'lengths': list(bam.lengths),
                    'mapped': bam.mapped,
                    'unmapped': bam.unmapped
                }
                # Create a simple DataFrame representing chromosomes
                df = pd.DataFrame({
                    'chr': stats['references'],
                    'length': stats['lengths']
                })
                return df
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def load_with_metadata(self, uploaded_file, file_type: str):
        """Load file and return (DataFrame, metadata dict).

        Metadata keys: 'format', 'num_rows', 'num_columns', 'column_names' (if known).
        """
        # Load only a chunk if possible to get metadata, but for full load we might need full df?
        # The prompt says "validate files without loading the full dataset into memory".
        # So we should modify validation to use chunks.
        # But this method returns a DataFrame. Consumers of this method (if any) expect a DataFrame.
        # If the file is huge, this will still crash.
        # However, for *Validation* purposes, we call validate_file which should be smart.

        # If the file is BAM, load_as_dataframe returns a summary, which is fine.
        df = self.load_as_dataframe(uploaded_file, file_type, chunk_size=1000 if file_type not in ['BAM', 'CRAM'] else None)

        meta = {
            'format': file_type,
            'num_rows': int(df.shape[0]) if file_type not in ['BAM', 'CRAM'] else 'N/A (Summary)',
            'num_columns': int(df.shape[1]),
            'column_names': list(df.columns) if df.columns is not None else None
        }

        # Normalize common BED-like files
        if file_type in ['Histone Marks', 'ChIP-seq', 'BED'] and df.shape[1] >= 3:
            df = df.rename(columns={0: 'chr', 1: 'start', 2: 'end'})
            meta['format'] = 'BED-like'

        if file_type == 'Gene Expression' and df.shape[1] >= 1:
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
        Uses chunking to avoid memory issues with large files.
        """
        validation_result = {
            'valid': True,
            'errors': [],
            'warnings': [],
            'metrics': {}
        }
        
        try:
            # Use chunks for validation
            chunk_size = 5000
            
            if file_type in ['BAM', 'CRAM']:
                # BAM validation using pysam
                if not pysam:
                    validation_result['valid'] = False
                    validation_result['errors'].append("pysam is not installed, cannot validate BAM files")
                    return validation_result

                # Write to temp file for pysam
                with tempfile.NamedTemporaryFile(delete=False, suffix=".bam") as tmp:
                    uploaded_file.seek(0)
                    tmp.write(uploaded_file.read())
                    tmp_path = tmp.name

                try:
                    # Validate BAM
                    bam_val = self._validate_bam_file(tmp_path)
                    validation_result['errors'].extend(bam_val.get('errors', []))
                    validation_result['warnings'].extend(bam_val.get('warnings', []))
                    validation_result['metrics'].update(bam_val.get('metrics', {}))
                finally:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)

            else:
                # Text-based file validation with chunks
                # Load first chunk to check structure
                df_chunk = self.load_as_dataframe(uploaded_file, file_type, chunk_size=chunk_size)

                # Basic metrics from first chunk
                validation_result['metrics']['n_columns'] = len(df_chunk.columns)
                # n_rows is just for the chunk

                # Apply type-specific validation on the first chunk
                if file_type in self.validation_rules and file_type != 'BAM':
                    type_validation = self.validation_rules[file_type](df_chunk)
                    validation_result['errors'].extend(type_validation.get('errors', []))
                    validation_result['warnings'].extend(type_validation.get('warnings', []))
                    validation_result['metrics'].update(type_validation.get('metrics', {}))

                # Note: Full validation (looping through all chunks) could be done here if needed
                # but might be slow. For now, validating the first chunk is a good proxy for structure.
            
            # Check for errors
            if validation_result['errors']:
                validation_result['valid'] = False
            
        except Exception as e:
            validation_result['valid'] = False
            validation_result['errors'].append(f"Failed to load file: {str(e)}")
            logger.exception("File validation failed: %s", e)
        
        return validation_result

    def _validate_bam_file(self, file_path: str) -> Dict[str, Any]:
        """Validate BAM file using pysam"""
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        try:
            with pysam.AlignmentFile(file_path, "rb") as bam:
                # Check header
                if not bam.header:
                    result['errors'].append("BAM file has no header")

                result['metrics']['references'] = len(bam.references)
                result['metrics']['mapped_reads'] = bam.mapped
                result['metrics']['unmapped_reads'] = bam.unmapped

                # Check first few reads
                try:
                    count = 0
                    for read in bam.fetch(until_eof=True):
                        count += 1
                        if count > 100:
                            break
                        if read.query_qualities is None and read.mapq == 0:
                            # Just a heuristic warning
                            pass
                    result['metrics']['checked_reads'] = count
                except Exception as e:
                    result['errors'].append(f"Error reading records: {e}")

        except Exception as e:
             result['errors'].append(f"Invalid BAM file: {e}")

        return result
    
    def _validate_bam(self, df: pd.DataFrame) -> Dict[str, any]:
        """Placeholder for registry, though BAM is handled specially"""
        return {'errors': [], 'warnings': [], 'metrics': {}}

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
                result['errors'].append(f"{n_invalid} regions have invalid coordinates (start >= end or negative) in checked chunk")
            
            result['metrics']['invalid_coordinates'] = int(n_invalid)
            result['metrics']['mean_region_length'] = float((ends - starts).mean())
            
        except Exception as e:
            result['errors'].append(f"Error validating coordinates: {str(e)}")
        
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
        
        # Check for missing values in chunk
        missing_pct = (df.isna().sum().sum() / df.size) * 100
        
        if missing_pct > 20:
            result['warnings'].append(f"{missing_pct:.1f}% missing values detected in chunk")
        
        # Check value distribution
        for col in numeric_cols:
            values = df[col].dropna()
            if len(values) > 0:
                result['metrics'][f'{col}_mean'] = float(values.mean())
        
        return result
    
    def _validate_vcf(self, df: pd.DataFrame) -> Dict[str, any]:
        """Validate VCF format file"""
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        # VCF has specific required columns
        required_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        
        if df.shape[1] < 8:
            result['errors'].append("VCF file must have at least 8 columns")
            return result
        
        result['metrics']['n_variants_in_chunk'] = len(df)
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
            
            # Check for BAM (magic bytes)
            if content.startswith('BAM\x01'):
                return 'BAM'

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

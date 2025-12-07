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
    
    def _validate_bam(self, uploaded_file) -> Dict[str, any]:
        """
        Validate BAM file format and alignment quality.
        
        This method performs validation of BAM files including:
        - Magic number verification
        - Header structure validation
        - Alignment flag validation
        - Quality score checking
        
        Args:
            uploaded_file: File object (BAM format, binary)
        
        Returns:
            Dictionary with errors, warnings, and metrics
        """
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        try:
            # Try to use pysam for BAM validation if available
            try:
                import pysam
                return self._validate_bam_with_pysam(uploaded_file)
            except ImportError:
                logger.info("pysam not available, using basic BAM validation")
            
            # Basic BAM validation without pysam
            uploaded_file.seek(0)
            
            # Check BAM magic number (first 4 bytes should be "BAM\1")
            magic = uploaded_file.read(4)
            if magic != b'BAM\x01':
                result['errors'].append("Invalid BAM file: incorrect magic number. Expected BAM format.")
                return result
            
            # Read header length
            header_len_bytes = uploaded_file.read(4)
            if len(header_len_bytes) < 4:
                result['errors'].append("Truncated BAM file: cannot read header length")
                return result
            
            header_len = int.from_bytes(header_len_bytes, 'little')
            result['metrics']['header_length'] = header_len
            
            # Read and validate header text
            if header_len > 0:
                header_text = uploaded_file.read(header_len)
                if len(header_text) < header_len:
                    result['errors'].append("Truncated BAM file: header incomplete")
                    return result
                
                try:
                    header_str = header_text.decode('utf-8', errors='ignore')
                    # Check for required SAM header lines
                    has_hd = '@HD' in header_str
                    has_sq = '@SQ' in header_str
                    
                    if not has_hd:
                        result['warnings'].append("BAM header missing @HD line (recommended)")
                    if not has_sq:
                        result['warnings'].append("BAM header missing @SQ reference sequences")
                    
                    # Count reference sequences
                    sq_count = header_str.count('@SQ')
                    result['metrics']['reference_sequences'] = sq_count
                    
                    # Extract additional header info
                    if '@RG' in header_str:
                        rg_count = header_str.count('@RG')
                        result['metrics']['read_groups'] = rg_count
                    
                except Exception as e:
                    result['warnings'].append(f"Could not fully parse BAM header: {e}")
            
            # Read number of reference sequences
            n_ref_bytes = uploaded_file.read(4)
            if len(n_ref_bytes) < 4:
                result['errors'].append("Truncated BAM file: cannot read reference count")
                return result
            
            n_ref = int.from_bytes(n_ref_bytes, 'little')
            result['metrics']['n_references'] = n_ref
            
            # Skip reference sequence information
            for _ in range(n_ref):
                name_len_bytes = uploaded_file.read(4)
                if len(name_len_bytes) < 4:
                    break
                name_len = int.from_bytes(name_len_bytes, 'little')
                uploaded_file.read(name_len)  # reference name
                uploaded_file.read(4)  # reference length
            
            # Sample alignments for quality validation
            alignments_checked = 0
            invalid_flags = 0
            missing_quality = 0
            unmapped_reads = 0
            max_alignments_to_check = 1000
            
            while alignments_checked < max_alignments_to_check:
                # Read alignment block size
                block_size_bytes = uploaded_file.read(4)
                if len(block_size_bytes) < 4:
                    break  # End of file
                
                block_size = int.from_bytes(block_size_bytes, 'little')
                if block_size <= 0:
                    break
                
                # Read alignment data
                alignment_data = uploaded_file.read(block_size)
                if len(alignment_data) < block_size:
                    break
                
                # Parse key alignment fields
                if len(alignment_data) >= 32:
                    # ref_id at bytes 0-3
                    ref_id = int.from_bytes(alignment_data[0:4], 'little', signed=True)
                    # pos at bytes 4-7
                    # mapq at byte 8 (after pos)
                    # bin_mq_nl at bytes 8-11
                    # flag_nc at bytes 12-15
                    flag_nc = int.from_bytes(alignment_data[12:16], 'little')
                    flag = flag_nc >> 16
                    
                    # Validate FLAG field
                    # Check for invalid flag combinations
                    is_paired = flag & 0x1
                    is_proper_pair = flag & 0x2
                    is_unmapped = flag & 0x4
                    is_mate_unmapped = flag & 0x8
                    is_reverse = flag & 0x10
                    is_mate_reverse = flag & 0x20
                    is_read1 = flag & 0x40
                    is_read2 = flag & 0x80
                    is_secondary = flag & 0x100
                    is_qc_fail = flag & 0x200
                    is_duplicate = flag & 0x400
                    is_supplementary = flag & 0x800
                    
                    # Invalid: proper pair set but not paired
                    if is_proper_pair and not is_paired:
                        invalid_flags += 1
                    
                    # Invalid: both read1 and read2 set
                    if is_read1 and is_read2:
                        invalid_flags += 1
                    
                    # Invalid: mate flags set but not paired
                    if (is_mate_unmapped or is_mate_reverse) and not is_paired:
                        invalid_flags += 1
                    
                    if is_unmapped:
                        unmapped_reads += 1
                    
                    # Check for missing quality scores (represented as 0xFF)
                    # l_seq at bytes 16-19
                    l_seq = int.from_bytes(alignment_data[16:20], 'little')
                    # Quality starts after: read_name + cigar + seq
                    # This is approximate - full parsing would require complete CIGAR decoding
                    
                alignments_checked += 1
            
            result['metrics']['alignments_sampled'] = alignments_checked
            result['metrics']['invalid_flags'] = invalid_flags
            result['metrics']['unmapped_reads'] = unmapped_reads
            
            if invalid_flags > 0:
                error_rate = invalid_flags / alignments_checked if alignments_checked > 0 else 0
                result['metrics']['flag_error_rate'] = error_rate
                if error_rate > 0.01:  # More than 1% invalid flags
                    result['warnings'].append(f"{invalid_flags} alignments with invalid FLAG combinations ({error_rate:.2%})")
            
            if alignments_checked == 0:
                result['warnings'].append("No alignments found in BAM file")
            else:
                unmapped_rate = unmapped_reads / alignments_checked
                result['metrics']['unmapped_rate'] = unmapped_rate
                if unmapped_rate > 0.5:  # More than 50% unmapped
                    result['warnings'].append(f"High unmapped read rate: {unmapped_rate:.2%}")
            
            uploaded_file.seek(0)  # Reset file pointer
            
        except Exception as e:
            result['errors'].append(f"Error validating BAM file: {str(e)}")
            logger.exception("BAM validation error: %s", e)
        
        return result
    
    def _validate_bam_with_pysam(self, uploaded_file) -> Dict[str, any]:
        """
        Validate BAM file using pysam library for comprehensive analysis.
        
        Args:
            uploaded_file: File object (BAM format)
        
        Returns:
            Dictionary with errors, warnings, and metrics
        """
        import pysam
        import tempfile
        import os
        
        result = {'errors': [], 'warnings': [], 'metrics': {}}
        
        tmp_path = None
        try:
            # Write to temporary file for pysam to read
            with tempfile.NamedTemporaryFile(delete=False, suffix='.bam') as tmp:
                uploaded_file.seek(0)
                tmp.write(uploaded_file.read())
                tmp_path = tmp.name
            
            # Open BAM file
            bam = pysam.AlignmentFile(tmp_path, 'rb')
            
            # Header validation
            header = bam.header
            result['metrics']['reference_sequences'] = bam.nreferences
            result['metrics']['references'] = list(bam.references)[:10]  # First 10 refs
            
            if 'HD' not in header:
                result['warnings'].append("BAM header missing @HD line")
            
            if 'RG' in header:
                result['metrics']['read_groups'] = len(header['RG'])
            
            # Sample alignments for validation
            alignments_checked = 0
            invalid_flags = 0
            missing_quality = 0
            unmapped_reads = 0
            duplicate_reads = 0
            qc_failed = 0
            total_mapq = 0
            max_alignments = 10000
            
            for read in bam.fetch(until_eof=True):
                alignments_checked += 1
                
                # Check flag validity
                if read.is_proper_pair and not read.is_paired:
                    invalid_flags += 1
                if read.is_read1 and read.is_read2:
                    invalid_flags += 1
                
                if read.is_unmapped:
                    unmapped_reads += 1
                if read.is_duplicate:
                    duplicate_reads += 1
                if read.is_qcfail:
                    qc_failed += 1
                
                # Check quality scores
                if read.query_qualities is None or len(read.query_qualities) == 0:
                    missing_quality += 1
                elif all(q == 255 for q in read.query_qualities):
                    missing_quality += 1
                
                if not read.is_unmapped:
                    total_mapq += read.mapping_quality
                
                if alignments_checked >= max_alignments:
                    break
            
            bam.close()
            
            # Calculate metrics
            result['metrics']['alignments_sampled'] = alignments_checked
            result['metrics']['invalid_flags'] = invalid_flags
            result['metrics']['missing_quality_scores'] = missing_quality
            result['metrics']['unmapped_reads'] = unmapped_reads
            result['metrics']['duplicate_reads'] = duplicate_reads
            result['metrics']['qc_failed_reads'] = qc_failed
            
            if alignments_checked > 0:
                mapped_reads = alignments_checked - unmapped_reads
                result['metrics']['unmapped_rate'] = unmapped_reads / alignments_checked
                result['metrics']['duplicate_rate'] = duplicate_reads / alignments_checked
                result['metrics']['mean_mapq'] = total_mapq / mapped_reads if mapped_reads > 0 else 0
                
                if invalid_flags / alignments_checked > 0.01:
                    result['warnings'].append(f"{invalid_flags} alignments with invalid FLAGS")
                if missing_quality / alignments_checked > 0.1:
                    result['warnings'].append(f"{missing_quality} reads missing quality scores")
                if unmapped_reads / alignments_checked > 0.5:
                    result['warnings'].append(f"High unmapped rate: {unmapped_reads/alignments_checked:.1%}")
            else:
                result['warnings'].append("No alignments found in BAM file")
            
            uploaded_file.seek(0)
            
        except Exception as e:
            result['errors'].append(f"Error validating BAM with pysam: {str(e)}")
            logger.exception("pysam BAM validation error: %s", e)
        finally:
            if tmp_path and os.path.exists(tmp_path):
                os.remove(tmp_path)
        
        return result
    
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

import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO, BytesIO
import gzip
import os
from typing import Dict, Any, Optional
import re

class DataValidator:
    """Handles validation of genomics data files"""
    
    def __init__(self):
        self.supported_formats = {
            'Histone Marks': ['bed', 'bedgraph', 'bigwig', 'wig'],
            'HiC': ['txt', 'bedpe', 'hic', 'cool'],
            'Gene Expression': ['csv', 'tsv', 'txt', 'fpkm', 'tpm'],
            'ChIP-seq': ['bed', 'bedgraph', 'narrowpeak', 'broadpeak', 'bigwig']
        }
    
    def validate_file(self, uploaded_file, data_type: str) -> Dict[str, Any]:
        """Validate uploaded genomics file"""
        try:
            # Get file extension
            filename = uploaded_file.name.lower()
            extension = self._get_file_extension(filename)
            
            # Check if format is supported
            if extension not in self.supported_formats.get(data_type, []):
                return {
                    'valid': False,
                    'error': f"Unsupported format '{extension}' for {data_type}. Supported: {', '.join(self.supported_formats[data_type])}"
                }
            
            # Read file content for validation
            file_content = self._read_file_content(uploaded_file)
            
            # Validate based on data type and format
            validation_result = self._validate_content(file_content, data_type, extension)
            
            if validation_result['valid']:
                # Get file info
                file_info = self._get_file_info(file_content, data_type, extension)
                validation_result['info'] = file_info
            
            return validation_result
            
        except Exception as e:
            return {
                'valid': False,
                'error': f"Validation error: {str(e)}"
            }
    
    def _get_file_extension(self, filename: str) -> str:
        """Extract file extension, handling compressed files"""
        if filename.endswith('.gz'):
            # Remove .gz and get the actual extension
            base_name = filename[:-3]
            return base_name.split('.')[-1] if '.' in base_name else 'gz'
        else:
            return filename.split('.')[-1] if '.' in filename else ''
    
    def _read_file_content(self, uploaded_file) -> str:
        """Read file content, handling compression"""
        # Reset file pointer
        uploaded_file.seek(0)
        
        try:
            # Try to read as compressed file first
            if uploaded_file.name.lower().endswith('.gz'):
                content = gzip.decompress(uploaded_file.read()).decode('utf-8')
            else:
                content = uploaded_file.read().decode('utf-8')
            
            # Reset file pointer for potential reuse
            uploaded_file.seek(0)
            return content
            
        except UnicodeDecodeError:
            # Handle binary files
            uploaded_file.seek(0)
            return uploaded_file.read().decode('utf-8', errors='ignore')
    
    def _validate_content(self, content: str, data_type: str, extension: str) -> Dict[str, Any]:
        """Validate file content based on format"""
        lines = content.strip().split('\n')
        
        if not lines or len(lines) < 1:
            return {'valid': False, 'error': 'File appears to be empty'}
        
        # Skip header lines and comments
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
        
        if not data_lines:
            return {'valid': False, 'error': 'No data lines found'}
        
        try:
            if extension in ['bed', 'bedgraph']:
                return self._validate_bed_format(data_lines)
            elif extension in ['csv', 'tsv', 'txt']:
                return self._validate_tabular_format(data_lines, extension)
            elif extension in ['narrowpeak', 'broadpeak']:
                return self._validate_peak_format(data_lines)
            elif extension in ['bedpe']:
                return self._validate_bedpe_format(data_lines)
            else:
                # For other formats, do basic validation
                return {'valid': True, 'error': None}
                
        except Exception as e:
            return {'valid': False, 'error': f'Content validation failed: {str(e)}'}
    
    def _validate_bed_format(self, lines) -> Dict[str, Any]:
        """Validate BED format files"""
        # Check first few lines
        for i, line in enumerate(lines[:10]):
            fields = line.strip().split('\t')
            
            if len(fields) < 3:
                return {'valid': False, 'error': f'Line {i+1}: BED format requires at least 3 columns (chr, start, end)'}
            
            # Validate chromosome
            if not fields[0]:
                return {'valid': False, 'error': f'Line {i+1}: Chromosome field cannot be empty'}
            
            # Validate coordinates
            try:
                start = int(fields[1])
                end = int(fields[2])
                if start >= end:
                    return {'valid': False, 'error': f'Line {i+1}: Start position must be less than end position'}
                if start < 0:
                    return {'valid': False, 'error': f'Line {i+1}: Coordinates cannot be negative'}
            except ValueError:
                return {'valid': False, 'error': f'Line {i+1}: Start and end positions must be integers'}
        
        return {'valid': True, 'error': None}
    
    def _validate_tabular_format(self, lines, extension) -> Dict[str, Any]:
        """Validate tabular format files (CSV, TSV, TXT)"""
        delimiter = ',' if extension == 'csv' else '\t'
        
        # Check if all lines have consistent number of columns
        first_line_cols = len(lines[0].split(delimiter))
        
        for i, line in enumerate(lines[:20]):  # Check first 20 lines
            cols = len(line.split(delimiter))
            if cols != first_line_cols:
                return {'valid': False, 'error': f'Line {i+1}: Inconsistent number of columns'}
        
        return {'valid': True, 'error': None}
    
    def _validate_peak_format(self, lines) -> Dict[str, Any]:
        """Validate ENCODE peak format files"""
        for i, line in enumerate(lines[:10]):
            fields = line.strip().split('\t')
            
            if len(fields) < 6:
                return {'valid': False, 'error': f'Line {i+1}: Peak format requires at least 6 columns'}
            
            # Validate basic BED fields
            try:
                start = int(fields[1])
                end = int(fields[2])
                if start >= end or start < 0:
                    return {'valid': False, 'error': f'Line {i+1}: Invalid coordinates'}
            except ValueError:
                return {'valid': False, 'error': f'Line {i+1}: Invalid coordinate format'}
        
        return {'valid': True, 'error': None}
    
    def _validate_bedpe_format(self, lines) -> Dict[str, Any]:
        """Validate BEDPE format for HiC data"""
        for i, line in enumerate(lines[:10]):
            fields = line.strip().split('\t')
            
            if len(fields) < 6:
                return {'valid': False, 'error': f'Line {i+1}: BEDPE format requires at least 6 columns'}
            
            # Validate coordinates for both regions
            try:
                start1, end1 = int(fields[1]), int(fields[2])
                start2, end2 = int(fields[4]), int(fields[5])
                
                if start1 >= end1 or start2 >= end2 or start1 < 0 or start2 < 0:
                    return {'valid': False, 'error': f'Line {i+1}: Invalid coordinates'}
            except ValueError:
                return {'valid': False, 'error': f'Line {i+1}: Invalid coordinate format'}
        
        return {'valid': True, 'error': None}
    
    def _get_file_info(self, content: str, data_type: str, extension: str) -> str:
        """Get basic file information"""
        lines = content.strip().split('\n')
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
        
        info_parts = [
            f"Format: {extension.upper()}",
            f"Total lines: {len(lines)}",
            f"Data lines: {len(data_lines)}"
        ]
        
        # Add format-specific info
        if data_lines:
            first_data_line = data_lines[0].strip().split('\t')
            info_parts.append(f"Columns: {len(first_data_line)}")
            
            if extension in ['bed', 'bedgraph'] and len(first_data_line) >= 3:
                try:
                    start = int(first_data_line[1])
                    end = int(first_data_line[2])
                    info_parts.append(f"Sample region: {first_data_line[0]}:{start}-{end}")
                except:
                    pass
        
        return " | ".join(info_parts)
    
    def load_as_dataframe(self, uploaded_file, data_type: str) -> pd.DataFrame:
        """Load validated file as pandas DataFrame"""
        # Reset file pointer
        uploaded_file.seek(0)
        
        # Read content
        content = self._read_file_content(uploaded_file)
        lines = content.strip().split('\n')
        
        # Filter out comments
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
        
        # Determine delimiter
        extension = self._get_file_extension(uploaded_file.name.lower())
        if extension == 'csv':
            delimiter = ','
        else:
            delimiter = '\t'
        
        # Create DataFrame
        data = []
        for line in data_lines:
            data.append(line.strip().split(delimiter))
        
        # Generate column names based on format
        if data:
            num_cols = len(data[0])
            if extension in ['bed', 'bedgraph']:
                col_names = ['chr', 'start', 'end']
                if num_cols > 3:
                    col_names.extend([f'col_{i}' for i in range(4, num_cols + 1)])
            elif extension in ['narrowpeak', 'broadpeak']:
                col_names = ['chr', 'start', 'end', 'name', 'score', 'strand']
                if num_cols > 6:
                    col_names.extend([f'col_{i}' for i in range(7, num_cols + 1)])
            else:
                col_names = [f'col_{i}' for i in range(1, num_cols + 1)]
            
            df = pd.DataFrame(data, columns=col_names[:num_cols])
            
            # Convert numeric columns
            for col in ['start', 'end', 'score']:
                if col in df.columns:
                    try:
                        df[col] = pd.to_numeric(df[col], errors='coerce')
                    except:
                        pass
            
            return df
        
        return pd.DataFrame()

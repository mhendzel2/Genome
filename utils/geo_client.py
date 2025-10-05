"""
Client for accessing Gene Expression Omnibus (GEO) database.
Provides methods to search, download, and parse GEO datasets.
"""

import requests
from typing import Dict, List, Any, Optional
import os
import gzip
import logging
import time
from xml.etree import ElementTree as ET


logger = logging.getLogger(__name__)


class GEOClient:
    """
    A client to interact with NCBI's Gene Expression Omnibus (GEO).
    """

    def __init__(self):
        self.base_url = "https://www.ncbi.nlm.nih.gov/geo"
        self.eutils_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.ftp_base = "https://ftp.ncbi.nlm.nih.gov/geo"
        
    def search_datasets(self, search_term: str, 
                       dataset_type: str = "gse",
                       max_results: int = 20) -> List[Dict[str, Any]]:
        """
        Search for GEO datasets using NCBI E-utilities.
        
        Args:
            search_term: Search query string
            dataset_type: 'gse' for series, 'gds' for datasets, 'gpl' for platforms
            max_results: Maximum number of results to return
        
        Returns:
            List of dataset dictionaries with metadata
        """
        try:
            # First, search for IDs
            search_url = f"{self.eutils_base}/esearch.fcgi"
            search_params = {
                'db': 'gds',
                'term': f"{search_term} AND {dataset_type.upper()}[Entry Type]",
                'retmax': max_results,
                'retmode': 'json'
            }
            
            response = requests.get(search_url, params=search_params, timeout=30)
            response.raise_for_status()
            search_results = response.json()
            
            id_list = search_results.get('esearchresult', {}).get('idlist', [])
            
            if not id_list:
                return []
            
            # Fetch summaries for the IDs
            summary_url = f"{self.eutils_base}/esummary.fcgi"
            summary_params = {
                'db': 'gds',
                'id': ','.join(id_list),
                'retmode': 'json'
            }
            
            time.sleep(0.34)  # NCBI rate limit: max 3 requests per second
            summary_response = requests.get(summary_url, params=summary_params, timeout=30)
            summary_response.raise_for_status()
            summaries = summary_response.json()
            
            # Format results
            formatted_results = []
            result_dict = summaries.get('result', {})
            
            for geo_id in id_list:
                if geo_id in result_dict:
                    record = result_dict[geo_id]
                    
                    # Extract sample count
                    n_samples = record.get('n_samples', 0)
                    
                    formatted_results.append({
                        'id': record.get('accession', geo_id),
                        'accession': record.get('accession', geo_id),
                        'title': record.get('title', 'No title'),
                        'description': record.get('summary', 'No description'),
                        'organism': record.get('taxon', 'Unknown'),
                        'platform': record.get('gpl', 'Unknown'),
                        'n_samples': n_samples,
                        'pubmed_id': record.get('pubmedids', [''])[0] if record.get('pubmedids') else '',
                        'submission_date': record.get('pdat', 'Unknown'),
                        'entry_type': record.get('entrytype', dataset_type.upper()),
                        'database': 'GEO'
                    })
            
            return formatted_results
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error querying GEO: %s", e)
            return []
        except Exception as e:
            logger.exception("Unexpected error in GEO search: %s", e)
            return []
    
    def get_dataset_details(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a specific GEO dataset.
        
        Args:
            accession: GEO accession (e.g., GSE12345)
        
        Returns:
            Dictionary with detailed metadata
        """
        try:
            # Use NCBI E-utilities to fetch detailed record
            fetch_url = f"{self.eutils_base}/efetch.fcgi"
            fetch_params = {
                'db': 'gds',
                'id': accession,
                'retmode': 'xml'
            }
            
            response = requests.get(fetch_url, params=fetch_params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            root = ET.fromstring(response.content)
            
            # Extract relevant information
            details = {
                'accession': accession,
                'title': '',
                'summary': '',
                'organism': '',
                'platform': '',
                'samples': [],
                'supplementary_files': []
            }
            
            # Parse based on entry type (this is simplified - full parsing would be more complex)
            for item in root.findall('.//Item'):
                name = item.get('Name', '')
                if name == 'title':
                    details['title'] = item.text or ''
                elif name == 'summary':
                    details['summary'] = item.text or ''
            
            return details
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error fetching GEO dataset details: %s", e)
            return None
        except Exception as e:
            logger.exception("Unexpected error in get_dataset_details: %s", e)
            return None
    
    def download_series_matrix(self, accession: str, 
                              download_dir: str = "downloads/geo") -> str:
        """
        Download series matrix file for a GEO Series.
        
        Args:
            accession: GEO Series accession (e.g., GSE12345)
            download_dir: Directory to save the file
        
        Returns:
            Path to downloaded file
        """
        os.makedirs(download_dir, exist_ok=True)
        
        try:
            # Construct FTP URL for series matrix file
            # Format: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE12345/matrix/
            series_stub = accession[:-3] + 'nnn'  # GSE12345 -> GSE12nnn
            matrix_url = f"{self.ftp_base}/series/{series_stub}/{accession}/matrix/"
            
            # Try to get the series matrix file
            matrix_filename = f"{accession}_series_matrix.txt.gz"
            file_url = matrix_url + matrix_filename
            
            logger.info(f"Downloading from: {file_url}")
            
            response = requests.get(file_url, stream=True, timeout=300)
            response.raise_for_status()
            
            local_path = os.path.join(download_dir, matrix_filename)
            
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            logger.info(f"Successfully downloaded to: {local_path}")
            return local_path
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error downloading series matrix: %s", e)
            raise Exception(f"Failed to download {accession}: {e}")
    
    def download_supplementary_files(self, accession: str,
                                    download_dir: str = "downloads/geo") -> List[str]:
        """
        Download supplementary files for a GEO dataset.
        
        Args:
            accession: GEO accession
            download_dir: Directory to save files
        
        Returns:
            List of paths to downloaded files
        """
        os.makedirs(download_dir, exist_ok=True)
        downloaded_files = []
        
        try:
            # Construct supplementary files URL
            series_stub = accession[:-3] + 'nnn'
            suppl_url = f"{self.ftp_base}/series/{series_stub}/{accession}/suppl/"
            
            logger.info(f"Accessing supplementary files at: {suppl_url}")
            
            # Get directory listing (this is simplified - full implementation would parse HTML)
            response = requests.get(suppl_url, timeout=30)
            response.raise_for_status()
            
            # This is a simplified approach - in practice, you'd parse the directory listing
            # For now, we'll try common supplementary file patterns
            common_patterns = [
                f"{accession}_RAW.tar",
                f"filelist.txt"
            ]
            
            for filename in common_patterns:
                try:
                    file_url = suppl_url + filename
                    file_response = requests.get(file_url, stream=True, timeout=300)
                    file_response.raise_for_status()
                    
                    local_path = os.path.join(download_dir, filename)
                    with open(local_path, 'wb') as f:
                        for chunk in file_response.iter_content(chunk_size=8192):
                            f.write(chunk)
                    
                    downloaded_files.append(local_path)
                    logger.info(f"Downloaded: {local_path}")
                except:
                    continue
            
            return downloaded_files
            
        except Exception as e:
            logger.exception("Error downloading supplementary files: %s", e)
            return downloaded_files
    
    def parse_series_matrix(self, file_path: str) -> Dict[str, Any]:
        """
        Parse a GEO series matrix file.
        
        Args:
            file_path: Path to series matrix file (can be gzipped)
        
        Returns:
            Dictionary with parsed data including expression matrix and metadata
        """
        try:
            # Open file (handle gzip if needed)
            if file_path.endswith('.gz'):
                with gzip.open(file_path, 'rt', encoding='utf-8') as f:
                    lines = f.readlines()
            else:
                with open(file_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
            
            metadata = {}
            data_start_idx = None
            
            # Parse metadata and find data start
            for i, line in enumerate(lines):
                line = line.strip()
                
                if line.startswith('!Series_'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0].replace('!Series_', '').replace('"', '')
                        value = '\t'.join(parts[1:]).replace('"', '')
                        metadata[key] = value
                
                elif line.startswith('!dataset_table_begin'):
                    data_start_idx = i + 1
                    break
            
            # Extract expression data if present
            expression_data = None
            if data_start_idx:
                # Find data end
                data_end_idx = None
                for i in range(data_start_idx, len(lines)):
                    if lines[i].strip().startswith('!dataset_table_end'):
                        data_end_idx = i
                        break
                
                if data_end_idx:
                    data_lines = [line.strip() for line in lines[data_start_idx:data_end_idx]]
                    # Parse as tab-delimited data
                    # First line is typically headers
                    if data_lines:
                        headers = data_lines[0].split('\t')
                        rows = [line.split('\t') for line in data_lines[1:]]
                        
                        expression_data = {
                            'headers': headers,
                            'n_rows': len(rows),
                            'n_columns': len(headers)
                        }
            
            return {
                'metadata': metadata,
                'expression_data': expression_data,
                'file_path': file_path
            }
            
        except Exception as e:
            logger.exception("Error parsing series matrix: %s", e)
            raise ValueError(f"Failed to parse series matrix: {e}")
    
    def get_sample_data(self, gse_accession: str, 
                       gsm_accession: str) -> Optional[Dict[str, Any]]:
        """
        Get data for a specific sample (GSM) within a series (GSE).
        
        Args:
            gse_accession: Series accession (e.g., GSE12345)
            gsm_accession: Sample accession (e.g., GSM123456)
        
        Returns:
            Dictionary with sample metadata and data
        """
        try:
            # Construct sample data URL
            sample_url = f"{self.base_url}/query/acc.cgi"
            params = {
                'acc': gsm_accession,
                'targ': 'self',
                'form': 'text',
                'view': 'quick'
            }
            
            response = requests.get(sample_url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse response (simplified - full parser would be more complex)
            content = response.text
            
            sample_data = {
                'accession': gsm_accession,
                'series': gse_accession,
                'raw_content': content[:1000]  # First 1000 chars for preview
            }
            
            return sample_data
            
        except Exception as e:
            logger.exception("Error getting sample data: %s", e)
            return None

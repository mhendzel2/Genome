"""
Client for accessing The Cancer Genome Atlas (TCGA) data via Genomic Data Commons (GDC).
Provides methods to search, filter, and download cancer genomics data.
"""

import requests
import json
import os
import logging
from typing import Dict, List, Any, Optional


logger = logging.getLogger(__name__)


class TCGAClient:
    """
    A client to interact with the TCGA/GDC data portal.
    """

    def __init__(self):
        self.api_base = "https://api.gdc.cancer.gov"
        self.data_endpoint = f"{self.api_base}/data"
        self.files_endpoint = f"{self.api_base}/files"
        self.cases_endpoint = f"{self.api_base}/cases"
        self.projects_endpoint = f"{self.api_base}/projects"
        
    def list_projects(self, program: str = "TCGA") -> List[Dict[str, Any]]:
        """
        List available TCGA/GDC projects.
        
        Args:
            program: Program name ('TCGA', 'TARGET', 'FM', etc.)
        
        Returns:
            List of project dictionaries
        """
        try:
            params = {
                'size': 100,
                'sort': 'project.project_id:asc'
            }
            
            if program:
                params['filters'] = json.dumps({
                    'op': '=',
                    'content': {
                        'field': 'program.name',
                        'value': program
                    }
                })
            
            response = requests.get(self.projects_endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            projects = []
            for hit in data.get('data', {}).get('hits', []):
                projects.append({
                    'project_id': hit.get('project_id', ''),
                    'name': hit.get('name', ''),
                    'primary_site': hit.get('primary_site', ['Unknown'])[0] if hit.get('primary_site') else 'Unknown',
                    'disease_type': hit.get('disease_type', ['Unknown'])[0] if hit.get('disease_type') else 'Unknown',
                    'program': hit.get('program', {}).get('name', ''),
                    'project_name': hit.get('project', {}).get('name', ''),
                    'case_count': hit.get('summary', {}).get('case_count', 0),
                    'file_count': hit.get('summary', {}).get('file_count', 0)
                })
            
            return projects
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error listing TCGA projects: %s", e)
            return []
        except Exception as e:
            logger.exception("Unexpected error in list_projects: %s", e)
            return []
    
    def search_files(self, 
                    project_id: Optional[str] = None,
                    data_category: Optional[str] = None,
                    data_type: Optional[str] = None,
                    experimental_strategy: Optional[str] = None,
                    workflow_type: Optional[str] = None,
                    access: str = "open",
                    size: int = 20) -> List[Dict[str, Any]]:
        """
        Search for files in GDC.
        
        Args:
            project_id: TCGA project ID (e.g., 'TCGA-BRCA')
            data_category: E.g., 'Transcriptome Profiling', 'Simple Nucleotide Variation'
            data_type: E.g., 'Gene Expression Quantification', 'Masked Somatic Mutation'
            experimental_strategy: E.g., 'RNA-Seq', 'WXS', 'WGS'
            workflow_type: E.g., 'HTSeq - Counts', 'STAR - Counts'
            access: 'open' or 'controlled'
            size: Maximum number of results
        
        Returns:
            List of file dictionaries
        """
        try:
            # Build filter
            filters = {'op': 'and', 'content': []}
            
            if project_id:
                filters['content'].append({
                    'op': '=',
                    'content': {
                        'field': 'cases.project.project_id',
                        'value': project_id
                    }
                })
            
            if data_category:
                filters['content'].append({
                    'op': '=',
                    'content': {
                        'field': 'data_category',
                        'value': data_category
                    }
                })
            
            if data_type:
                filters['content'].append({
                    'op': '=',
                    'content': {
                        'field': 'data_type',
                        'value': data_type
                    }
                })
            
            if experimental_strategy:
                filters['content'].append({
                    'op': '=',
                    'content': {
                        'field': 'experimental_strategy',
                        'value': experimental_strategy
                    }
                })
            
            if workflow_type:
                filters['content'].append({
                    'op': '=',
                    'content': {
                        'field': 'analysis.workflow_type',
                        'value': workflow_type
                    }
                })
            
            filters['content'].append({
                'op': '=',
                'content': {
                    'field': 'access',
                    'value': access
                }
            })
            
            params = {
                'filters': json.dumps(filters),
                'fields': 'file_id,file_name,data_category,data_type,file_size,cases.project.project_id,experimental_strategy,analysis.workflow_type',
                'size': size,
                'format': 'JSON'
            }
            
            response = requests.get(self.files_endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            files = []
            for hit in data.get('data', {}).get('hits', []):
                files.append({
                    'file_id': hit.get('file_id', ''),
                    'file_name': hit.get('file_name', ''),
                    'data_category': hit.get('data_category', ''),
                    'data_type': hit.get('data_type', ''),
                    'file_size': hit.get('file_size', 0),
                    'file_size_mb': round(hit.get('file_size', 0) / 1048576, 2),
                    'experimental_strategy': hit.get('experimental_strategy', ''),
                    'workflow_type': hit.get('analysis', {}).get('workflow_type', ''),
                    'project_id': hit.get('cases', [{}])[0].get('project', {}).get('project_id', '') if hit.get('cases') else '',
                    'access': hit.get('access', ''),
                    'database': 'TCGA/GDC'
                })
            
            return files
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error searching TCGA files: %s", e)
            return []
        except Exception as e:
            logger.exception("Unexpected error in search_files: %s", e)
            return []
    
    def get_file_metadata(self, file_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed metadata for a specific file.
        
        Args:
            file_id: GDC file UUID
        
        Returns:
            Dictionary with file metadata
        """
        try:
            params = {
                'expand': 'cases,cases.project,cases.samples,analysis',
                'format': 'JSON'
            }
            
            url = f"{self.files_endpoint}/{file_id}"
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            file_data = data.get('data', {})
            
            return {
                'file_id': file_data.get('file_id', ''),
                'file_name': file_data.get('file_name', ''),
                'data_category': file_data.get('data_category', ''),
                'data_type': file_data.get('data_type', ''),
                'file_size': file_data.get('file_size', 0),
                'md5sum': file_data.get('md5sum', ''),
                'access': file_data.get('access', ''),
                'created_datetime': file_data.get('created_datetime', ''),
                'updated_datetime': file_data.get('updated_datetime', ''),
                'cases': file_data.get('cases', []),
                'analysis': file_data.get('analysis', {})
            }
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error getting file metadata: %s", e)
            return None
        except Exception as e:
            logger.exception("Unexpected error in get_file_metadata: %s", e)
            return None
    
    def download_file(self, file_id: str, 
                     download_dir: str = "downloads/tcga",
                     file_name: Optional[str] = None) -> Optional[str]:
        """
        Download a file from GDC.
        
        Args:
            file_id: GDC file UUID
            download_dir: Directory to save the file
            file_name: Optional custom filename
        
        Returns:
            Path to downloaded file, or None if failed
        """
        os.makedirs(download_dir, exist_ok=True)
        
        try:
            # Get file metadata to get the filename if not provided
            if not file_name:
                metadata = self.get_file_metadata(file_id)
                if metadata:
                    file_name = metadata.get('file_name', file_id)
                else:
                    file_name = file_id
            
            url = f"{self.data_endpoint}/{file_id}"
            
            logger.info(f"Downloading file {file_id} from GDC...")
            
            response = requests.get(url, stream=True, timeout=300)
            response.raise_for_status()
            
            local_path = os.path.join(download_dir, file_name)
            
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            logger.info(f"Successfully downloaded to: {local_path}")
            return local_path
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error downloading file: %s", e)
            return None
        except Exception as e:
            logger.exception("Unexpected error in download_file: %s", e)
            return None
    
    def download_manifest(self, file_ids: List[str], 
                         manifest_path: str = "downloads/tcga/manifest.txt") -> str:
        """
        Create a GDC download manifest file.
        
        Args:
            file_ids: List of GDC file UUIDs
            manifest_path: Path to save the manifest
        
        Returns:
            Path to manifest file
        """
        try:
            os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
            
            # Get metadata for all files
            manifest_lines = ['id\tfilename\tmd5\tsize\tstate\n']
            
            for file_id in file_ids:
                metadata = self.get_file_metadata(file_id)
                if metadata:
                    line = f"{metadata['file_id']}\t{metadata['file_name']}\t{metadata['md5sum']}\t{metadata['file_size']}\tlive\n"
                    manifest_lines.append(line)
            
            with open(manifest_path, 'w') as f:
                f.writelines(manifest_lines)
            
            logger.info(f"Manifest created: {manifest_path}")
            return manifest_path
            
        except Exception as e:
            logger.exception("Error creating manifest: %s", e)
            raise
    
    def get_clinical_data(self, project_id: str) -> List[Dict[str, Any]]:
        """
        Get clinical data for a TCGA project.
        
        Args:
            project_id: TCGA project ID (e.g., 'TCGA-BRCA')
        
        Returns:
            List of case dictionaries with clinical information
        """
        try:
            filters = {
                'op': '=',
                'content': {
                    'field': 'project.project_id',
                    'value': project_id
                }
            }
            
            params = {
                'filters': json.dumps(filters),
                'fields': 'case_id,submitter_id,demographic,diagnoses,exposures',
                'size': 1000,
                'format': 'JSON',
                'expand': 'demographic,diagnoses,exposures'
            }
            
            response = requests.get(self.cases_endpoint, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            cases = []
            for hit in data.get('data', {}).get('hits', []):
                case_data = {
                    'case_id': hit.get('case_id', ''),
                    'submitter_id': hit.get('submitter_id', ''),
                    'demographic': hit.get('demographic', {}),
                    'diagnoses': hit.get('diagnoses', []),
                    'exposures': hit.get('exposures', [])
                }
                cases.append(case_data)
            
            return cases
            
        except requests.exceptions.RequestException as e:
            logger.exception("Error getting clinical data: %s", e)
            return []
        except Exception as e:
            logger.exception("Unexpected error in get_clinical_data: %s", e)
            return []
    
    def get_mutation_data(self, gene_symbol: str, 
                         project_id: Optional[str] = None,
                         size: int = 100) -> List[Dict[str, Any]]:
        """
        Get mutation data for a specific gene.
        
        Args:
            gene_symbol: Gene symbol (e.g., 'TP53')
            project_id: Optional TCGA project to filter by
            size: Maximum number of results
        
        Returns:
            List of mutation dictionaries
        """
        try:
            # Note: This uses the mutations endpoint which requires more complex filtering
            # This is a simplified implementation
            
            filters = {'op': 'and', 'content': [
                {
                    'op': '=',
                    'content': {
                        'field': 'ssm.consequence.transcript.gene.symbol',
                        'value': gene_symbol
                    }
                }
            ]}
            
            if project_id:
                filters['content'].append({
                    'op': '=',
                    'content': {
                        'field': 'cases.project.project_id',
                        'value': project_id
                    }
                })
            
            # This is a placeholder - full implementation would use the proper mutations endpoint
            logger.info(f"Mutation query for {gene_symbol} in project {project_id}")
            
            return []
            
        except Exception as e:
            logger.exception("Error getting mutation data: %s", e)
            return []

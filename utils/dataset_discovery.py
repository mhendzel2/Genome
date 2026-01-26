import requests
from typing import Dict, List, Any, Optional
import os
import logging

class DatasetDiscovery:
    """A class to discover datasets from public genomics databases."""

    def __init__(self):
        self.base_url = "https://www.encodeproject.org"
        self.ols_base_url = "https://www.ebi.ac.uk/ols/api"

    def semantic_search(self, query: str) -> List[str]:
        """
        Expand user query using EBI OLS (Ontology Lookup Service).
        Example: "brain cancer" -> ["glioblastoma", "glioma", "astrocytoma"]
        """
        logger = logging.getLogger(__name__)
        expanded_terms = [query]

        try:
            # Search EBI OLS
            params = {
                'q': query,
                'ontology': 'efo,mondo,ncit', # Disease ontologies
                'rows': 5
            }
            response = requests.get(f"{self.ols_base_url}/search", params=params, headers={'accept': 'application/json'})

            if response.status_code == 200:
                data = response.json()
                if 'response' in data and 'docs' in data['response']:
                    for doc in data['response']['docs']:
                        if 'label' in doc:
                            expanded_terms.append(doc['label'])
                        # Also check synonyms
                        if 'synonym' in doc and doc['synonym']:
                            expanded_terms.extend(doc['synonym'][:3]) # Take top 3

            # De-duplicate
            expanded_terms = list(set(expanded_terms))
            logger.info(f"Semantic search expanded '{query}' to: {expanded_terms}")

        except requests.exceptions.RequestException as e:
            logger.warning(f"Semantic search failed: {e}")

        return expanded_terms

    def search_datasets(self,
                       search_term: str = "",
                       data_types: Optional[List[str]] = None,
                       tissues: Optional[List[str]] = None,
                       organism: str = "Homo sapiens",
                       use_semantic_search: bool = True) -> List[Dict[str, Any]]:
        """Search for datasets on the ENCODE platform."""

        # Expand search term if requested
        search_terms = [search_term]
        if use_semantic_search and search_term:
            expanded = self.semantic_search(search_term)
            if expanded:
                search_terms = expanded

        all_results = []
        seen_accessions = set()

        for term in search_terms:
            search_params = {
                "type": "Experiment",
                "status": "released",
                "limit": "25", # Limit per term to avoid overload
                "format": "json"
            }

            if term:
                search_params["searchTerm"] = term

            if organism:
                search_params["replicates.library.biosample.donor.organism.scientific_name"] = organism

            if data_types:
                search_params["assay_title"] = data_types

            if tissues:
                search_params["biosample_ontology.term_name"] = tissues

            import logging
            logger = logging.getLogger(__name__)
            try:
                response = requests.get(f"{self.base_url}/search/", params=search_params, headers={'accept': 'application/json'})
                response.raise_for_status()
                search_results = response.json()

                formatted = self._format_results(search_results.get("@graph", []))

                for res in formatted:
                    if res['accession'] not in seen_accessions:
                        all_results.append(res)
                        seen_accessions.add(res['accession'])

            except requests.exceptions.RequestException as e:
                logger.exception("Error querying ENCODE API: %s", e)
                continue

        return all_results

    def _format_results(self, results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Format ENCODE results into the application's data structure."""
        formatted_results = []
        for result in results:
            biosample = result.get('biosample_summary', '')
            
            formatted_results.append({
                'id': result.get('accession'),
                'title': result.get('description', 'No title available'),
                'description': result.get('summary', 'No description available'),
                'data_type': result.get('assay_title', 'N/A'),
                'tissue': biosample,
                'organism': result.get('replicates', [{}])[0].get('library', {}).get('biosample', {}).get('organism', {}).get('scientific_name', 'N/A'),
                'database': 'ENCODE',
                'accession': result.get('accession'),
            })
        return formatted_results

    def get_dataset_details(self, dataset_id: str) -> Optional[Dict[str, Any]]:
        """Get detailed information about a specific dataset from ENCODE."""
        import logging
        logger = logging.getLogger(__name__)
        try:
            response = requests.get(f"{self.base_url}/experiments/{dataset_id}/?frame=embedded", headers={'accept': 'application/json'})
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.exception("Error fetching dataset details from ENCODE: %s", e)
            return None

    def get_sra_download_links(self, bioproject_id: str) -> List[str]:
        """
        Fetch fastq-dump commands or download links for a given BioProject ID via ENA or SRA.
        Simplified implementation using ENA API (often easier than NCBI E-utils for direct links).
        """
        logger = logging.getLogger(__name__)
        links = []

        ena_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
        params = {
            'accession': bioproject_id,
            'result': 'read_run',
            'fields': 'run_accession,fastq_ftp',
            'format': 'json'
        }

        try:
            response = requests.get(ena_url, params=params)
            if response.status_code == 200:
                data = response.json()
                for run in data:
                    run_acc = run.get('run_accession')
                    ftps = run.get('fastq_ftp', '').split(';')
                    if ftps:
                        # Construct fastq-dump command as alternative
                        cmd = f"fastq-dump --split-files {run_acc}"
                        links.append({'run': run_acc, 'ftp': ftps, 'command': cmd})
            else:
                logger.warning(f"ENA API returned {response.status_code}")

        except Exception as e:
             logger.warning(f"Failed to fetch SRA links: {e}")

        return links

    def download_dataset(self, dataset_id: str, file_types: Optional[List[str]] = None,
                        max_files: int = 10) -> str:
        """Download files for a given dataset from ENCODE.
        
        Args:
            dataset_id: ENCODE experiment accession (e.g., ENCSR675FLJ)
            file_types: Optional list of file types to download (e.g., ['bigWig', 'bed'])
                       If None, downloads processed files (bigWig, bed, bigBed)
            max_files: Maximum number of files to download (default 10)
        
        Returns:
            Path to download directory
        """
        dataset_details = self.get_dataset_details(dataset_id)
        if not dataset_details:
            raise Exception("Could not retrieve dataset details.")

        download_dir = os.path.join("downloads", "encode", dataset_id)
        os.makedirs(download_dir, exist_ok=True)
        import logging
        logger = logging.getLogger(__name__)

        # Default to processed file types if not specified
        if file_types is None:
            file_types = ['bigWig', 'bed', 'bigBed', 'bam']
        
        files_to_download = []
        for file_info in dataset_details.get('files', []):
            file_type = file_info.get('file_type', '')
            file_format = file_info.get('file_format', '')
            output_type = file_info.get('output_type', '')
            status = file_info.get('status', '')
            
            # Skip archived or revoked files
            if status not in ['released', 'in progress']:
                continue
            
            # Filter by file type
            if file_type in file_types or file_format in file_types:
                files_to_download.append(file_info)
        
        # Limit number of files
        files_to_download = files_to_download[:max_files]
        
        if not files_to_download:
            logger.warning(f"No files matching types {file_types} found for {dataset_id}")
            # Fallback: get first few files of any type
            files_to_download = dataset_details.get('files', [])[:3]

        logger.info(f"Downloading {len(files_to_download)} files for {dataset_id}")
        
        downloaded_count = 0
        for file_info in files_to_download:
            href = file_info.get('href', '')
            if not href:
                continue
                
            file_url = self.base_url + href
            accession = file_info.get('accession', 'file')
            
            # Extract filename from href (e.g., /files/ENCFF871EJI/@@download/ENCFF871EJI.fastq.gz)
            filename = href.split('/')[-1] if '/' in href else f"{accession}.dat"
            local_filepath = os.path.join(download_dir, filename)

            try:
                logger.info(f"Downloading {filename} from {file_url}")
                with requests.get(file_url, stream=True, timeout=600) as r:
                    r.raise_for_status()
                    total_size = int(r.headers.get('content-length', 0))
                    
                    with open(local_filepath, 'wb') as f:
                        downloaded = 0
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
                            downloaded += len(chunk)
                    
                    logger.info(f"Downloaded {filename} ({downloaded / 1e6:.1f} MB)")
                    downloaded_count += 1
                    
            except requests.exceptions.RequestException as e:
                logger.exception("Error downloading %s: %s", file_url, e)
                continue

        logger.info(f"Successfully downloaded {downloaded_count} files to {download_dir}")
        return download_dir

    def _generate_mock_datasets(self) -> list:
        """Generate a small set of mock dataset entries for local/dev use.

        This provides a stable structure expected by other modules (e.g. DatabaseManager.populate_sample_datasets).
        It's intentionally small and safe for unit tests and local development.
        """
        return [
            {
                'title': 'Mock H3K27ac ChIP-seq in K562',
                'description': 'Synthetic dataset for testing',
                'data_type': 'ChIP-seq',
                'tissue': 'blood',
                'organism': 'Human',
                'database': 'Mock',
                'accession': 'MOCK0001',
                'cell_type': 'K562',
                'cell_line': 'K562',
                'file_count': 2,
                'file_size_gb': 0.05,
                'publication_date': '2020-01-01',
                'treatment': None,
            },
            {
                'title': 'Mock RNA-seq in HepG2',
                'description': 'Synthetic gene expression dataset',
                'data_type': 'Gene Expression',
                'tissue': 'liver',
                'organism': 'Human',
                'database': 'Mock',
                'accession': 'MOCK0002',
                'cell_type': 'HepG2',
                'cell_line': 'HepG2',
                'file_count': 1,
                'file_size_gb': 0.01,
                'publication_date': '2021-06-01',
                'treatment': None,
            },
            {
                'title': 'Mock Hi-C sample',
                'description': 'Synthetic Hi-C interactions',
                'data_type': 'HiC',
                'tissue': 'brain',
                'organism': 'Human',
                'database': 'Mock',
                'accession': 'MOCK0003',
                'cell_type': None,
                'cell_line': None,
                'file_count': 1,
                'file_size_gb': 0.5,
                'publication_date': '2019-11-11',
                'treatment': None,
            }
        ]

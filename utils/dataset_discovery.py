import requests
from typing import Dict, List, Any, Optional
import os

class DatasetDiscovery:
    """A class to discover datasets from public genomics databases."""

    def __init__(self):
        self.base_url = "https://www.encodeproject.org"

    def search_datasets(self,
                       search_term: str = "",
                       data_types: Optional[List[str]] = None,
                       tissues: Optional[List[str]] = None,
                       organism: str = "Human") -> List[Dict[str, Any]]:
        """Search for datasets on the ENCODE platform."""
        search_params = {
            "type": "Experiment",
            "status": "released",
            "limit": "all",
            "frame": "object"
        }

        if search_term:
            search_params["searchTerm"] = search_term

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
            return self._format_results(search_results.get("@graph", []))
        except requests.exceptions.RequestException as e:
            logger.exception("Error querying ENCODE API: %s", e)
            return []

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

    def download_dataset(self, dataset_id: str) -> str:
        """Download all files for a given dataset from ENCODE."""
        dataset_details = self.get_dataset_details(dataset_id)
        if not dataset_details:
            raise Exception("Could not retrieve dataset details.")

        download_dir = os.path.join("downloads", dataset_id)
        os.makedirs(download_dir, exist_ok=True)
        import logging
        logger = logging.getLogger(__name__)

        def _sanitize_name(s: str) -> str:
            return ''.join(c if c.isalnum() or c in ('_', '-', '.') else '_' for c in (s or ''))[:200]

        for file_info in dataset_details.get('files', []):
            file_url = self.base_url + file_info.get('href')
            accession = file_info.get('accession', 'file')
            title = file_info.get('title', '')
            safe_name = f"{_sanitize_name(accession)}_{_sanitize_name(title)}"
            local_filename = os.path.join(download_dir, safe_name)

            try:
                with requests.get(file_url, stream=True) as r:
                    r.raise_for_status()
                    with open(local_filename, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
            except requests.exceptions.RequestException as e:
                logger.exception("Error downloading %s: %s", file_url, e)
                continue

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

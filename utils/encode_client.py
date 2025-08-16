import requests
from typing import Dict, List, Any, Optional
import os
from .database_client import DatabaseClient

class ENCODEClient(DatabaseClient):
    """A client for accessing the ENCODE database."""

    def __init__(self):
        self.base_url = "https://www.encodeproject.org"

    def search(self, search_term: str, **kwargs) -> List[Dict[str, Any]]:
        """Search for datasets on the ENCODE platform."""

        search_params = {
            "type": "Experiment",
            "status": "released",
            "limit": "all",
            "frame": "object"
        }

        if search_term:
            search_params["searchTerm"] = search_term

        # Add any additional kwargs to the search parameters
        search_params.update(kwargs)

        try:
            response = requests.get(f"{self.base_url}/search/", params=search_params, headers={'accept': 'application/json'})
            response.raise_for_status()
            search_results = response.json()
            return self._format_results(search_results.get("@graph", []))
        except requests.exceptions.RequestException as e:
            print(f"Error querying ENCODE API: {e}")
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
        try:
            response = requests.get(f"{self.base_url}/experiments/{dataset_id}/?frame=embedded", headers={'accept': 'application/json'})
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching dataset details from ENCODE: {e}")
            return None

    def download(self, dataset_id: str) -> str:
        """Download all files for a given dataset from ENCODE."""
        dataset_details = self.get_dataset_details(dataset_id)
        if not dataset_details:
            raise Exception("Could not retrieve dataset details.")

        download_dir = os.path.join("downloads", dataset_id)
        os.makedirs(download_dir, exist_ok=True)

        for file_info in dataset_details.get('files', []):
            file_url = self.base_url + file_info.get('href')
            local_filename = os.path.join(download_dir, file_info['accession'] + '_' + file_info.get('title', ''))

            try:
                with requests.get(file_url, stream=True) as r:
                    r.raise_for_status()
                    with open(local_filename, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
            except requests.exceptions.RequestException as e:
                print(f"Error downloading {file_url}: {e}")
                continue

        return download_dir

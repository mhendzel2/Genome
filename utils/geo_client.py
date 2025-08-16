import requests
import xml.etree.ElementTree as ET
from typing import Dict, List, Any
from .database_client import DatabaseClient

class GEOClient(DatabaseClient):
    """A client for accessing the GEO database."""

    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def search(self, search_term: str, **kwargs) -> List[Dict[str, Any]]:
        """Search for datasets on the GEO platform."""
        # E-utils require a two-step process: esearch to get UIDs, then esummary to get details.

        # 1. esearch
        esearch_params = {
            "db": "gds",
            "term": search_term,
            "retmax": "100", # Limit to 100 results for now
            "usehistory": "y"
        }
        try:
            esearch_response = requests.get(f"{self.base_url}/esearch.fcgi", params=esearch_params)
            esearch_response.raise_for_status()
            esearch_root = ET.fromstring(esearch_response.content)

            id_list = [id_node.text for id_node in esearch_root.findall(".//Id")]
            if not id_list:
                return []

            # 2. esummary
            esummary_params = {
                "db": "gds",
                "id": ",".join(id_list),
                "retmode": "json" # Request JSON for easier parsing
            }
            esummary_response = requests.get(f"{self.base_url}/esummary.fcgi", params=esummary_params)
            esummary_response.raise_for_status()
            esummary_data = esummary_response.json()

            return self._format_results(esummary_data['result'])

        except requests.exceptions.RequestException as e:
            print(f"Error querying GEO API: {e}")
            return []
        except ET.ParseError as e:
            print(f"Error parsing GEO XML response: {e}")
            return []

    def _format_results(self, results: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Format GEO results into the application's data structure."""
        formatted_results = []
        # The result is a dictionary of UIDs, so we iterate through its values
        for uid, result in results.items():
            if uid == 'uids': continue # Skip the list of uids

            formatted_results.append({
                'id': result.get('accession'),
                'title': result.get('title', 'No title available'),
                'description': result.get('summary', 'No description available'),
                'data_type': result.get('gdstype', 'N/A'),
                'tissue': 'N/A', # GEO metadata is complex, hard to extract tissue consistently
                'organism': result.get('taxon', 'N/A'),
                'database': 'GEO',
                'accession': result.get('accession'),
            })
        return formatted_results

    def download(self, dataset_id: str) -> str:
        """Download a dataset from GEO."""
        # Downloading from GEO requires constructing an FTP URL.
        # This is a placeholder for the more complex FTP download logic.
        return f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{dataset_id[:5]}nnn/{dataset_id}/matrix/{dataset_id}_series_matrix.txt.gz"

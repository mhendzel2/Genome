from pridepy import PridePy
from typing import Dict, List, Any
from .database_client import DatabaseClient

class PRIDEClient(DatabaseClient):
    """A client for accessing the PRIDE database."""

    def __init__(self):
        self.client = PridePy()

    def search(self, search_term: str, **kwargs) -> List[Dict[str, Any]]:
        """Search for projects on the PRIDE platform."""
        try:
            results = self.client.search_projects(keyword=search_term, page_size=100)
            return self._format_results(results)
        except Exception as e:
            print(f"Error querying PRIDE API: {e}")
            return []

    def _format_results(self, results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Format PRIDE results into the application's data structure."""
        formatted_results = []
        for result in results:
            formatted_results.append({
                'id': result.get('accession'),
                'title': result.get('title', 'No title available'),
                'description': result.get('projectDescription', 'No description available'),
                'data_type': 'Proteomics',
                'tissue': ', '.join(result.get('tissues', [])),
                'organism': ', '.join(result.get('species', [])),
                'database': 'PRIDE',
                'accession': result.get('accession'),
            })
        return formatted_results

    def download(self, dataset_id: str) -> str:
        """Download a dataset from PRIDE."""
        # This is a placeholder for the more complex download logic.
        # PRIDE datasets are typically downloaded via FTP or Aspera.
        return f"PRIDE download for {dataset_id} is not yet implemented."

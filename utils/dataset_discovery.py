from typing import Dict, List, Any
from .encode_client import ENCODEClient
from .geo_client import GEOClient
from .pride_client import PRIDEClient

class DatasetDiscovery:
    """A manager for different database clients."""

    def __init__(self):
        self.clients = {
            "ENCODE": ENCODEClient(),
            "GEO": GEOClient(),
            "PRIDE": PRIDEClient(),
        }

    def search(self, databases: List[str], search_term: str, **kwargs) -> List[Dict[str, Any]]:
        """Search across multiple databases."""
        all_results = []
        for db_name in databases:
            if db_name in self.clients:
                client = self.clients[db_name]
                results = client.search(search_term, **kwargs)
                all_results.extend(results)
        return all_results

    def download(self, database: str, dataset_id: str) -> str:
        """Download a dataset from a specific database."""
        if database in self.clients:
            client = self.clients[database]
            return client.download(dataset_id)
        else:
            raise ValueError(f"Database '{database}' not supported.")

from typing import Dict, List, Any, Optional

class DatabaseClient:
    """A base class for database clients."""

    def search(self, search_term: str) -> List[Dict[str, Any]]:
        """
        Search the database for a given term.

        This method should be implemented by subclasses.
        """
        raise NotImplementedError

    def download(self, dataset_id: str) -> str:
        """
        Download a dataset from the database.

        This method should be implemented by subclasses.
        """
        raise NotImplementedError

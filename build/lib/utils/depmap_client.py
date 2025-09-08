import requests
from typing import List, Dict, Any

class DepMapClient:
    """
    A client to interact with the DepMap data portal (via Figshare).
    """

    def __init__(self):
        self.figshare_api_base = "https://api.figshare.com/v2"

    def list_files_in_release(self, release_id: str = "27993248") -> List[Dict[str, Any]]:
        """
        Lists all files in a given DepMap release on Figshare.

        :param release_id: The Figshare article ID for the release.
        :return: A list of dictionaries, where each dictionary represents a file.
        """
        import logging
        logger = logging.getLogger(__name__)
        try:
            response = requests.get(f"{self.figshare_api_base}/articles/{release_id}/files")
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.exception("Error fetching file list from Figshare: %s", e)
            return []

    def download_file(self, download_url: str, local_filename: str) -> bool:
        """
        Downloads a file from a given URL.

        :param download_url: The URL to download the file from.
        :param local_filename: The local path to save the file to.
        :return: True if the download was successful, False otherwise.
        """
        import logging
        logger = logging.getLogger(__name__)
        try:
            with requests.get(download_url, stream=True) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            return True
        except requests.exceptions.RequestException as e:
            logger.exception("Error downloading file: %s", e)
            return False

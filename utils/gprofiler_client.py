from gprofiler import GProfiler
import pandas as pd

class GProfilerClient:
    """
    A client to interact with the g:Profiler API.
    """

    def __init__(self):
        """
        Initializes the GProfilerClient.
        """
        self.gp = GProfiler(user_agent='Jules/1.0')

    def perform_enrichment(self, query, organism='hsapiens'):
        """
        Performs enrichment analysis using g:Profiler.

        :param query: A list of gene IDs.
        :param organism: The organism to perform the analysis on (e.g., 'hsapiens').
        :return: A pandas DataFrame with the enrichment results.
        """
        enrichment_results = self.gp.profile(
            organism=organism,
            query=query
        )
        return pd.DataFrame(enrichment_results)

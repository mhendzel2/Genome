import requests
from typing import Dict, List, Any, Optional
import os
import logging
import json
import hashlib
from pathlib import Path

logger = logging.getLogger(__name__)


class DatasetDiscovery:
    """A class to discover datasets from public genomics databases."""

    def __init__(self, cache_dir: Optional[str] = None):
        self.base_url = "https://www.encodeproject.org"
        self.ols_base_url = "https://www.ebi.ac.uk/ols4/api"
        
        # Set up caching
        if cache_dir:
            self.cache_dir = Path(cache_dir)
        else:
            self.cache_dir = Path.home() / '.genome_cache' / 'ontology'
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Local ontology fallback dictionary
        self._local_ontology = self._load_local_ontology()

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
    
    def _load_local_ontology(self) -> Dict[str, List[str]]:
        """Load local ontology dictionary for fallback.
        
        Returns:
            Dictionary mapping terms to their synonyms
        """
        # Common genomics/biology term synonyms
        return {
            # Tissues
            'brain': ['cerebral cortex', 'cerebellum', 'hippocampus', 'neural tissue', 'CNS'],
            'liver': ['hepatic', 'hepatocyte'],
            'kidney': ['renal', 'nephron'],
            'heart': ['cardiac', 'myocardium', 'cardiomyocyte'],
            'lung': ['pulmonary', 'respiratory'],
            'blood': ['hematopoietic', 'peripheral blood', 'PBMC', 'leukocyte'],
            'skin': ['dermis', 'epidermis', 'keratinocyte'],
            'muscle': ['skeletal muscle', 'myocyte', 'muscular'],
            'bone': ['skeletal', 'osteocyte', 'osseous'],
            'breast': ['mammary', 'mammary gland'],
            'colon': ['colorectal', 'intestinal', 'colonic'],
            'pancreas': ['pancreatic', 'islet'],
            'prostate': ['prostatic'],
            'ovary': ['ovarian'],
            'testis': ['testicular'],
            
            # Cell types
            'stem cell': ['pluripotent', 'progenitor', 'iPSC', 'ESC', 'hESC'],
            'T cell': ['T lymphocyte', 'CD4+', 'CD8+', 'T-cell'],
            'B cell': ['B lymphocyte', 'B-cell'],
            'macrophage': ['monocyte', 'phagocyte'],
            'fibroblast': ['stromal cell'],
            'epithelial': ['epithelium'],
            'endothelial': ['vascular', 'HUVEC'],
            'neuron': ['neural', 'neuronal'],
            'astrocyte': ['glial', 'glia'],
            
            # Assays
            'ChIP-seq': ['chromatin immunoprecipitation', 'ChIP sequencing'],
            'RNA-seq': ['transcriptome', 'gene expression', 'mRNA sequencing'],
            'ATAC-seq': ['chromatin accessibility', 'open chromatin'],
            'DNase-seq': ['DNase hypersensitivity', 'DHS'],
            'Hi-C': ['chromosome conformation', '3D genome', 'chromatin interaction'],
            'methylation': ['bisulfite sequencing', 'WGBS', 'RRBS', 'DNA methylation'],
            'proteomics': ['mass spectrometry', 'protein expression'],
            'WGS': ['whole genome sequencing'],
            'WES': ['whole exome sequencing', 'exome'],
            
            # Histone marks
            'H3K4me3': ['active promoter', 'promoter mark'],
            'H3K27ac': ['active enhancer', 'enhancer mark'],
            'H3K27me3': ['repressive', 'polycomb'],
            'H3K9me3': ['heterochromatin', 'repressive mark'],
            'H3K36me3': ['gene body', 'transcribed region'],
            
            # Cancer types
            'breast cancer': ['BRCA', 'breast carcinoma', 'mammary tumor'],
            'lung cancer': ['LUAD', 'LUSC', 'lung adenocarcinoma'],
            'colorectal cancer': ['COAD', 'colon cancer', 'rectal cancer'],
            'leukemia': ['AML', 'ALL', 'CML', 'CLL', 'blood cancer'],
            'lymphoma': ['Hodgkin', 'non-Hodgkin', 'NHL'],
            'melanoma': ['skin cancer', 'SKCM'],
            'glioblastoma': ['GBM', 'brain tumor', 'glioma'],
        }
    
    def semantic_search(self, query: str, 
                       n_synonyms: int = 10,
                       ontologies: Optional[List[str]] = None,
                       use_cache: bool = True) -> Dict[str, Any]:
        """
        Perform semantic search to expand query terms using ontology services.
        
        Queries the EBI Ontology Lookup Service (OLS) to find synonyms and
        related terms, enabling more comprehensive dataset discovery.
        
        Args:
            query: Search term to expand
            n_synonyms: Maximum number of synonyms to return
            ontologies: List of ontology IDs to search (e.g., ['EFO', 'UBERON', 'CL'])
                       If None, searches across all ontologies
            use_cache: Whether to use cached results
        
        Returns:
            Dictionary with:
            - original_term: The input query
            - expanded_terms: List of synonyms and related terms
            - ontology_matches: Detailed ontology term matches
            - source: Where the terms came from ('ols', 'cache', 'local')
        """
        result = {
            'original_term': query,
            'expanded_terms': [query],  # Always include original
            'ontology_matches': [],
            'source': 'local'
        }
        
        # Check cache first
        cache_key = self._get_cache_key(query, ontologies)
        if use_cache:
            cached = self._load_from_cache(cache_key)
            if cached:
                cached['source'] = 'cache'
                return cached
        
        # Try OLS API
        ols_result = self._query_ols(query, n_synonyms, ontologies)
        
        if ols_result['expanded_terms']:
            result = ols_result
            result['source'] = 'ols'
            
            # Cache successful results
            if use_cache:
                self._save_to_cache(cache_key, result)
        else:
            # Fallback to local dictionary
            result = self._local_search(query, n_synonyms)
            result['source'] = 'local'
        
        return result
    
    def _query_ols(self, query: str, n_synonyms: int,
                  ontologies: Optional[List[str]] = None) -> Dict[str, Any]:
        """Query the EBI Ontology Lookup Service."""
        result = {
            'original_term': query,
            'expanded_terms': [query],
            'ontology_matches': []
        }
        
        try:
            # Build search URL
            search_url = f"{self.ols_base_url}/search"
            params = {
                'q': query,
                'rows': n_synonyms * 2,  # Request more to filter
                'exact': 'false',
                'local': 'false'
            }
            
            if ontologies:
                params['ontology'] = ','.join(ontologies)
            else:
                # Default to biology-relevant ontologies
                params['ontology'] = 'efo,uberon,cl,doid,go,hp'
            
            response = requests.get(search_url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()
            
            docs = data.get('response', {}).get('docs', [])
            
            seen_terms = {query.lower()}
            
            for doc in docs:
                # Get the term label
                label = doc.get('label', '')
                if label and label.lower() not in seen_terms:
                    result['expanded_terms'].append(label)
                    seen_terms.add(label.lower())
                
                # Get synonyms
                synonyms = doc.get('synonym', [])
                if isinstance(synonyms, str):
                    synonyms = [synonyms]
                
                for syn in synonyms:
                    if syn and syn.lower() not in seen_terms:
                        result['expanded_terms'].append(syn)
                        seen_terms.add(syn.lower())
                
                # Store ontology match details
                result['ontology_matches'].append({
                    'term': label,
                    'ontology': doc.get('ontology_name', ''),
                    'iri': doc.get('iri', ''),
                    'description': doc.get('description', [''])[0] if doc.get('description') else '',
                    'synonyms': synonyms[:5]
                })
                
                if len(result['expanded_terms']) >= n_synonyms:
                    break
            
            # Limit to requested number
            result['expanded_terms'] = result['expanded_terms'][:n_synonyms]
            result['ontology_matches'] = result['ontology_matches'][:5]
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"OLS API request failed: {e}")
        except Exception as e:
            logger.warning(f"Error querying OLS: {e}")
        
        return result
    
    def _local_search(self, query: str, n_synonyms: int) -> Dict[str, Any]:
        """Search local ontology dictionary for synonyms."""
        result = {
            'original_term': query,
            'expanded_terms': [query],
            'ontology_matches': []
        }
        
        query_lower = query.lower()
        
        # Check for exact matches
        if query_lower in self._local_ontology:
            synonyms = self._local_ontology[query_lower]
            result['expanded_terms'].extend(synonyms[:n_synonyms - 1])
        
        # Check for partial matches
        for term, synonyms in self._local_ontology.items():
            if query_lower in term or any(query_lower in s.lower() for s in synonyms):
                result['expanded_terms'].append(term)
                for syn in synonyms:
                    if len(result['expanded_terms']) < n_synonyms:
                        if syn.lower() not in [t.lower() for t in result['expanded_terms']]:
                            result['expanded_terms'].append(syn)
            
            if len(result['expanded_terms']) >= n_synonyms:
                break
        
        # Remove duplicates while preserving order
        seen = set()
        unique_terms = []
        for term in result['expanded_terms']:
            if term.lower() not in seen:
                seen.add(term.lower())
                unique_terms.append(term)
        
        result['expanded_terms'] = unique_terms[:n_synonyms]
        
        return result
    
    def _get_cache_key(self, query: str, ontologies: Optional[List[str]]) -> str:
        """Generate cache key for query."""
        key_data = f"{query.lower()}:{sorted(ontologies) if ontologies else 'all'}"
        return hashlib.md5(key_data.encode()).hexdigest()
    
    def _load_from_cache(self, cache_key: str) -> Optional[Dict[str, Any]]:
        """Load cached ontology expansion result."""
        cache_file = self.cache_dir / f"{cache_key}.json"
        
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.debug(f"Failed to load cache: {e}")
        
        return None
    
    def _save_to_cache(self, cache_key: str, data: Dict[str, Any]) -> None:
        """Save ontology expansion result to cache."""
        cache_file = self.cache_dir / f"{cache_key}.json"
        
        try:
            # Remove source field from cached data (will be set on retrieval)
            cache_data = {k: v for k, v in data.items() if k != 'source'}
            
            with open(cache_file, 'w') as f:
                json.dump(cache_data, f)
        except Exception as e:
            logger.debug(f"Failed to save cache: {e}")
    
    def clear_cache(self) -> int:
        """Clear all cached ontology expansions.
        
        Returns:
            Number of cache files removed
        """
        count = 0
        for cache_file in self.cache_dir.glob('*.json'):
            try:
                cache_file.unlink()
                count += 1
            except Exception:
                pass
        
        logger.info(f"Cleared {count} cached ontology expansions")
        return count
    
    def search_datasets_semantic(self, query: str,
                                n_synonyms: int = 5,
                                **kwargs) -> List[Dict[str, Any]]:
        """
        Search datasets using semantic expansion of query terms.
        
        Expands the search query using ontology synonyms and searches
        across multiple terms for broader results.
        
        Args:
            query: Search term
            n_synonyms: Number of synonyms to use for expansion
            **kwargs: Additional parameters passed to search_datasets
        
        Returns:
            List of dataset results from expanded search
        """
        # Expand query semantically
        expansion = self.semantic_search(query, n_synonyms=n_synonyms)
        expanded_terms = expansion['expanded_terms']
        
        logger.info(f"Expanded '{query}' to {len(expanded_terms)} terms: {expanded_terms}")
        
        # Search with each expanded term
        all_results = []
        seen_ids = set()
        
        for term in expanded_terms:
            try:
                results = self.search_datasets(search_term=term, **kwargs)
                
                for r in results:
                    result_id = r.get('id') or r.get('accession')
                    if result_id and result_id not in seen_ids:
                        seen_ids.add(result_id)
                        r['matched_term'] = term
                        all_results.append(r)
            except Exception as e:
                logger.debug(f"Search failed for term '{term}': {e}")
        
        return all_results

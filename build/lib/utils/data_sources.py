"""
Registry of public genomics data sources and helper utilities.
This module provides a lightweight list of public data portals with
metadata and sample endpoints that can be used to pre-program dataset
lookups or automated downloads.

Add entries as needed; this is intentionally read-only metadata—actual
API clients should live elsewhere (e.g., dataset_discovery, depmap_client).
"""
from typing import List, Dict

DATA_SOURCES = [
    {
        'name': 'ENCODE',
        'homepage': 'https://www.encodeproject.org',
        'api_search': 'https://www.encodeproject.org/search/',
        'notes': 'Large catalog of functional genomics experiments (ChIP-seq, RNA-seq, ATAC-seq, etc).'
    },
    {
        'name': 'GEO',
        'homepage': 'https://www.ncbi.nlm.nih.gov/geo/',
        'api_search': 'https://www.ncbi.nlm.nih.gov/gds',
        'notes': 'Gene Expression Omnibus; many processed and raw gene expression datasets.'
    },
    {
        'name': 'SRA',
        'homepage': 'https://www.ncbi.nlm.nih.gov/sra',
        'api_search': 'https://trace.ncbi.nlm.nih.gov/Traces/sra/',
        'notes': 'Sequence Read Archive — raw sequencing reads. Use fasterq-dump / sra-tools for downloads.'
    },
    {
        'name': 'ENA',
        'homepage': 'https://www.ebi.ac.uk/ena',
        'api_search': 'https://www.ebi.ac.uk/ena/browser/api/',
        'notes': 'European Nucleotide Archive; alternative to SRA with convenient ftp/http endpoints.'
    },
    {
        'name': 'UCSC',
        'homepage': 'https://genome.ucsc.edu',
        'api_search': 'https://api.genome.ucsc.edu/',
        'notes': 'Genome assemblies and many annotation tracks; useful for prebuilt bigWig/bigBed files.'
    },
    {
        'name': 'TCGA (GDC)',
        'homepage': 'https://portal.gdc.cancer.gov',
        'api_search': 'https://api.gdc.cancer.gov',
        'notes': 'Cancer genomics; requires controlled access for many raw files.'
    },
    {
        'name': 'DepMap/Figshare',
        'homepage': 'https://depmap.org',
        'api_search': 'https://api.figshare.com/v2/articles/<release_id>/files',
        'notes': 'DepMap data releases published via Figshare — large CSV/TSV tables useful for cell-line analyses.'
    }
]


def list_sources() -> List[Dict]:
    """Return the list of available data sources (readonly)."""
    return DATA_SOURCES.copy()


def find_by_name(name: str) -> Dict:
    """Find a source by case-insensitive name; raises KeyError if not found."""
    nl = name.lower()
    for s in DATA_SOURCES:
        if s['name'].lower() == nl:
            return s
    raise KeyError(f"Data source not found: {name}")

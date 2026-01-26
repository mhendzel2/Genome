#!/usr/bin/env python
"""Test download functionality"""
import sys
sys.path.insert(0, '.')

import logging
logging.basicConfig(level=logging.DEBUG)

from utils.dataset_discovery import DatasetDiscovery
from utils.geo_client import GEOClient

print("=" * 60)
print("Testing Download Functions")
print("=" * 60)

# Initialize clients
discovery = DatasetDiscovery()
geo_client = GEOClient()

# Test 1: Get ENCODE dataset details
print("\n--- Test 1: Get ENCODE dataset details for ENCSR675FLJ ---")
details = discovery.get_dataset_details("ENCSR675FLJ")
if details:
    print(f"  Title: {details.get('description', 'N/A')}")
    print(f"  Status: {details.get('status', 'N/A')}")
    files = details.get('files', [])
    print(f"  Files count: {len(files)}")
    if files:
        print("  First 3 files:")
        for f in files[:3]:
            print(f"    - {f.get('accession')}: {f.get('file_type')} | href: {f.get('href', 'N/A')[:50]}")
else:
    print("  FAILED: Could not get details")

# Test 2: Check if files have download URLs
print("\n--- Test 2: Check ENCODE file download structure ---")
if details and details.get('files'):
    file_info = details['files'][0]
    print(f"  File keys: {list(file_info.keys())[:10]}")
    print(f"  href: {file_info.get('href')}")
    print(f"  accession: {file_info.get('accession')}")
    print(f"  file_type: {file_info.get('file_type')}")
    print(f"  output_type: {file_info.get('output_type')}")
    
    # Test download URL construction
    base_url = "https://www.encodeproject.org"
    download_url = base_url + file_info.get('href', '')
    print(f"  Full download URL: {download_url}")
    
    # Check if URL is accessible
    import requests
    try:
        r = requests.head(download_url, allow_redirects=True, timeout=10)
        print(f"  URL status: {r.status_code}")
        print(f"  Content-Type: {r.headers.get('Content-Type', 'N/A')}")
        print(f"  Content-Length: {r.headers.get('Content-Length', 'N/A')}")
    except Exception as e:
        print(f"  Error checking URL: {e}")

# Test 3: GEO download URLs
print("\n--- Test 3: Check GEO download URL for GSE135771 ---")
accession = "GSE135771"
series_stub = accession[:-3] + 'nnn'  # GSE135nnn
ftp_base = "https://ftp.ncbi.nlm.nih.gov/geo"
matrix_url = f"{ftp_base}/series/{series_stub}/{accession}/matrix/{accession}_series_matrix.txt.gz"
print(f"  Matrix URL: {matrix_url}")

import requests
try:
    r = requests.head(matrix_url, allow_redirects=True, timeout=10)
    print(f"  URL status: {r.status_code}")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "=" * 60)
print("Tests complete")
print("=" * 60)

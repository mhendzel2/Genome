#!/usr/bin/env python
"""Test download functionality after fixes"""
import sys
sys.path.insert(0, '.')

import logging
logging.basicConfig(level=logging.INFO)

from utils.dataset_discovery import DatasetDiscovery
from utils.geo_client import GEOClient

print("=" * 60)
print("Testing Download Functions (After Fixes)")
print("=" * 60)

# Initialize clients
discovery = DatasetDiscovery()
geo_client = GEOClient()

# Test 1: Download GEO series matrix for GSE135771
print("\n--- Test 1: Download GEO series matrix for GSE135771 ---")
try:
    files = geo_client.download_series_matrix("GSE135771", download_dir="downloads/geo/GSE135771")
    print(f"  SUCCESS: Downloaded {len(files)} file(s)")
    for f in files:
        print(f"    - {f}")
except Exception as e:
    print(f"  FAILED: {e}")

# Test 2: Check GEO supplementary files for GSE135771
print("\n--- Test 2: List GEO supplementary files for GSE135771 ---")
try:
    import requests
    stub = geo_client._get_series_stub("GSE135771")
    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/GSE135771/suppl/"
    r = requests.get(suppl_url, timeout=30)
    if r.status_code == 200:
        import re
        file_pattern = re.compile(r'href="([^"]+\.(txt|csv|gz|tar|bed|wig|bw|bam|fastq|CEL|RAW)[^"]*)"', re.IGNORECASE)
        files_found = file_pattern.findall(r.text)
        filenames = [f[0] for f in files_found if not f[0].startswith('/')]
        print(f"  Found {len(filenames)} supplementary files:")
        for f in filenames[:5]:
            print(f"    - {f}")
    else:
        print(f"  Status: {r.status_code}")
except Exception as e:
    print(f"  Error: {e}")

# Test 3: List ENCODE files for ENCSR675FLJ
print("\n--- Test 3: List ENCODE files for ENCSR675FLJ ---")
details = discovery.get_dataset_details("ENCSR675FLJ")
if details:
    files = details.get('files', [])
    print(f"  Total files: {len(files)}")
    
    # Count file types
    file_types = {}
    for f in files:
        ft = f.get('file_type', 'unknown')
        file_types[ft] = file_types.get(ft, 0) + 1
    
    print("  File types:")
    for ft, count in sorted(file_types.items()):
        print(f"    - {ft}: {count}")
    
    # Show a few bigWig/bed files
    processed = [f for f in files if f.get('file_type') in ['bigWig', 'bed', 'bigBed']]
    print(f"\n  Processed files (bigWig/bed): {len(processed)}")
    for f in processed[:3]:
        print(f"    - {f.get('accession')}: {f.get('file_type')} | {f.get('output_type')}")

print("\n" + "=" * 60)
print("Tests complete")
print("=" * 60)

#!/usr/bin/env python
"""Test dataset discovery functionality"""
import sys
sys.path.insert(0, '.')

import logging
logging.basicConfig(level=logging.INFO)

from utils.dataset_discovery import DatasetDiscovery
from utils.geo_client import GEOClient

print("=" * 60)
print("Testing Dataset Discovery")
print("=" * 60)

# Initialize clients
discovery = DatasetDiscovery()
geo_client = GEOClient()

# Test 1: ENCODE search for H4K5ac in H1-hESC (ENCSR675FLJ)
print("\n--- Test 1: ENCODE search for 'ENCSR675FLJ' ---")
results = discovery.search_datasets(search_term="ENCSR675FLJ", use_semantic_search=False)
print(f"Results: {len(results)}")
for r in results[:5]:
    print(f"  - {r['accession']}: {r['title'][:60]}...")

# Test 2: ENCODE search for H4K5ac
print("\n--- Test 2: ENCODE search for 'H4K5ac' ---")
results = discovery.search_datasets(search_term="H4K5ac", use_semantic_search=False)
print(f"Results: {len(results)}")
for r in results[:5]:
    print(f"  - {r['accession']}: {r['title'][:60]}...")

# Test 3: ENCODE search for H4K5ac H1-hESC
print("\n--- Test 3: ENCODE search for 'H4K5ac H1-hESC' ---")
results = discovery.search_datasets(search_term="H4K5ac H1-hESC", use_semantic_search=False)
print(f"Results: {len(results)}")
for r in results[:5]:
    print(f"  - {r['accession']}: {r['title'][:60]}...")

# Test 4: GEO search for GSE135771 (NSL complex)
print("\n--- Test 4: GEO search for 'GSE135771' ---")
results = geo_client.search_datasets(search_term="GSE135771", max_results=10)
print(f"Results: {len(results)}")
for r in results[:5]:
    print(f"  - {r['accession']}: {r['title'][:60]}...")

# Test 5: GEO search for KANSL1
print("\n--- Test 5: GEO search for 'KANSL1' ---")
results = geo_client.search_datasets(search_term="KANSL1", max_results=10)
print(f"Results: {len(results)}")
for r in results[:5]:
    print(f"  - {r['accession']}: {r['title'][:60]}...")

# Test 6: GEO search for MSL1 ChIP-seq
print("\n--- Test 6: GEO search for 'MSL1 ChIP-seq' ---")
results = geo_client.search_datasets(search_term="MSL1 ChIP-seq", max_results=10)
print(f"Results: {len(results)}")
for r in results[:5]:
    print(f"  - {r['accession']}: {r['title'][:60]}...")

# Test 7: Direct ENCODE API test for specific accession
print("\n--- Test 7: Direct ENCODE API fetch for ENCSR675FLJ ---")
import requests
try:
    resp = requests.get("https://www.encodeproject.org/experiments/ENCSR675FLJ/?frame=object", 
                       headers={'accept': 'application/json'}, timeout=30)
    if resp.status_code == 200:
        data = resp.json()
        print(f"  Accession: {data.get('accession')}")
        print(f"  Description: {data.get('description', 'N/A')[:80]}")
        print(f"  Assay: {data.get('assay_title')}")
        print(f"  Target: {data.get('target', {}).get('label', 'N/A') if isinstance(data.get('target'), dict) else data.get('target', 'N/A')}")
        print(f"  Status: {data.get('status')}")
    else:
        print(f"  HTTP {resp.status_code}: {resp.text[:200]}")
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "=" * 60)
print("Tests complete")
print("=" * 60)

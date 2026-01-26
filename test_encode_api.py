#!/usr/bin/env python
"""Test ENCODE API formats"""
import requests

# Test different URL formats
tests = [
    # Original format with frame=object (current code)
    "https://www.encodeproject.org/search/?type=Experiment&searchTerm=H4K5ac&frame=object",
    # With format=json 
    "https://www.encodeproject.org/search/?type=Experiment&searchTerm=H4K5ac&format=json",
    # Simpler format
    "https://www.encodeproject.org/search/?type=Experiment&searchTerm=H4K5ac",
    # With limit
    "https://www.encodeproject.org/search/?type=Experiment&searchTerm=H4K5ac&limit=10&format=json",
    # Direct accession lookup
    "https://www.encodeproject.org/search/?type=Experiment&accession=ENCSR675FLJ&format=json",
    # With organism filter (the problematic one)
    "https://www.encodeproject.org/search/?type=Experiment&searchTerm=H4K5ac&replicates.library.biosample.donor.organism.scientific_name=Human",
    # Try Homo sapiens instead of Human
    "https://www.encodeproject.org/search/?type=Experiment&searchTerm=H4K5ac&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens",
]

for url in tests:
    try:
        r = requests.get(url, headers={'accept': 'application/json'}, timeout=30)
        print(f"Status {r.status_code} | {url[:90]}...")
        if r.status_code == 200:
            data = r.json()
            if '@graph' in data:
                print(f"   -> Found {len(data['@graph'])} results")
            elif 'notification' in data:
                print(f"   -> Notification: {data.get('notification')}")
            else:
                print(f"   -> Keys: {list(data.keys())[:5]}")
    except Exception as e:
        print(f"Error: {e}")
    print()

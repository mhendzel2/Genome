#!/usr/bin/env python
"""Test correct GEO URL format"""
import requests

# GSE135771 -> series_stub should be GSE135xxx for 6-digit, GSE13xnnn for 5-digit etc
# Actually the rule is: take all but last 3 digits, add nnn
# GSE135771 -> GSE135 + nnn -> GSE135nnn  (WRONG - this is what we had)
# Actual format: GSEnnnnn -> take first 3 chars after GSE, then "nnn"

accession = "GSE135771"
# Correct format: GSE + first digits (varies) + nnn
# For 6+ digit accessions: GSE135nnn (first 3 digits of number)
# For smaller: different

# Let's test different URL patterns
ftp_base = "https://ftp.ncbi.nlm.nih.gov/geo"

# Pattern 1: GSE135nnn (our current - WRONG)
url1 = f"{ftp_base}/series/GSE135nnn/{accession}/matrix/"

# Pattern 2: GSE13nnn  
url2 = f"{ftp_base}/series/GSE13nnn/{accession}/matrix/"

# Pattern 3: Calculate properly - for GSE135771, it's GSE135nnn
# The series number is 135771, so we need to bucket by first N-3 digits
# 135771 has 6 digits, so bucket is 135 -> GSE135nnn

# Actually let me check the NCBI format more carefully
# For a 6-digit number like 135771, the format is GSE + 135 + nnn = GSE135nnn
# But for GSE12345 (5 digits), it's GSE + 12 + nnn = GSE12nnn

def get_geo_series_stub(accession):
    """Get correct GEO FTP path bucket for a series accession."""
    # Extract number part
    num_str = accession.replace("GSE", "")
    # Take all but last 3 digits, padded with zeros if needed
    if len(num_str) <= 3:
        return "GSEnnn"
    else:
        prefix = num_str[:-3]
        return f"GSE{prefix}nnn"

print("Testing GEO URL patterns...")
for acc in ["GSE135771", "GSE12345", "GSE1234", "GSE123", "GSE283747"]:
    stub = get_geo_series_stub(acc)
    print(f"  {acc} -> {stub}")

# Test the correct URL for GSE135771
stub = get_geo_series_stub("GSE135771")
matrix_url = f"{ftp_base}/series/{stub}/{accession}/matrix/{accession}_series_matrix.txt.gz"
print(f"\nTesting URL: {matrix_url}")
r = requests.head(matrix_url, timeout=10)
print(f"Status: {r.status_code}")

# Check the matrix directory
print(f"\nListing matrix directory...")
matrix_dir_url = f"{ftp_base}/series/{stub}/{accession}/matrix/"
r3 = requests.get(matrix_dir_url, timeout=10)
print(f"Matrix dir status: {r3.status_code}")
if r3.status_code == 200:
    print("Matrix directory contents:")
    print(r3.text[:1500])

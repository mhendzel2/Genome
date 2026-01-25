#!/usr/bin/env python
"""Quick import test script"""
import sys
import traceback

print("Testing imports...")

# Test 1: genomics_analysis
print("\n1. Testing analysis.genomics_analysis...")
try:
    from analysis.genomics_analysis import GenomicsAnalyzer
    print("   OK: GenomicsAnalyzer imported")
except Exception as e:
    print(f"   FAILED: {e}")
    traceback.print_exc()

# Test 2: analyzer
print("\n2. Testing analysis.analyzer...")
try:
    from analysis.analyzer import GenomicsAnalyzer as GA2
    print("   OK: analyzer.GenomicsAnalyzer imported")
except Exception as e:
    print(f"   FAILED: {e}")
    traceback.print_exc()

# Test 3: r_integration
print("\n3. Testing utils.r_integration...")
try:
    from utils.r_integration import RIntegration
    print("   OK: RIntegration imported")
except Exception as e:
    print(f"   FAILED: {e}")
    traceback.print_exc()

# Test 4: analyzer_adapter
print("\n4. Testing analysis.analyzer_adapter...")
try:
    from utils.r_integration import RIntegration
    r_int = RIntegration()
    from analysis.analyzer_adapter import AnalyzerAdapter
    adapter = AnalyzerAdapter(r_int)
    print("   OK: AnalyzerAdapter created")
except Exception as e:
    print(f"   FAILED: {e}")
    traceback.print_exc()

# Test 5: cooler
print("\n5. Testing cooler...")
try:
    import cooler
    print("   OK: cooler imported")
except Exception as e:
    print(f"   FAILED: {e}")

# Test 6: cooltools
print("\n6. Testing cooltools...")
try:
    import cooltools
    print("   OK: cooltools imported")
except Exception as e:
    print(f"   FAILED (optional on Windows): {e}")

print("\nDone.")

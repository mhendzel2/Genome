import sys
import os
import pytest
import pandas as pd

# Ensure project root is importable
proj_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)

from analysis.analyzer_adapter import AnalyzerAdapter


class MockAnalyzer:
    def basic_statistics(self, uploaded_files):
        return {'total_regions': 10}

    def differential_expression(self, uploaded_files, **params):
        return pd.DataFrame({'gene': ['A', 'B'], 'log2FC': [1.0, -1.2], 'pvalue': [0.001, 0.02]})

    def enrichment_analysis(self, uploaded_files, **params):
        return [{'name': 'pathway1', 'p_value': 0.0001}]

    def quality_control(self, uploaded_files):
        return {'sample1.bed': {'quality_score': 'PASS'}}


def test_adapter_wraps_dict_result():
    adapter = AnalyzerAdapter.__new__(AnalyzerAdapter)
    adapter._analyzer = MockAnalyzer()

    res = adapter.basic_statistics({'f': {'file': None}})
    assert isinstance(res, dict)
    assert res['status'] == 'ok'
    assert isinstance(res['result'], dict)


def test_adapter_wraps_dataframe_result():
    adapter = AnalyzerAdapter.__new__(AnalyzerAdapter)
    adapter._analyzer = MockAnalyzer()

    res = adapter.differential_expression({'f': {'file': None}})
    assert isinstance(res, dict)
    assert res['status'] == 'ok'
    assert isinstance(res['result'], pd.DataFrame)


def test_adapter_wraps_list_result():
    adapter = AnalyzerAdapter.__new__(AnalyzerAdapter)
    adapter._analyzer = MockAnalyzer()

    res = adapter.enrichment_analysis({'f': {'file': None}})
    assert res['status'] == 'ok'
    assert isinstance(res['result'], list)


def test_adapter_qc_result():
    adapter = AnalyzerAdapter.__new__(AnalyzerAdapter)
    adapter._analyzer = MockAnalyzer()

    res = adapter.quality_control({'f': {'file': None}})
    assert res['status'] == 'ok'
    assert isinstance(res['result'], dict)

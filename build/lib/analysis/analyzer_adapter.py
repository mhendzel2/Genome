"""Adapter to provide a stable analyzer API for the app.

This adapter wraps whichever GenomicsAnalyzer implementation is present and
normalizes outputs to predictable dictionaries with documented keys.
"""
from typing import Any, Dict, Optional
import importlib
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class AnalyzerAdapter:
    """Adapter for genomics analyzers.

    Methods return a dict with at least these keys where applicable:
    - 'status': 'ok' or 'error'
    - 'result': DataFrame or dict/list describing results
    - 'message': optional human message
    """

    def __init__(self, r_integration=None):
        # Try to import the main analyzer implementations; prefer the one that
        # accepts r_integration in the constructor (analysis/genomics_analysis.py)
        analyzer = None
        try:
            mod = importlib.import_module('analysis.genomics_analysis')
            analyzer_cls = getattr(mod, 'GenomicsAnalyzer', None)
            if analyzer_cls:
                analyzer = analyzer_cls(r_integration) if r_integration is not None else analyzer_cls(r_integration)
        except Exception:
            logger.debug('analysis.genomics_analysis not available or failed to import')

        if analyzer is None:
            try:
                mod2 = importlib.import_module('analysis.analyzer')
                analyzer_cls2 = getattr(mod2, 'GenomicsAnalyzer', None)
                if analyzer_cls2:
                    analyzer = analyzer_cls2()
            except Exception:
                logger.debug('analysis.analyzer not available')

        if analyzer is None:
            raise ImportError('No GenomicsAnalyzer implementation available')

        self._analyzer = analyzer

    def _wrap(self, func_name: str, *args, **kwargs) -> Dict[str, Any]:
        try:
            func = getattr(self._analyzer, func_name)
            res = func(*args, **kwargs)

            # Normalize some common return types
            if isinstance(res, pd.DataFrame):
                return {'status': 'ok', 'result': res, 'message': ''}
            if isinstance(res, dict):
                return {'status': 'ok', 'result': res, 'message': ''}
            if isinstance(res, list):
                return {'status': 'ok', 'result': res, 'message': ''}

            # Fallback: wrap whatever is returned
            return {'status': 'ok', 'result': res, 'message': ''}

        except Exception as e:
            logger.exception('Analyzer %s failed', func_name)
            return {'status': 'error', 'result': None, 'message': str(e)}

    def basic_statistics(self, uploaded_files: Dict[str, Any]) -> Dict[str, Any]:
        return self._wrap('basic_statistics', uploaded_files)

    def differential_expression(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        return self._wrap('differential_expression', uploaded_files, **params)

    def enrichment_analysis(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        return self._wrap('enrichment_analysis', uploaded_files, **params)

    def quality_control(self, uploaded_files: Dict[str, Any]) -> Dict[str, Any]:
        return self._wrap('quality_control', uploaded_files)

    def peak_calling(self, uploaded_files: Dict[str, Any], **params) -> Dict[str, Any]:
        return self._wrap('peak_calling', uploaded_files, **params)

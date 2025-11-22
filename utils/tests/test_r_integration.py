import os
import json
import tempfile
import subprocess
from types import SimpleNamespace
import io

import pytest

from utils.r_integration import RIntegration


class DummyCompletedProcess:
    def __init__(self, returncode=0, stdout='{}', stderr=''):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_execute_r_script_writes_files_and_replaces_placeholders(monkeypatch, tmp_path):
    ri = RIntegration()

    # Force r_available True for testing path (we'll mock subprocess)
    ri.r_available = True

    # Prepare a small script with placeholders
    script = 'cat("PLACEHOLDER: ", "{file:myfile}", "\\n")\n'

    # Data files payload
    data_files = {'myfile': 'chr1\t100\t200\nchr1\t300\t400\n'}

    captured = {}

    def fake_run(cmd, cwd, capture_output, text, timeout, env):
        # Validate that cwd contains the temp dir and an analysis.R exists
        assert os.path.isdir(cwd)
        script_path = os.path.join(cwd, 'analysis.R')
        assert os.path.exists(script_path)

        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # The placeholder should have been replaced with a path inside cwd
        assert '{file:myfile}' not in content
        assert 'myfile' in content

        # Simulate R stdout returning the path for verification
        return DummyCompletedProcess(returncode=0, stdout='RESPONSE')

    monkeypatch.setattr(subprocess, 'run', fake_run)

    out = ri._execute_r_script(script, data_files=data_files)
    assert 'RESPONSE' in out


def test_basic_statistics_fallback_json_parsing(monkeypatch):
    ri = RIntegration()
    ri.r_available = True

    # Mock _execute_r_script to return JSON
    def fake_exec(script, data_files=None, args_files=None, extra_args=None):
        return json.dumps({'total_regions': 2, 'mean_length': 100})

    monkeypatch.setattr(ri, '_execute_r_script', fake_exec)

    fake_file = SimpleNamespace()
    fake_file.read = lambda: b'chr1\t100\t200\nchr1\t300\t400\n'
    fake_file.seek = lambda x: None

    data_files = {'test.bed': {'file': fake_file}}

    res = ri.basic_statistics(data_files)
    assert 'test.bed' in res

    # If there was an error, print it
    if 'error' in res['test.bed']:
        pytest.fail(f"Basic statistics failed with error: {res['test.bed']['error']}")

    assert res['test.bed']['total_regions'] == 2


def test_peak_calling_parses_json_and_attaches_metadata(monkeypatch):
    ri = RIntegration()
    ri.r_available = True

    # Make _execute_r_script return a JSON array
    def fake_exec(script, data_files=None, args_files=None, extra_args=None):
        return json.dumps([{'chr': 'chr1', 'start': 100, 'end': 200, 'score': 10}])

    monkeypatch.setattr(ri, '_execute_r_script', fake_exec)

    fake_file = io.StringIO('chr1\t100\t200\t10\n')

    data_files = {'chip_sample.bed': {'file': fake_file}}

    peaks = ri.peak_calling(data_files)
    assert isinstance(peaks, list)
    assert len(peaks) > 0
    assert peaks[0]['file'] == 'chip_sample.bed'
    assert peaks[0]['length'] == 100

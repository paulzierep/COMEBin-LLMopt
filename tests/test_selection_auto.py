import logging
import os
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'COMEBin'))

import get_final_result  # noqa: E402


def _args(mode, tmp_path):
    return SimpleNamespace(final_selection=mode, checkm2_python=None, output_path=str(tmp_path))


def test_auto_without_checkm2_keeps_partition(monkeypatch, tmp_path):
    monkeypatch.setattr(get_final_result, 'checkm2_python', lambda explicit=None: None)
    out = get_final_result.checkm2_bin_selection(
        logging.getLogger('t'), _args('auto', tmp_path), tmp_path, [], [], 'base', 1)
    assert out is None


def test_explicit_checkm2_without_checkm2_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(get_final_result, 'checkm2_python', lambda explicit=None: None)
    with pytest.raises(RuntimeError, match='needs CheckM2'):
        get_final_result.checkm2_bin_selection(
            logging.getLogger('t'), _args('checkm2', tmp_path), tmp_path, [], [], 'base', 1)


def test_auto_with_failing_checkm2_keeps_partition(monkeypatch, tmp_path):
    monkeypatch.setattr(get_final_result, 'checkm2_python', lambda explicit=None: [sys.executable, '-c', 'raise SystemExit(3)', '--'])
    args = _args('auto', tmp_path)
    args.contig_file, args.selection_weight = 'x.fa', 4.0
    out = get_final_result.checkm2_bin_selection(logging.getLogger('t'), args, tmp_path, [], [], 'base', 1)
    assert out is None

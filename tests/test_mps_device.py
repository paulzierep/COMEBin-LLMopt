"""Apple MPS device option (upstream #24); only the error paths can run on Linux."""
import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from train_CLmodel import _resolve_training_device  # noqa: E402


def test_mps_without_apple_silicon_fails_clearly():
    if getattr(torch.backends, 'mps', None) is not None and torch.backends.mps.is_available():
        assert _resolve_training_device('mps').type == 'mps'
    else:
        with pytest.raises(RuntimeError, match='MPS'):
            _resolve_training_device('MPS')


def test_unknown_device_lists_mps():
    with pytest.raises(ValueError, match='mps'):
        _resolve_training_device('tpu')


def test_cpu_still_works():
    assert _resolve_training_device('cpu').type == 'cpu'

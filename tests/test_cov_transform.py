import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'COMEBin'))

from get_augfeature import _log1p_checked, _normalise_sample  # noqa: E402


def test_log_transform_compresses_range_before_normalization():
    import logging
    cov = np.array([0.0, 1.0, 10.0, 1000.0])
    plain = cov.copy()
    _normalise_sample(plain, 'max', 'x', logging.getLogger('t'))
    logged = cov.copy()
    _log1p_checked(logged, 'x')
    _normalise_sample(logged, 'max', 'x', logging.getLogger('t'))
    assert np.isclose(logged.max(), 1.0) and np.isclose(plain.max(), 1.0)
    assert logged[2] > 0.3 > plain[2]  # coverage 10 is no longer squashed to ~0.01
    assert np.all(np.diff(logged) > 0)


def test_log_transform_rejects_negative_values():
    with pytest.raises(ValueError):
        _log1p_checked(np.array([1.0, -1.0]), 'x')

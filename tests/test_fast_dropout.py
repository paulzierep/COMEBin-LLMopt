"""FastDropout must behave like nn.Dropout (identity in eval, unbiased in train)."""
import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from models.fast_dropout import FastDropout  # noqa: E402


def test_eval_is_identity():
    layer = FastDropout(0.2).eval()
    x = torch.randn(64, 32)
    assert torch.equal(layer(x), x)


@pytest.mark.parametrize('p', [0.1, 0.2, 0.5])
def test_train_drop_rate_and_scale(p):
    torch.manual_seed(0)
    layer = FastDropout(p).train()
    x = torch.ones(2000, 500)
    y = layer(x)
    kept = y != 0
    assert abs(1 - kept.float().mean().item() - p) < 0.005
    assert torch.allclose(y[kept], torch.full_like(y[kept], 1 / (1 - p)))
    assert abs(y.mean().item() - 1.0) < 0.01


def test_gradient_matches_mask():
    torch.manual_seed(1)
    layer = FastDropout(0.3).train()
    x = torch.randn(100, 20, requires_grad=True)
    y = layer(x)
    y.sum().backward()
    expected = (y != 0).float() / 0.7
    assert torch.allclose(x.grad, expected)


def test_invalid_p():
    with pytest.raises(ValueError):
        FastDropout(1.0)

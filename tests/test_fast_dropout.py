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
    keep_rate = 1 - round(p * 256) / 256  # CPU masks quantize p to 1/256
    assert abs(kept.float().mean().item() - keep_rate) < 0.005
    assert abs(1 - kept.float().mean().item() - p) < 0.005
    assert torch.allclose(y[kept], torch.full_like(y[kept], 1 / keep_rate))
    assert abs(y.mean().item() - 1.0) < 0.01


def test_gradient_matches_mask():
    torch.manual_seed(1)
    layer = FastDropout(0.3).train()
    x = torch.randn(100, 20, requires_grad=True)
    y = layer(x)
    y.sum().backward()
    expected = (y != 0).float() / (1 - round(0.3 * 256) / 256)
    assert torch.allclose(x.grad, expected)


def test_invalid_p():
    with pytest.raises(ValueError):
        FastDropout(1.0)


@pytest.mark.parametrize('shape', [(1,), (7,), (13, 3), (64, 33, 2)])
def test_odd_shapes(shape):
    torch.manual_seed(2)
    y = FastDropout(0.2).train()(torch.ones(shape))
    assert y.shape == shape


def test_seeded_masks_are_reproducible():
    layer = FastDropout(0.2).train()
    x = torch.ones(300, 70)
    torch.manual_seed(5)
    a = layer(x)
    torch.manual_seed(5)
    b = layer(x)
    assert torch.equal(a, b)


def test_mask_is_not_periodic():
    """Bytes of one draw must not repeat across neighbouring elements."""
    torch.manual_seed(3)
    keep = (FastDropout(0.5).train()(torch.ones(700000)) != 0).float()
    for lag in (1, 7, 8):
        corr = torch.corrcoef(torch.stack([keep[:-lag], keep[lag:]]))[0, 1]
        assert abs(corr.item()) < 0.01

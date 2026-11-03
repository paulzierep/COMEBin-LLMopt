"""The memory-efficient InfoNCE loss must equal the original SimCLR formulation."""
import sys
import types
from pathlib import Path

import pytest
import torch
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from simclr import SimCLR  # noqa: E402
from utils import accuracy  # noqa: E402


def reference_loss(features, batch_size, n_views, temperature):
    """Original COMEBin v1.1.0 implementation (materialised logits)."""
    sample_ids = torch.arange(batch_size).repeat(n_views)
    positive_mask = sample_ids.unsqueeze(0).eq(sample_ids.unsqueeze(1))
    features = F.normalize(features, dim=1)
    similarity_matrix = features @ features.T
    diagonal_mask = torch.eye(positive_mask.shape[0], dtype=torch.bool)
    positive_mask = positive_mask[~diagonal_mask].view(positive_mask.shape[0], -1)
    similarity_matrix = similarity_matrix[~diagonal_mask].view(similarity_matrix.shape[0], -1)
    positives = similarity_matrix[positive_mask].view(-1, 1)
    negatives = similarity_matrix[~positive_mask].view(similarity_matrix.shape[0], -1)
    negatives = negatives[:, None].expand(-1, n_views - 1, -1).flatten(0, 1)
    logits = torch.cat([positives, negatives], dim=1) / temperature
    labels = torch.zeros(logits.shape[0], dtype=torch.long)
    loss = torch.nn.functional.cross_entropy(logits, labels)
    top1, top5 = accuracy(logits, labels, topk=(1, min(5, logits.shape[1])))
    return loss, top1, top5


@pytest.mark.parametrize('batch_size,n_views,temperature', [(32, 6, 0.07), (7, 2, 0.15), (16, 3, 0.15)])
def test_loss_and_gradients_match(batch_size, n_views, temperature):
    torch.manual_seed(0)
    model = types.SimpleNamespace(args=types.SimpleNamespace(n_views=n_views, temperature=temperature))
    base = torch.randn(batch_size, 16)
    features = (base.repeat(n_views, 1) + 0.8 * torch.randn(batch_size * n_views, 16)).double()

    f1 = features.clone().requires_grad_(True)
    loss, top1, top5 = SimCLR.info_nce_loss(model, f1, with_accuracy=True)
    loss.backward()

    f2 = features.clone().requires_grad_(True)
    ref_loss, ref_top1, ref_top5 = reference_loss(f2, batch_size, n_views, temperature)
    ref_loss.backward()

    assert torch.allclose(loss, ref_loss, rtol=1e-10, atol=1e-12)
    assert torch.allclose(f1.grad, f2.grad, rtol=1e-8, atol=1e-12)
    assert top1.item() == pytest.approx(ref_top1.item())
    assert top5.item() == pytest.approx(ref_top5.item())

"""InMemoryBatchLoader must yield exactly the batches torch's DataLoader yields."""
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin'))
from train_CLmodel import InMemoryBatchLoader  # noqa: E402


def test_matches_seeded_dataloader():
    views = [torch.randn(103, 7) + i for i in range(6)]
    reference = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(*views), batch_size=16, shuffle=True,
        num_workers=0, drop_last=True, generator=torch.Generator().manual_seed(5))
    loader = InMemoryBatchLoader(views, batch_size=16, generator=torch.Generator().manual_seed(5))
    assert len(loader) == len(reference) == 6
    for _ in range(3):  # several epochs: the generator state must advance identically
        for expected, actual in zip(reference, loader, strict=True):
            for e, a in zip(expected, actual, strict=True):
                assert torch.equal(e, a)


def test_unseeded_covers_each_sample_once():
    views = [torch.arange(40).unsqueeze(1).float()]
    loader = InMemoryBatchLoader(views, batch_size=8)
    seen = torch.cat([batch[0] for batch in loader]).flatten()
    assert sorted(seen.tolist()) == list(range(40))

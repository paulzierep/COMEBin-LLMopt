"""The minimum-length prefilter keeps exactly the records >= MIN_LENGTH."""
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'COMEBin' / 'scripts'))
from filter_min_length import filter_fasta  # noqa: E402

FASTA = '>a desc\nACGT\nAC\n>b\nACGTACGTAC\n>c\nA\n>d\nACGTA'


@pytest.mark.parametrize('compress', [False, True])
def test_filter(tmp_path, compress):
    src = tmp_path / ('in.fa.gz' if compress else 'in.fa')
    if compress:
        with gzip.open(src, 'wt') as handle:
            handle.write(FASTA)
    else:
        src.write_text(FASTA)
    out = tmp_path / 'out.fa'
    assert filter_fasta(src, out, 5) == (3, 4)
    assert out.read_text() == '>a desc\nACGT\nAC\n>b\nACGTACGTAC\n>d\nACGTA\n'


def test_nothing_kept(tmp_path):
    src = tmp_path / 'in.fa'
    src.write_text(FASTA)
    with pytest.raises(SystemExit):
        filter_fasta(src, tmp_path / 'out.fa', 100)
    assert not (tmp_path / 'out.fa').exists()

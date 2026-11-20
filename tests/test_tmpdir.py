import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'COMEBin', 'scripts'))
sys.path.insert(0, os.path.join(HERE, '..', 'COMEBin'))

unitem_profile = pytest.importorskip('unitem_profile')
RUN = os.path.join(HERE, '..', 'COMEBin', 'run_comebin.sh')


def test_checkm_gets_tmpdir(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setenv('TMPDIR', str(tmp_path))
    monkeypatch.setattr(unitem_profile.subprocess, 'run',
                        lambda cmd, check: seen.append(cmd) or subprocess.CompletedProcess(cmd, 0))
    unitem_profile._run_checkm(['checkm', 'analyze', '-t', '2', 'ms', 'bins', 'out'])
    assert seen[0][:4] == ['checkm', 'analyze', '--tmpdir', str(tmp_path)]


def test_checkm_without_tmpdir(monkeypatch):
    monkeypatch.delenv('TMPDIR', raising=False)
    assert unitem_profile._tmpdir_args() == []


def _run(args, env):
    # run_comebin.sh calls `python`: use the interpreter running the tests.
    env = dict(env, PATH=os.path.dirname(sys.executable) + os.pathsep + env.get('PATH', ''))
    return subprocess.run(['bash', RUN] + args, capture_output=True, text=True, env=env)


def test_T_creates_and_reports_dir(tmp_path):
    target = tmp_path / 'scratch' / 'tmp'
    # No -a/-o/-p: the script stops at the parameter check, after setting TMPDIR.
    out = _run(['-T', str(target)], dict(os.environ))
    assert target.is_dir()
    assert f'Temporary files: {target}' in out.stdout


def test_unwritable_tmpdir_fails(tmp_path):
    env = dict(os.environ, TMPDIR=str(tmp_path / 'missing'))
    out = _run([], env)
    assert out.returncode == 1
    assert 'not writable' in out.stderr


def test_gzip_assembly_is_decompressed(tmp_path):
    import gzip
    asm = tmp_path / 'asm.fasta.gz'
    with gzip.open(asm, 'wt') as handle:
        handle.write('>c1\n' + 'A' * 3000 + '\n>c2\n' + 'C' * 1500 + '\n')
    (tmp_path / 'bams').mkdir()
    out = tmp_path / 'out'
    # No BAMs: the run stops in data augmentation, after decompressing and filtering.
    result = _run(['-a', str(asm), '-o', str(out), '-p', str(tmp_path / 'bams'),
                   '-S', 'partition', '-L', '2500', '-d', 'cpu'], dict(os.environ))
    assert result.returncode != 0
    plain = out / 'input_asm.fasta'
    assert plain.read_text().count('>') == 2
    assert (out / 'contigs_min2500.fa').read_text().count('>') == 1


@pytest.mark.skipif(shutil.which('busybox') is None, reason='busybox not installed')
def test_runs_with_busybox_coreutils(tmp_path):
    """BioContainers use BusyBox: its realpath/od/head accept fewer options than GNU coreutils."""
    shim = tmp_path / 'shim'
    shim.mkdir()
    for tool in ['realpath', 'od', 'head', 'tr', 'dirname', 'basename', 'mkdir', 'mv', 'readlink']:
        (shim / tool).symlink_to(shutil.which('busybox'))
    env = dict(os.environ, PATH=str(shim) + os.pathsep + os.environ.get('PATH', ''))
    asm = tmp_path / 'asm.fasta'
    asm.write_text('>c1\n' + 'A' * 3000 + '\n')
    (tmp_path / 'bams').mkdir()
    result = _run(['-a', str(asm), '-o', str(tmp_path / 'out'), '-p', str(tmp_path / 'bams'),
                   '-T', str(tmp_path / 't'), '-S', 'partition', '-d', 'cpu'], env)
    assert 'No such file or directory' not in result.stderr, result.stderr
    assert 'Running data augmentation' in result.stdout

import os
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

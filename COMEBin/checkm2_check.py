"""Check before training whether CheckM2 can run the final-bin selection (run_comebin.sh -S).

usage: checkm2_check.py [--db PATH] [--checkm2_python PYTHON]
Prints one line (what was found, or what is missing) and exits 0 when CheckM2 and its DIAMOND
database are usable, 1 otherwise.
"""
import argparse
import subprocess
import sys
from pathlib import Path

from get_final_result import checkm2_python

DB_PROBE = 'from checkm2.fileManager import DiamondDB; print(DiamondDB().get_DB_location())'


def checkm2_status(db=None, explicit_python=None):
    """Return (usable, message)."""
    python = checkm2_python(explicit_python)
    if python is None:
        return False, 'CheckM2 was not found (no checkm2 on PATH, and this Python cannot import it)'
    if not db:
        probe = subprocess.run(python + ['-c', DB_PROBE], capture_output=True, text=True)
        lines = probe.stdout.strip().splitlines() if probe.returncode == 0 else []
        db = lines[-1].strip() if lines else ''
        if not db or db == 'None':
            return False, ('CheckM2 was found, but no CheckM2 database is configured '
                           '(pass -D, or run "checkm2 database --download")')
    if not Path(db).is_file():
        return False, f'the CheckM2 database does not exist: {db}'
    return True, f'CheckM2 ({python[0]}) with database {db}'


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--db')
    p.add_argument('--checkm2_python')
    a = p.parse_args(argv)
    ok, message = checkm2_status(a.db, a.checkm2_python)
    print(message)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

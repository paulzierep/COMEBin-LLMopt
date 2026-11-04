#!/usr/bin/env python
"""Stream a FASTA (optionally gzipped) and keep records of at least MIN_LENGTH bp.

Headers and sequence lines are copied unchanged. The output is written to a
temporary file and renamed, so an interrupted run never leaves a partial file.

usage: filter_min_length.py INPUT OUTPUT MIN_LENGTH
"""
import gzip
import os
import sys


def _open(path):
    with open(path, 'rb') as handle:
        magic = handle.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'rt', encoding='utf-8')


def filter_fasta(input_path, output_path, min_length):
    kept = total = 0
    tmp_path = f'{output_path}.tmp{os.getpid()}'

    with _open(input_path) as src, open(tmp_path, 'w', encoding='utf-8') as dst:
        header, lines, length = None, [], 0

        def flush():
            nonlocal kept, total
            if header is None:
                return
            total += 1
            if length >= min_length:
                kept += 1
                dst.write(header)
                dst.writelines(lines)

        for line in src:
            if line.startswith('>'):
                flush()
                header, lines, length = line, [], 0
            else:
                if not line.endswith('\n'):
                    line += '\n'
                lines.append(line)
                length += len(line.strip())
        flush()

    if kept == 0:
        os.unlink(tmp_path)
        raise SystemExit(f'No contigs of at least {min_length} bp in {input_path}.')
    os.replace(tmp_path, output_path)
    return kept, total


def main(argv):
    if len(argv) != 4:
        raise SystemExit(__doc__)
    kept, total = filter_fasta(argv[1], argv[2], int(argv[3]))
    print(f'Kept {kept} of {total} contigs with length >= {argv[3]} bp.')


if __name__ == '__main__':
    main(sys.argv)

#!/usr/bin/env python3
"""Regenerate the synthetic test fixtures (deterministic; no real data).

Run from the repository root:  python3 test/fixtures/make-fixtures.py
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(20261005)


def qual(n, lo=2, hi=40):
    return "".join(chr(33 + rng.randint(lo, hi)) for _ in range(n))


def read(i, mate, n=60):
    seq = "".join(rng.choice("ACGT") for _ in range(n))
    if i == 3:  # a no-call and a stretch of low quality
        seq = seq[:10] + "N" + seq[11:]
    q = qual(n, 2 if i in (3, 7) else 25, 40 if i != 7 else 15)
    head = f"@SYN01:7:FC123ABXX:1:1101:{1000 + i}:{2000 + i} {mate}:N:0:ACGTACGT+TTGCAACC"
    return f"{head}\n{seq}\n+\n{q}\n"


def write(name, text, gz=False):
    path = os.path.join(HERE, name)
    data = text.encode("ascii")
    if gz:
        # mtime=0 and no file name: identical bytes on every run
        with open(path, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
                f.write(data)
    else:
        with open(path, "wb") as f:
            f.write(data)


r1 = "".join(read(i, 1) for i in range(12))
rng = random.Random(20261006)
r2 = "".join(read(i, 2) for i in range(12))
write("sample_R1_001.fastq", r1)
write("sample_R2_001.fastq", r2)
write("sample_R1_001.fastq.gz", r1, gz=True)

# Illumina 1.3-1.7 Phred+64, pre-1.8 headers
old = "".join(
    f"@HWUSI-EAS100R:6:73:941:{1973 + i}#0/1\n{'ACGTN' * 8}\n+\n"
    + "".join(chr(64 + 10 + (j % 30)) for j in range(40)) + "\n"
    for i in range(4)
)
write("phred64.fastq", old)

# record 2 has a length mismatch
bad = read(0, 1) + "@SYN01:7:FC123ABXX:1:1101:1:1 1:N:0:A\nACGT\n+\nIII\n"
write("invalid.fastq", bad)

# no trailing newline
write("no-final-newline.fastq", r1.rstrip("\n"))

#!/usr/bin/env python3
"""Write examples/demo_R1_001.fastq.gz and demo_R2_001.fastq.gz.

Synthetic paired-end reads that look like a modern Illumina run (150 bp,
CASAVA 1.8 headers, NovaSeq-style binned qualities with a few weak tails
and no-calls).  Deterministic; contains no real data.
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
BINS = "#-9I"  # Q2, Q12, Q24, Q40 as written by binned instruments


def make(mate, seed, n=2000):
    rng = random.Random(seed)
    out = []
    for i in range(n):
        length = 150 if rng.random() > 0.1 else rng.randint(60, 149)
        seq = [rng.choice("ACGT") for _ in range(length)]
        qual = ["I" if rng.random() > 0.08 else rng.choice("9-") for _ in range(length)]
        if rng.random() < 0.15:  # a fading tail
            start = rng.randint(length - 40, length - 10)
            for j in range(start, length):
                qual[j] = rng.choice("#-9")
        if rng.random() < 0.05:
            j = rng.randrange(length)
            seq[j], qual[j] = "N", "#"
        tile, x, y = 1101 + i // 400, 1000 + (i * 37) % 30000, 1000 + (i * 53) % 9000
        out.append(f"@LH00999:12:22DEMOLT4:1:{tile}:{x}:{y} {mate}:N:0:ACGTGACTAG+TTAGCCATGA\n"
                   f"{''.join(seq)}\n+\n{''.join(qual)}\n")
    path = os.path.join(HERE, f"demo_R{mate}_001.fastq.gz")
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
            f.write("".join(out).encode("ascii"))


make(1, 1)
make(2, 2)

#!/usr/bin/env python3
"""Generate the synthetic FASTQ files under examples/ (no real data).

Deterministic: every file comes from a fixed seed, and gzip members are
written with mtime 0, so re-running produces byte-identical output.

  demo_R{1,2}_001.fastq.gz         small paired set used by `make screenshots`
  samples/NA00001_S1_L001_R{1,2}_001.fastq.gz
                                   NovaSeq-style paired-end run: 151 bp,
                                   CASAVA 1.8 headers, 2-colour binned
                                   qualities that decay along the read,
                                   occasional no-calls, adapter read-through
                                   ending in poly-G on short inserts
  samples/lowqual_S2_L001_R1_001.fastq.gz
                                   a failing older-chemistry run: full-range
                                   Phred, steep 3' decay, many N calls
  samples/legacy_phred64.fastq     Illumina 1.5 era: pre-1.8 headers
                                   (@HWI-...#0/1) and Phred+64 qualities

Reads are sampled from one random "genome" so R1 and R2 are true mates
(R2 is the reverse complement of the other end of the same fragment).
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = os.path.join(HERE, "samples")

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # TruSeq read 1 read-through
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # TruSeq read 2 read-through
NOVASEQ_BINS = {2: "#", 12: "-", 23: "8", 37: "F"}  # Q -> Phred+33 char
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def genome(rng, length, gc=0.41):
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


def write_gz(path, records):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
            f.write("".join(records).encode("ascii"))


def record(header, seq, qual):
    return f"{header}\n{seq}\n+\n{qual}\n"


def novaseq_quality(rng, length, tail_start):
    """Binned qualities: mostly Q37, more Q23/Q12 towards the 3' end."""
    out = []
    for pos in range(length):
        drift = pos / length
        r = rng.random()
        if pos >= tail_start:
            q = 2
        elif r < 0.80 - 0.30 * drift:
            q = 37
        elif r < 0.95 - 0.10 * drift:
            q = 23
        else:
            q = 12
        out.append(q)
    return out


def call(rng, seq, quals, nchance):
    """Turn quals into chars and put an N wherever the call failed."""
    seq = list(seq)
    for i in range(len(seq)):
        if rng.random() < nchance or quals[i] == 2 and rng.random() < 0.3:
            seq[i], quals[i] = "N", 2
    return "".join(seq), "".join(NOVASEQ_BINS[q] for q in quals)


def mate_read(rng, fragment, adapter, read_len):
    """First READ_LEN bases of FRAGMENT, running into adapter then poly-G."""
    seq = (fragment + adapter + "G" * read_len)[:read_len]
    # 2-colour chemistry reads "no signal" as high-confidence G, so the
    # poly-G tail keeps decent qualities until the very end.
    quals = novaseq_quality(rng, read_len, read_len)
    if len(fragment) + len(adapter) < read_len:
        for i in range(len(fragment) + len(adapter), read_len):
            quals[i] = rng.choice([12, 23, 23, 37])
    return call(rng, seq, quals, 0.002)


def novaseq_pair(n, seed=11, read_len=151):
    rng = random.Random(seed)
    ref = genome(rng, 200_000)
    r1, r2 = [], []
    for i in range(n):
        insert = int(rng.gauss(320, 70))
        if rng.random() < 0.08:
            insert = rng.randint(60, read_len - 10)  # short insert: read-through
        insert = max(40, insert)
        start = rng.randrange(len(ref) - insert)
        fragment = ref[start:start + insert]
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        tile = 1101 + (i // 600) * 100 + (i // 150) % 4
        x, y = 1000 + (i * 7919) % 32000, 1000 + (i * 104729) % 30000
        name = f"@LH00123:47:22FJKLLT3:1:{tile}:{x}:{y}"
        index = "ACAGTGGTCA+TGCTCGATAC"
        filt = "Y" if rng.random() < 0.01 else "N"
        s1, q1 = mate_read(rng, fragment, ADAPTER_R1, read_len)
        s2, q2 = mate_read(rng, revcomp(fragment), ADAPTER_R2, read_len)
        r1.append(record(f"{name} 1:{filt}:0:{index}", s1, q1))
        r2.append(record(f"{name} 2:{filt}:0:{index}", s2, q2))
    return r1, r2


def lowqual(n, seed=23, read_len=101):
    """Older four-colour run going wrong: full Phred range, steep decay."""
    rng = random.Random(seed)
    ref = genome(rng, 50_000, gc=0.46)
    out = []
    for i in range(n):
        start = rng.randrange(len(ref) - read_len)
        seq = list(ref[start:start + read_len])
        qual = []
        for pos in range(read_len):
            mean = 36 - 30 * (pos / read_len) ** 1.6 - rng.choice([0, 0, 4, 10])
            q = max(2, min(41, int(rng.gauss(mean, 4))))
            if q <= 5 and rng.random() < 0.35:
                seq[pos], q = "N", 2
            qual.append(chr(q + 33))
        tile, x, y = 2101 + i // 250, 2000 + (i * 31) % 20000, 3000 + (i * 17) % 19000
        out.append(record(f"@M00777:212:000000000-BX3K9:1:{tile}:{x}:{y} 1:N:0:2",
                          "".join(seq), "".join(qual)))
    return out


def legacy_phred64(n, seed=37, read_len=76):
    """Illumina 1.5 (GA II) era: pre-1.8 headers, Phred+64, B = Q2 tails."""
    rng = random.Random(seed)
    ref = genome(rng, 20_000)
    out = []
    for i in range(n):
        start = rng.randrange(len(ref) - read_len)
        seq = ref[start:start + read_len]
        qual = []
        tail = read_len - rng.randint(0, 20) if rng.random() < 0.4 else read_len
        for pos in range(read_len):
            q = 2 if pos >= tail else max(3, min(40, int(rng.gauss(36 - 12 * pos / read_len, 3))))
            qual.append(chr(q + 64))
        x, y = 1000 + (i * 13) % 1900, 1000 + (i * 29) % 1900
        out.append(record(f"@HWI-EAS209_0006:6:{1 + i // 100}:{x}:{y}#0/1", seq, "".join(qual)))
    return out


def demo_pair():
    """The small set `make screenshots` renders (kept stable on purpose)."""
    for mate, seed in ((1, 1), (2, 2)):
        rng = random.Random(seed)
        out = []
        for i in range(2000):
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
            out.append(record(f"@LH00999:12:22DEMOLT4:1:{tile}:{x}:{y} {mate}:N:0:ACGTGACTAG+TTAGCCATGA",
                              "".join(seq), "".join(qual)))
        write_gz(os.path.join(HERE, f"demo_R{mate}_001.fastq.gz"), out)


def main():
    os.makedirs(SAMPLES, exist_ok=True)
    demo_pair()
    r1, r2 = novaseq_pair(3000)
    write_gz(os.path.join(SAMPLES, "NA00001_S1_L001_R1_001.fastq.gz"), r1)
    write_gz(os.path.join(SAMPLES, "NA00001_S1_L001_R2_001.fastq.gz"), r2)
    write_gz(os.path.join(SAMPLES, "lowqual_S2_L001_R1_001.fastq.gz"), lowqual(1500))
    with open(os.path.join(SAMPLES, "legacy_phred64.fastq"), "w", newline="\n") as f:
        f.write("".join(legacy_phred64(400)))


if __name__ == "__main__":
    main()

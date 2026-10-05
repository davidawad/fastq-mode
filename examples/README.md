# Sample FASTQ files

Synthetic reads for trying fastq-mode and for the screenshots in
[`docs/README.md`](../docs/README.md). Nothing here is real sequencing data.
All of it is written by [`make-demo.py`](make-demo.py) from fixed seeds, so
re-running it (`python3 examples/make-demo.py`, or `make fixtures`) produces
byte-identical files.

| File | What it shows |
|---|---|
| `samples/NA00001_S1_L001_R{1,2}_001.fastq.gz` | A NovaSeq-style paired-end run: 3,000 pairs of 151 bp reads, CASAVA 1.8 headers with dual indexes, two-colour binned qualities (Q2/12/23/37) that decay along the read, occasional N calls, and ~8% short inserts that read through the TruSeq adapter into a poly-G tail. R2 is the reverse complement of the other end of the same fragment, so the pairs are true mates. |
| `samples/lowqual_S2_L001_R1_001.fastq.gz` | A failing older-chemistry run: 1,500 reads of 101 bp with full-range Phred scores, a steep 3' quality drop and many N calls. |
| `samples/legacy_phred64.fastq` | An Illumina 1.5-era file: 400 reads of 76 bp, pre-1.8 `@HWI-...#0/1` headers and Phred+64 qualities with `B` (Q2) tails. Plain text, to show uncompressed files too. |
| `demo_R{1,2}_001.fastq.gz` | The small set `make screenshots` renders headlessly. |

Open any of them with `M-x fastq-view`, or just visit the file once
`fastq-auto-view-mode` is on.

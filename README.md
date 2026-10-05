# fastq-mode

An Emacs major mode and streaming viewer for FASTQ sequencing reads —
the raw output of DNA sequencers. Open a 3 GB `.fastq.gz` in a few
milliseconds: the file is never loaded, only the page of reads on screen.

[![fastq-mode viewing a NovaSeq-style paired-end FASTQ file](docs/screenshots/neodoom-viewer.png)](docs/screenshots/neodoom-viewer.png)

Bases coloured by nucleotide, qualities by Phred bin, weak bases (below
Q20) underlined, and eldoc giving the Phred score and error probability of
the base at point. The read with the long `GGGG...` tail ran through its
adapter into poly-G.

More screenshots (failing runs, header decoding, R1/R2 mates, statistics,
Phred+64): [docs/README.md](docs/README.md). All of them use synthetic
sample data from [`examples/`](examples/).

## Why

As of 2026-10 there is no FASTQ mode for Emacs — see [Prior art](#prior-art).
Emacs alone can't cope with these files either: `find-file` on a
`.fastq.gz` decompresses the whole thing into memory (tens of gigabytes
for one whole-genome lane), and plain font-lock can't tell a quality line
that starts with `@` from a header.

## What it does

- **Streams files of any size**, plain or gzipped. Paging reads one page.
  Jumping to read N in a plain file seeks from a remembered byte offset
  (kept every 10,000 reads). A gzip file can't seek, so the jump
  decompresses up to N: about 1.3 s for read 1,000,000 on a NovaSeq X lane.
- **Highlighting by record role**: header (with the comment dimmed), bases
  (A/C/G/T/N), the `+` line, and qualities in four Phred bins
  (`fastq-quality-bins`, default 10/20/30). Bases whose call is below
  `fastq-low-quality-threshold` (Q20) are marked in the sequence itself,
  read from the quality line under it.
- **Eldoc**: on a base or quality character, `Read 2, base 17/150: G  Q12
  (error 6.31%, poor)`; on a header, the decoded fields.
- **Header decoding** for Illumina CASAVA 1.8+ (instrument, run, flowcell,
  lane, tile, x/y, UMI, read, filter flag, index) and the pre-1.8 form;
  `h` shows every field.
- **Paired-end mates**: `m` opens the R2 file for an R1 (and back) at the
  same read number, and warns if the read names differ.
- **Phred+33 and +64**, detected per file (`fastq-phred-offset`).
- **Statistics** (`s`): reads, bases, length distribution, GC%, N%, mean Q,
  % bases at or above Q20 and Q30. When [seqkit](https://bioinf.shenwei.me/seqkit/)
  is installed, an exact whole-file `seqkit stats -a` summary is added
  asynchronously.
- **Validation** (`v`): first record with a broken header, separator, or a
  sequence/quality length mismatch.

## Install

Requires Emacs 29.1+. Optional: `gzip` (any OS, for `.gz` files of any
size; without it, `.gz` files under 64 MB use Emacs' built-in zlib) and
`seqkit`.

```elisp
(use-package fastq-mode
  :straight (:host gitlab :repo "davidawad/fastq-mode")
  :config
  ;; open .fastq.gz and large .fastq files in the streaming viewer
  (fastq-auto-view-mode 1))
```

`.fastq` and `.fq` files open in `fastq-mode` automatically. Without
`fastq-auto-view-mode`, use `M-x fastq-open` (or `M-x fastq-view`) for
compressed and large files.

## Use

`M-x fastq-open RET reads_R1_001.fastq.gz`. In the viewer:

| Key | Action |
|---|---|
| `n` / `p` | next / previous read (crosses pages) |
| `SPC` / `DEL` | next / previous page (`fastq-view-page-records`, 200) |
| `g` | go to read number |
| `<` / `>` | first / last page (`>` counts the reads once) |
| `#` | count reads |
| `m` | mate file (R1 ↔ R2) at this read |
| `s` | statistics (`C-u s`: whole file) |
| `h` | decode the header of this read |
| `v` | validate the reads on this page |
| `w` | copy this read |
| `q` | quit (`M-x revert-buffer` reloads the page) |

In `fastq-mode` (small, editable files) the same commands are on `C-c C-n`,
`C-c C-p`, `C-c C-g`, `C-c C-m`, `C-c C-s`, `C-c C-h` and `C-c C-v`.

## Limits

- Four-line FASTQ only (what every current instrument and converter
  writes). Wrapped, multi-line FASTQ is reported by `fastq-validate`.
- gzip has no random access, so in a `.gz` file a jump to read N, and `>`,
  cost a decompression scan of everything before it. Plain files seek.
- Statistics in Emacs Lisp process roughly 60,000 reads a second
  (byte-compiled), so the default is a 200,000-read sample. For exact
  whole-file numbers on large files, install seqkit.

## Platforms

Linux, macOS and native Windows. External programs are found with
`executable-find` and run with argument lists (no shell). CI:
`.github/workflows/test.yml` is a fast Linux check on every push;
`ci-full.yml` runs the full Linux/macOS/Windows × Emacs 29.1/30.1 matrix
on demand, on tags and weekly.

## Development

```sh
make test       # ERT suite
make compile    # byte-compile, warnings are errors
make checkdoc
make fixtures   # regenerate the synthetic fixtures and demo files
make screenshots
```

Without make: `emacs -Q --batch -l test/run-tests.el [test|compile|checkdoc|all]`.

## Prior art

Checked 2026-10-05:

- **Package archives**: no package names or describes FASTQ in
  [MELPA](https://melpa.org/archive.json) (6,331 packages),
  [GNU ELPA](https://elpa.gnu.org/packages/archive-contents) or
  [NonGNU ELPA](https://elpa.nongnu.org/nongnu/archive-contents).
- **GitHub**: no repository named `fastq-mode`. A code search for `fastq`
  in Emacs Lisp finds only user configs. One config maps `.fastq` to a
  `fastq-mode` from an unpublished private file.
- **EmacsWiki** blocks scripted search, and a web search for an Emacs FASTQ
  mode returns nothing.
- **Nearest Emacs packages**, all FASTA-only:
  - [sequed](https://github.com/brannala/sequed)
  - [fasta.el](https://github.com/davep/fasta.el)
  - [SEQEL](https://github.com/rnaer/seqel)
  - [emacs-fasta-mode](https://github.com/vaiteaopuu/emacs-fasta-mode): FASTQ support was on its to-do list in 2017 and never written.
  - [bioseq-mode](https://github.com/mnbram/bioseq-mode)
  - [dna-mode](https://github.com/mikpom/dna-mode)
- **Other editors**: [bioSyntax](https://github.com/bioSyntax/bioSyntax)
  provides FASTQ syntax colouring for Vim, Sublime, gedit and `less`.

Format references:

- Cock et al. 2010, [The Sanger FASTQ file format for sequences with
  quality scores, and the Solexa/Illumina FASTQ
  variants](https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/), Nucleic
  Acids Research 38(6).
- [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) on Wikipedia,
  which covers the Illumina header layouts.

## Related

- [genetics.el](https://gitlab.com/davidawad/genetics-el) (variant-level
  data) and [genome-cli](https://gitlab.com/davidawad/genome-cli) (its
  `genome pipeline` aligns FASTQ reads into a VCF).

## Licence

MIT.

## utag: barbell demux for custom CDI/UDI amplicon libraries (ONT)

Demultiplexes combinatorial dual-index (CDI) or unique dual-index (UDI) amplicon libraries sequenced on Oxford Nanopore with
[barbell](https://github.com/rickbeeloo/barbell). For CDI, each amplicon carries a utag + primer at both ends; the **pair** of utags
(e.g. a 32 x 12 grid = 384 wells) identifies the sample. Works for NBD (native barcoding) and LSK (ligation) libraries and
any primer pair (default 16S 27F / 1492R, pass `--primer-f <SEQ> --primer-r <SEQ>` for other PCR primers). See below for UDI with the MAB kit.

### Quick start

**Requirements:** Python >= 3.8 and `barbell` on your PATH
(`RUSTFLAGS="-C target-cpu=native" cargo install barbell` or `conda install -c bioconda barbell`). Optional, to call it as `utag`:

```bash
chmod +x utag.py && ln -s "$(pwd)/utag.py" /usr/local/bin/utag
```

```bash
# all 384 wells, output named Plate1-A1.fastq ...
utag demux -r 'fastq_pass/**/*.fastq.gz' -c nbd -o demux_out -t 16

# your own sample names / a subset of wells
utag demux -r 'fastq_pass/**/*.fastq.gz' -c nbd --samplesheet samples.csv -o demux_out -t 16

# ONT barcodes also split the pool: ONT barcode first, then utags inside each bin
utag twopass -r 'fastq_pass/**/*.fastq.gz' -k SQK-NBD114-96 --samplesheet samples.csv -o twopass_out -t 16

# MAB unique dual indexes: 24 AB barcodes x 12 utags = 288 samples (see "MAB: unique dual indexes")
# (-k SQK-MAB* defaults to --plates data/mab-udi-platemap.csv --tags data/utags-10bp-sequence.csv)
utag twopass -r 'fastq_pass/**/*.fastq.gz' -k SQK-MAB114-24 --samplesheet samples.csv -o mab_out -t 16
```

`-r` takes files, globs or directories (searched recursively for
`.fastq/.fq[.gz]`). `-h` on any subcommand lists its options and the input file formats.

### Choosing the command and context
Diagrams of amplicon ONT libraries, utag contexts, and use cases: https://angelovangel.github.io/utag/

```
5'-[ONT context][utag-fwd][forward primer] ...insert... [reverse primer][utag-rev][ONT context]-3'
```

The "context" is the constant sequence directly 5' of the utag, which barbell uses as part of the flank. `demux` defaults to
`bare`; for untrimmed NBD (reads in `fastq_pass/barcode01`) pass `-c nbd --ont-barcode NB01`, for LSK reads pass `-c lsk`:

| `-c` | Use for | Sequence 5' of the utag |
|---|---|---|
| `nbd` | SQK-NBD114 | `CAGCACCT` (after every NB barcode); with `--ont-barcode NB07`: `ATTGCTAAGGTTAA` + NB07 + `CAGCACCT` |
| `lsk` | SQK-LSK114, no ONT barcode | `CTTCGTTCAGTTACGTATTGCT` (end of the Ligation Adapter) |
| `mab` | SQK-MAB114-24 reads not split by `barbell kit` (`twopass` uses `bare` in pass 2) | `GCCCTT` (end of the TOPO site after the AB barcode) |
| `bare` (default) | Reads whose ONT adapter/barcode is already trimmed (`twopass` uses this in pass 2) | nothing |

* **`demux`**: utags are the only sample identifier (no ONT barcodes, or after demux on ONT barcodes with e.g. dorado or barbell).
* **`twopass`**: several ONT barcodes in one run (e.g. NB01-NB04). Pass 1 `barbell kit` splits by ONT barcode and trims it;
  pass 2 runs `demux -c bare` inside each bin, so the same utag pair can be reused in different bins.
* Don't use `-c nbd` for LSK libraries. MAB kits can't use a tag *pair* (the utag at the
  read end is low quality/trimmed); use them as unique dual indexes, see below.
* If unsure what is really upstream of your utag, run `utag context` (see Troubleshooting) and pass its suggestion as `--prefix`.

### NBD: combinatorial dual indexes

The native barcode is ligated to both amplicon ends, so every read carries a utag at each end and the **pair** (fwd + rev)
names the well. Reads come from either strand (`Ftag ... Rtag(rc)` or `Rtag ... Ftag(rc)`):

```
                      ┌─── twopass pass 1: barbell kit ───┐┌─────────────┐ pass 2: -c bare
5'-[native adapter ..][ATTGCTAAGGTTAA][NB01-NB96][CAGCACCT][utag-fwd][27F] ...insert... -3'
                      └──── -c nbd --ont-barcode NBxx ────┘└─────────────┘ utag + primer
                                                 └ -c nbd ┘
```

The far read end is the same layout, reverse complemented, and is searched the same way.

`CAGCACCT` follows every NB barcode, so `demux -c nbd` (no `--ont-barcode`) searches `CAGCACCT` + utag + primer in reads of
any NB barcode. With one NB barcode per library, `--ont-barcode NB07` adds `ATTGCTAAGGTTAA` + NB07 to the flank (longer, more
specific). With several NB barcodes in a run, `twopass -k SQK-NBD114-96` splits by NB barcode first, then searches utag +
primer (`bare`) inside each bin, so the same utag pair can be reused across NB barcodes.

### LSK: combinatorial dual indexes

No ONT barcode: the Ligation Adapter is ligated straight onto the dA-tailed amplicon, at
both ends. As for NBD, the utag **pair** determines the well and reads can come from either strand:

```
   ┌────── Ligation Adapter ──────┐
5'-[CCTGTA][CTTCGTTCAGTTACGTATTGCT][utag-fwd][27F] ...insert... -3'
           └ -c lsk (last 22 nt) ─┘└─────────────┘ utag + primer
```

The far read end is the same layout, reverse complemented, and is searched the same way.

Read starts are often miscalled, so `-c lsk` uses only the last 22 nt of the adapter. Don't use `-c nbd` here.

### MAB: unique dual indexes

With SQK-MAB114-24 the AB barcode attaches by chance to either amplicon end, and the read always starts at it, so a read starts
with the fwd **or** the rev utag; with unique dual indexes **either** tag determines the sample. The sample is **AB barcode + the utag
at the read start**; the far read end (low quality/truncated) is not used. Here `--single-end` is used - one utag at the read start (fwd or rev) determines which sample.
The utag follows the TOPO site directly:

```
                     ┌── twopass pass 1: barbell kit ──┐        ┌── pass 2: -c bare --single-end ──┐
5'-[rapid adapter ..][TTTAACC][AB01-AB24][CCATATCCGTGTC][GCCCTT][utag-fwd or utag-rev][27F or 1492R] ...insert... -3'
                                                        └──────┘ -c mab
                                                                └───────── utag + primer ──────────┘
```

The AB barcode is attached by TOPO (`CCCTT` is the vaccinia TOPO site). Pass 1 uses barbell's built-in MAB kit (flanks
`TTTAACC` / `CCATATCCGTGTC`); its trim leaves `GCCCTT` and sometimes cuts into it, so pass 2 searches utag + primer with no
prefix (`bare`). `utag twopass -k SQK-MAB*` switches on `--single-end` in pass2 by itself. For untrimmed reads, `demux -c mab` puts `GCCCTT`
before the utag.

* **Plate map** `data/mab-udi-platemap.csv` (`plate,well,ont-barcode,utag-fwd,utag-rev`): 3 plates; rows A-H = AB01-AB08 on
  plate 1, AB09-AB16 on plate 2, AB17-AB24 on plate 3; column 1-12 = UDI 1-12 = utag-fwd `U001`-`U012` + utag-rev
  `U017 U018 U021 U023 U025 U028 U030 U033 U034 U040 U042 U043` from `data/utags-10bp-sequence.csv`. So plate 2 well C5 = AB11 +
  U005/U025. The 24 tags are >= 4 edits apart, fwd and rev together (both can be at the read start). Edit the rev column to
  match your oligos.
* **Few reads in AB bins?** MAB read starts are noisy, and pass 1 (`barbell kit`) allows only 3 errors in its 44 nt
  flank + AB barcode by default. On real reads `--kit-flank-max-errors 4` raised the binned reads from 52-84% to 70-91% per AB
  barcode with no reads in a wrong bin (5 binned more but put 0.2-0.4% in a wrong bin). If the reads are already split into
  `barcodeNN` folders by MinKNOW/Dorado (adapters not trimmed), you can skip pass 1 and run `demux -c mab --single-end` per
  folder: 81% assigned vs 46% with the default `twopass` on the same AB19 reads.
* **Sample sheet**: `sample,plate,well` as usual; the AB barcode comes from the plate map, so no `ont-barcode` column is needed.
* `--single-end` needs every tag (fwd and rev) to name only one sample per AB bin. Reads keep their sequenced orientation; the
  optional `--flip` (off by default) turns reads that start with the rev tag into forward-primer orientation.
  

Other subcommands: `queries` (write the query FASTAs only), `tagqc` (check tags, plate map and sample sheet), `context`.

### Input files

| Option | Default | Format | Example row |
|---|---|---|---|
| `--tags` | `data/utags-sequence.csv` | `utag,sequence` (header optional) | `U013,CAACAGTAG` |
| `--plates` | `data/utags-platemap.csv` | `plate,well,utag-fwd,utag-rev` | `1,A1,U013,U001` |
| `--samplesheet` | none | `sample,plate,well[,ont-barcode]`, csv or `.xlsx` | `Soil_A,1,A1` |

The tag file and the plate map describe the oligos you have, not the run. The tag file is just the sequence of each tag
ID. The plate map says which two tags (fwd and rev) were put into each well of the oligo plates, so it is the layout of the
whole library: 4 plates, 384 wells. Both come with the repo and you only change them if your oligos or their layout are
different.

The sample sheet is the part that changes from run to run. It says which of your samples went into which well, for example
`Soil_A,1,A1`. The script looks up plate 1 well A1 in the plate map, finds `U013` + `U001`, takes their sequences from the tag
file, and writes every read with that tag pair to `samples/Soil_A.fastq`. If you already know the tags of a sample you can
skip the plate map and write them directly (`sample,utag-fwd,utag-rev`).

Without a sample sheet every well in the plate map is treated as a sample and the files are named after the well
(`Plate1-A1.fastq`). With a sample sheet only your samples are written. The plate map is still read in that case: the tags of
the wells you did not list are searched as decoys, so reads from those wells end up in `unexpected.tsv` instead of in one of
your samples. For samples given by their tags, `summary.tsv` also looks up their plate and well in the plate map.

**Tags** (`--tags`): 96 tags of 9 nt (`U001`-`U096`), 5'->3' as in the oligo. Only tags used by the plate map / sample sheet go into
the queries (`--all-tags` uses all). Every tag ID you reference must exist here. Check your setup with `tagqc`: tag edit distances, the plate map's layout per plate (which tags run along rows / columns),
duplicate wells, missing tag sequences, and with `--samplesheet` the samples per plate, ONT bins and decoys.
`tagqc --print columnwise|rowwise|plate` also prints every plate-map well with its tags: listed column by column (A1, B1, ..),
row by row (A1, A2, ..), or drawn as a plate grid.

**Plate map** (`--plates`): the full library design, 4 plates x 96 wells. Column 1-12 = rev tags `U001`-`U012`; rows A-H = fwd tags
`U013`-`U020` on plate 1, `U021`-`U028` on plate 2, `U029`-`U036` on plate 3, `U037`-`U044` on plate 4 (so plate 2 well C7 = `U023` + `U007`).
You don't edit or pass it unless your physical tag layout differs. A plate map may add an `ont-barcode` column (the ONT/AB
barcode of each well, as in `data/mab-udi-platemap.csv`); a sample sheet with `plate,well` then inherits it.

**Sample sheet** (`--samplesheet`): your sample names and their plates and wells.

```
sample,plate,well
Soil_A,1,A1
Water_2,1,H12
```

* `plate` may be written `1`, `Plate1` or `P1`; `well` as `A1` or `A01` (or `Plate1-A1` without a plate column).
* Output is `samples/<sample>.fastq`; spaces and special characters in names become `_`. Only listed samples are written.
* Instead of `plate,well` you can give the tags directly (`utag-fwd,utag-rev`, e.g. `U013,U001`), e.g. for tags without a plate layout.
* Headers are case-insensitive; column order and extra columns don't matter; blank rows are skipped.
  Excel: first worksheet, or `--sheet-name`; old `.xls` is not supported.
* Duplicate names are an error, as are two samples in the same well (except in different `ont-barcode` bins).
* `ont-barcode` is used by `twopass` only: each bin is searched for its own samples, and bins without samples are skipped.
  Without the column every bin is searched for every sample.
* **Decoys:** all plate-map tags are searched even if your sheet doesn't use them, so reads from unlisted wells are reported
  as unexpected instead of being forced onto the most similar listed tag. The script prints how many decoys it added, and warns
  if the plate map is missing. `--no-decoys` turns this off (more mis-assignments).

### Output

Results and the script's own files are at the top; everything barbell writes is in `barbell/`.

```
demux_out/
  samples/<sample>.fastq[.gz]   the result: one file per sample, utag and primer trimmed off
  summary.tsv                   reads per sample with plate, well and tags (0-read samples too); input for plate-heatmap.py
  unexpected.tsv                tag pairs found on reads but not in the plate map / sample sheet, with read counts
  run_summary.csv               totals (also printed): reads in, with a hit, assigned, tag at both / one end, decoys,
                                unexpected (in decoy wells / other tag pairs)
  utag.log               everything the script printed, with date and command line
  barbell/
    inspect.txt                 most common read patterns
    trimmed/                    per-tag-pair files (below)
    queries/                    the query FASTAs used
    annotation.tsv, filtered.tsv, dropped.tsv, trim_failed.txt, filters.txt, *.log
```

barbell names `trimmed/` files `<tag at read start>__<tag at read end>`, so each well gives two files, one per strand:

```
barbell/trimmed/U013__U001.trimmed.fastq   forward-strand reads
barbell/trimmed/U001__U013.trimmed.fastq   reverse-strand reads    -> both merged into samples/Plate1-A1.fastq
```

Pairs not in the plate map / sample sheet (e.g. `U001__U037`) are not merged; they are listed in `unexpected.tsv`.
Reads keep their sequenced orientation; add `--flip` to reverse-complement reverse-strand reads (forward primer -> reverse primer).
Reads are trimmed (ONT adapter/barcode, utags and primers removed); add `--skip-trim` to write whole, untrimmed reads to `samples/`
(twopass: pass 1 still removes the ONT adapter/barcode, `--skip-trim` keeps the utags and primers).
Once you have checked the results, `barbell/` can be deleted.

```
twopass_out/
  samples/<ONT bin>/<sample>.fastq
  summary.tsv, unexpected.tsv   as above, with an ont_barcode column
  ont_bins.tsv                  reads per ONT bin: in bin, with hit, assigned
  run_summary.csv, utag.log
  barbell/
    pass1_kit/, pass1_kit.log   barbell kit output (reads per ONT bin)
    pass2_utag/<bin>/           the demux run of each bin (its own log, summaries and barbell/ folder)
```

### Troubleshooting: few reads assigned

`demux` prints what was found per read (from `barbell/inspect.txt`):

| Line | Meaning |
|---|---|
| utag at both ends | tag pair read: only these reads can be assigned |
| utag at one end, primer-flank only at the other | the other end's tag is damaged or missing |
| utag at ONE end only | no tag pair possible |
| primer-flank only | primer found, utag unreadable |

* **Mostly one end only:** the far end has no utag in this prep, is truncated, or was trimmed by the basecaller (use untrimmed
  reads), or the primer/context at that end is wrong.
* **Mostly primer-flank only / few hits:** wrong primers or wrong context.
* **Other primers** (ITS, COI, 18S, ...): give both 5'->3' as ordered, IUPAC allowed: `--primer-f <SEQ> --primer-r <SEQ>`
  (`demux`, `twopass`, `queries`, `context`). The primer is part of the flank and is trimmed from the output.
* **Wrong context:** `context` finds utag + primer hits (primer within 250 nt of a read end with <= 2 edits, utag <= 2 edits
  right before it; `--exact` for error-free utag + 10 nt of primer) and aligns the sequence upstream (indels allowed) to the
  known adapters (LSK Ligation Adapter, NB rear, MAB TOPO site), separately for read starts and read ends. It reports how
  much of each adapter is there (% within k edits), which `-c` context fits, the NB/AB barcodes for NBD/MAB, and the command
  to use. A per-position consensus is printed too, for unknown kits. On real LSK reads the whole adapter is at the read start,
  but the read end stops inside it (only the last ~6 nt), which is due to how ONT works (there is some distance between the motor and the reader).

  ```bash
  utag context -r reads.fastq.gz                           # also writes context_report.html (--html FILE to rename, --html '' to skip)
  utag demux -r reads.fastq.gz -c <suggested context> ...   # or -c bare --prefix <seq> for unknown kits
  ```
* Tuning: `--flank-max-errors N` (relax if reads are missed), `--min-score-diff X` (raise, e.g. 0.2, for fewer wrong-well calls),
  `--max-offset 250` (max distance of the utag from the read end).

### Reference
**Barbell paper:** https://doi.org/10.1093/bioinformatics/btag349  
**Workflow:** `demux` runs barbell `annotate` (find utag + flank at both ends, both strands) -> `inspect` (pattern summary; the
two top patterns should be `Ftag ... Rtag(rc)` and `Rtag ... Ftag(rc)`) -> `filter` (keep reads with a tag at both ends) ->
`trim` (cut tags and primers, one file per tag pair) -> rename tag pairs to samples. It adapts to both the release barbell
(comma-separated `-q`/`-b`, single input) and newer builds, decompresses `.gz` input first, and raises `ulimit -n` for `trim`
(use `ulimit -n 65535` if you run barbell by hand).

**Query FASTAs.** One per end; each entry is the full oligo `context + utag + primer`, 5'->3'. Ready-made files for the default primers are in `data/fasta/`; regenerate with `queries --primer-f ... --primer-r ...`.

**This tag set.**
* Rev tags `U001`-`U012` all start with `A`, so barbell folds that base into the flank and compares 8 nt instead of 9 (printed; works).
* Minimum edit distance is 2 for rev tags (`U005`/`U011`), 3 for fwd tags: expect 0-2 wrong-well calls per ~1,400 assigned reads
  in simulation. Raise `--min-score-diff` if that matters.
* The 384-well grid is complete, so with the full plate map chimeras/index hopping cannot be detected; they land in a wrong well.
  With a sample sheet covering part of the grid, they show up in `unexpected.tsv`.

**ONT sequences.** The 96 NB barcodes and the NBD/LSK flanks match ONT's Chemistry Technical Document
(CHTD_500_v1_revAV_24July2026). You can check unknown reads with `context`.


### Repository

| File | Description |
|---|---|
| `utag.py` | The tool (`demux`, `twopass`, `context`, `queries`, `tagqc`) |
| `data/utags-sequence.csv`, `data/utags-platemap.csv` | Default tags and plate map |
| `data/samplesheet-example.csv`, `data/samplesheet-example-twopass.csv` | Example sample sheets (the second with `ont-barcode`) |
| `data/utags-10bp-sequence.csv`, `data/utags-10bp-platemap.csv` | An alternative 10 nt tag set (`U001`-`U044`) and its plate map (rows = rev, columns = fwd tags); use with `--tags` and `--plates` |
| `data/mab-udi-platemap.csv` | MAB unique dual index plate map: 24 AB barcodes x 12 utags (10 nt set), for `twopass -k SQK-MAB114-24` |
| `data/fasta/` | Ready-made query FASTAs for the default primers |
| `simulate-reads.py` | Fake reads with known truth: `python3 simulate-reads.py nbd 3000 sim.fastq` |
| `plate-heatmap.py` | HTML heatmap of reads per well: `python3 plate-heatmap.py demux_out/summary.tsv -o plates.html` (samples are placed by the `plate`/`well` columns, so custom names need no extra file). Mark negative-control wells (e.g. an unused plate) and wells to exclude in the page to get sensitivity / specificity at a read cutoff, the cutoff for 100% specificity and the utags that collect reads in the negatives |
| `compare-runs.py` | Compare runs from their `utag.log`: `python3 compare-runs.py run1 run2 -o compare_out` writes into the folder `compare_out/` (the default): `run_comparison.csv` (metrics x runs, differences to `--baseline`), `samples.csv` (reads per sample, from each `summary.tsv`) and `run_comparison.html` (charts, metrics table, top patterns, ONT bins, per-sample scatter, and negative controls: mark empty wells on the plate grids to get sensitivity / specificity per run) |

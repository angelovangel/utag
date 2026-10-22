#!/usr/bin/env python3
"""
utag.py - barbell query FASTAs + automated demux for custom utag amplicon
libraries on Oxford Nanopore: combinatorial (CDI) or unique (UDI) dual indexes.
Any primer pair works (--primer-f / --primer-r); defaults are 16S 27F / 1492R.

Library model (5'->3' of each tagged PCR primer):

    [ONT context] [utag] [primer]
      left  : ctx + utag-fwd + forward primer   (default 27F    {PRIMER_F})
      right : ctx + utag-rev + reverse primer   (default 1492R  {PRIMER_R})

Both queries are written in oligo (5'->3') orientation, which is what barbell
expects for Ftag/Rtag FASTAs; barbell finds them on either read strand.

Only needs Python >= 3.8 (stdlib) and the `barbell` binary.
"""

import argparse, csv, glob, gzip, importlib.util, itertools, json, os, random, re, resource, shutil
import subprocess, sys, time
from collections import Counter, defaultdict

# ---- defaults ------------------------------------------------------------------------------------------------------
PRIMER_F = "AGRGTTYGATYMTGGCTCAG"            # --primer-f  16S 27F
PRIMER_R = "RGYTACCTTGTTACGACTT"             # --primer-r  16S 1492R
TAGS = "data/utags-sequence.csv"             # --tags      9 nt utags (paths relative to this script's folder)
PLATES = "data/utags-platemap.csv"           # --plates    CDI: 32 fwd x 12 rev utags = 4 x 96 wells
MAB_TAGS = "data/utags-10bp-sequence.csv"    # --tags      for twopass -k SQK-MAB* (the MAB plate map's tag IDs)
MAB_PLATES = "data/mab-udi-platemap.csv"     # --plates    for twopass -k SQK-MAB*: 24 AB barcodes x 12 utags (UDI)

EXAMPLES = """\
examples:
  utag demux   -r 'fastq_pass/**/*.fastq.gz' -c nbd -o demux_out -t 16
  utag demux   -r fastq_pass/ -c nbd --samplesheet samples.csv
  utag twopass -r fastq_pass/ -k SQK-NBD114-96 --samplesheet samples.csv -o twopass_out
  utag context -r fastq_pass/                     # which -c fits these reads (+ context_report.html)
  utag tagqc   --samplesheet samples.csv          # check tags, plate map and sample sheet
"""

INPUT_HELP = f"""\
input file formats (defaults live in data/ next to this script):
  --tags         {TAGS:<26} utag,sequence
                                            U013,CAACAGTAG
  --plates       {PLATES:<26} plate,well,utag-fwd,utag-rev
                                            1,A1,U013,U001
                 the full 384-well design: 4 plates; column 1-12 -> U001-U012,
                 rows A-H -> U013-U020 on plate 1, U021-U028 on plate 2, ...
  --samplesheet  your .csv or .xlsx         sample,plate,well[,ont-barcode]
                                            Soil_A,1,A1[,NB01]
                 (or sample,utag-fwd,utag-rev to give the tags directly)

  --tags and --plates describe the oligos and stay the same; the sample sheet is
  the only file you write. Soil_A,1,A1 -> plate 1 well A1 = U013 + U001 (plate map)
  -> tag sequences (--tags) -> reads with that pair go to samples/Soil_A.fastq.

  No --samplesheet: every plate-map well is a sample (samples/Plate1-A1.fastq).
  With --samplesheet: only the listed samples are written (samples/Soil_A.fastq);
  the plate map's tags are still searched as decoys, so reads from unlisted wells
  are reported as unexpected (--no-decoys turns that off).
  ont-barcode (e.g. NB01) is only used by twopass. See README "Input files".
"""

ANSI = {
    "reset": "\033[0m", "bold": "\033[1m", "dim": "\033[2m",
    "red": "\033[31m", "green": "\033[32m", "yellow": "\033[33m", "blue": "\033[34m",
    "magenta": "\033[35m", "cyan": "\033[36m", "gray": "\033[90m",
}


def _color_default(stream):
    if os.environ.get("NO_COLOR") is not None:
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    return stream.isatty()


ENABLE_COLOR = _color_default(sys.stdout)
ERR_COLOR = _color_default(sys.stderr)


def set_color(mode):
    """--color auto|always|never"""
    global ENABLE_COLOR, ERR_COLOR
    if mode == "always":
        ENABLE_COLOR = ERR_COLOR = True
    elif mode == "never":
        ENABLE_COLOR = ERR_COLOR = False
    else:
        ENABLE_COLOR, ERR_COLOR = _color_default(sys.stdout), _color_default(sys.stderr)


UNI = (getattr(sys.stdout, "encoding", None) or "").lower().replace("-", "").startswith("utf")
SYM = dict(rule="\u2500" if UNI else "-", dot="\u2022" if UNI else "*", ok="\u2713" if UNI else "ok",
           warn="!", full="\u2588" if UNI else "#", empty="\u2591" if UNI else ".", arrow="\u2192" if UNI else "->",
           times="\u00d7" if UNI else "x", prompt="$")
_ANSI_RE = re.compile(r"\033\[[0-9;]*m")


def style(text, *codes, err=False):
    if not (ERR_COLOR if err else ENABLE_COLOR):
        return str(text)
    seq = "".join(ANSI[c] for c in codes if c in ANSI)
    return f"{seq}{text}{ANSI['reset']}" if seq else str(text)


def vlen(text):
    """visible length (ANSI escape codes do not count)"""
    return len(_ANSI_RE.sub("", str(text)))


def pad(text, width, align="<"):
    text = str(text)
    gap = max(0, width - vlen(text))
    return text + " " * gap if align == "<" else " " * gap + text if align == ">" else \
        " " * (gap // 2) + text + " " * (gap - gap // 2)


def cpad(text, width, align="<", *codes):
    """pad the PLAIN text first, then colour it: keeps columns aligned with colour on or off"""
    return style(pad(text, width, align), *codes)


def term_width():
    return max(64, min(shutil.get_terminal_size((100, 24)).columns, 110))


def fmt_int(n):
    return f"{int(n):,}"


def pct_color(p):
    return "green" if p >= 80 else "yellow" if p >= 50 else "red"


def bar(frac, width=22, code="green"):
    n = max(0, min(width, int(round(frac * width))))
    return style(SYM["full"] * n, code) + style(SYM["empty"] * (width - n), "gray")


def print_step(title):
    w = term_width()
    print("\n" + style(SYM["rule"] * 2, "gray") + " " + style(title, "bold", "cyan") + " "
          + style(SYM["rule"] * max(3, w - vlen(title) - 4), "gray"))


def _block(sym, msg, indent=2):
    lines = str(msg).split("\n")
    print(" " * indent + sym + " " + lines[0])
    for ln in lines[1:]:
        print(" " * (indent + 2) + ln.strip())


def print_info(msg):
    _block(style(SYM["dot"], "cyan"), msg)


def print_note(msg):
    _block(style(SYM["warn"], "yellow", "bold"), msg)


def print_ok(msg):
    _block(style(SYM["ok"], "green", "bold"), msg)


def print_error(msg, *, file=None):
    file = file or sys.stderr                # looked up at call time, so the log tee sees it
    print(f"{style('ERROR', 'bold', 'red', err=file is sys.stderr)}: {msg}", file=file)


def print_table(headers, rows, aligns=None, indent=4, gap=2, rule=True):
    """Aligned table. `rows` cells may contain ANSI colour; widths use visible length."""
    aligns = aligns or "<" * len(headers)
    cols = range(len(headers))
    widths = [max([vlen(headers[i])] + [vlen(r[i]) for r in rows]) for i in cols]
    sep = " " * gap
    n = len(widths)
    avail = term_width() - indent - sum(widths[:-1]) - gap * (n - 1)
    rule_w = list(widths)
    if aligns[-1] == "<" and widths[-1] > avail > 8:        # long free-text last column: do not run past the terminal
        rule_w[-1] = avail
    print(" " * indent + sep.join(cpad(headers[i], widths[i], aligns[i], "bold") for i in cols).rstrip())
    if rule:
        print(" " * indent + sep.join(style(SYM["rule"] * rule_w[i], "gray") for i in cols))
    for r in rows:
        print(" " * indent + sep.join(pad(r[i], widths[i], aligns[i]) for i in cols).rstrip())


def print_kv(rows, indent=4, gap=2, key_code="cyan"):
    """Aligned key / value list; value may be pre-coloured."""
    w = max(vlen(k) for k, _ in rows)
    for k, v in rows:
        print(" " * indent + cpad(k, w, "<", key_code) + " " * gap + str(v))


DISPLAY_BASE = [None]


def shorten(text):
    """paths relative to the output folder (and ~) for display only"""
    t = str(text)
    if DISPLAY_BASE[0]:
        t = t.replace(DISPLAY_BASE[0] + os.sep, "./")
    home = os.path.expanduser("~")
    return t.replace(home, "~") if home and home != "/" else t


def rel(path):
    try:
        r = os.path.relpath(path)
        return r if not r.startswith("..") else shorten(path)
    except ValueError:
        return shorten(path)


KIND_LABEL = {"tag_both_ends": "utag at both ends (can be assigned)",
              "tag_plus_flank_other_end": "utag at one end, primer-flank only at the other",
              "tag_one_end_only": "utag at ONE end only, nothing at the other",
              "flank_only": "primer-flank only, no readable utag"}
KIND_COLOR = {"tag_both_ends": "green", "tag_plus_flank_other_end": "yellow",
              "tag_one_end_only": "red", "flank_only": "yellow"}
BASE_COLOR = {"A": "green", "C": "blue", "G": "yellow", "T": "red"}


def color_seq(seq):
    return "".join(style(c, BASE_COLOR.get(c, "gray")) for c in seq)


# --- ONT context, taken from barbell src/kits/kits.rs (Dorado barcode_kits) ---
NB_FRONT = "ATTGCTAAGGTTAA"   # native adapter side, precedes the NB barcode
NB_REAR = "CAGCACCT"         # constant 8 nt directly after the NB barcode
# LSK (SQK-LSK114, no ONT barcode): the Ligation Adapter (LA) is joined straight onto the dA-tailed amplicon ends.
# ONT CHTD p.93: LA top strand 5'-CCTGTACTTCGTTCAGTTACGTATTGCT-3' (the final T is the 3' T overhang that pairs with
# the dA tail), so the insert (here: utag + primer) follows directly. Read starts are often miscalled, so only the
# last 22 nt are used.
LSK_ADAPTER = "CCTGTACTTCGTTCAGTTACGTATTGCT"
LSK_ADAPTER_TAIL = LSK_ADAPTER[-22:]
# MAB (SQK-MAB114-24): read start = rapid adapter ... GCTTGGGTGTTTAACC + AB barcode + CCATATCCGTGTCGCCCTT (Dorado),
# then directly utag + primer (seen in real reads: demo/mab). The AB barcode attaches to either amplicon end, so a read
# starts with utag-fwd + fwd primer or utag-rev + rev primer. CCCTT is the vaccinia TOPO site; barbell's kit uses only
# TTTAACC / CCATATCCGTGTC and its trim sometimes cuts into GCCCTT, so twopass pass 2 uses no prefix (bare).
MAB_TOPO = "CCATATCCGTGTCGCCCTT"
MAB_FRONT = "GCTTGGGTGTTTAACC"                 # end of the rapid adapter, directly 5' of the AB barcode
MAB_PREFIX = MAB_TOPO[-6:]                     # GCCCTT: directly 5' of the utag on untrimmed reads (-c mab)
# RAD / RBK (SQK-RAD114, SQK-RBK114): rapid adapter ... GCTTGGGTGTTTAACC [+ RB barcode] + transposase flank, then a
# random genomic/amplicon position (Dorado RBK114 rear). RAD has no barcode between the two.
RAPID_FLANK = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
MAB_SEQS = {    # Dorado TP01-TP24 = barbell AB01-AB24
    "AB01": "GCACCTGGAACTTGTGCCTTCCAC", "AB02": "CCGAAATAGGTTATCTGTTGTTGT",
    "AB03": "ATCAATCGCTGGACGATGGATTAG", "AB04": "CCACCCGCTCCTGCCGGTGGGCGT",
    "AB05": "AGACTCTTGGGCTCGCCACGTCCC", "AB06": "TCTGTATCCGGAGACGGGATGGAC",
    "AB07": "TTTCGGATCAATCGACCGCAAACG", "AB08": "ACTCAAACATTCTGTTAGATCGCG",
    "AB09": "AAATGGAACCCGGATATGTTTACT", "AB10": "TAAATCGACCTATGATGAACACAG",
    "AB11": "ACATGTTGGAGTGAAAGTCGGGTA", "AB12": "CCTGGACCACGATCATTGTAACAT",
    "AB13": "TATGGTGGATCTCCCTCTATCTTC", "AB14": "AAGTAAATGGGACGCCCACTCCGA",
    "AB15": "TGTTCGCGGCTTGATCTAATATTA", "AB16": "AGAGAGCTTCCCGGGAGGGTGGTC",
    "AB17": "TTGTGAATATCTGTCACAAACACC", "AB18": "CAATCGTACCAGGGAACATAAAGT",
    "AB19": "CACACCCAAACAATATGGACCCGT", "AB20": "AATAACCACATCCGCCCTCCGCAC",
    "AB21": "TCCTAATAATGTGTAGATCGGTCC", "AB22": "AGTCGATGGAACAAGAGAAGTTAT",
    "AB23": "AAACTCACTGTATGTCGTTTCTAT", "AB24": "TGACATCACTGATCGAGGAAGATC",
}


NB_SEQS = {
    "NB01": "CACAAAGACACCGACAACTTTCTT", "NB02": "ACAGACGACTACAAACGGAATCGA",
    "NB03": "CCTGGTAACTGGGACACAAGACTC", "NB04": "TAGGGAAACACGATAGAATCCGAA",
    "NB05": "AAGGTTACACAAACCCTGGACAAG", "NB06": "GACTACTTTCTGCCTTTGCGAGAA",
    "NB07": "AAGGATTCATTCCCACGGTAACAC", "NB08": "ACGTAACTTGGTTTGTTCCCTGAA",
    "NB09": "AACCAAGACTCGCTGTGCCTAGTT", "NB10": "GAGAGGACAAAGGTTTCAACGCTT",
    "NB11": "TCCATTCCCTCCGATAGATGAAAC", "NB12": "TCCGATTCTGCTTCTTTCTACCTG",
    "NB13": "AGAACGACTTCCATACTCGTGTGA", "NB14": "AACGAGTCTCTTGGGACCCATAGA",
    "NB15": "AGGTCTACCTCGCTAACACCACTG", "NB16": "CGTCAACTGACAGTGGTTCGTACT",
    "NB17": "ACCCTCCAGGAAAGTACCTCTGAT", "NB18": "CCAAACCCAACAACCTAGATAGGC",
    "NB19": "GTTCCTCGTGCAGTGTCAAGAGAT", "NB20": "TTGCGTCCTGTTACGAGAACTCAT",
    "NB21": "GAGCCTCTCATTGTCCGTTCTCTA", "NB22": "ACCACTGCCATGTATCAAAGTACG",
    "NB23": "CTTACTACCCAGTGAACCTCCTCG", "NB24": "GCATAGTTCTGCATGATGGGTTAG",
    "NB25": "GTAAGTTGGGTATGCAACGCAATG", "NB26": "CATACAGCGACTACGCATTCTCAT",
    "NB27": "CGACGGTTAGATTCACCTCTTACA", "NB28": "TGAAACCTAAGAAGGCACCGTATC",
    "NB29": "CTAGACACCTTGGGTTGACAGACC", "NB30": "TCAGTGAGGATCTACTTCGACCCA",
    "NB31": "TGCGTACAGCAATCAGTTACATTG", "NB32": "CCAGTAGAAGTCCGACAACGTCAT",
    "NB33": "CAGACTTGGTACGGTTGGGTAACT", "NB34": "GGACGAAGAACTCAAGTCAAAGGC",
    "NB35": "CTACTTACGAAGCTGAGGGACTGC", "NB36": "ATGTCCCAGTTAGAGGAGGAAACA",
    "NB37": "GCTTGCGATTGATGCTTAGTATCA", "NB38": "ACCACAGGAGGACGATACAGAGAA",
    "NB39": "CCACAGTGTCAACTAGAGCCTCTC", "NB40": "TAGTTTGGATGACCAAGGATAGCC",
    "NB41": "GGAGTTCGTCCAGAGAAGTACACG", "NB42": "CTACGTGTAAGGCATACCTGCCAG",
    "NB43": "CTTTCGTTGTTGACTCGACGGTAG", "NB44": "AGTAGAAAGGGTTCCTTCCCACTC",
    "NB45": "GATCCAACAGAGATGCCTTCAGTG", "NB46": "GCTGTGTTCCACTTCATTCTCCTG",
    "NB47": "GTGCAACTTTCCCACAGGTAGTTC", "NB48": "CATCTGGAACGTGGTACACCTGTA",
    "NB49": "ACTGGTGCAGCTTTGAACATCTAG", "NB50": "ATGGACTTTGGTAACTTCCTGCGT",
    "NB51": "GTTGAATGAGCCTACTGGGTCCTC", "NB52": "TGAGAGACAAGATTGTTCGTGGAC",
    "NB53": "AGATTCAGACCGTCTCATGCAAAG", "NB54": "CAAGAGCTTTGACTAAGGAGCATG",
    "NB55": "TGGAAGATGAGACCCTGATCTACG", "NB56": "TCACTACTCAACAGGTGGCATGAA",
    "NB57": "GCTAGGTCAATCTCCTTCGGAAGT", "NB58": "CAGGTTACTCCTCCGTGAGTCTGA",
    "NB59": "TCAATCAAGAAGGGAAAGCAAGGT", "NB60": "CATGTTCAACCAAGGCTTCTATGG",
    "NB61": "AGAGGGTACTATGTGCCTCAGCAC", "NB62": "CACCCACACTTACTTCAGGACGTA",
    "NB63": "TTCTGAAGTTCCTGGGTCTTGAAC", "NB64": "GACAGACACCGTTCATCGACTTTC",
    "NB65": "TTCTCAGTCTTCCTCCAGACAAGG", "NB66": "CCGATCCTTGTGGCTTCTAACTTC",
    "NB67": "GTTTGTCATACTCGTGTGCTCACC", "NB68": "GAATCTAAGCAAACACGAAGGTGG",
    "NB69": "TACAGTCCGAGCCTCATGTGATCT", "NB70": "ACCGAGATCCTACGAATGGAGTGT",
    "NB71": "CCTGGGAGCATCAGGTAGTAACAG", "NB72": "TAGCTGACTGTCTTCCATACCGAC",
    "NB73": "AAGAAACAGGATGACAGAACCCTC", "NB74": "TACAAGCATCCCAACACTTCCACT",
    "NB75": "GACCATTGTGATGAACCCTGTTGT", "NB76": "ATGCTTGTTACATCAACCCTGGAC",
    "NB77": "CGACCTGTTTCTCAGGGATACAAC", "NB78": "AACAACCGAACCTTTGAATCAGAA",
    "NB79": "TCTCGGAGATAGTTCTCACTGCTG", "NB80": "CGGATGAACATAGGATAGCGATTC",
    "NB81": "CCTCATCTTGTGAAGTTGTTTCGG", "NB82": "ACGGTATGTCGAGTTCCAGGACTA",
    "NB83": "TGGCTTGATCTAGGTAAGGTCGAA", "NB84": "GTAGTGGACCTAGAACCTGTGCCA",
    "NB85": "AACGGAGGAGTTAGTTGGATGATC", "NB86": "AGGTGATCCCAACAAGCGTAAGTA",
    "NB87": "TACATGCTCCTGTTGTTAGGGAGG", "NB88": "TCTTCTACTACCGATCCGAAGCAG",
    "NB89": "ACAGCATCAATGTTTGGCTAGTTG", "NB90": "GATGTAGAGGGTACGGTTTGAGGC",
    "NB91": "GGCTCCATAGGAACTCACGCTACT", "NB92": "TTGTGAGTGGAAAGATACAGGACC",
    "NB93": "AGTTTCCATCACTTCAGACTTGGG", "NB94": "GATTGTCCTCAAACTGCCACCTAC",
    "NB95": "CCTGTCTGGAAGAAGAATGGACTT", "NB96": "CTGAACGGTCATAGAGTCCACCAT",
}

CONTEXTS = ["bare", "nbd", "lsk", "mab"]

COMP = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")


def rc(s):
    return s.translate(COMP)[::-1]


# RBK114 (barbell RBK01-RBK96): RB01-RB12 are the reverse complements of NB01-NB12, RB13-RB96 the NB sequences,
# except six with their own sequence (barbell's table; seen in real RBK114-96 reads)
RB_SEQS = {"RB" + k[2:]: rc(v) if int(k[2:]) <= 12 else v for k, v in NB_SEQS.items()}
RB_SEQS.update({"RB26": "ACTATGCCTTTCCGTGAAACAGTT", "RB39": "TCTGCCACACACTCGTAAGTCCTT",
                "RB40": "GTCGATACTGGACCTATCCCTTGG", "RB48": "GAGTCCGTGACAACTTCTGAAAGC",
                "RB54": "GGGTGCCAACTACATACCAAACCT", "RB60": "GAACCCTACTTTGGACAGACACCT"})


# ----------------------------------------------------------------- inputs --
def read_tags(path):
    tags = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if len(row) >= 2 and row[0].strip() and not row[0].lower().startswith("utag"):
                tags[row[0].strip()] = row[1].strip().upper()
    return tags


def read_xlsx(path, sheet=None):
    """Minimal stdlib .xlsx reader -> list of rows (lists of str). First sheet unless `sheet` is a name."""
    import zipfile
    import xml.etree.ElementTree as ET
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        shared = []
        if "xl/sharedStrings.xml" in names:
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns):
                shared.append("".join(t.text or "" for t in si.iter("{%s}t" % ns["m"])))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = {r.get("Id"): r.get("Target") for r in
                ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
        sheets = [(sh.get("name"), rels[sh.get("{%s}id" % ns["r"])]) for sh in wb.find("m:sheets", ns)]
        if sheet:
            hit = [t for n, t in sheets if n == sheet]
            if not hit:
                sys.exit(f"{path}: no sheet named '{sheet}' (have: {', '.join(n for n, _ in sheets)})")
            target = hit[0]
        else:
            target = sheets[0][1]
        target = target.lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        rows = []
        for row in ET.fromstring(z.read(target)).iter("{%s}row" % ns["m"]):
            cells = {}
            for c in row.findall("m:c", ns):
                col = 0
                for ch in re.match(r"[A-Z]+", c.get("r")).group():
                    col = col * 26 + ord(ch) - 64
                t, v = c.get("t"), c.find("m:v", ns)
                if t == "s" and v is not None:
                    val = shared[int(v.text)]
                elif t == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter("{%s}t" % ns["m"]))
                elif v is not None and v.text is not None:
                    val = v.text
                    if t not in ("str", "b", "e"):
                        try:                       # numbers: 101.0 -> "101"
                            f = float(val)
                            val = str(int(f)) if f == int(f) else val
                        except ValueError:
                            pass
                else:
                    val = ""
                cells[col - 1] = val.strip()
            if cells:
                rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
    return rows


def read_table(path, sheet=None):
    """CSV or XLSX -> list of dicts keyed by lower-case column name."""
    if path.lower().endswith((".xlsx", ".xlsm")):
        raw = read_xlsx(path, sheet)
    elif path.lower().endswith(".xls"):
        sys.exit(f"{path}: old .xls is not supported; save as .xlsx or .csv")
    else:
        if not os.path.exists(path):
            sys.exit(f"{path}: file not found")
        with open(path, newline="", encoding="utf-8-sig") as fh:
            raw = list(csv.reader(fh))
    raw = [r for r in raw if any(c.strip() for c in r)]
    if not raw:
        return []
    head = [h.strip().lower() for h in raw[0]]
    return [{h: (r[i].strip() if i < len(r) else "") for i, h in enumerate(head) if h} for r in raw[1:]]


PLATE_WELL_RE = re.compile(r"^(?:(?P<plate>.+?)[-_ :/.])?(?P<row>[A-Pa-p])0*(?P<col>\d{1,2})$")


def plate_key(plate):
    """'Plate1', 'P01', '1' -> '1'; other names -> lower case. Used to match a sheet's plate to the plate map."""
    p = str(plate).strip()
    m = re.fullmatch(r"(?i)(?:plate|p)?\s*0*(\d+)", p)
    return m.group(1) if m else p.lower()


def split_plate_well(text):
    """'Plate1-A01' -> ('Plate1', 'A1'); 'B7' -> ('', 'B7'); anything else -> None."""
    m = PLATE_WELL_RE.match(str(text).strip())
    if not m or not 1 <= int(m.group("col")) <= 24:
        return None
    return m.group("plate") or "", f"{m.group('row').upper()}{int(m.group('col'))}"


def plate_label(plate):
    """'1', 'P01', 'Plate1' -> 'Plate1'; other names unchanged."""
    k = plate_key(plate)
    return f"Plate{k}" if k.isdigit() else str(plate).strip()


def design_index(design_rows):
    """Plate map rows -> {(plate key, well): row} and {(fwd, rev): (plate label, well)}."""
    by_pos, by_pair = {}, {}
    for r in design_rows:
        if r["well"]:
            by_pos[(plate_key(r["plate"]), r["well"])] = r
            pair = (r["utag-fwd"], r["utag-rev"])
            # a pair in several wells (MAB plate map: the same UDI in every AB bin) does not name one well
            by_pair[pair] = ("", "") if pair in by_pair else (r["plate"], r["well"])
    return by_pos, by_pair


def read_plates(path, sheet=None, design=None):
    """Plate map or sample sheet -> rows with keys 'plate-well' (the output name), 'utag-fwd', 'utag-rev',
    'plate', 'well', 'ont'.
      plate map:     plate,well,utag-fwd,utag-rev            (output named Plate1-A1)
      sample sheet:  sample,plate,well[,ont-barcode]         (tags looked up in the plate map `design`)
                     sample,utag-fwd,utag-rev[,ont-barcode]  (tags given directly)
    The old plate map format `plate-well,utag-fwd,utag-rev` is still read."""
    rows = read_table(path, sheet)
    if not rows:
        sys.exit(f"{path}: empty")
    cols = set(rows[0])
    name_col = "sample" if "sample" in cols else "plate-well" if "plate-well" in cols else None
    has_tags = {"utag-fwd", "utag-rev"} <= cols
    has_well = "well" in cols
    if not (has_tags or has_well) or not (name_col or has_well):
        sys.exit(f"{path}: need columns sample,plate,well (e.g. Soil_A,1,A1), sample,utag-fwd,utag-rev, "
                 f"or for a plate map plate,well,utag-fwd,utag-rev; found {sorted(cols)}")
    by_pos, by_pair = {}, {}
    if design and os.path.exists(design) and os.path.abspath(design) != os.path.abspath(path):
        by_pos, by_pair = design_index(read_plates(design))
    elif not has_tags:
        sys.exit(f"{path}: plate,well columns need the plate map to look up the tags, but {design} was not found (--plates)")
    plates_in_design = sorted({k for k, _ in by_pos}, key=lambda k: (len(k), k))
    rows = [r for r in rows if r[name_col]] if name_col else [r for r in rows if r["well"]]
    seen_name, seen_pair = {}, {}
    for r in rows:
        pw = split_plate_well(r.get("well", "")) if has_well else None
        if has_well and not pw:
            sys.exit(f"{path}: bad well '{r.get('well', '')}' (use e.g. A1 or A01)")
        plate = (r.get("plate") or (pw[0] if pw else "")).strip()
        if name_col and r[name_col]:
            name = r[name_col]
        else:                                 # plate map: the output is named Plate1-A1
            name = f"{plate_label(plate)}-{pw[1]}" if plate else pw[1]
        r["plate-well"] = re.sub(r"[^A-Za-z0-9._-]+", "_", name)   # safe file name
        r["ont"] = r.get("ont-barcode", "")
        if not has_tags:                      # plate + well -> tags from the plate map
            if not plate:
                if len(plates_in_design) != 1:
                    sys.exit(f"{path}: sample '{name}' has no plate; add a plate column "
                             f"(the plate map has plates {', '.join(plates_in_design)})")
                plate = plates_in_design[0]
            hit = by_pos.get((plate_key(plate), pw[1]))
            if not hit:
                sys.exit(f"{path}: sample '{name}': plate {plate} well {pw[1]} is not in the plate map ({design})")
            r["utag-fwd"], r["utag-rev"] = hit["utag-fwd"], hit["utag-rev"]
            r["ont"] = r["ont"] or hit.get("ont", "")   # e.g. the AB barcode of the well in a MAB plate map
        # where the sample sits, for summary.tsv and plate-heatmap.py. Its own plate/well first: in a MAB
        # plate map the same tag pair is in every AB bin, so the pair alone does not name one well
        loc = (plate_label(plate), pw[1]) if pw and plate else by_pair.get((r["utag-fwd"], r["utag-rev"]))
        if not loc and pw:
            loc = ("", pw[1])
        if not loc and name_col == "plate-well":
            loc = split_plate_well(name)
            loc = loc and (plate_label(loc[0]) if loc[0] else "", loc[1])
        r["plate"], r["well"] = loc or ("", "")
        pair = (r["ont"], r["utag-fwd"], r["utag-rev"])
        if r["plate-well"] in seen_name:
            sys.exit(f"{path}: duplicate sample name '{name}'")
        if pair in seen_pair:
            sys.exit(f"{path}: samples '{seen_pair[pair]}' and '{name}' share "
                     f"{pair[1]}/{pair[2]}" + (f" in {pair[0]}" if pair[0] else ""))
        seen_name[r["plate-well"]], seen_pair[pair] = 1, name
    return rows


def load_layout(a):
    if getattr(a, "samplesheet", None):
        return read_plates(a.samplesheet, getattr(a, "sheet_name", None), design=a.plates)
    return read_plates(a.plates)


def tag_rows(a, layout):
    """Rows that decide which tags go into the query FASTAs. With a sample sheet, the tags of the
    library design (--plates) are added as decoys, so reads carrying a tag that is not in the sheet
    are reported as unexpected instead of being forced onto the nearest listed tag."""
    if not getattr(a, "samplesheet", None) or getattr(a, "no_decoys", False):
        return layout
    if not os.path.exists(a.plates):
        print_note(f"plate map {a.plates} not found: NO decoy tags are searched, so reads from wells that are\n"
                   "not in the sample sheet are forced onto the nearest listed tag instead of being reported\n"
                   "as unexpected. Pass --plates <your plate map>, or --no-decoys to silence this.")
        return layout
    try:
        design = read_plates(a.plates)
    except SystemExit as e:
        print_note(f"plate map {a.plates} unreadable ({e}): no decoy tags are searched")
        return layout
    used = {(r["utag-fwd"], r["utag-rev"]) for r in layout}
    fwd = {r["utag-fwd"] for r in design} - {r["utag-fwd"] for r in layout}
    rev = {r["utag-rev"] for r in design} - {r["utag-rev"] for r in layout}
    pairs = {(r["utag-fwd"], r["utag-rev"]) for r in design} - used
    a.decoys = dict(fwd=fwd, rev=rev, pairs=pairs)      # for run_summary.csv
    print_info(f"decoys from {os.path.basename(a.plates)}: {len(fwd)} fwd + {len(rev)} rev tags not in the sample sheet; "
               f"{len(pairs)} tag pairs can be reported as unexpected")
    return layout + design


def udi_problem(rows):
    """Single-end (MAB) demux sees one utag per read, so every tag must name exactly one sample. Returns an error text or None."""
    n_tag = Counter(t for r in rows for t in {r["utag-fwd"], r["utag-rev"]})
    dup = sorted((t for t, n in n_tag.items() if n > 1), key=lambda t: -n_tag[t])
    if not dup:
        return None
    ex = ", ".join(f"{t} ({n_tag[t]} samples)" for t in dup[:3])
    return (f"{len(dup)} tags are shared by several samples, e.g. {ex}. Only one utag is read per read, so this needs\n"
            f"unique dual indexes (each fwd and rev tag in one sample per ONT/AB barcode), e.g. --plates {MAB_PLATES}.\n"
            f"A combinatorial (CDI) plate map like {PLATES} cannot be demultiplexed from one tag.")


def edit_distance(a, b):
    d = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        p, d[0] = d[:], i
        for j, y in enumerate(b, 1):
            d[j] = min(p[j] + 1, d[j - 1] + 1, p[j - 1] + (x != y))
    return d[-1]


# ---------------------------------------------------------------- context --
def context_prefix(context, ont_barcode=None, prefix=None):
    """Constant sequence placed immediately 5' of the utag."""
    if prefix is not None:
        return prefix.upper()
    if context == "bare":
        return ""
    if context == "lsk":
        if ont_barcode:
            sys.exit("the lsk context has no ONT barcode (SQK-LSK is unbarcoded); drop --ont-barcode")
        return LSK_ADAPTER_TAIL
    if context == "mab":
        if ont_barcode:
            sys.exit("the mab context does not put the AB barcode in the utag query (the TOPO join is variable); "
                     "drop --ont-barcode and split by AB barcode with twopass -k SQK-MAB114-24")
        return MAB_PREFIX
    if context == "nbd":
        if ont_barcode:
            if ont_barcode not in NB_SEQS:
                sys.exit(f"unknown NBD barcode {ont_barcode} (use NB01..NB96)")
            return NB_FRONT + NB_SEQS[ont_barcode] + NB_REAR
        return NB_REAR
    sys.exit(f"unknown context {context}")


def used_tags(plates, tags, all_tags=False):
    if all_tags:
        ids = sorted(tags)
        return ids, ids
    fwd = sorted({r["utag-fwd"] for r in plates})
    rev = sorted({r["utag-rev"] for r in plates})
    for t in fwd + rev:
        if t not in tags:
            sys.exit(f"tag {t} is in the sample sheet/plate map but not in the sequence file")
    # barbell needs >= 2 sequences per query FASTA: pad with unused tags (reads hitting them are
    # reported as unexpected, never written to a sample)
    spare = [t for t in sorted(tags) if t not in fwd and t not in rev]
    for side in (fwd, rev):
        while len(side) < 2 and spare:
            side.append(spare.pop(0))
    return fwd, rev


def effective_barcode_region(seqs):
    """Mimic barbell: flanks = longest common prefix/suffix of the whole file."""
    p = len(os.path.commonprefix(seqs))
    s = len(os.path.commonprefix([x[::-1] for x in seqs]))
    return p, s, len(seqs[0]) - p - s


def write_fasta(path, ids, tags, prefix, primer, label=""):
    seqs = [prefix + tags[i] + primer for i in ids]
    if len(ids) < 2:
        sys.exit("barbell needs >= 2 sequences per query FASTA")
    with open(path, "w") as fh:
        for i, s in zip(ids, seqs):
            fh.write(f">{label}{i}\n{s}\n")
    return seqs


def build_queries(outdir, plates, tags, context, ont_barcode=None, prefix=None,
                  all_tags=False, name=None, quiet=False, label="",
                  primer_f=None, primer_r=None):
    primer_f = (primer_f or PRIMER_F).upper()
    primer_r = (primer_r or PRIMER_R).upper()
    os.makedirs(outdir, exist_ok=True)
    name = name or context + (f"_{ont_barcode}" if ont_barcode else "")
    pre = context_prefix(context, ont_barcode, prefix)
    fwd, rev = used_tags(plates, tags, all_tags)
    lp, rp = (os.path.join(outdir, f"{name}_{s}.fasta") for s in ("left", "right"))
    ls = write_fasta(lp, fwd, tags, pre, primer_f, label)
    rs = write_fasta(rp, rev, tags, pre, primer_r, label)
    if not quiet:
        rows, folded = [], 0
        for side, seqs in (("left", ls), ("right", rs)):
            p, s_, m = effective_barcode_region(seqs)
            sh = len(tags[fwd[0]]) - m
            folded = max(folded, sh)
            rows.append([style(side, "cyan"), fmt_int(len(seqs)), f"{p} + {s_} nt",
                         style(f"{m} nt", "yellow" if sh else "green") + (style(f"  (-{sh} folded)", "gray") if sh else "")])
        print_table(["query", "seqs", "flank 5'+3'", "barcode region"], rows, "<>><", indent=4)
        if folded:
            print_note(f"tags share {folded} nt at the tag edge: barbell folds it into the flank (barcode region shrinks by {folded} nt)")
    return lp, rp


# -------------------------------------------------------------- barbell ---
def raise_nofile():
    """Best effort: trim opens one file per tag pair. Never fails (macOS rejects very large values)."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        for want in (65535, 20480, 10240, 4096, 2048):
            want = want if hard == resource.RLIM_INFINITY else min(want, hard)
            if soft >= want:
                return
            try:
                resource.setrlimit(resource.RLIMIT_NOFILE, (want, hard))
                return
            except (ValueError, OSError):
                continue
    except Exception:
        pass


def run(cmd, log=None):
    shown = shorten(" ".join(cmd)).split(" ")
    shown[0] = os.path.basename(shown[0]) if os.sep in shown[0] else shown[0]
    head = style(" ".join(shown[:2]), "bold", "cyan")
    print(f"  {style(SYM['prompt'], 'gray')} {head} {style(' '.join(shown[2:]), 'gray')}", flush=True)
    kw = dict(preexec_fn=raise_nofile)
    t0 = time.perf_counter()
    if log:
        with open(log, "w") as fh:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, **kw)
    else:
        r = subprocess.run(cmd, **kw)
    if r.returncode:
        if log and os.path.exists(log):
            tail = [ln.rstrip() for ln in open(log, errors="replace").readlines()
                    if ln.strip() and not re.search(r"[█╚╔║═╗╝]|\[===", ln)][-15:]
            print_error(f"command failed ({r.returncode}); last lines of {rel(log)}:", file=sys.stderr)
            for ln in tail:
                print("    " + style(ln, "gray", err=True), file=sys.stderr)
        sys.exit(f"command failed ({r.returncode}); see {rel(log) if log else 'output above'}")
    print(f"    {style(SYM['ok'], 'green')} {style(f'{time.perf_counter() - t0:.1f}s', 'gray')}")


def count_fastq_reads(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rb") as fh:
        return sum(1 for _ in fh) // 4


def help_text(barbell, sub):
    r = subprocess.run([barbell, sub, "--help"], capture_output=True, text=True)
    return r.stdout + r.stderr


def merge_reads(reads, dest):
    """Concatenate (possibly gzipped) FASTQs into one plain FASTQ for single-input barbell builds."""
    with open(dest, "wb") as o:
        for f in reads:
            with (gzip.open(f, "rb") if f.endswith(".gz") else open(f, "rb")) as i:
                shutil.copyfileobj(i, o)
    return dest


def rc_fastq_to(src, dst, gz):
    """Append src FASTQ to dst with every read reverse-complemented."""
    op = gzip.open if src.endswith(".gz") else open
    out = gzip.open(dst, "ab") if gz else open(dst, "ab")
    with op(src, "rt") as i, out as o:
        while True:
            h = i.readline()
            if not h:
                break
            seq, plus, q = i.readline().strip(), i.readline(), i.readline().strip()
            o.write(f"{h}{rc(seq)}\n+\n{q[::-1]}\n".encode())


# ----------------------------------------------------------------- demux --
FILTERS = [
    # forward-strand read: left utag-fwd ... right utag-rev (reverse complemented)
    "Ftag[fw, *, @left(0..{m}), >>]__Rtag[<<, rc, *, @right(0..{m})]",
    # reverse-strand read: left utag-rev ... right utag-fwd (reverse complemented)
    "Rtag[fw, *, @left(0..{m}), >>]__Ftag[<<, rc, *, @right(0..{m})]",
]
# --single-end (MAB): one utag at the read start, fwd or rev (the AB barcode attaches to either amplicon end);
# whatever is at the far end is not used (as barbell's own single-label kits)
FILTERS_SINGLE = ["Ftag[fw, *, @left(0..{m}), >>]", "Rtag[fw, *, @left(0..{m}), >>]"]


def label_pair(fname, tag_ids, n=2):
    """'U013__U001.trimmed.fastq[.gz]' (or with _fw/_rc suffixes) -> ('U013','U001'); n=1: 'U013.trimmed.fastq' -> ('U013',)"""
    base = re.sub(r"\.trimmed\.fastq(\.gz)?$", "", os.path.basename(fname))
    pat = "|".join(sorted(map(re.escape, tag_ids), key=len, reverse=True))
    ids = re.findall(r"(?<![A-Za-z0-9])(" + pat + r")(?![A-Za-z0-9])", base)
    return tuple(ids) if len(ids) == n else None


FASTQ_EXT = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def expand_reads(patterns):
    """-r accepts files, shell-style globs (** allowed), or directories (searched recursively for fastq/fq[.gz])."""
    found, problems = [], []
    for pat in patterns:
        p = os.path.expanduser(pat)
        if os.path.isfile(p):
            found.append(p)
        elif os.path.isdir(p):
            hits = sorted(os.path.join(d, f) for d, _, fs in os.walk(p) for f in fs if f.lower().endswith(FASTQ_EXT))
            if hits:
                found += hits
            else:
                problems.append(f"  {pat!r}: directory contains no {'/'.join(FASTQ_EXT)} files")
        else:
            hits = sorted(h for h in glob.glob(p, recursive=True) if os.path.isfile(h))
            if hits:
                found += hits
                continue
            why = "no such file or glob matches nothing"
            if os.path.islink(p):
                why = f"broken symlink -> {os.readlink(p)}"
            d = os.path.dirname(p) or "."
            near = []
            if os.path.isdir(d):
                near = sorted(f for f in os.listdir(d) if f.lower().endswith(FASTQ_EXT) or os.path.isdir(os.path.join(d, f)))[:8]
                why += f"; '{d}' exists and contains: {', '.join(near) if near else '(no fastq files or folders)'}"
            else:
                why += f"; directory '{d}' does not exist"
            problems.append(f"  {pat!r}: {why}")
    if problems:
        sys.exit("cannot find input reads (current directory: " + os.getcwd() + ")\n" + "\n".join(problems))
    seen, out = set(), []
    for f in found:
        k = os.path.abspath(f)
        if k not in seen:
            seen.add(k)
            out.append(f)
    if not out:
        sys.exit("no input reads found")
    return out


def pattern_kind(pat):
    kinds = []
    for half in pat.split("__"):
        h = half.strip()
        kinds.append("tag" if re.match(r"[FR]tag", h) else "flank" if re.match(r"[FR]flank", h) else "other")
    if kinds.count("tag") == 2:
        return "tag_both_ends"
    if kinds.count("tag") == 1 and len(kinds) == 2:
        return "tag_plus_flank_other_end"
    if kinds.count("tag") == 1:
        return "tag_one_end_only"
    return "flank_only"


def parse_patterns(path):
    # barbell colours its output when CLICOLOR_FORCE is set (e.g. RStudio/Shiny), even into a file
    txt = _ANSI_RE.sub("", open(path, errors="replace").read())
    m = re.search(r"Found (\d+) unique patterns", txt)
    pats = [(int(n), pt) for n, pt in re.findall(r"Pattern \d+: (\d+) occurrences\s*\n\s*(\S[^\n]*)", txt)]
    return (int(m.group(1)) if m else len(pats)), pats


def summarise_patterns(path):
    """Group barbell inspect patterns by what was found: tag/flank at one or both read ends."""
    cats = Counter()
    for n, pt in parse_patterns(path)[1]:
        cats[pattern_kind(pt)] += n
    return cats


def color_pattern(pat, code):
    parts = re.split(r"([FR](?:tag|flank))", pat)
    return "".join(style(x, "bold", code) if re.fullmatch(r"[FR](?:tag|flank)", x) else style(x, code) for x in parts)


def show_patterns(path, top=8):
    uniq, pats = parse_patterns(path)
    tot = sum(n for n, _ in pats) or 1
    rows = [[style(str(i), "gray"), fmt_int(n), f"{100 * n / tot:5.1f}%",
             color_pattern(pt, KIND_COLOR[pattern_kind(pt)])] for i, (n, pt) in enumerate(pats[:top], 1)]
    print_table(["#", "reads", "share", f"pattern (top {min(top, len(pats))} of {uniq})"], rows, ">>><", indent=4)


def print_found_per_read(pc):
    tot_p = sum(pc.values()) or 1
    rows = []
    for k in ("tag_both_ends", "tag_plus_flank_other_end", "tag_one_end_only", "flank_only"):
        frac = pc[k] / tot_p
        rows.append([style(KIND_LABEL[k], KIND_COLOR[k]), fmt_int(pc[k]), f"{100 * frac:5.1f}%", bar(frac, 22, KIND_COLOR[k])])
    print_table(["found", "reads", "share", ""], rows, "<>>< ", indent=4)


def print_summary(stats, n_total, extra_zero=()):
    """aligned run summary; percentages coloured (green >= 80%, yellow >= 50%, else red)"""
    rows = []
    for k, v in stats:
        pct, pc_code = "", "gray"
        if n_total and (k.startswith(("reads", "pattern_reads", "unexpected_reads")) or k == "total_reads"):
            p = 100 * v / n_total
            pct = f"{p:5.1f}%"
            if k != "total_reads" and not k.startswith("unexpected"):
                pc_code = pct_color(100 - p) if k.endswith("one_end_only") else pct_color(p)
        vcode = ("yellow" if v else "gray") if k.startswith("unexpected") else "bold"
        rows.append([style(k, "cyan"), style(fmt_int(v), vcode), style(pct, pc_code)])
    print_table(["metric", "value", "of total"], rows, "<>>", indent=4)


def cmd_demux(a):
    reads = expand_reads(a.reads)
    barbell = shutil.which(a.barbell) or sys.exit(f"barbell not found ({a.barbell}); install it first")
    tags, plates = read_tags(a.tags), load_layout(a)
    single = a.single_end
    if single:              # unique dual indexes: the utag at the read start (fwd OR rev) names the sample on its own
        err = udi_problem(plates)
        if err:
            sys.exit(f"--single-end: {err}")
        well_of = {(r[s],): r["plate-well"] for r in plates for s in ("utag-fwd", "utag-rev")}
    else:
        well_of = {(r["utag-fwd"], r["utag-rev"]): r["plate-well"] for r in plates}
    trows = tag_rows(a, plates)
    fwd_ids = {r["utag-fwd"] for r in trows}
    rev_ids = {r["utag-rev"] for r in trows}

    out = os.path.abspath(a.outdir)
    DISPLAY_BASE[0] = out
    bb = os.path.join(out, "barbell")         # barbell's own files; samples/ and the summaries stay in out
    for d in (os.path.join(bb, "queries"), os.path.join(bb, "trimmed"), os.path.join(out, "samples")):
        os.makedirs(d, exist_ok=True)
    print_step(f"[1/6] query FASTAs ({a.context}{' ' + a.ont_barcode if a.ont_barcode else ''})")
    lq, rq = build_queries(os.path.join(bb, "queries"), trows, tags, a.context,
                           a.ont_barcode, a.prefix, a.all_tags, name="cdi",
                           primer_f=a.primer_f, primer_r=a.primer_r)

    ann_help, filt_help = help_text(barbell, "annotate"), help_text(barbell, "filter")
    trim_help = help_text(barbell, "trim")
    if a.skip_trim and "--skip-trim" not in trim_help:     # fail before the long annotate step
        sys.exit("--skip-trim: this barbell build has no 'barbell trim --skip-trim'; update barbell")
    multi = "--input <INPUT>..." in ann_help          # newer builds: space-separated, many inputs
    print_info(f"barbell CLI style: {'new (space-separated, multi-input)' if multi else 'release (comma-separated, single input)'}")
    merged = None
    if any(f.endswith(".gz") for f in reads) or (not multi and len(reads) > 1):
        # gzipped input is decompressed here (works with every barbell build); old builds also need one input file
        print_info("decompressing/merging input -> reads_merged.fastq")
        merged = merge_reads(reads, os.path.join(bb, "reads_merged.fastq"))
        in_reads = [merged]
    else:
        in_reads = reads
    anno = os.path.join(bb, "annotation.tsv")
    print_step("[2/6] annotate")
    qargs = (["-q", lq, rq, "-b", "Ftag", "Rtag"] if multi
             else ["-q", f"{lq},{rq}", "-b", "Ftag,Rtag"])
    cmd = [barbell, "annotate", *qargs,
           "-i", *in_reads, "-o", anno, "-t", str(a.threads)]
    if a.flank_max_errors is not None:
        cmd += ["--flank-max-errors", str(a.flank_max_errors)]
    if a.min_score_diff is not None:
        cmd += ["--min-score-diff", str(a.min_score_diff)]
    run(cmd, os.path.join(bb, "annotate.log"))

    print_step("[3/6] inspect (top patterns -> inspect.txt)")
    run([barbell, "inspect", "-i", anno, "-n", "500"], os.path.join(bb, "inspect.txt"))
    show_patterns(os.path.join(bb, "inspect.txt"))
    pc = summarise_patterns(os.path.join(bb, "inspect.txt"))
    tot_p = sum(pc.values()) or 1
    print()
    print_found_per_read(pc)
    if not single and pc["tag_both_ends"] < 0.5 * tot_p and (pc["tag_one_end_only"] + pc["tag_plus_flank_other_end"]) > 0.3 * tot_p:
        print()
        print_note("most reads carry a utag at one end only, so no tag PAIR can be formed and few reads can be assigned.\n"
                   "Usual causes: the far end has no utag in this library/prep,\n"
                   "the far end is truncated or was trimmed by the basecaller (use untrimmed reads),\n"
                   "or a primer/context mismatch at that end (check --primer-r and run `context`).")

    print_step("[4/6] filter")
    filt = os.path.join(bb, "filters.txt")
    with open(filt, "w") as fh:
        fh.write("\n".join(f.format(m=a.max_offset) for f in (FILTERS_SINGLE if single else FILTERS)) + "\n")
    filtered = os.path.join(bb, "filtered.tsv")
    cmd = [barbell, "filter", "-i", anno, "-f", filt, "-o", filtered]
    if "--dropped" in filt_help:
        cmd += ["--dropped", os.path.join(bb, "dropped.tsv")]
    run(cmd, os.path.join(bb, "filter.log"))

    print_step("[5/6] trim" + (" (--skip-trim: reads are sorted, not trimmed)" if a.skip_trim else ""))
    tdir = os.path.join(bb, "trimmed")
    for f in glob.glob(os.path.join(tdir, "*")):
        os.remove(f)
    cmd = [barbell, "trim", "-i", filtered, "-r", *in_reads, "-o", tdir]
    for opt in ("--no-orientation", "--no-flanks"):
        if opt in trim_help:
            cmd.append(opt)
    if "--failed-out" in trim_help:
        cmd += ["--failed-out", os.path.join(bb, "trim_failed.txt")]
    native_flip = "--flip" in trim_help and not single   # barbell flips on an rc Ftag; single-end reads are flipped here
    if native_flip and not a.no_flip:
        cmd.append("--flip")          # every read -> forward-primer -> reverse-primer orientation
    do_flip_here = (not native_flip) and (not a.no_flip)
    gz_native = "--gzip" in trim_help
    if a.gzip and gz_native:
        cmd.append("--gzip")
    if a.skip_trim:
        cmd.append("--skip-trim")     # samples/ get whole reads: ONT adapter/barcode, utags and primers kept
    run(cmd, os.path.join(bb, "trim.log"))
    if merged and not a.keep_merged and os.path.exists(merged):
        os.remove(merged)

    print_step(f"[6/6] rename {'utags' if single else 'tag pairs'} -> sample names")
    ext = ".fastq.gz" if a.gzip else ".fastq"
    count_open = gzip.open if a.gzip else open
    per_well = defaultdict(list)           # sample -> [(file, needs_rc)]
    unexpected = []
    for f in sorted(glob.glob(os.path.join(tdir, "*.trimmed.fastq*"))):
        pair = label_pair(f, tags, 1 if single else 2)
        if not pair:
            unexpected.append(f)
            continue
        if single:                       # read starts with utag-rev: reverse-primer orientation (flipped only with --flip)
            rev_strand, key = pair[0] in rev_ids and pair[0] not in fwd_ids, pair
        else:
            x, y = pair
            rev_strand = x in rev_ids and y in fwd_ids       # file is Rtag__Ftag
            key = (y, x) if rev_strand else (x, y)
        if key in well_of:
            per_well[well_of[key]].append((f, rev_strand and do_flip_here))
        else:
            unexpected.append(f)
    print_info(f"{fmt_int(sum(len(v) for v in per_well.values()))} tag-pair files merged into {fmt_int(len(per_well))} sample files"
               + (style(f", {len(unexpected)} not in the sample list", "yellow") if unexpected else ""))
    counts = Counter()
    for well, files in per_well.items():
        dest = os.path.join(out, "samples", well + ext)
        if os.path.exists(dest):
            os.remove(dest)
        for f, needs_rc in files:
            if needs_rc or (a.gzip and not f.endswith(".gz")):
                if needs_rc:
                    rc_fastq_to(f, dest, a.gzip)
                else:
                    with open(f, "rb") as i, gzip.open(dest, "ab") as o:
                        shutil.copyfileobj(i, o)
            else:
                with open(f, "rb") as i, open(dest, "ab") as o:   # raw concat (gzip members concatenate validly)
                    shutil.copyfileobj(i, o)
        with count_open(dest, "rb") as i:
            counts[well] = sum(1 for _ in i) // 4
    summary = os.path.join(out, "summary.tsv")
    with open(summary, "w") as fh:
        fh.write("sample\tplate\twell\tutag_fwd\tutag_rev\treads\n")
        for r in plates:
            w = r["plate-well"]
            fh.write(f"{w}\t{r['plate']}\t{r['well']}\t{r['utag-fwd']}\t{r['utag-rev']}\t{counts.get(w, 0)}\n")
    with open(anno) as fh:
        next(fh, None)            # empty when barbell found nothing (e.g. a near-empty ONT bin)
        n_ann = len({ln.split("\t", 1)[0] for ln in fh})
    n_total = sum(count_fastq_reads(f) for f in reads)
    tot = sum(counts.values())
    stats = [("total_reads", n_total),
             ("reads_with_ge1_tag_or_flank_hit", n_ann),
             ("reads_assigned_to_wells", tot)]
    if single:
        stats.append(("pattern_reads_with_utag", sum(n for k, n in pc.items() if k != "flank_only")))
    else:
        stats += [("pattern_reads_tag_both_ends", pc["tag_both_ends"]),
                  ("pattern_reads_tag_one_end_only", pc["tag_one_end_only"] + pc["tag_plus_flank_other_end"])]
    stats += [("wells_with_reads", len(counts)),
             ("wells_expected", len(plates)),
             ("unexpected_or_unparsed_files", len(unexpected))]
    unexp_tsv = os.path.join(out, "unexpected.tsv")
    n_unexp = n_decoy = 0
    dec = getattr(a, "decoys", None) or dict(fwd=set(), rev=set(), pairs=set())
    with open(unexp_tsv, "w") as fh:          # tag pairs found on reads but not in the sample list
        fh.write("file\tutag_fwd\tutag_rev\treads\n")
        for f in unexpected:
            pair = label_pair(f, tags, 1 if single else 2)
            if pair and not single and pair[0] in rev_ids and pair[1] in fwd_ids:
                pair = pair[::-1]
            if pair and single:
                pair = ("", pair[0]) if pair[0] in rev_ids and pair[0] not in fwd_ids else (pair[0], "")
            pair = pair or ("", "")
            n = count_fastq_reads(f)
            n_unexp += n
            if (pair[0] in dec["fwd"] or pair[1] in dec["rev"]) if single else pair in dec["pairs"]:
                n_decoy += n                  # a plate-map well that is not in the sample sheet
            fh.write(f"barbell/trimmed/{os.path.basename(f)}\t{pair[0]}\t{pair[1]}\t{n}\n")
    stats += [("decoy_tags_fwd", len(dec["fwd"])), ("decoy_tags_rev", len(dec["rev"])), ("decoy_wells", len(dec["pairs"])),
              ("unexpected_reads", n_unexp),
              ("unexpected_reads_decoy_wells", n_decoy),
              ("unexpected_reads_other", n_unexp - n_decoy)]
    run_csv = os.path.join(out, "run_summary.csv")
    with open(run_csv, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["metric", "value", "pct_of_total_reads"])
        for k, v in stats:
            pct = f"{100 * v / n_total:.2f}" if n_total and (k.startswith(("reads", "pattern_reads", "unexpected_reads")) or k == "total_reads") else ""
            w.writerow([k, v, pct])
    print_step("run summary")
    print_summary(stats, n_total)
    print()
    report = write_run_report(out)
    print_kv([("samples", f"{rel(os.path.join(out, 'samples'))}/<sample>{ext}"),
              ("summary", rel(summary)), ("unexpected", rel(unexp_tsv)), ("run summary", rel(run_csv)), ("patterns", rel(os.path.join(bb, 'inspect.txt'))),
              ("barbell files", rel(bb))] + ([("report", rel(report))] if report else []))
    print()
    print_note("check barbell/inspect.txt - the two top patterns should be the two filter patterns.")


def write_run_report(out):
    """run_summary.html: compare-runs.py's page for this run alone (from utag.log and summary.tsv in out)"""
    src, log = os.path.join(os.path.dirname(os.path.realpath(__file__)), "compare-runs.py"), os.path.join(out, "utag.log")
    if os.environ.get("UTAG_NO_REPORT") or not os.path.exists(src) or not os.path.exists(log):
        return None                         # UTAG_NO_REPORT: twopass pass 2, the bins go into the twopass report
    spec = importlib.util.spec_from_file_location("compare_runs", src)
    cr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cr)
    sys.stdout.flush()
    r = cr.parse_log(log)
    r["label"] = os.path.basename(os.path.normpath(out))
    r["samples"], r["meta"] = cr.read_samples(log, r["subcommand"])
    path = os.path.join(out, "run_summary.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(cr.render_html([r], r["label"], f"utag run summary: {r['label']}", 8))
    return path


BIN_DIR_RE = re.compile(r"^(barcode\d+)$", re.I)
# the barcode Dorado / MinKNOW wrote into a read header: al:Z:barcode01, RG:Z:..._barcode01, or barcode=barcode01
HEADER_BC_RE = re.compile(r"(?:\bal:Z:|_|\bbarcode=)(barcode\d+|unclassified)\b")


def bin_num(name):
    """barcode01, NB01, RB01, 01 -> 1 (the ONT bin a sample sheet ont-barcode means); other names unchanged"""
    m = re.search(r"(\d+)$", name or "")
    return int(m.group(1)) if m else name


def dorado_bins(reads, workdir):
    """Reads already demultiplexed by Dorado / MinKNOW -> ({bin: [fastq files]}, [unclassified files]).
    Files under a barcodeNN/ folder belong to that bin; other files are split by the barcode in each read header
    into workdir/dorado_bins/<bin>.fastq."""
    bins, unclassified, other = defaultdict(list), [], []
    for f in reads:
        parts = os.path.normpath(os.path.abspath(f)).split(os.sep)[:-1]
        b = next((p for p in reversed(parts) if BIN_DIR_RE.match(p)), None)
        if b:
            bins[b.lower()].append(f)
        elif "unclassified" in parts:
            unclassified.append(f)
        else:
            other.append(f)
    if other:
        os.makedirs(os.path.join(workdir, "dorado_bins"), exist_ok=True)
        outs, no_bc = {}, 0
        for f in other:
            with (gzip.open if f.endswith(".gz") else open)(f, "rt") as fh:
                while True:
                    rec = [fh.readline() for _ in range(4)]
                    if not rec[0]:
                        break
                    m = HEADER_BC_RE.search(rec[0])
                    b = m.group(1).lower() if m else "unclassified"
                    no_bc += not m
                    if b not in outs:
                        outs[b] = open(os.path.join(workdir, "dorado_bins", b + ".fastq"), "w")
                    outs[b].writelines(rec)
        for b, fh in outs.items():
            fh.close()
            (unclassified if b == "unclassified" else bins[b]).append(fh.name)
        if no_bc:
            print_note(f"{fmt_int(no_bc)} reads in {len(other)} file(s) outside barcodeNN/ folders have no barcode in the header: "
                       "counted as unclassified")
    return dict(sorted(bins.items())), unclassified


def detect_context(files, n=2000):
    """Pass-2 context of Dorado-demultiplexed reads: is the ONT barcode still on them (untrimmed) and which kit.
    Looks for the kit sequence 5' of the utag in the first 150 nt of both read ends. Returns (context, shares)."""
    per = max(1, n // max(1, len(files)))
    seqs = []
    for f in files:
        with (gzip.open if f.endswith(".gz") else open)(f, "rt") as fh:
            for i, ln in enumerate(fh):
                if i // 4 >= per or len(seqs) >= n:
                    break
                if i % 4 == 1:
                    seqs.append(ln.strip().upper())
    probes = {"mab": lambda w: locate(MAB_TOPO, w, present_edits(len(MAB_TOPO))),
              "nbd": lambda w: NB_REAR in w and locate(NB_FRONT, w, present_edits(len(NB_FRONT))),
              "lsk": lambda w: locate(LSK_ADAPTER_TAIL, w, present_edits(len(LSK_ADAPTER_TAIL)))}
    share = {c: sum(1 for x in seqs if any(f(w) for w in (x[:150], rc(x[-150:])))) / max(1, len(seqs)) for c, f in probes.items()}
    ctx = next((c for c in ("mab", "nbd", "lsk") if share[c] >= 0.3), "bare")
    return ctx, share


def barbell_pass1(a, reads, barbell, bb, kitdir):
    """twopass pass 1: barbell kit splits the reads by ONT barcode (and trims it) -> {bin: [fastq]}"""
    print_step(f"[pass 1] barbell kit -k {a.kit}")
    kit_help = help_text(barbell, "kit")
    in_reads = reads
    merged = None
    if any(f.endswith(".gz") for f in reads) or ("--input <INPUT>..." not in kit_help and len(reads) > 1):
        print_info("decompressing/merging input -> reads_merged.fastq")
        merged = merge_reads(reads, os.path.join(bb, "reads_merged.fastq"))
        in_reads = [merged]
    cmd = [barbell, "kit", "-k", a.kit, "-i", *in_reads, "-o", kitdir, "-t", str(a.threads)]
    if "--maximize" in kit_help and not a.no_maximize:
        cmd.append("--maximize")
    if a.kit_flank_max_errors is not None:
        cmd += ["--flank-max-errors", str(a.kit_flank_max_errors)]
    if a.gzip and "--gzip" in kit_help:
        cmd.append("--gzip")
    run(cmd, os.path.join(bb, "pass1_kit.log"))
    if merged:
        os.remove(merged)
    bins = {}
    for f in sorted(glob.glob(os.path.join(kitdir, "*.trimmed.fastq*"))):
        bins[re.sub(r"\.trimmed\.fastq(\.gz)?$", "", os.path.basename(f))] = [f]
    if not bins:
        sys.exit("pass 1 produced no ONT barcode bins; check pass1_kit.log and --kit")
    print_info(f"{len(bins)} ONT bins: {', '.join(bins)}")
    return bins


def dorado_pass1(a, reads, bb):
    """twopass --dorado: the bins Dorado / MinKNOW already made -> {bin: [fastq]}; unclassified reads are left out"""
    print_step("[pass 1] Dorado bins (barbell kit skipped)")
    bins, unclassified = dorado_bins(reads, bb)
    if not bins:
        sys.exit("--dorado: no barcodeNN/ folders and no barcode in the read headers (al:Z:, RG:Z:..._barcodeNN, barcode=); "
                 "use -k to demultiplex with barbell kit")
    print_info(f"{len(bins)} ONT bins: {', '.join(bins)}")
    if unclassified:
        n = sum(count_fastq_reads(f) for f in unclassified)
        print_note(f"unclassified: {fmt_int(n)} reads in {len(unclassified)} file(s), not demultiplexed "
                   "(Dorado found no ONT barcode; twopass -k would try them again with barbell kit)")
    return bins


def pass1_only(reads, bins, out, barbell, kitdir):
    """twopass --pass1-only: the ONT bins of barbell kit are the result. Move them to ont_bins/, count reads, summarise
    (with barbell's top read patterns, as demux does)."""
    dest = os.path.join(out, "ont_bins")
    os.makedirs(dest, exist_ok=True)
    counts = {}
    for b, fs in bins.items():
        counts[b] = sum(count_fastq_reads(f) for f in fs)
        for f in fs:
            shutil.move(f, os.path.join(dest, b + (".fastq.gz" if f.endswith(".gz") else ".fastq")))
    n_total = sum(count_fastq_reads(f) for f in reads)
    n_bins = sum(counts.values())
    totals = [("total_reads", n_total), ("reads_in_ont_bins", n_bins), ("ont_bins_found", len(bins))]
    with open(os.path.join(out, "ont_bins.tsv"), "w") as fh:
        fh.write("ont_barcode\treads_in_bin\tpct_of_total_reads\n")
        fh.writelines(f"{b}\t{n}\t{100 * n / max(1, n_total):.2f}\n" for b, n in counts.items())
    with open(os.path.join(out, "run_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["metric", "value", "pct_of_total_reads"])
        for k, v in totals:
            w.writerow([k, v, f"{100 * v / n_total:.2f}" if n_total and k != "ont_bins_found" else ""])
    anno, insp = os.path.join(kitdir, "annotation.tsv"), os.path.join(out, "barbell", "pass1_inspect.txt")
    if os.path.exists(anno):
        print_step("pass 1 inspect (top patterns -> pass1_inspect.txt)")
        run([barbell, "inspect", "-i", anno, "-n", "500"], insp)
        show_patterns(insp)
        print_info(f"{fmt_int(sum(n for n, _ in parse_patterns(insp)[1]))} reads with a pattern (share = of these)")
    print_step("pass 1 summary")
    print_summary(totals, n_total)
    print()
    rows = [[style(b, "cyan"), fmt_int(n), f"{100 * n / max(1, n_total):5.1f}%", bar(n / max(1, n_total), 18)]
            for b, n in sorted(counts.items(), key=lambda x: -x[1])]
    print_table(["ONT bin", "reads", "of total", ""], rows, "<>><", indent=4)
    print()
    report = write_run_report(out)
    print_kv([("ONT bins", f"{rel(dest)}/<ONT bin>.fastq"), ("tables", f"{rel(os.path.join(out, 'ont_bins.tsv'))}, run_summary.csv"),
              ("barbell files", rel(os.path.join(out, "barbell")))] + ([("report", rel(report))] if report else []))
    print_note("pass 1 only: no utag demux; run utag demux -c bare on a bin, or twopass without --pass1-only")


def cmd_twopass(a):
    """Pass 1: barbell kit (ONT barcode), or the bins Dorado already made (--dorado). Pass 2: barbell utag demux
    inside every ONT bin."""
    if not a.kit and not a.dorado:
        sys.exit("-k/--kit is required (or --dorado for reads already demultiplexed by Dorado / MinKNOW)")
    if a.pass1_only and a.dorado:
        sys.exit("--pass1-only splits the reads with barbell kit; with --dorado they are split already (barcodeNN/ folders)")
    reads = expand_reads(a.reads)
    barbell = shutil.which(a.barbell) or sys.exit(f"barbell not found ({a.barbell}); install it first")
    mab = (a.kit or "").upper().startswith("SQK-MAB")
    if a.mab_defaults:
        print_info(f"MAB kit: using the unique dual index design --plates {rel(a.plates)} --tags {rel(a.tags)}")
    rows = load_layout(a)
    ctx2 = a.context or "bare"                          # pass 1 trimmed the ONT barcode and its flanks (--dorado: detected below)
    single = a.single_end or mab                        # MAB reads carry barcode + utag at the read start only
    if single:                                          # check before the long pass 1, per ONT/AB bin
        by_bin = defaultdict(list)
        for r in rows:
            by_bin[r["ont"]].append(r)
        for b, rs in sorted(by_bin.items()):
            err = udi_problem(rs)
            if err:
                sys.exit(f"--single-end{' (' + b + ')' if b else ''}: {err}")
    if mab:
        print_info("MAB kit: pass 2 searches one utag (fwd or rev) + primer at the read start; the far read end is not used")
    out = os.path.abspath(a.outdir)
    DISPLAY_BASE[0] = out
    bb = os.path.join(out, "barbell")         # barbell's own files; samples/ and the summaries stay in out
    kitdir, p2dir = os.path.join(bb, "pass1_kit"), os.path.join(bb, "pass2_utag")
    for d in (kitdir, p2dir):
        os.makedirs(d, exist_ok=True)

    # ---- pass 1: ONT barcodes (bins as lists of fastq files)
    if a.dorado:
        bins = dorado_pass1(a, reads, bb)
        if not a.context:                       # Dorado may or may not have trimmed the barcode: look at the reads
            big = max(bins, key=lambda b: sum(os.path.getsize(f) for f in bins[b]))
            ctx2, share = detect_context(bins[big])
            print_info(f"pass 2 context: -c {ctx2} (detected in {big}: " +
                       ", ".join(f"{c} adapter in {100 * v:.0f}%" for c, v in share.items()) +
                       "; >= 30% -> that context, else bare = barcode trimmed; set -c to override)")
            if ctx2 == "mab" and not single:
                single = True
                print_info("MAB reads: pass 2 searches one utag (fwd or rev) + primer at the read start (--single-end)")
                by_bin = defaultdict(list)
                for r in rows:
                    by_bin[r["ont"]].append(r)
                for b, rs in sorted(by_bin.items()):
                    err = udi_problem(rs)
                    if err:
                        sys.exit(f"--single-end{' (' + b + ')' if b else ''}: {err}")
    else:
        bins = barbell_pass1(a, reads, barbell, bb, kitdir)
        if a.pass1_only:
            return pass1_only(reads, bins, out, barbell, kitdir)

    # ---- pass 2: utags inside each bin
    has_ont_col = any(r["ont"] for r in rows)
    bin_rows, skipped = {}, []
    for b in bins:
        bin_rows[b] = [r for r in rows if (not has_ont_col) or bin_num(r["ont"]) == bin_num(b)]   # NB01 = barcode01
        if not bin_rows[b]:
            skipped.append(b)
    if skipped:
        print_note(f"no samples listed for bin(s) {', '.join(skipped)}: skipped in pass 2")
    stats, per_sample, per_unexp = [], [], []
    unexp_tot = Counter()
    for b, fs in bins.items():
        if not bin_rows[b]:
            continue
        print_step(f"[pass 2] {b}: {len(bin_rows[b])} samples")
        bdir = os.path.join(p2dir, b)
        os.makedirs(bdir, exist_ok=True)
        sheet = os.path.join(bdir, "sheet.csv")
        with open(sheet, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["sample", "plate", "well", "utag-fwd", "utag-rev"])
            for r in bin_rows[b]:
                w.writerow([r["plate-well"], r["plate"], r["well"], r["utag-fwd"], r["utag-rev"]])
        cmd = [sys.executable, os.path.abspath(__file__), "demux", "-r", *fs, "-c", ctx2,
               "-o", bdir, "-t", str(a.threads), "--barbell", barbell,
               "--samplesheet", sheet, "--plates", a.plates, "--tags", a.tags,
               "--max-offset", str(a.max_offset)]
        for opt, val in (("--primer-f", a.primer_f), ("--primer-r", a.primer_r),
                         ("--flank-max-errors", a.flank_max_errors),
                         ("--min-score-diff", a.min_score_diff)):
            if val is not None:
                cmd += [opt, str(val)]
        cmd += [x for x, on in (("--flip", not a.no_flip), ("--gzip", a.gzip), ("--skip-trim", a.skip_trim),
                                ("--no-decoys", a.no_decoys), ("--single-end", single)) if on]
        r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "NO_COLOR": "1", "UTAG_NO_REPORT": "1"})
        if r.returncode:
            for ln in (r.stdout + r.stderr).strip().splitlines()[-10:]:
                print("    " + style(ln, "gray", err=True), file=sys.stderr)
            sys.exit(f"pass 2 failed for {b}; see {rel(os.path.join(bdir, 'utag.log'))}")
        m = {k: v for k, v, *_ in csv.reader(open(os.path.join(bdir, "run_summary.csv")))}
        stats.append((b, int(m["total_reads"]), int(m["reads_with_ge1_tag_or_flank_hit"]),
                      int(m["reads_assigned_to_wells"]), int(m["wells_with_reads"]), len(bin_rows[b])))
        for k in ("unexpected_reads", "unexpected_reads_decoy_wells", "unexpected_reads_other"):
            unexp_tot[k] += int(m.get(k, 0))
        # final layout: samples/<ONT bin>/<sample>.fastq
        dest = os.path.join(out, "samples", b)
        os.makedirs(dest, exist_ok=True)
        for fq in glob.glob(os.path.join(bdir, "samples", "*")):
            shutil.move(fq, os.path.join(dest, os.path.basename(fq)))
        shutil.rmtree(os.path.join(bdir, "samples"), ignore_errors=True)   # emptied: moved to samples/<bin>/
        with open(os.path.join(bdir, "summary.tsv")) as fh:
            next(fh)
            per_sample += [[b] + ln.rstrip("\n").split("\t") for ln in fh]
        ux = os.path.join(bdir, "unexpected.tsv")
        if os.path.exists(ux):
            with open(ux) as fh:
                next(fh)
                for ln in fh:
                    f, *rest = ln.rstrip("\n").split("\t")
                    per_unexp.append([b, os.path.relpath(os.path.join(bdir, f), out)] + rest)
        st = stats[-1]
        print_info(f"{fmt_int(st[1])} reads in bin, {style(fmt_int(st[3]), 'bold', pct_color(100 * st[3] / max(1, st[1])))} assigned "
                   + style(f"({100 * st[3] / max(1, st[1]):.1f}%)", "gray"))

    with open(os.path.join(out, "ont_bins.tsv"), "w") as fh:
        fh.write("ont_barcode\treads_in_bin\treads_with_ge1_tag_or_flank_hit\treads_assigned\tsamples_with_reads\tsamples_expected\n")
        for st in stats:
            fh.write("\t".join(map(str, st)) + "\n")
    with open(os.path.join(out, "unexpected.tsv"), "w") as fh:
        fh.write("ont_barcode\tfile\tutag_fwd\tutag_rev\treads\n")
        fh.writelines("\t".join(x) + "\n" for x in per_unexp)
    with open(os.path.join(out, "summary.tsv"), "w") as fh:
        fh.write("ont_barcode\tsample\tplate\twell\tutag_fwd\tutag_rev\treads\n")
        fh.write("\n".join("\t".join(x) for x in per_sample) + "\n")
    n_total = sum(count_fastq_reads(f) for f in reads)
    n_bins = sum(st[1] for st in stats) + sum(count_fastq_reads(f) for b in skipped for f in bins[b])
    totals = [("total_reads", n_total), ("reads_in_ont_bins", n_bins),
              ("reads_with_ge1_tag_or_flank_hit", sum(st[2] for st in stats)),
              ("reads_assigned_to_samples", sum(st[3] for st in stats)),
              ("ont_bins_found", len(bins)), ("ont_bins_processed", len(stats)),
              ("samples_with_reads", sum(st[4] for st in stats)),
              ("samples_expected", sum(st[5] for st in stats))] + \
             [(k, unexp_tot[k]) for k in ("unexpected_reads", "unexpected_reads_decoy_wells", "unexpected_reads_other")]
    with open(os.path.join(out, "run_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["metric", "value", "pct_of_total_reads"])
        for k, v in totals:
            w.writerow([k, v, f"{100 * v / n_total:.2f}" if n_total and (k.startswith(("reads", "unexpected_reads")) or k == "total_reads") else ""])
    print_step("two-pass summary")
    print_summary(totals, n_total)
    print()
    rows = []
    for b, nb, hit, asg, w_with, w_exp in stats:
        frac = asg / max(1, nb)
        rows.append([style(b, "cyan"), fmt_int(nb), fmt_int(hit), style(fmt_int(asg), "bold"),
                     style(f"{100 * frac:5.1f}%", pct_color(100 * frac)), bar(frac, 18, pct_color(100 * frac)),
                     f"{w_with}/{w_exp}"])
    print_table(["ONT bin", "reads", "with hit", "assigned", "assigned %", "", "samples with reads"], rows, "<>>>>><>", indent=4)
    print()
    report = write_run_report(out)
    print_kv([("samples", f"{rel(os.path.join(out, 'samples'))}/<ONT bin>/<sample>.fastq"),
              ("tables", f"{rel(os.path.join(out, 'ont_bins.tsv'))}, summary.tsv, unexpected.tsv, run_summary.csv"),
              ("barbell files", rel(bb))] + ([("report", rel(report))] if report else []))


IUPAC_RE = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]",
            "S": "[CG]", "W": "[AT]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]", "I": "[ACGT]"}


def suffix_edits(adapter, up):
    """d[k] = fewest edits (indels allowed) to align adapter[-k:] so that it ends exactly where `up` ends (at the
    utag); it may start anywhere in `up`. One DP over the reversed strings gives every suffix length at once."""
    q, t = adapter[::-1], up[::-1][:len(adapter) + 12]
    row = list(range(len(t) + 1))           # 0 adapter bases used: skipping read bases next to the utag costs
    d = [0]
    for x in q:
        new = [row[0] + 1]
        for j, y in enumerate(t, 1):
            new.append(min(row[j] + 1, new[j - 1] + 1, row[j - 1] + (x != y)))
        row = new
        d.append(min(row))
    return d


HIT_PRIMER_NT, HIT_EDITS = 18, 2       # context hit: primer start (<= 2 edits) + utag (<= 2 edits) right before it
IUPAC_SET = {k: set(v.strip("[]")) for k, v in IUPAC_RE.items()}


def primer_start(p, text, k):
    """(start, edits) of the whole primer p (IUPAC) in text within k edits, or None. DP on the reversed strings:
    the primer end is free, so the best row-end gives where it starts."""
    q, t = p[::-1], text[::-1]
    prev = [0] * (len(t) + 1)
    for i, a in enumerate(q, 1):
        ok, cur = IUPAC_SET.get(a, {a}), [i]
        for j, b in enumerate(t, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (b not in ok)))
        prev = cur
    e = min(prev)
    return (len(t) - prev.index(e), e) if e <= k else None


def tag_before(tag, region):
    """(edits, read bases used) for tag aligned so that it ends where region ends (indels allowed)"""
    q, t = tag[::-1], region[::-1]
    row = list(range(len(t) + 1))
    for x in q:
        new = [row[0] + 1]
        for j, y in enumerate(t, 1):
            new.append(min(row[j] + 1, new[j - 1] + 1, row[j - 1] + (x != y)))
        row = new
    e = min(row)
    return e, row.index(e)


def primer_seeds(p):
    """the primer's first HIT_PRIMER_NT nt in 3 pieces: with <= 2 edits one piece is exact (pigeonhole)"""
    p = p.upper()[:HIT_PRIMER_NT]
    b = [0, len(p) // 3, 2 * len(p) // 3, len(p)]
    return p, [(re.compile("".join(IUPAC_RE[c] for c in p[b[i]:b[i + 1]])), b[i]) for i in range(3)]


# known sequence directly 5' of the utag: (context, name, sequence, -c query length)
KNOWN_CONTEXTS = [("lsk", "Ligation Adapter (SQK-LSK114)", LSK_ADAPTER, len(LSK_ADAPTER_TAIL)),
                  ("nbd", "NB rear (SQK-NBD114)", NB_REAR, len(NB_REAR)),
                  ("mab", "TOPO site (SQK-MAB114)", MAB_TOPO, len(MAB_PREFIX))]


BARCODE_SAMPLE, BARCODE_PER, BARCODE_MAX = 150, 30, 3000  # hits checked for the ONT barcode (96 alignments each):
                                                           # more when many barcodes show up, ~30 per barcode
# what sits 5' of the utag for kits with an ONT barcode, 5' -> 3': [label, sequence or None (= the barcode), length]
BARCODE_LAYOUT = {"lsk": [["Ligation Adapter", LSK_ADAPTER[:-1], len(LSK_ADAPTER) - 1], ["T", "T", 1]],   # T: 3' overhang on the dA tail
                  "nbd": [["NB front", NB_FRONT, len(NB_FRONT)], ["NB barcode", None, 24], ["NB rear", NB_REAR, len(NB_REAR)]],
                  "mab": [["rapid adapter", MAB_FRONT, len(MAB_FRONT)], ["AB barcode", None, 24], ["TOPO site", MAB_TOPO, len(MAB_TOPO)]]}


def present_edits(n):
    """edits allowed for a known adapter to count as present (~15%, at least 1)"""
    return max(1, round(0.15 * n))


def best_barcode(seqs, up_list, flank_len):
    """Which ONT barcode sits just 5' of a known flank of flank_len nt: Counter of barcode -> reads (or '?')."""
    c = Counter()
    for up in up_list:
        region = up[:-flank_len] if flank_len else up
        scores = sorted((suffix_edits(s, region)[-1], b) for b, s in seqs.items())
        best, second = scores[0], scores[1]
        c[best[1] if best[0] <= 6 and second[0] - best[0] >= 2 else "?"] += 1
    return c


KIT_WINDOW, KIT_SAMPLE = 160, 1000       # read-start nt and reads checked for the kit adapters (no utag hits)
# adapters looked for at the read start: key -> sequence
KIT_ADAPTERS = {"la": LSK_ADAPTER_TAIL, "nb_front": NB_FRONT, "nb_rear": NB_REAR,
                "rapid_flank": RAPID_FLANK, "topo": MAB_TOPO}
# kit -> (layout, barcode set or None); checked in this order, the first that fits a read wins
KITS = {"SQK-NBD114": ("ligation adapter + NB front + NB barcode + NB rear", NB_SEQS),
        "SQK-MAB114": ("(rapid adapter + AB barcode +) TOPO site", MAB_SEQS),
        "SQK-RAD114": ("rapid adapter + transposase flank, no barcode", None),
        "SQK-RBK114": ("(rapid adapter + RB barcode +) transposase flank", RB_SEQS),
        "SQK-LSK114": ("ligation adapter only", None)}
# what sits 5' of the insert, 5' -> 3' (as BARCODE_LAYOUT): the kit alignment in the report
KIT_LAYOUT = {"SQK-NBD114": BARCODE_LAYOUT["nbd"], "SQK-MAB114": BARCODE_LAYOUT["mab"], "SQK-LSK114": BARCODE_LAYOUT["lsk"],
              "SQK-RBK114": [["rapid adapter", MAB_FRONT, len(MAB_FRONT)], ["RB barcode", None, 24],
                             ["transposase flank", RAPID_FLANK, len(RAPID_FLANK)]],
              "SQK-RAD114": [["rapid adapter", MAB_FRONT, len(MAB_FRONT)], ["transposase flank", RAPID_FLANK, len(RAPID_FLANK)]]}


def locate(p, text, k):
    """(start, end, edits) of p in text within k edits (indels allowed), or None. Exact seeds (k+1 pieces: one is
    exact, pigeonhole) pick the candidates; a semi-global DP around each gives the end, primer_start the start."""
    b = [i * len(p) // (k + 1) for i in range(k + 2)]
    cands, best = sorted({m.start() - b[i] for i in range(k + 1)
                          for m in re.finditer("(?=" + p[b[i]:b[i + 1]] + ")", text)}), None
    done = -1
    for c in cands:
        if c <= done:
            continue
        lo, hi = max(0, c - k), min(len(text), c + len(p) + k)
        done = c + k
        t = text[lo:hi]
        row = [0] * (len(t) + 1)
        for i, x in enumerate(p, 1):
            new = [i]
            for j, y in enumerate(t, 1):
                new.append(min(row[j] + 1, new[j - 1] + 1, row[j - 1] + (x != y)))
            row = new
        e = min(row)
        if e <= k and (best is None or e < best[2]):
            end = lo + row.index(e)
            best = (lo + primer_start(p, t[:end - lo], k)[0], end, e)
    return best


def read_kit(x):
    """(kit, read up to the end of its ONT barcode, read up to the insert) for the start of read x;
    (None, None, None) if no kit fits"""
    f = {key: locate(seq, x, present_edits(len(seq))) for key, seq in KIT_ADAPTERS.items()}
    if f["nb_front"] and f["nb_rear"] and 18 <= f["nb_rear"][0] - f["nb_front"][1] <= 30:   # a 24 nt barcode fits
        return "SQK-NBD114", x[:f["nb_rear"][0]], x[:f["nb_rear"][1]]
    if f["topo"]:                               # distinctive; the rapid adapter before the AB barcode is often garbled
        return "SQK-MAB114", x[:f["topo"][0]], x[:f["topo"][1]]
    if f["rapid_flank"]:                        # rapid adapter end right before the flank (RAD) or ~24 nt before (RBK)
        r = x[:f["rapid_flank"][0]]
        rad = suffix_edits(MAB_FRONT, r)[-1]
        rbk = min(suffix_edits(MAB_FRONT, r[:-n])[-1] for n in range(20, 29))
        return ("SQK-RAD114", None, x[:f["rapid_flank"][1]]) if rad <= rbk else ("SQK-RBK114", r, x[:f["rapid_flank"][1]])
    if f["la"]:
        return "SQK-LSK114", None, x[:f["la"][1]]
    return None, None, None


def base_profile(L):
    """per position of the equal-length strings L ('-' = no base): ([A, C, G, T] share of all strings, consensus
    base or '.', % of all strings with the consensus base)"""
    cons, supp, logo = [], [], []
    for pos in range(len(L[0]) if L else 0):
        c = Counter(x[pos] for x in L if x[pos] != "-")
        logo.append([c[b] / len(L) for b in "ACGT"])
        n = sum(c.values())
        b, k = c.most_common(1)[0] if c else ("-", 0)
        cons.append(b if n and k / len(L) >= 0.5 else ".")
        supp.append(100 * k / len(L) if n else 0)
    return logo, "".join(cons), supp


def kit_alignment(kit, ups, barcodes):
    """The kit's reads aligned at the start of the insert: consensus over the kit layout (terminal), and the base
    frequencies for the report"""
    lay = KIT_LAYOUT[kit]
    U = min(KIT_WINDOW, sum(n for _, _, n in lay) + 10)
    L = [u[-U:].rjust(U, "-") for u in ups]
    logo, cons, supp = base_profile(L)
    one = barcodes[0] if barcodes and barcodes[0][1] >= 0.8 else None     # one barcode: show its sequence
    ref = "".join(sq or (one[2] if one else "·" * n) for _, sq, n in lay)[-U:].rjust(U)
    print()
    print_kv([("consensus", color_seq(cons) + "  " + style("<insert>", "bold", "magenta")),
              ("support", "".join(style(str(min(9, int(v // 10))), pct_color(v)) for v in supp)),
              (kit, color_seq(ref)),
              ("", style("  ".join(f"{lab} {n}" for lab, _, n in lay) + f"; {fmt_int(len(L))} reads aligned at the insert", "gray"))],
             indent=4)
    return {"logo": logo, "consensus": cons, "support": supp, "layout": lay, "reads": len(L)}


def kit_from_starts(starts):
    """No utag hits: say which ONT kit the read starts look like (the first KIT_SAMPLE starts) and, for barcoded
    kits, which barcodes (more starts when many barcodes show up). Returns a dict for the HTML report."""
    print_step("context: kit from the read starts (no utag hits)")
    calls = [read_kit(x) for x in starts[:KIT_SAMPLE]]
    c = Counter(k for k, _, _ in calls)
    n = max(1, len(calls))
    rows = [[k, lay, style(f"{100 * c[k] / n:5.1f}%", pct_color(100 * c[k] / n))] for k, (lay, _) in KITS.items()]
    rows.append(["none", "no known adapter in the first " + str(KIT_WINDOW) + " nt", f"{100 * c[None] / n:5.1f}%"])
    print_info(f"first {KIT_WINDOW} nt of {fmt_int(len(calls))} reads, known adapters within ~15% edits")
    print_table(["kit", "read start", "reads"], rows, "<<>", indent=4)
    K = {"n": len(calls), "window": KIT_WINDOW,
         "rows": [[k, lay, c[k] / n] for k, (lay, _) in KITS.items()] + [["none", "", c[None] / n]]}
    kit, m = max(((k, c[k]) for k in KITS), key=lambda kv: kv[1])
    if m < 0.2 * n:
        K["verdict"] = "no known ONT adapter at the read starts: the reads look trimmed or the kit is unknown"
        print_note(K["verdict"])
        return K
    K["kit"] = kit
    print_ok(f"read starts look like {style(kit, 'bold')} ({100 * m / n:.0f}% of reads)")
    seqs = KITS[kit][1]
    if seqs:                                    # which ONT barcodes sit 5' of the flank (a subsample)
        more = (read_kit(x) for x in starts[KIT_SAMPLE:])
        bc, nc = call_barcodes(seqs, (r for k, r, _ in itertools.chain(calls, more) if k == kit), 0)
        top = [(b, v) for b, v in bc.most_common() if b != "?"]
        K.update(barcodes=[[b, v / nc, seqs[b]] for b, v in top], checked=nc, unclear=bc["?"] / nc, set=list(seqs.items()))
    K["align"] = kit_alignment(kit, [u for k, _, u in calls if k == kit], K.get("barcodes"))
    if seqs:
        print()
        print_barcodes(bc, nc, "reads")
    if kit in ("SQK-RAD114", "SQK-RBK114"):
        K["verdict"] = ("rapid kit: the transposase cuts the amplicons at random, so the utags are rarely at a read end "
                        "and utag cannot demultiplex these reads; use a ligation kit (SQK-LSK114 / SQK-NBD114) on the "
                        "intact amplicons")
    else:
        K["verdict"] = "the kit adapter is there but no utag + primer: check --primer-f/--primer-r, the tag file and the sample sheet"
    print_note(K["verdict"])
    return K

def call_barcodes(seqs, regions, flank_len):
    """best_barcode over a sample of regions (any iterable, read lazily): BARCODE_SAMPLE first, then more until
    there are ~BARCODE_PER per barcode found (up to BARCODE_MAX). Returns (Counter, regions checked)."""
    it = iter(regions)
    bc = best_barcode(seqs, list(itertools.islice(it, BARCODE_SAMPLE)), flank_len)
    while True:
        n = sum(bc.values())
        want = min(BARCODE_MAX, BARCODE_PER * sum(1 for k in bc if k != "?"))
        more = list(itertools.islice(it, want - n)) if want > n else []
        if not more:
            return bc, n
        bc += best_barcode(seqs, more, flank_len)


def print_barcodes(bc, n, what):
    """all barcodes found, most reads first, wrapped to the terminal; returns [(barcode, reads)]"""
    top = [(k, v) for k, v in bc.most_common() if k != "?"]
    items = [f"{k} {100 * v / n:.1f}%" for k, v in top] or ["none identified"]
    lines, width = [""], max(40, term_width() - 22)
    for it in items:
        if lines[-1] and len(lines[-1]) + len(it) + 2 > width:
            lines[-1] += ","
            lines.append("")
        lines[-1] += (", " if lines[-1] else "") + it
    print_kv([("ONT barcodes", lines[0])] + [("", ln) for ln in lines[1:]]
             + [("", style(f"{len(top)} barcodes in {fmt_int(n)} {what} checked; "
                           f"{100 * bc['?'] / max(1, n):.0f}% unclear", "gray"))], indent=4)
    return top


def cmd_context(a):
    """Read the constant sequence upstream of the utags off the reads: align it against the known kit adapters
    (indels allowed) and say which -c context fits. Read starts and read ends are reported apart: the read end
    often stops inside the adapter. The per-position consensus is shown for unknown kits."""
    reads = expand_reads(a.reads)
    tags, plates = read_tags(a.tags), load_layout(a)
    fwd, rev = used_tags(plates, tags)
    pf = re.compile("".join(IUPAC_RE[c] for c in (a.primer_f or PRIMER_F).upper()[:10]))
    pr = re.compile("".join(IUPAC_RE[c] for c in (a.primer_r or PRIMER_R).upper()[:10]))
    U = a.upstream
    ups = {"start": [], "end": []}           # utag + primer read 5'->3' on the read (start) or on its rc (end)
    sides = {"start": Counter(), "end": Counter()}
    with_hit = Counter()                     # reads with >= 1 hit: "start", "end", "any"
    starts = []                              # read starts, for the kit check when there are no hits
    def scan_exact(where, x):
        for side, ids, pat in (("fwd", fwd, pf), ("rev", rev, pr)):
            for t in ids:
                tg = tags[t]
                i = x.find(tg)
                while i >= 0:
                    if pat.match(x, i + len(tg)):
                        ups[where].append(x[max(0, i - U):i].rjust(U, "-"))
                        sides[where][side] += 1
                    i = x.find(tg, i + 1)

    K, TL = HIT_EDITS, max(len(tags[t]) for t in fwd + rev)
    seeds = (("fwd", [tags[t] for t in fwd], *primer_seeds(a.primer_f or PRIMER_F)),
             ("rev", [tags[t] for t in rev], *primer_seeds(a.primer_r or PRIMER_R)))

    def scan_relaxed(where, x):
        """one hit per read end: the primer within --max-offset of the read end, the utag right before it"""
        w, best = x[:a.max_offset + TL + HIT_PRIMER_NT + K], None
        for side, tg, ps, pieces in seeds:
            tried = []
            for c in sorted({m.start() - off for rx, off in pieces for m in rx.finditer(w)}):
                if c < 6 or any(abs(c - d) <= K for d in tried):
                    continue
                tried.append(c)
                lo = c - K
                hit = primer_start(ps, w[lo:c + len(ps) + K], K)
                if not hit:
                    continue
                st, pe = lo + hit[0], hit[1]
                region = w[max(0, st - TL - 3):st]
                ex = next((t for t in tg if region.endswith(t)), None)
                te, used = (0, len(ex)) if ex else min(tag_before(t, region) for t in tg)
                if te <= K and (best is None or pe + te < best[0]):
                    best = (pe + te, side, st - used)
        if best:
            ups[where].append(x[max(0, best[2] - U):best[2]].rjust(U, "-"))
            sides[where][best[1]] += 1

    def scan(seq):
        starts.append(seq[:KIT_WINDOW])
        before = {w: len(ups[w]) for w in ups}
        for where, x in (("start", seq), ("end", rc(seq))):
            (scan_exact if a.exact else scan_relaxed)(where, x)
        got = [w for w in ups if len(ups[w]) > before[w]]
        with_hit.update(got + (["any"] if got else []) + (["both"] if len(got) == 2 else []))

    def read_from(path, skip, take):
        """scan reads skip+1 .. skip+take of a file; returns (reads scanned, file ended)"""
        op = gzip.open if path.endswith(".gz") else open
        n = 0
        with op(path, "rt") as fh:
            for _ in range(4 * skip):
                if not fh.readline():
                    return 0, True
            while n < take:
                if not fh.readline():
                    return n, True
                seq = fh.readline().strip().upper()
                fh.readline(); fh.readline()
                scan(seq)
                n += 1
        return n, False

    # spread -n over all files (fastq_pass/ = every barcode folder and unclassified), not the first files in path
    # order; files that run out leave their share to the others in the next round
    done, live, seen = [0] * len(reads), list(range(len(reads))), 0
    while seen < a.n_reads and live:
        quota, keep = -(-(a.n_reads - seen) // len(live)), []
        for i in live:
            take = min(quota, a.n_reads - seen)
            if take <= 0:
                keep.append(i)
                continue
            n, ended = read_from(reads[i], done[i], take)
            done[i] += n
            seen += n
            if not ended:
                keep.append(i)
        live = keep
    n_used = sum(1 for d in done if d)
    rnd = random.Random(1)                   # hits in file order -> the --n-align / barcode subsets would come from
    for w in ups:                            # the first files again; shuffle (fixed seed: same result every run)
        rnd.shuffle(ups[w])
    hit_rule = ("exact match of utag + first 10 nt of the primer, anywhere in the read (--exact)" if a.exact else
                f"primer (first {HIT_PRIMER_NT} nt, <= {HIT_EDITS} edits) within {a.max_offset} nt of a read end + "
                f"utag (<= {HIT_EDITS} edits) right before it")
    R = {"hit_rule": hit_rule, "command": " ".join(sys.argv), "started": time.strftime("%Y-%m-%d %H:%M"), "cwd": os.getcwd(),
         "reads": [rel(f) for f in reads[:5]], "input": " ".join(a.reads), "n_files": len(reads), "files_used": n_used, "scanned": seen, "with_hit": dict(with_hit),
         "n_align": a.n_align, "upstream": U, "min_support": a.min_support, "tags": rel(a.tags),
         "primer_f": a.primer_f or PRIMER_F, "primer_r": a.primer_r or PRIMER_R,
         "known": [{"ctx": c, "name": nm, "seq": sq, "qlen": ql, "layout": BARCODE_LAYOUT.get(c)} for c, nm, sq, ql in KNOWN_CONTEXTS],
         "ends": {}}
    print_step("context: sequence upstream of the utags")
    print_info(f"scanned {fmt_int(seen)} reads from {fmt_int(n_used)} of {fmt_int(len(reads))} files "
               f"(up to {fmt_int(-(-a.n_reads // max(1, len(reads))))} per file); hit = {hit_rule}")
    print(style("    read start: the hit is in the read as sequenced, so the adapter before the utag is the start of the read\n"
                "    read end:   the hit is in the reverse complement of the read, so the utag is near the end of the read;\n"
                "                the adapter after it is shown reverse-complemented, in the same orientation as at the start", "gray"))
    print_info(f"reads with a hit: {fmt_int(with_hit['any'])} ({100 * with_hit['any'] / max(1, seen):.0f}%); "
               f"at the read start {fmt_int(with_hit['start'])}, at the read end {fmt_int(with_hit['end'])}, at both {fmt_int(with_hit['both'])}")
    if a.exact:
        print(style("    exact hits need ~20 error-free bases; at Q15 near the read start many reads have an error there", "gray"))
    print_info(f"the sequence before each hit is aligned (indels allowed) to the known adapters, ending at the utag; "
               f"all % below are of the first {fmt_int(a.n_align)} hits at that read end")
    dists = {}                   # where -> context -> [suffix_edits per hit]
    for where, label in (("start", "read start  (utag near the 5' end of the read)"),
                         ("end", "read end  (utag near the 3' end of the read, shown reverse-complemented)")):
        L = ups[where]
        n_side = ", ".join(f"{k} {fmt_int(v)}" for k, v in sorted(sides[where].items()))
        E = R["ends"][where] = {"hits": len(L), "sides": dict(sides[where])}
        print(f"\n  {style(label, 'bold')}  {style(fmt_int(len(L)) + ' hits', 'green' if len(L) >= 10 else 'red')}"
              + style(f"  ({n_side})" if n_side else "", "gray"))
        if len(L) < 10:
            print_note("too few hits: check --primer-f/--primer-r, the tag file and the sample sheet")
            E["too_few"] = True
            continue
        S = L[:a.n_align]
        E["used"] = len(S)
        print(style(f"    % = of {fmt_int(len(S))} hits" + (f" (first {fmt_int(len(S))} of {fmt_int(len(L))})"
                                                             if len(S) < len(L) else ""), "gray"))
        pc = lambda n: style(f"{100 * n / len(S):5.1f}%", pct_color(100 * n / len(S)))
        D = dists[where] = {ctx: [suffix_edits(seq, u) for u in S] for ctx, _, seq, _ in KNOWN_CONTEXTS}
        # 1. which known adapter is there: the whole adapter within ~15% edits
        rows = []
        for ctx, name, seq, _ in KNOWN_CONTEXTS:
            full = [d[-1] for d in D[ctx]]
            k = present_edits(len(seq))
            rows.append([style(ctx, "cyan"), name, str(len(seq)), pc(sum(e == 0 for e in full)),
                         *(pc(sum(e <= m for e in full)) for m in (1, 2, 3)), f"{pc(sum(e <= k for e in full))} (<= {k})"])
            E.setdefault("known", []).append({"ctx": ctx, "name": name, "len": len(seq), "k": k,
                                              "within": [sum(e <= m for e in full) / len(S) for m in (0, 1, 2, 3)],
                                              "present": sum(e <= k for e in full) / len(S)})
        print_table(["-c", "whole adapter 5' of the utag", "nt", "exact", "<=1", "<=2", "<=3 edits", "present"],
                    rows, "<<>>>>>>", indent=4)
        print(style("    hits with the whole adapter right before the utag, within that many edits; "
                    "present = within ~15% of its length", "gray"))
        # 2. for the adapter that fits best: how much of it is readable (suffixes ending at the utag)
        ctx, name, seq, qlen = max(KNOWN_CONTEXTS, key=lambda c: sum(context_ok(D[c[0]], len(c[2]))))
        if sum(context_ok(D[ctx], len(seq))) >= 0.3 * len(S) or sum(context_ok(D[ctx], 6)) >= 0.5 * len(S):
            print(f"\n    {style(ctx, 'cyan')} adapter: % of hits with its last n nt right before the utag, "
                  "within that many edits")
            rows = []
            E["suffix"] = {"ctx": ctx, "rows": []}
            for n in sorted({n for n in (6, 8, 12, 16, qlen, len(seq)) if n <= len(seq)}):
                note = "  (-c " + ctx + ")" if n == qlen else "  (whole)" if n == len(seq) else ""
                rows.append([f"last {n}{note}", color_seq(seq[-n:]),
                             *(pc(sum(d[n] <= m for d in D[ctx])) for m in (0, 1, 2, 3))])
                E["suffix"]["rows"].append({"n": n, "note": note.strip(" ()"), "seq": seq[-n:],
                                            "within": [sum(d[n] <= m for d in D[ctx]) / len(S) for m in (0, 1, 2, 3)]})
            print_table(["suffix", "sequence  <utag>", "exact", "<=1", "<=2", "<=3 edits"], rows, "<>>>>>", indent=4)
        # 3. per-position consensus: shows sequence beyond the known adapters (unknown kits)
        logo, cs, supp = base_profile(L)
        j = U
        while j > 0 and supp[j - 1] >= a.min_support:
            j -= 1
        top = Counter(x[-12:] for x in L).most_common(3)
        E.update(consensus=cs, support=supp, prefix=cs[j:], logo=logo, top12=top)
        print()
        print_kv([
            ("consensus", color_seq(cs) + "  " + style("<utag>", "bold", "magenta")),
            ("support", "".join(style(str(min(9, int(v // 10))), pct_color(v)) for v in supp)),
            ("consensus prefix", f"{color_seq(cs[j:]) if cs[j:] else '(none)'}  "
                                 + style(f"({U - j} nt >= {a.min_support:.0f}% at the same position; "
                                         "indels and short read ends shorten it)", "gray")),
            ("most common 12 nt", "   ".join(f"{color_seq(k)} {style(SYM['times'] + str(v), 'gray')}" for k, v in top)),
        ], indent=4)
    print(style("\n  consensus: majority base at each position before the utag ('.' = none); support digit = "
                "tens of % of all hits at that read end with that base (9 = 90-100%)", "gray"))
    if dists:
        R["verdict"] = context_verdict(dists, ups, a.n_align)
    else:
        rnd.shuffle(starts)
        R["kit"] = kit_from_starts(starts)
    if a.html:
        with open(a.html, "w", encoding="utf-8") as fh:
            fh.write(CONTEXT_HTML.replace("__DATA__", json.dumps(R, separators=(",", ":")).replace("</", "<\\/")))
        print()
        print_kv([("report", rel(a.html))])


def context_ok(D, n):
    """per hit: the last n nt of the adapter are there (within ~15% edits)"""
    k = present_edits(n)
    return [d[n] <= k for d in D]


def context_verdict(dists, ups, n_align):
    """Say which -c context the reads look like, whether its query fits at both read ends, and which ONT barcodes
    are there for nbd / mab. Returns the same as a dict for the HTML report."""
    print_step("context: suggestion")
    frac = lambda where, ctx, n: sum(context_ok(dists[where][ctx], n)) / len(dists[where][ctx])
    # the kit is decided at the read start (full adapter there); fall back to the read end
    where0 = "start" if "start" in dists else "end"
    ctx, name, seq, qlen = max(KNOWN_CONTEXTS, key=lambda c: frac(where0, c[0], len(c[2])))
    f0 = frac(where0, ctx, len(seq))
    if f0 < 0.3:
        print_info(f"no known adapter directly 5' of the utag (best: {ctx}, {100 * f0:.0f}% at the read {where0}): "
                   "the reads look trimmed (adapter/barcode removed) or the kit is unknown")
        print_kv([("use", style("-c bare", "bold") + "  (or --prefix from the consensus above)")], indent=4)
        return {"status": "none", "ctx": ctx, "frac": f0, "where": where0, "use": "utag demux -c bare",
                "alt": "utag demux -c bare --prefix <consensus prefix>"}
    n0 = len(dists[where0][ctx])
    print_ok(f"{name}: the whole adapter is right before the utag in {100 * f0:.0f}% of {fmt_int(n0)} hits at the "
             f"read {where0}  ->  kit context " + style("-c " + ctx, "bold"))
    # can the -c query (last qlen nt) be found at both read ends? the read end often stops inside the adapter
    V = {"status": "ok", "ctx": ctx, "name": name, "frac": f0, "where": where0, "hits": n0, "rows": [], "notes": []}
    rows, longest = [], None
    for n in sorted({n for n in (6, 8, 12, 16, qlen, len(seq)) if n <= len(seq)}):
        fr = {w: frac(w, ctx, n) for w in dists}
        V["rows"].append({"n": n, "query": n == qlen, "seq": seq[-n:], "edits": present_edits(n), **fr})
        rows.append([f"last {n}" + ("  (-c " + ctx + ")" if n == qlen else ""), color_seq(seq[-n:]), f"<= {present_edits(n)}",
                     *(style(f"{100 * fr[w]:5.1f}%", pct_color(100 * fr[w])) if w in fr else "-" for w in ("start", "end"))])
        if all(v >= 0.7 for v in fr.values()):
            longest = n
    print()
    hdr = [f"read {w} ({fmt_int(len(dists[w][ctx]))} hits)" if w in dists else f"read {w}" for w in ("start", "end")]
    print_table(["adapter suffix", "sequence  <utag>", "edits", *hdr], rows, "<>>>>", indent=4)
    print(style("    % of hits at that read end with that adapter suffix right before the utag (within the edits shown); "
                "the -c query should be >= 70% at both ends", "gray"))
    q_ok = all(frac(w, ctx, qlen) >= 0.7 for w in dists)
    use = f"utag demux -c {ctx}"
    alt = None
    if not q_ok and longest:
        V["notes"].append(f"the far end of most reads stops inside the adapter: only the last {longest} nt are there, not the "
                          f"{qlen} nt of the -c {ctx} query. barbell tolerates a partial flank, so -c {ctx} still works; "
                          f"--prefix {seq[-longest:]} matches both ends and may assign slightly more reads")
        print_note(V["notes"][-1])
        alt = f"utag demux --prefix {seq[-longest:]}"
    elif not q_ok:
        V["notes"].append("even a short adapter suffix is missing at one read end: try -c bare (utag + primer only)")
        print_note(V["notes"][-1])
        use = "utag demux -c bare"
    V.update(q_ok=q_ok, use=use, alt=alt)
    if ctx == "lsk":
        print_kv([("use", style(use, "bold"))] + ([("or", alt)] if alt else []), indent=4)
        return V
    # nbd / mab: which ONT barcodes are 5' of the flank (a subsample of the read starts where the flank is present)
    seqs, flank, kit = (NB_SEQS, NB_REAR, "SQK-NBD114-96") if ctx == "nbd" else (MAB_SEQS, MAB_TOPO, "SQK-MAB114-24")
    cand = (u for u in ups[where0] if suffix_edits(flank, u)[-1] <= present_edits(len(flank)))
    bc, nc = call_barcodes(seqs, cand, len(flank))
    called = sum(v for k, v in bc.items() if k != "?")
    top = print_barcodes(bc, nc, "hits")
    main = top[0] if top else None
    if ctx == "mab":
        use = f"utag twopass -k {kit}   (or utag demux -c mab --single-end)"
    elif len(top) > 1 and main[1] < 0.8 * called:
        use = f"utag twopass -k {kit}   (several NB barcodes)"
    elif main and q_ok:
        use = f"utag demux -c nbd --ont-barcode {main[0]}   (or plain -c nbd)"
    print_kv([("use", style(use, "bold"))], indent=4)
    V.update(use=use, barcodes=[[k, v / max(1, nc), seqs[k]] for k, v in top], checked=nc, set=list(seqs.items()),
             unclear=bc["?"] / max(1, nc))
    return V


def cmd_queries(a):
    tags, plates = read_tags(a.tags), load_layout(a)
    ctxs = a.context or CONTEXTS
    for c in ctxs:
        print_step(f"query FASTAs: {c}")
        build_queries(a.outdir, tag_rows(a, plates), tags, c, a.ont_barcode, a.prefix, a.all_tags,
                      primer_f=a.primer_f, primer_r=a.primer_r)


def cmd_tagqc(a):
    tags, plates = read_tags(a.tags), load_layout(a)
    fwd, rev = used_tags(plates, tags)
    print_step(f"tag set: {rel(a.tags)}")
    rows = []
    for nm, ids in (("fwd", fwd), ("rev", rev)):
        d = sorted((edit_distance(tags[x], tags[y]), x, y) for x, y in itertools.combinations(ids, 2))
        seqs = [tags[i] for i in ids]
        shared = len(os.path.commonprefix(seqs)), len(os.path.commonprefix([s_[::-1] for s_ in seqs]))
        n_min = sum(1 for x in d if x[0] == d[0][0])
        rows.append([style(nm, "cyan"), fmt_int(len(ids)), style(str(d[0][0]), "green" if d[0][0] >= 3 else "yellow" if d[0][0] == 2 else "red"),
                     f"{n_min} (e.g. {d[0][1]}/{d[0][2]})", f"{shared[0]} / {shared[1]}"])
    print_table(["set", "tags", "min edit dist", "pairs at min", "shared 5' / 3' nt"], rows, "<>>>>", indent=4)
    if any(len(os.path.commonprefix([tags[i] for i in ids])) for ids in (fwd, rev)):
        print()
        print_note("tags share a leading base: barbell folds it into the flank (the barcode region is 1 nt shorter)")
    plate_map_report(a, tags, plates)
    if a.print_map and os.path.exists(a.plates):
        print_plate_map(a.plates, a.print_map)


def id_ranges(ids):
    """['U001','U002','U003','U007'] -> 'U001-U003, U007'"""
    runs = []
    for i in sorted(set(ids), key=lambda x: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", x)]):
        m = re.fullmatch(r"(.*?)(\d+)", i)
        if runs and m and runs[-1][0] == m.group(1) and runs[-1][1] + 1 == int(m.group(2)):
            runs[-1][1], runs[-1][3] = int(m.group(2)), i
        else:
            runs.append([m.group(1) if m else None, int(m.group(2)) if m else -9, i, i])
    return ", ".join(r[2] if r[2] == r[3] else f"{r[2]}-{r[3]}" for r in runs)


def plate_map_report(a, tags, layout):
    """tagqc: layout of the plate map (which tag runs along rows / columns), consistency checks, sample sheet coverage."""
    print_step(f"plate map: {rel(a.plates)}")
    if not os.path.exists(a.plates):
        print_note("plate map not found (--plates)")
        return
    design = read_plates(a.plates)
    well_rc = lambda w: (w[0], int(w[1:]))
    by_plate = defaultdict(list)
    for r in design:
        by_plate[r["plate"]].append(r)
    rows = []
    for plate in sorted(by_plate, key=lambda p: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", p)]):
        prs = [r for r in by_plate[plate] if r["well"]]
        letters = sorted({well_rc(r["well"])[0] for r in prs})
        cols = sorted({well_rc(r["well"])[1] for r in prs})

        def along(axis):                    # which tag is constant along each row (axis 0) / column (axis 1)
            for side in ("utag-fwd", "utag-rev", "ont"):
                g = defaultdict(set)
                for r in prs:
                    g[well_rc(r["well"])[axis]].add(r[side])
                if all(len(v) == 1 and "" not in v for v in g.values()):
                    return f"{side[5:] or 'ONT'} {id_ranges(next(iter(v)) for v in g.values())}"
            return style("mixed", "yellow")
        rows.append([style(plate or "?", "cyan"), fmt_int(len(prs)),
                     f"{letters[0]}-{letters[-1]}" if letters else "", f"{cols[0]}-{cols[-1]}" if cols else "",
                     along(0), along(1)])
    print_table(["plate", "wells", "rows", "cols", "per row", "per column"], rows, "<>>><<", indent=4)

    d_fwd = {r["utag-fwd"] for r in design}
    d_rev = {r["utag-rev"] for r in design}
    pairs = {(r["utag-fwd"], r["utag-rev"]) for r in design}
    print()
    print_kv([("plates", fmt_int(len(by_plate))), ("wells", fmt_int(len(design))),
              ("fwd tags", f"{len(d_fwd)}  ({id_ranges(d_fwd)})"), ("rev tags", f"{len(d_rev)}  ({id_ranges(d_rev)})"),
              ("tag pairs", f"{fmt_int(len(pairs))} of {fmt_int(len(d_fwd) * len(d_rev))} possible fwd x rev combinations")])
    issues = []
    seen = Counter((r["plate"], r["well"]) for r in design if r["well"])
    dup = [f"{p}-{w}" for (p, w), n in seen.items() if n > 1]
    if dup:
        issues.append(f"wells listed more than once: {', '.join(dup[:10])}")
    if d_fwd & d_rev:
        issues.append(f"tags used as both fwd and rev: {id_ranges(d_fwd & d_rev)}")
    missing = (d_fwd | d_rev) - set(tags)
    if missing:
        issues.append(f"tags not in {os.path.basename(a.tags)}: {id_ranges(missing)}")
    print()
    for msg in issues:
        print_note(msg)
    if not issues:
        print_info(style(SYM["ok"], "green") + " no duplicate wells; fwd and rev tags are separate sets; all tags have a sequence")
    if len(pairs) == len(d_fwd) * len(d_rev):
        print_info("the grid is complete: every fwd x rev combination is a well, so chimeras land in a wrong well "
                   "(with a sample sheet, unlisted wells are reported as unexpected)")
    else:
        print_info(f"{fmt_int(len(d_fwd) * len(d_rev) - len(pairs))} fwd x rev combinations are not wells: "
                   "reads with them are reported as unexpected")
    unused = set(tags) - d_fwd - d_rev
    if unused:
        print_info(f"{len(unused)} tags in {os.path.basename(a.tags)} are not in the plate map: {id_ranges(unused)}")

    if getattr(a, "samplesheet", None):
        print_step(f"sample sheet: {rel(a.samplesheet)}")
        rows = []
        for plate in sorted(by_plate, key=lambda p: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", p)]):
            n = sum(1 for r in layout if r["plate"] == plate)
            if n:
                rows.append([style(plate, "cyan"), fmt_int(n), fmt_int(len(by_plate[plate]))])
        off = sum(1 for r in layout if not r["well"])
        if off:
            rows.append([style("not in plate map", "yellow"), fmt_int(off), ""])
        print_table(["plate", "samples", "wells in plate map"], rows, "<>>", indent=4)
        bins = Counter(r["ont"] for r in layout if r["ont"])
        print()
        kv = [("samples", fmt_int(len(layout)))]
        if bins:
            kv.append(("ONT bins", ", ".join(f"{b} ({n})" for b, n in sorted(bins.items()))))
        used = {(r["utag-fwd"], r["utag-rev"]) for r in layout}
        kv.append(("decoys", f"{len(d_fwd - {r['utag-fwd'] for r in layout})} fwd + {len(d_rev - {r['utag-rev'] for r in layout})} rev tags; "
                             f"{fmt_int(len(pairs - used))} plate-map wells not in the sheet are reported as unexpected"))
        print_kv(kv)


TAG_PALETTE = [196, 34, 33, 208, 129, 37, 201, 70, 63, 172, 31, 163,     # 256-colour codes, readable on dark and light
               160, 28, 25, 166, 91, 30, 125, 64, 99, 136, 39, 168]


def print_plate_map(path, mode):
    """tagqc --print: every plate-map well with its tags, listed column by column (A1, B1, ..), row by row
    (A1, A2, ..) or drawn as a plate grid (one line per tag in each plate row)."""
    design = [r for r in read_plates(path) if r["well"]]
    nat = lambda p: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", p)]
    well_rc = lambda w: (w[0], int(w[1:]))
    has_ont = any(r["ont"] for r in design)
    sides = [("fwd", "utag-fwd"), ("rev", "utag-rev")] + ([("ONT", "ont")] if has_ont else [])
    # one colour per tag ID (U001 is the same colour as fwd or rev, in the list and the grid); neighbouring IDs differ
    ids = sorted({r[k] for r in design for _, k in sides if r[k]}, key=nat)
    for c in TAG_PALETTE:
        ANSI.setdefault(f"c{c}", f"\033[38;5;{c}m")
    tag_code = {t: f"c{TAG_PALETTE[i % len(TAG_PALETTE)]}" for i, t in enumerate(ids)}
    tc = lambda t: style(t, tag_code.get(t, ""))
    print_step(f"plate map wells ({mode}): {rel(path)}")
    if mode != "plate":
        order = (lambda r: (nat(r["plate"]), well_rc(r["well"])[1], well_rc(r["well"])[0])) if mode == "columnwise" else \
                (lambda r: (nat(r["plate"]), well_rc(r["well"])))
        print_table(["plate", "well"] + [k if k != "ont" else "ont-barcode" for _, k in sides],
                    [[style(r["plate"], "cyan"), style(r["well"], "bold")] + [tc(r[k]) for _, k in sides]
                     for r in sorted(design, key=order)], indent=4)
        return
    by_plate = defaultdict(dict)
    for r in design:
        by_plate[r["plate"]][well_rc(r["well"])] = r
    for plate in sorted(by_plate, key=nat):
        wells = by_plate[plate]
        letters = [chr(c) for c in range(ord("A"), ord(max(l for l, _ in wells)) + 1)]
        cols = range(1, max(c for _, c in wells) + 1)
        w = max(len(r[k]) for r in wells.values() for _, k in sides)
        # box lines around the wells: one cell per well, one box row per plate row
        hz, vt = ("─", "│") if UNI else ("-", "|")
        rule = lambda l, m, r: "    " + " " * 9 + style(l + m.join(hz * (w + 2) for _ in cols) + r, "gray")
        corners = [("┌", "┬", "┐"), ("├", "┼", "┤"), ("└", "┴", "┘")] \
            if UNI else [("+", "+", "+")] * 3
        bar_ = style(vt, "gray")
        print()
        print("    " + cpad(plate, 9, "<", "bold", "cyan") + " " + " ".join(cpad(c, w + 2, "^", "bold") for c in cols))
        print(rule(*corners[0]))
        for j, l in enumerate(letters):
            for i, (lab, k) in enumerate(sides):
                cells = [pad(tc(wells[(l, c)][k]), w, ">") if (l, c) in wells else cpad(SYM["dot"], w, ">", "gray")
                         for c in cols]
                print("    " + cpad(l if i == 0 else "", 4, "<", "bold") + cpad(lab, 5, "<", "gray")
                      + bar_ + bar_.join(f" {x} " for x in cells) + bar_)
            print(rule(*corners[1 if j < len(letters) - 1 else 2]))


class Tee:
    """Write to the terminal and, without colour codes, to a log file."""
    def __init__(self, stream, log):
        self.stream, self.log = stream, log

    def write(self, text):
        self.stream.write(text)
        self.log.write(_ANSI_RE.sub("", text))
        return len(text)

    def flush(self):
        self.stream.flush()
        self.log.flush()

    def __getattr__(self, name):            # isatty, encoding, fileno, ... -> the terminal stream
        return getattr(self.stream, name)


def start_log(outdir, name="utag.log"):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    log = open(path, "w", buffering=1)
    log.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}  cwd: {os.getcwd()}\n"
              f"# {' '.join(sys.argv)}\n\n")
    sys.stdout, sys.stderr = Tee(sys.stdout, log), Tee(sys.stderr, log)
    return path


CONTEXT_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ONT/utag context report</title>
<style>
:root{color-scheme:light;--bg:#f9f9f7;--panel:#fcfcfb;--panel2:#f1f0ec;--ink:#0b0b0b;--dim:#52514e;--muted:#898781;--grid:#e1e0d9;
  --ring:rgba(11,11,11,.10);--s1:#2a78d6;--good:#006300;--warn:#a86400;--bad:#d03b3b;
  --q0:#f1f0ec;--q1:#cde2fb;--q2:#9ec5f4;--q3:#6da7ec;--q4:#3987e5;--q5:#256abf;
  --A:#1a9850;--C:#2a78d6;--G:#e6a100;--T:#d03b3b;--hb:#2a78d6;--hs:#6da7ec;--he:#1baf7a;--hn:#e1e0d9}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{padding:18px 24px 12px;border-bottom:1px solid var(--grid)}h1{margin:0;font-size:19px;font-weight:650}
.sub{color:var(--dim);font-size:12px;margin-top:8px;padding:6px 12px;background:#e7eef8;border:1px solid var(--ring);border-radius:8px;word-break:break-all}.sub a{color:inherit;text-decoration:underline}
main{max-width:1280px;margin:0 auto;padding:16px 24px 48px;display:grid;gap:16px}
.panel{background:var(--panel);border:1px solid var(--ring);border-radius:12px;padding:14px 16px 16px;min-width:0}
h2{margin:-14px -16px 12px;padding:9px 16px;background:#e7eef8;border-bottom:1px solid var(--ring);border-radius:12px 12px 0 0;
  font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:#1c3d66;font-weight:650}
h3{margin:14px 0 6px;font-size:12.5px;font-weight:600}
.hitsum{font-size:13px;color:var(--dim);margin-bottom:8px}.hitsum b{color:var(--ink);font-size:15px}
#hitcols{display:grid;grid-template-columns:max-content minmax(0,1fr);gap:16px 40px;align-items:start}
#bcs table,#kitbcs table{width:auto}#bcs .bcbar,#kitbcs .bcbar{width:140px}#bcs td,#kitbcs td{vertical-align:middle}
th.sort{cursor:pointer;user-select:none}th.sort:hover{color:var(--ink)}
.barsw{display:inline-block;margin-left:12px;font:12px system-ui;color:var(--dim);font-weight:normal}
.barsw button{font:12px system-ui;border:1px solid var(--ring);background:var(--panel);color:var(--dim);padding:1px 8px;cursor:pointer;margin-left:-1px}
.barsw button:first-child{margin-left:6px;border-radius:6px 0 0 6px}.barsw button:last-child{border-radius:0 6px 6px 0}
.barsw button.on{color:var(--ink);background:var(--panel2);font-weight:600}.pbar.rel i{background:#d4831a}
@media (max-width:1000px){#hitcols{grid-template-columns:1fr}}
.reads{display:grid;grid-template-columns:max-content 72px 320px max-content;gap:7px 12px;align-items:center;max-width:100%}
.reads .ends{display:flex;justify-content:space-between;font-size:11px;color:var(--muted);margin-bottom:-4px}
.rk{font-size:12.5px;color:var(--dim)}
.read{position:relative;height:12px;background:linear-gradient(#a9a89f,#a9a89f) center/100% 4px no-repeat}
.read .tag{position:absolute;top:0;bottom:0;width:12px;background:#d03b3b}
.pbar{position:relative;height:14px;background:var(--panel2);border:1px solid var(--grid)}
.pbar i{position:absolute;left:0;top:0;bottom:0;background:#2a78d6}.read .tag.l{left:0}.read .tag.r{right:0}
.tagkey{display:inline-block;width:14px;height:9px;border-radius:0;background:#d03b3b;vertical-align:-1px}
.rv{font-size:12.5px;white-space:nowrap}.rv b{font-weight:650}.rv em{font-style:normal;color:var(--muted);margin-left:4px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px}
.tile{background:var(--panel2);border-radius:10px;padding:8px 12px}.tile .v{font-size:20px;font-weight:650}.tile .k{color:var(--dim);font-size:12px}
.use{font:600 15px/1.4 ui-monospace,Menlo,monospace;background:var(--panel2);border-radius:8px;padding:8px 12px;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
.use button{margin-left:auto;font:12px system-ui;border:1px solid var(--ring);background:var(--panel);border-radius:6px;padding:2px 9px;cursor:pointer}
.verdict{font-size:14px;margin:0 0 10px}.verdict b{font-weight:650}
.ok{color:var(--good)}.no{color:var(--warn)}
.two{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:16px}
@media (max-width:1000px){.two{grid-template-columns:1fr}header,main{padding-left:16px;padding-right:16px}}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
th,td{padding:4px 8px;border-bottom:1px solid var(--grid);text-align:right;white-space:nowrap}
th{color:var(--dim);font-weight:600}th:first-child,td:first-child{text-align:left}
td.l,th.l{text-align:left}.seq{font-family:ui-monospace,Menlo,monospace;font-size:12px;letter-spacing:.02em}
.cell{display:inline-block;min-width:52px;border-radius:4px;padding:1px 6px;text-align:right}
.bar{display:inline-block;width:90px;height:8px;border-radius:4px;background:var(--q0);vertical-align:middle;margin-left:6px;overflow:hidden}
.bar i{display:block;height:100%;background:var(--s1);border-radius:4px}
.note{color:var(--dim);font-size:12px;line-height:1.5;margin:6px 0}
.warn{color:var(--dim)}.warn::before{content:"! ";color:var(--warn);font-weight:700}
.kv{display:grid;grid-template-columns:auto 1fr;gap:3px 14px;font-size:12.5px;margin-top:8px}.kv dt{color:var(--dim)}.kv dd{margin:0;word-break:break-all}
svg text{fill:var(--dim);font-size:10px}
.legend{display:flex;flex-wrap:wrap;gap:4px 12px;color:var(--dim);font-size:11.5px;margin:4px 0}.legend span{white-space:nowrap}
#ends{grid-auto-rows:auto;row-gap:0}
#ends>.panel{display:grid;grid-row:span 6;grid-template-rows:subgrid;row-gap:0;align-content:start}
#ends>.panel>div{min-width:0}
#ends table.tk,#ends table.ts{table-layout:fixed}#ends .tk td:nth-child(2){overflow:hidden;text-overflow:ellipsis}
#ends th,#ends td{padding:3px 4px;font-size:12px}#ends .ts .seq{font-size:11px;letter-spacing:0}#ends .tk td:nth-child(2){font-size:11.5px}
#ends .kv dd span{word-break:normal}
tr.pick td{background:#e7eef8}tr.pick td:first-child{border-radius:4px 0 0 4px}tr.pick td:last-child{border-radius:0 4px 4px 0}th.r{text-align:right!important}th.grp{text-align:center!important;border-bottom:1px solid var(--axis,#c3c2b7);padding-bottom:2px}th[rowspan]{vertical-align:bottom}
#ends .tk td,#ends .tk th{overflow:hidden;font-size:11.5px;padding:3px 3px}#ends .cell{min-width:44px;padding:1px 4px}#ends .bar{width:50px}.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:-1px}
details summary{cursor:pointer;color:var(--dim)}
</style></head><body>
<header><h1>ONT/utag context report</h1><div class="sub" id="sub"></div></header>
<main>
  <section class="panel"><h2>Hits</h2><div id="hitcols"><div id="tiles"></div><div id="bcs"><h3 style="margin-top:0">ONT barcodes</h3><div class="note">No ONT barcodes with utags found</div></div></div>
    <div class="note">hit: <span id="rule"></span></div></section>
  <div class="two" id="ends"></div>
  <section class="panel"><h2>Suggestion</h2><div id="verdict"></div></section>
  <section class="panel"><h2>Run</h2><dl class="kv" id="run"></dl></section>
</main>
<script>
const R = __DATA__;
const $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const fmt = n => Number(n).toLocaleString("en-US");
const pf = f => (100 * f).toFixed(1) + "%";
const qv = ["--q0", "--q1", "--q2", "--q3", "--q4", "--q5"];
const heat = f => { const i = f <= 0 ? 0 : Math.min(5, 1 + Math.floor(f * 5)); return `<span class="cell" style="background:var(${qv[i]});${i >= 4 ? "color:#fff" : ""}">${pf(f)}</span>`; };
const bar = f => `<span class="bar"><i style="width:${(100 * f).toFixed(1)}%"></i></span>`;
const cseq = s => [...s].map(c => `<span style="color:var(--${"ACGT".includes(c) ? c : "muted"})">${esc(c)}</span>`).join("");
$("rule").textContent = R.hit_rule;
$("sub").innerHTML = `Report generated on <b>${esc(R.started)}</b> by <a href="https://github.com/angelovangel/utag">utag</a>` +
  `<br><span class="seq">${esc(R.input)}</span>, ${fmt(R.scanned)} reads scanned from ${fmt(R.files_used)} of ${fmt(R.n_files)} files, ` +
  `${R.upstream} nt upstream, first ${fmt(R.n_align)} hits per read end aligned`;

// ---- suggestion
const V = R.verdict, E = R.ends;
let vh = "";
const K = R.kit;
if (!V && K) {
  vh = `<p class="verdict no">no utag hits${K.kit ? ` · read starts look like <b>${esc(K.kit)}</b>` : ""}</p><div class="note warn">${esc(K.verdict)}</div>` +
    `<h3>kit adapters at the read start</h3><div class="scroll"><table><tr><th class="l">kit</th><th class="l">read start</th><th>reads</th></tr>` +
    K.rows.map(([k, lay, f]) => `<tr class="${k === K.kit ? "pick" : ""}"><td class="l"><b>${esc(k)}</b></td><td class="l">${esc(lay)}</td><td>${heat(f)}</td></tr>`).join("") +
    `</table></div><div class="note">first ${K.window} nt of ${fmt(K.n)} reads, known adapters within ~15% edits</div>`;
  if (K.align) vh += `<h3>${esc(K.kit)} reads aligned at the insert</h3><div class="legend">${["A", "C", "G", "T"].map(b => `<span><i style="background:var(--${b})"></i>${b}</span>`).join("")}` +
    `<span style="white-space:normal">letter height = fraction of the ${fmt(K.align.reads)} ${esc(K.kit)} reads with that base</span></div>` +
    logo(K.align, {label: "kit", ctx: K.kit, layout: K.align.layout, qlen: 0}, K.barcodes || [], "insert");
  if (K.barcodes) vh += `<div id="kitbcs"></div>`;
}
else if (!V) vh = `<p class="verdict no">too few hits at both read ends: check --primer-f/--primer-r, the tag file and the sample sheet</p>`;
else if (V.status === "none") vh = `<p class="verdict no">no known adapter directly 5' of the utag (best: <b>${esc(V.ctx)}</b>, ${pf(V.frac)} at the read ${esc(V.where)}): the reads look trimmed (adapter/barcode removed) or the kit is unknown</p>`;
else vh = `<p class="verdict"><span class="ok">✓</span> <b>${esc(V.name)}</b>: the whole adapter is right before the utag in <b>${pf(V.frac)}</b> of ${fmt(V.hits)} hits at the read ${esc(V.where)} → kit context <b>-c ${esc(V.ctx)}</b></p>`;
if (V) {
  vh += `<div class="use"><span>${esc(V.use)}</span><button onclick="navigator.clipboard&&navigator.clipboard.writeText(this.previousElementSibling.textContent);this.textContent='copied'">copy</button></div>`;
  if (V.alt) vh += `<div class="note">or <span class="seq">${esc(V.alt)}</span></div>`;
  (V.notes || []).forEach(n => vh += `<div class="note warn">${esc(n)}</div>`);
  if (V.rows) {
    const ends = ["start", "end"].filter(w => E[w] && !E[w].too_few);
    vh += `<h3>${esc(V.ctx)} adapter suffix at both read ends</h3><div class="scroll"><table><tr><th>adapter suffix</th><th>sequence &lt;utag&gt;</th><th>edits</th>${ends.map(w => `<th>read ${w}</th>`).join("")}</tr>` +
      V.rows.map(r => `<tr><td>last ${r.n}${r.query ? ` <b>(-c ${esc(V.ctx)})</b>` : ""}</td><td class="seq">${cseq(r.seq)}</td><td>≤ ${r.edits}</td>${ends.map(w => `<td>${r[w] === undefined ? "–" : heat(r[w])}</td>`).join("")}</tr>`).join("") +
      `</table></div><div class="note">% of hits at that read end with the suffix right before the utag, within the edits shown; the -c query should be ≥ 70% at both ends</div>`;
  }
}
$("verdict").innerHTML = vh;
// ONT barcode tables, sortable by barcode or count (click the header again to reverse); top BC_TOP rows, rest on expand
const BT = {}, BC_TOP = 6;
// rows shown: "top" (BC_TOP most reads), "found" (all barcodes with reads) or "set" (the kit's whole set, 0 included)
function bcTable(id, rows, checked, unit, note, set) {
  const got = new Set(rows.map(r => r[0]));
  BT[id] = {rows, checked, unit, note, key: "n", dir: -1, view: "top", rel: false,
            set: rows.concat((set || []).filter(([b]) => !got.has(b)).map(([b, sq]) => [b, 0, sq]))};
  drawBc(id);
}
function viewBc(id, view) { BT[id].view = view; drawBc(id); }
function barBc(id, rel) { BT[id].rel = rel; drawBc(id); }       // bar = share of all (false) or of the top barcode (true)
function sortBc(id, key) { const t = BT[id]; t.dir = t.key === key ? -t.dir : key === "b" ? 1 : -1; t.key = key; drawBc(id); }
function drawBc(id) {
  const t = BT[id], top = Math.max(...t.rows.map(r => r[1]));
  const rows = [...(t.view === "set" ? t.set : t.rows)].sort((a, b) => t.dir * (t.key === "b" ? a[0].localeCompare(b[0], "en", {numeric: true}) : a[1] - b[1]));
  const th = (key, label, attr) => `<th class="sort${key === "n" ? " l" : ""}"${attr} onclick="sortBc('${id}','${key}')">${label}${t.key === key ? (t.dir > 0 ? " ▲" : " ▼") : ""}</th>`;
  const vw = (v, label) => `<button class="${t.view === v ? "on" : ""}" onclick="viewBc('${id}','${v}')">${label}</button>`;
  const sw = (rel, label) => `<button class="${t.rel === rel ? "on" : ""}" onclick="barBc('${id}',${rel})">${label}</button>`;
  $(id).innerHTML = `<h3 style="margin-top:${id === "bcs" ? 0 : ""}">ONT barcodes · ${rows.length}` +
    `<span class="barsw">show ${vw("top", "top " + BC_TOP)}${vw("found", "found (" + t.rows.length + ")")}${vw("set", "whole set (" + t.set.length + ")")}</span>` +
    `<span class="barsw">bar ${sw(false, "share of all")}${sw(true, "share of top")}</span></h3><table><tr>${th("b", "barcode", "")}<th class="l">sequence</th>${th("n", t.unit, ' colspan="2"')}</tr>` +
    rows.slice(0, t.view === "top" ? BC_TOP : rows.length).map(([b, f, sq]) => `<tr><td>${esc(b)}</td><td class="l seq">${cseq(sq || "")}</td>` +
      `<td class="l"><div class="pbar bcbar${t.rel ? " rel" : ""}"><i style="width:${(100 * (t.rel ? f / top : f)).toFixed(2)}%"></i></div></td>` +
      `<td class="l rv"><b>${pf(f)}</b> <em>${fmt(Math.round(f * t.checked))}</em></td></tr>`).join("") +
    `</table>` +
    `<div class="note">${t.note}</div>`;
}
if (V && V.barcodes && V.barcodes.length)      // shown next to the read hits
  bcTable("bcs", V.barcodes, V.checked, "hits", `of ${fmt(V.checked)} hits at the read start checked; ${pf(V.unclear)} unclear`, V.set);
if (!V && K && K.barcodes)
  bcTable("kitbcs", K.barcodes, K.checked, "reads", `of ${fmt(K.checked)} ${esc(K.kit)} reads checked; ${pf(K.unclear)} unclear`, K.set);

// ---- hits
const W = R.with_hit, sc = Math.max(1, R.scanned);
// reads in four groups that add up to 100%, each drawn as a read: utag markers at the ends with a hit,
// filled to the share of reads in that group
{ const both = W.both || 0, st = W.start || 0, en = W.end || 0, any = W.any || 0;
  const rows = [["both ends", both, 1, 1], ["read start only", st - both, 1, 0], ["read end only", en - both, 0, 1], ["no hit", sc - any, 0, 0]];
  $("tiles").innerHTML = `<div class="hitsum"><b>${fmt(R.scanned)}</b> reads scanned · <b>${fmt(any)}</b> with a hit (${pf(any / sc)})</div>` +
    `<div class="reads"><div></div><div class="ends"><span>5'</span><span>3'</span></div><div></div><div></div>` +
    rows.map(([k, n, l, r]) => { const f = n / sc;
      return `<div class="rk">${k}</div><div class="read">${l ? '<i class="tag l"></i>' : ""}${r ? '<i class="tag r"></i>' : ""}</div>` +
        `<div class="pbar" title="${k}: ${fmt(n)} reads (${pf(f)})"><i style="width:${(100 * f).toFixed(2)}%"></i></div>` +
        `<div class="rv"><b>${pf(f)}</b> <em>${fmt(n)}</em></div>`; }).join("") +
    `</div>`; }

// ---- per read end
function logo(e, known, bcs, end = "utag") {    // base frequencies per position before the utag (share of all hits)
  const lay = known && known.layout, U = e.logo.length, cw = 12, lw = 44, H = 96, Wd = lw + U * cw + 40, top = lay ? 24 : 8, base = top + H;
  const xi = p => lw + (U - p) * cw;      // left edge of position -p (p = 1 is the base right before the utag)
  let h = `<svg viewBox="0 0 ${Wd} ${base + 58}" width="100%" role="img" aria-label="base frequencies before the utag">`;
  // kit layout (nbd / mab): constant parts and the ONT barcode as bands behind the letters
  let bcSeq = "";
  if (lay) {
    bcs = bcs || (V && V.ctx === known.ctx && V.barcodes ? V.barcodes : []);
    const bc = bcs.length && bcs[0][1] >= 0.8 ? bcs[0] : null;           // one barcode: show its sequence
    const bcLab = bc ? ` (${bc[0]})` : bcs.length > 1 ? ` (mixed: ${bcs.slice(0, 3).map(b => b[0]).join(", ")})` : "";
    let p = 0, alt = 0;                       // alt: neighbouring constant parts in two shades
    [...lay].reverse().forEach(([lab, sq, n]) => {
      const p0 = p + 1, p1 = p + n; p = p1;                  // positions -p1 .. -p0
      bcSeq = (sq || (bc ? bc[2] : "·".repeat(n))) + bcSeq;
      if (p0 > U) return;
      const x0 = xi(Math.min(p1, U)), x1 = xi(p0) + cw, isBc = !sq;
      const fill = isBc ? "#fdf0d5" : alt++ % 2 ? "#e2e9f2" : "#eef2f7";
      h += `<rect x="${x0 + 1}" y="${top - 16}" width="${x1 - x0 - 2}" height="${H + 16}" fill="${fill}" rx="3"/>`;
      h += `<text x="${(x0 + x1) / 2}" y="${top - 5}" text-anchor="middle" style="fill:${isBc ? "#8a5a00" : "var(--dim)"};font-weight:${isBc ? 650 : 500}">${esc(lab)}${isBc ? esc(bcLab) : ""}${p1 > U ? " …" : ""}</text>`;
    });
  }
  for (const f of [0, .5, 1]) h += `<line x1="${lw}" x2="${lw + U * cw}" y1="${base - f * H}" y2="${base - f * H}" stroke="var(--grid)"/><text x="${lw - 4}" y="${base - f * H + 3}" text-anchor="end">${f * 100}%</text>`;
  e.logo.forEach((fr, i) => { let y = base;
    ["A", "C", "G", "T"].map((b, k) => [b, fr[k]]).sort((a, b) => a[1] - b[1]).forEach(([b, f]) => { const hh = f * H; if (hh < 1) return;
      h += `<text transform="translate(${lw + i * cw + 0.5} ${y.toFixed(2)}) scale(${(cw - 2) / 7.2} ${(hh / 9.4).toFixed(3)})" style="font:400 13px ui-monospace,Menlo,monospace;fill:var(--${b})">${b}<title>-${U - i}: ${b} ${pf(f)}</title></text>`; y -= hh; }); });
  h += `<text x="${lw + U * cw + 4}" y="${base - 2}" style="fill:var(--ink);font-weight:600">${end}</text>`;
  for (let i = 0; i < U; i += 10) h += `<text x="${lw + i * cw + cw / 2}" y="${base + 12}" text-anchor="middle">-${U - i}</text>`;
  const row = (y, lab, s, title) => { let t = `<text x="${lw - 4}" y="${y}" text-anchor="end">${lab}</text>`;
    [...s].forEach((c, i) => { const x = lw + (U - s.length + i) * cw + cw / 2; if (x < lw) return;
      t += `<text x="${x}" y="${y}" text-anchor="middle" style="font:600 11px ui-monospace,Menlo,monospace;fill:var(--${"ACGT".includes(c) ? c : "muted"})">${esc(c)}</text>`; });
    return t + (title ? `<text x="${lw + U * cw + 4}" y="${y}">${esc(title)}</text>` : ""); };
  h += row(base + 30, "cons.", e.consensus, "");
  if (known) {
    const q = Math.min(known.qlen || 0, U);
    if (q) h += `<rect x="${xi(q)}" y="${base + 37}" width="${q * cw}" height="15" rx="3" fill="#e7eef8"><title>-c ${esc(known.ctx)} query: last ${known.qlen} nt</title></rect>`;
    h += row(base + 48, known.label || "-c " + known.ctx, lay ? bcSeq : known.seq, q ? "query" : "");
  }
  return h + "</svg>";
}
$("ends").innerHTML = ["start", "end"].map(w => { const e = E[w]; if (!e) return "";
  const title = w === "start" ? "Read start <span style='text-transform:none;letter-spacing:0;font-weight:400'>(utag near the 5' end of the read)</span>" : "Read end <span style='text-transform:none;letter-spacing:0;font-weight:400'>(utag near the 3' end, shown reverse-complemented)</span>";
  const slot = [`<h2>${title}</h2>`], out = () => `<section class="panel">${[...slot, "", "", "", "", ""].slice(0, 6).map(x => x.startsWith("<h2") ? x : `<div>${x}</div>`).join("")}</section>`;
  slot.push(`<div class="note" style="margin-top:0">${fmt(e.hits)} hits (${Object.entries(e.sides).map(([k, v]) => k + " " + fmt(v)).join(", ")})${e.used ? `; % of the first ${fmt(e.used)} hits` : ""}</div>`);
  if (e.too_few) { slot.push(`<div class="note warn">too few hits: check --primer-f/--primer-r, the tag file and the sample sheet</div>`); return out(); }
  slot.push(`<h3>whole adapter 5' of the utag</h3><div><table class="tk"><colgroup><col style="width:7%"><col style="width:35%"><col style="width:6%"><col style="width:10.4%"><col style="width:10.4%"><col style="width:10.4%"><col style="width:10.4%"><col style="width:10.4%"></colgroup><tr><th rowspan="2">-c</th><th rowspan="2" class="l">adapter</th><th rowspan="2">nt</th><th rowspan="2">exact</th><th colspan="3" class="grp">edits</th><th rowspan="2">present</th></tr><tr><th class="r">≤1</th><th>≤2</th><th>≤3</th></tr>` +
    e.known.map(k => `<tr class="${V && V.status === "ok" && k.ctx === V.ctx ? "pick" : ""}"><td><b>${esc(k.ctx)}</b></td><td class="l" title="${esc(k.name)}">${esc(k.name)}</td><td>${k.len}</td>${k.within.map(f => `<td>${pf(f)}</td>`).join("")}<td title="within ${k.k} edits (~15% of its length)"><b>${pf(k.present)}</b></td></tr>`).join("") +
    `</table></div><div class="note">hits with the whole adapter right before the utag, within 0-3 edits; present = within ~15% of its length</div>`);
  slot.push(!e.suffix ? `<div class="note">no known adapter in enough hits for a suffix table</div>` : `<h3>${esc(e.suffix.ctx)} adapter: last n nt right before the utag</h3><div class="scroll"><table class="ts"><colgroup><col style="width:19%"><col style="width:41%"><col style="width:10%"><col style="width:10%"><col style="width:10%"><col style="width:10%"></colgroup><tr><th rowspan="2">suffix</th><th rowspan="2">sequence &lt;utag&gt;</th><th rowspan="2">exact</th><th colspan="3" class="grp">edits</th></tr><tr><th class="r">≤1</th><th>≤2</th><th>≤3</th></tr>` +
    e.suffix.rows.map(r => `<tr><td>last ${r.n}${r.note ? ` <span style="color:var(--dim)">(${esc(r.note)})</span>` : ""}</td><td class="seq">${cseq(r.seq)}</td>${r.within.map(f => `<td>${heat(f)}</td>`).join("")}</tr>`).join("") + "</table></div>");
  const known = e.suffix ? R.known.find(k => k.ctx === e.suffix.ctx) : null;
  slot.push(`<h3>base frequency, ${e.logo.length} nt 5' of the utag</h3><div class="legend">${["A", "C", "G", "T"].map(b => `<span><i style="background:var(--${b})"></i>${b}</span>`).join("")}<span style="white-space:normal">letter height = fraction of all hits with that base</span></div>${logo(e, known)}`);
  slot.push(`<dl class="kv"><dt>consensus prefix</dt><dd class="seq">${e.prefix ? cseq(e.prefix) : "(none)"} <span style="color:var(--dim);font-family:system-ui">${e.prefix.length} nt ≥ ${R.min_support}% at the same position; indels and short read ends shorten it</span></dd>` +
    `<dt>most common 12 nt</dt><dd class="seq">${e.top12.map(([s, n]) => `${cseq(s)} <span style="color:var(--dim)">×${fmt(n)}</span>`).join(" &nbsp; ")}</dd></dl>`);
  return out(); }).join("");

// ---- run
$("run").innerHTML = [["command", `<span class="seq">${esc(R.command)}</span>`], ["cwd", esc(R.cwd)], ["input", esc(R.reads.join(", ")) + (R.n_files > R.reads.length ? ` (+${R.n_files - R.reads.length} more)` : "")],
  ["tags", esc(R.tags)], ["primers", `<span class="seq">${esc(R.primer_f)}</span> / <span class="seq">${esc(R.primer_r)}</span>`],
  ["settings", `${fmt(R.scanned)} reads scanned from ${fmt(R.files_used)} of ${fmt(R.n_files)} files, ${R.upstream} nt upstream, first ${fmt(R.n_align)} hits per read end aligned`]]
  .map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");
</script></body></html>
"""


_OPT_RE = re.compile(r"(?<![\w-])(--?[a-zA-Z][\w-]*)")
_SEQ_RE = re.compile(r"\b[ACGTRYMKSWBDHVN]{12,}\b")
_DEF_RE = re.compile(r"\((default[^()\033]*)\)")


def _help_text(t):
    """options green, sequences by base, (default ...) gray"""
    t = _SEQ_RE.sub(lambda m: color_seq(m.group(0)), t)
    t = _OPT_RE.sub(lambda m: style(m.group(1), "green"), t)
    return _DEF_RE.sub(lambda m: style("(" + m.group(1) + ")", "gray"), t)


def color_help(text):
    """colour argparse help after layout (ANSI codes would break argparse's column widths)"""
    out, sec = [], ""
    for ln in text.split("\n"):
        if ln.startswith("usage:"):
            head, _, rest = ln.partition(" ")
            prog, _, rest = rest.partition(" [")
            out.append(style(head, "bold") + " " + style(prog, "bold") + (" [" + _OPT_RE.sub(lambda m: style(m.group(1), "green"), rest) if _ else ""))
            sec = "usage"
            continue
        m = re.match(r"^(\S.*):\s*$", ln)
        if m:
            sec = m.group(1)
            out.append(style(ln, "bold", "cyan"))
            continue
        if sec == "usage" and ln.startswith(" "):
            out.append(_OPT_RE.sub(lambda m: style(m.group(1), "green"), ln))
            continue
        if sec == "subcommands":
            m = re.match(r"^(    )(\w+)(\s{2,})(.*)$", ln)
            if m:
                out.append(m.group(1) + style(m.group(2), "bold", "green") + m.group(3) + m.group(4))
                continue
        if sec == "examples":
            m = re.match(r"^(  )(utag)( +)(\w+)(.*?)(\s+#.*)?$", ln)
            if m:
                out.append(m.group(1) + style(m.group(2), "bold") + m.group(3) + style(m.group(4), "bold", "green")
                           + m.group(5) + (style(m.group(6), "gray") if m.group(6) else ""))
                continue
        if sec.startswith("input file formats"):
            m = re.match(r"^(  )(--[\w-]+)(\s+)(\S.*?)(\s{2,})(\S.*)$", ln)
            if m:
                out.append(m.group(1) + style(m.group(2), "green") + m.group(3) + style(m.group(4), "yellow")
                           + m.group(5) + m.group(6))
                continue
            if re.match(r"^ {40,}\S", ln):                      # example rows under a format
                out.append(style(ln, "gray"))
                continue
        m = re.match(r"^(  )(-\S.*?)(\s{2,}|$)(.*)$", ln)      # option line: invocation, then help
        if m and sec not in ("examples", "subcommands") and not sec.startswith("input file formats"):
            inv = " ".join(style(t, "green") if t.startswith("-") else style(t, "yellow")
                           for t in m.group(2).split(" "))
            out.append(m.group(1) + inv + m.group(3) + _help_text(m.group(4)))
            continue
        out.append(_help_text(ln))
    return "\n".join(out)


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    def format_help(self):
        text = super().format_help()
        return color_help(text) if ENABLE_COLOR else text


def main():
    here = os.path.dirname(os.path.realpath(__file__))
    width = min(112, shutil.get_terminal_size((100, 24)).columns)
    for i, x in enumerate(sys.argv[1:], 1):          # --color also for -h (printed before the args are parsed)
        if x == "--color" and i + 1 < len(sys.argv):
            set_color(sys.argv[i + 1])
        elif x.startswith("--color="):
            set_color(x.split("=", 1)[1])
    raw = lambda prog: HelpFormatter(prog, max_help_position=32, width=width)
    ap = argparse.ArgumentParser(description=__doc__.format(PRIMER_F=PRIMER_F, PRIMER_R=PRIMER_R),
                                 epilog=EXAMPLES + "\n" + INPUT_HELP, formatter_class=raw)
    sp = ap.add_subparsers(dest="cmd", required=True, title="subcommands", metavar="<subcommand>",
                          description="run 'utag.py <subcommand> -h' for its options")

    def sub(name, help_, fn):
        p = sp.add_parser(name, help=help_, description=help_, epilog=INPUT_HELP, formatter_class=raw)
        p.set_defaults(fn=fn)
        return p

    def io(p, reads=True, outdir=None, kit=False):
        g = p.add_argument_group("input / output")
        if reads:
            g.add_argument("-r", "--reads", nargs="+", required=True, metavar="FASTQ",
                           help="FASTQ files, globs or folders (searched recursively; .fastq/.fq[.gz])")
        if kit:
            g.add_argument("-k", "--kit", help="barbell kit, e.g. SQK-NBD114-96, SQK-NBD114-24, SQK-MAB114-24 "
                                                "(required unless --dorado; with --dorado it only sets the SQK-MAB defaults)")
            g.add_argument("--pass1-only", action="store_true",
                           help="only split the reads by ONT barcode (barbell kit, barcode trimmed) into <outdir>/ont_bins/; "
                                "no utag demux")
            g.add_argument("--dorado", action="store_true",
                           help="reads already demultiplexed by Dorado / MinKNOW: the barcodeNN/ folders (or the barcode "
                                "in each read header) are the ONT bins and barbell kit (pass 1) is skipped; unclassified "
                                "reads are left out; the pass 2 context is detected (untrimmed NBD reads: -c nbd)")
        if outdir:
            g.add_argument("-o", "--outdir", default=outdir, metavar="DIR", help=f"output folder (default {outdir})")

    def library(p, context=None, single=None, prefix=True):
        g = p.add_argument_group("library")
        if context:
            g.add_argument("-c", "--context", **context)
        if single:
            g.add_argument("--single-end", action="store_true", help=single)
        if prefix:
            g.add_argument("--ont-barcode", metavar="NBxx", help="NB01..NB96 (-c nbd only): add that barcode to the flank")
            g.add_argument("--prefix", metavar="SEQ", help="constant sequence directly 5' of every utag (overrides -c)")
        g.add_argument("--primer-f", metavar="SEQ", help=f"forward primer 5'->3' (default 16S 27F {PRIMER_F})")
        g.add_argument("--primer-r", metavar="SEQ", help=f"reverse primer 5'->3' (default 16S 1492R {PRIMER_R})")

    def files(p):
        g = p.add_argument_group("input files (formats: end of this help)")
        g.add_argument("--tags", metavar="CSV", help=f"tag sequences (default {TAGS}; twopass -k SQK-MAB*: {MAB_TAGS})")
        g.add_argument("--plates", metavar="CSV", help=f"plate map, the tags in each well (default {PLATES}; twopass -k SQK-MAB*: {MAB_PLATES})")
        g.add_argument("--samplesheet", metavar="CSV|XLSX", help="sample names and wells (default: every plate-map well is a sample)")
        g.add_argument("--sheet-name", metavar="NAME", help="worksheet of an .xlsx sample sheet (default: first sheet)")
        g.add_argument("--no-decoys", action="store_true", help="with --samplesheet: do not search the other plate-map tags as decoys")
        g.add_argument("--all-tags", action="store_true", help="search all tags in --tags, not only those in the plate map")

    def output(p, skip_help):
        g = p.add_argument_group("output reads")
        g.add_argument("--flip", dest="no_flip", action="store_false", default=True,
                       help="reverse-complement reads into forward primer -> reverse primer orientation")
        g.add_argument("--gzip", action="store_true", help="write samples/*.fastq.gz")
        g.add_argument("--skip-trim", action="store_true", help=skip_help)
        return g

    def run(p, barbell=True):
        g = p.add_argument_group("run")
        if barbell:
            g.add_argument("-t", "--threads", type=int, default=10, metavar="N", help="barbell threads (default 10)")
            g.add_argument("--barbell", default="barbell", metavar="PATH", help="barbell binary (default: barbell on PATH)")
        g.add_argument("--color", choices=["auto", "always", "never"], default="auto",
                       help="coloured output (default auto: on a terminal only; NO_COLOR / FORCE_COLOR honoured)")

    # ---- demux
    d = sub("demux", "annotate -> inspect -> filter -> trim -> split reads into samples", cmd_demux)
    io(d, outdir="demux_out")
    library(d, context=dict(choices=CONTEXTS, default="bare",
                            help="sequence 5' of the utag: nbd, lsk, mab (end of the MAB TOPO site), or bare for reads "
                                 "already trimmed of ONT adapters (default bare)"),
            single="one utag at the read start (fwd or rev) names the sample: unique dual indexes, e.g. MAB; "
                   "the far read end is not used")
    files(d)
    g = d.add_argument_group("matching")
    g.add_argument("--flank-max-errors", type=int, metavar="N", help="relax if too many reads are missed (default: barbell decides)")
    g.add_argument("--min-score-diff", type=float, metavar="F",
                   help="min score difference between the best and second-best tag (default barbell, 0.1); raise for fewer wrong-well calls")
    g.add_argument("--max-offset", type=int, default=250, metavar="NT", help="max distance of the utag from the read end (default 250)")
    output(d, "write untrimmed reads to samples/: ONT adapter/barcode, utags and primers kept (barbell trim --skip-trim)") \
        .add_argument("--keep-merged", action="store_true", help="keep the decompressed/merged input FASTQ (deleted by default)")
    run(d)

    # ---- twopass
    t2 = sub("twopass", "barbell kit (ONT barcode) or the Dorado bins (--dorado) first, then utag demux inside each ONT bin", cmd_twopass)
    io(t2, outdir="twopass_out", kit=True)
    library(t2, context=dict(choices=CONTEXTS, help="pass 2 context 5' of the utag (default bare: pass 1 trimmed the ONT barcode)"),
            single="pass 2 uses only the utag at the read start, fwd or rev (default for SQK-MAB kits)", prefix=False)
    files(t2)
    g = t2.add_argument_group("matching")
    g.add_argument("--no-maximize", action="store_true", help="pass 1 without barbell kit --maximize")
    g.add_argument("--kit-flank-max-errors", type=int, metavar="N",
                   help="pass 1 (barbell kit) flank max errors (default: barbell decides, often 3). Raise to 4-5 if many "
                        "reads end up in no ONT bin (real MAB reads: 3 -> 53%%, 4 -> 71%%, 5 -> 86%% binned)")
    g.add_argument("--flank-max-errors", type=int, metavar="N", help="pass 2 (utag + primer flank) max errors (default: barbell decides)")
    g.add_argument("--min-score-diff", type=float, metavar="F",
                   help="min score difference between the best and second-best tag (default barbell, 0.1); raise for fewer wrong-well calls")
    g.add_argument("--max-offset", type=int, default=250, metavar="NT", help="max distance of the utag from the read end (default 250)")
    output(t2, "pass 2 does not trim: utags and primers stay on the reads in samples/ "
               "(pass 1, barbell kit, always removes the ONT adapter/barcode)")
    run(t2)

    # ---- context
    cx = sub("context", "find the kit context: align the sequence upstream of the utags to the known adapters", cmd_context)
    io(cx)
    library(cx, prefix=False)
    files(cx)
    g = cx.add_argument_group("scan")
    g.add_argument("-n", "--n-reads", type=int, default=20000, metavar="N", help="reads to scan, spread evenly over the input files (default 20000)")
    g.add_argument("--max-offset", type=int, default=250, metavar="NT", help="max distance of the utag from the read end (default 250)")
    g.add_argument("--exact", action="store_true",
                   help="hit = exact utag + first 10 nt of the primer anywhere in the read (default: primer <= 2 edits "
                        "within --max-offset of a read end, utag <= 2 edits right before it)")
    g.add_argument("--upstream", type=int, default=50, metavar="NT", help="nt upstream of the utag to report (default 50)")
    g.add_argument("--n-align", type=int, default=2000, metavar="N", help="hits per read end aligned to the known adapters (default 2000)")
    g.add_argument("--min-support", type=float, default=60, metavar="PCT", help="%% of hits per base for the consensus prefix (default 60)")
    g = cx.add_argument_group("report")
    g.add_argument("--html", default="context_report.html", metavar="FILE",
                   help="HTML report (default context_report.html; --html '' to skip)")
    run(cx, barbell=False)

    # ---- queries
    q = sub("queries", "write left/right query FASTAs for one or more contexts", cmd_queries)
    io(q, reads=False, outdir="fasta")
    library(q, context=dict(action="append", choices=CONTEXTS, help="context(s) to write, repeatable (default: all)"))
    files(q)
    run(q, barbell=False)

    # ---- tagqc
    t = sub("tagqc", "check the tag set, plate map layout and sample sheet coverage", cmd_tagqc)
    files(t)
    g = t.add_argument_group("print")
    g.add_argument("--print", dest="print_map", choices=["columnwise", "rowwise", "plate"],
                   help="also print every plate-map well with its tags: column by column (A1, B1, ..), "
                        "row by row (A1, A2, ..), or as a plate grid")
    library(t, prefix=False)
    run(t, barbell=False)

    a = ap.parse_args()
    # defaults that depend on the library: a MAB kit is read from one tag, so it needs the unique dual index design
    mab = a.cmd == "twopass" and (a.kit or "").upper().startswith("SQK-MAB")
    a.mab_defaults = mab and a.plates is None
    if a.plates is None:
        a.plates = os.path.join(here, MAB_PLATES if mab else PLATES)
    if a.tags is None:
        udi = os.path.basename(a.plates) == os.path.basename(MAB_PLATES)    # its tag IDs refer to MAB_TAGS
        a.tags = os.path.join(here, MAB_TAGS if udi else TAGS)
    set_color(a.color)
    if getattr(a, "outdir", None) and a.cmd != "queries":
        start_log(a.outdir)           # everything printed is also saved to <outdir>/utag.log
    try:
        a.fn(a)
    except SystemExit as e:
        if isinstance(e.code, str):          # sys.exit("message") -> coloured error, exit code 1
            print_error(e.code)
            sys.exit(1)
        raise
    except KeyboardInterrupt:
        print_error("interrupted")
        sys.exit(130)


if __name__ == "__main__":
    main()

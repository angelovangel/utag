#!/usr/bin/env python3
"""
utag.py - barbell query FASTAs + automated demux for custom utag amplicon
libraries on Oxford Nanopore: combinatorial (CDI) or unique (UDI) dual indexes.
Any primer pair works (--primer-f / --primer-r); defaults are 16S 27F / 1492R.

Library model (5'->3' of each tagged PCR primer):

    [ONT context] [utag] [primer]
      left  : ctx + utag-fwd + forward primer   (default 27F    AGRGTTYGATYMTGGCTCAG)
      right : ctx + utag-rev + reverse primer   (default 1492R  RGYTACCTTGTTACGACTT)

Both queries are written in oligo (5'->3') orientation, which is what barbell
expects for Ftag/Rtag FASTAs; barbell finds them on either read strand.

Only needs Python >= 3.8 (stdlib) and the `barbell` binary.
"""

INPUT_HELP = """\
input file formats (defaults live in data/ next to this script):
  --tags         data/utags-sequence.csv    no header: tag ID,sequence
                                            U013,CAACAGTAG
  --plates       data/utags-platemap.csv    plate,well,utag-fwd,utag-rev
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
import argparse, csv, glob, gzip, itertools, os, re, resource, shutil
import subprocess, sys, time
from collections import Counter, defaultdict

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


PRIMER_F = "AGRGTTYGATYMTGGCTCAG"   # 27F
PRIMER_R = "RGYTACCTTGTTACGACTT"    # 1492R

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
MAB_PREFIX = MAB_TOPO[-6:]                     # GCCCTT: directly 5' of the utag on untrimmed reads (-c mab)
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
            "unique dual indexes (each fwd and rev tag in one sample per ONT/AB barcode), e.g. --plates data/mab-udi-platemap.csv.\n"
            "A combinatorial (CDI) plate map like data/utags-platemap.csv cannot be demultiplexed from one tag.")


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
    txt = open(path, errors="replace").read()
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
    print_kv([("samples", f"{rel(os.path.join(out, 'samples'))}/<sample>{ext}"),
              ("summary", rel(summary)), ("unexpected", rel(unexp_tsv)), ("run summary", rel(run_csv)), ("patterns", rel(os.path.join(bb, 'inspect.txt'))),
              ("barbell files", rel(bb))])
    print()
    print_note("check barbell/inspect.txt - the two top patterns should be the two filter patterns.")


def cmd_twopass(a):
    """Pass 1: barbell kit (ONT barcode). Pass 2: barbell utag demux inside every ONT bin."""
    reads = expand_reads(a.reads)
    barbell = shutil.which(a.barbell) or sys.exit(f"barbell not found ({a.barbell}); install it first")
    mab = a.kit.upper().startswith("SQK-MAB")
    if a.mab_defaults:
        print_info(f"MAB kit: using the unique dual index design --plates {rel(a.plates)} --tags {rel(a.tags)}")
    rows = load_layout(a)
    ctx2 = a.context or "bare"                          # pass 1 trimmed the ONT barcode and its flanks
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

    # ---- pass 1: ONT barcodes
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
        bins[re.sub(r"\.trimmed\.fastq(\.gz)?$", "", os.path.basename(f))] = f
    if not bins:
        sys.exit("pass 1 produced no ONT barcode bins; check pass1_kit.log and --kit")
    print_info(f"{len(bins)} ONT bins: {', '.join(bins)}")

    # ---- pass 2: utags inside each bin
    has_ont_col = any(r["ont"] for r in rows)
    bin_rows, skipped = {}, []
    for b in bins:
        bin_rows[b] = [r for r in rows if (not has_ont_col) or r["ont"] == b]
        if not bin_rows[b]:
            skipped.append(b)
    if skipped:
        print_note(f"no samples listed for bin(s) {', '.join(skipped)}: skipped in pass 2")
    stats, per_sample, per_unexp = [], [], []
    unexp_tot = Counter()
    for b, f in bins.items():
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
        cmd = [sys.executable, os.path.abspath(__file__), "demux", "-r", f, "-c", ctx2,
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
        r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "NO_COLOR": "1"})
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
    n_bins = sum(st[1] for st in stats) + sum(count_fastq_reads(bins[b]) for b in skipped)
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
    print_kv([("samples", f"{rel(os.path.join(out, 'samples'))}/<ONT bin>/<sample>.fastq"),
              ("tables", f"{rel(os.path.join(out, 'ont_bins.tsv'))}, summary.tsv, unexpected.tsv, run_summary.csv"),
              ("barbell files", rel(bb))])


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


# known sequence directly 5' of the utag: (context, name, sequence, -c query length)
KNOWN_CONTEXTS = [("lsk", "Ligation Adapter (SQK-LSK114)", LSK_ADAPTER, len(LSK_ADAPTER_TAIL)),
                  ("nbd", "NB rear (SQK-NBD114)", NB_REAR, len(NB_REAR)),
                  ("mab", "TOPO site (SQK-MAB114)", MAB_TOPO, len(MAB_PREFIX))]


BARCODE_SAMPLE = 150     # hits checked for the NB / AB barcode (96 alignments each)


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
    seen = 0
    for path in reads:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            while seen < a.n_reads:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip().upper()
                fh.readline(); fh.readline()
                seen += 1
                before = {w: len(ups[w]) for w in ups}
                for where, x in (("start", seq), ("end", rc(seq))):
                    for side, ids, pat in (("fwd", fwd, pf), ("rev", rev, pr)):
                        for t in ids:
                            tg = tags[t]
                            i = x.find(tg)
                            while i >= 0:
                                if pat.match(x, i + len(tg)):
                                    ups[where].append(x[max(0, i - U):i].rjust(U, "-"))
                                    sides[where][side] += 1
                                i = x.find(tg, i + 1)
                got = [w for w in ups if len(ups[w]) > before[w]]
                with_hit.update(got + (["any"] if got else []))
        if seen >= a.n_reads:
            break
    print_step("context: sequence upstream of the utags")
    print_info(f"scanned {fmt_int(seen)} reads; a hit = an exact match of utag + first 10 nt of the primer")
    print(style("    read start: the hit is in the read as sequenced, so the adapter before the utag is the start of the read\n"
                "    read end:   the hit is in the reverse complement of the read, so the utag is near the end of the read;\n"
                "                the adapter after it is shown reverse-complemented, in the same orientation as at the start", "gray"))
    print_info(f"reads with a hit: {fmt_int(with_hit['any'])} ({100 * with_hit['any'] / max(1, seen):.0f}%); "
               f"at the read start {fmt_int(with_hit['start'])}, at the read end {fmt_int(with_hit['end'])}")
    print(style("    hits need ~20 error-free bases (utag + 10 nt of primer), so ~50% of reads is normal at Q15 near the "
                "read start; demux allows errors and assigns more", "gray"))
    print_info(f"the sequence before each hit is aligned (indels allowed) to the known adapters, ending at the utag; "
               f"all % below are of the first {fmt_int(a.n_align)} hits at that read end")
    dists = {}                   # where -> context -> [suffix_edits per hit]
    for where, label in (("start", "read start  (utag near the 5' end of the read)"),
                         ("end", "read end  (utag near the 3' end of the read, shown reverse-complemented)")):
        L = ups[where]
        n_side = ", ".join(f"{k} {fmt_int(v)}" for k, v in sorted(sides[where].items()))
        print(f"\n  {style(label, 'bold')}  {style(fmt_int(len(L)) + ' hits', 'green' if len(L) >= 10 else 'red')}"
              + style(f"  ({n_side})" if n_side else "", "gray"))
        if len(L) < 10:
            print_note("too few hits: check --primer-f/--primer-r, the tag file and the sample sheet")
            continue
        S = L[:a.n_align]
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
            for n in sorted({n for n in (6, 8, 12, 16, qlen, len(seq)) if n <= len(seq)}):
                note = "  (-c " + ctx + ")" if n == qlen else "  (whole)" if n == len(seq) else ""
                rows.append([f"last {n}{note}", color_seq(seq[-n:]),
                             *(pc(sum(d[n] <= m for d in D[ctx])) for m in (0, 1, 2, 3))])
            print_table(["suffix", "sequence  <utag>", "exact", "<=1", "<=2", "<=3 edits"], rows, "<>>>>>", indent=4)
        # 3. per-position consensus: shows sequence beyond the known adapters (unknown kits)
        cons, supp = [], []
        for pos in range(U):
            c = Counter(x[pos] for x in L if x[pos] != "-")
            n = sum(c.values())
            b, k = c.most_common(1)[0] if c else ("-", 0)
            cons.append(b if n and k / len(L) >= 0.5 else ".")
            supp.append(100 * k / len(L) if n else 0)
        cs = "".join(cons)
        j = U
        while j > 0 and supp[j - 1] >= a.min_support:
            j -= 1
        top = Counter(x[-12:] for x in L).most_common(3)
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
        context_verdict(dists, ups, a.n_align)


def context_ok(D, n):
    """per hit: the last n nt of the adapter are there (within ~15% edits)"""
    k = present_edits(n)
    return [d[n] <= k for d in D]


def context_verdict(dists, ups, n_align):
    """Say which -c context the reads look like, whether its query fits at both read ends, and which ONT barcodes
    are there for nbd / mab."""
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
        return
    n0 = len(dists[where0][ctx])
    print_ok(f"{name}: the whole adapter is right before the utag in {100 * f0:.0f}% of {fmt_int(n0)} hits at the "
             f"read {where0}  ->  kit context " + style("-c " + ctx, "bold"))
    # can the -c query (last qlen nt) be found at both read ends? the read end often stops inside the adapter
    rows, longest = [], None
    for n in sorted({n for n in (6, 8, 12, 16, qlen, len(seq)) if n <= len(seq)}):
        fr = {w: frac(w, ctx, n) for w in dists}
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
        print_note(f"the far end of most reads stops inside the adapter: only the last {longest} nt are there, not the "
                   f"{qlen} nt of the -c {ctx} query. barbell tolerates a partial flank, so -c {ctx} still works; "
                   f"--prefix {seq[-longest:]} matches both ends and may assign slightly more reads")
        alt = f"utag demux --prefix {seq[-longest:]}"
    elif not q_ok:
        print_note(f"even a short adapter suffix is missing at one read end: try -c bare (utag + primer only)")
        use = "utag demux -c bare"
    if ctx == "lsk":
        print_kv([("use", style(use, "bold"))] + ([("or", alt)] if alt else []), indent=4)
        return
    # nbd / mab: which ONT barcodes are 5' of the flank (a subsample of the read starts where the flank is present)
    seqs, flank, kit = (NB_SEQS, NB_REAR, "SQK-NBD114-96") if ctx == "nbd" else (MAB_SEQS, MAB_TOPO, "SQK-MAB114-24")
    cand = [u for u in ups[where0][:n_align] if suffix_edits(flank, u)[-1] <= present_edits(len(flank))][:BARCODE_SAMPLE]
    bc = best_barcode(seqs, cand, len(flank))
    called = sum(v for k, v in bc.items() if k != "?")
    top = [(k, v) for k, v in bc.most_common() if k != "?"][:6]
    print_kv([("ONT barcodes", ", ".join(f"{k} {100 * v / len(cand):.0f}%" for k, v in top) or "none identified"),
              ("", style(f"of {len(cand)} hits checked; {100 * bc['?'] / max(1, len(cand)):.0f}% unclear", "gray"))],
             indent=4)
    main = top[0] if top else None
    if ctx == "mab":
        use = f"utag twopass -k {kit}   (or utag demux -c mab --single-end)"
    elif len(top) > 1 and main[1] < 0.8 * called:
        use = f"utag twopass -k {kit}   (several NB barcodes)"
    elif main and q_ok:
        use = f"utag demux -c nbd --ont-barcode {main[0]}   (or plain -c nbd)"
    print_kv([("use", style(use, "bold"))], indent=4)


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


def main():
    here = os.path.dirname(os.path.realpath(__file__))
    raw = argparse.RawDescriptionHelpFormatter
    ap = argparse.ArgumentParser(description=__doc__, epilog=INPUT_HELP, formatter_class=raw)
    sp = ap.add_subparsers(dest="cmd", required=True, title="subcommands", metavar="<subcommand>",
                          description="run 'utag.py <subcommand> -h' for its options")

    def common(p):
        p.add_argument("--color", choices=["auto", "always", "never"], default="auto",
                       help="coloured output (default auto: only on a terminal; NO_COLOR / FORCE_COLOR are honoured)")
        g = p.add_argument_group("input files (formats and defaults: see the end of this help)")
        g.add_argument("--tags", metavar="CSV",
                       help="tag sequences (default: data/utags-sequence.csv; twopass -k SQK-MAB*: data/utags-10bp-sequence.csv)")
        g.add_argument("--plates", metavar="CSV",
                       help="plate map: the tags in each well (default: data/utags-platemap.csv; twopass -k SQK-MAB*: data/mab-udi-platemap.csv)")
        g.add_argument("--samplesheet", metavar="CSV|XLSX", help="your sample names and wells (without it every plate-map well is a sample)")
        g.add_argument("--sheet-name", help="worksheet name for an .xlsx sample sheet (default: first sheet)")
        g.add_argument("--no-decoys", action="store_true", help="with --samplesheet: do not search the plate map tags as decoys")
        p.add_argument("--primer-f", help=f"forward primer 5'->3' (default: 16S 27F {PRIMER_F})")
        p.add_argument("--primer-r", help=f"reverse primer 5'->3' (default: 16S 1492R {PRIMER_R})")
        p.add_argument("--all-tags", action="store_true", help="use all tags in --tags, not just those in the plate map")

    def ctx(p, default):
        p.add_argument("--ont-barcode", help="NB01..NB96 (nbd only): add that barcode to the flank")
        p.add_argument("--prefix", help="override: constant sequence placed directly 5' of every utag")

    q = sp.add_parser("queries", epilog=INPUT_HELP, formatter_class=raw, help="write left/right query FASTAs for one or more contexts")
    common(q); ctx(q, None)
    q.add_argument("-o", "--outdir", default="fasta")
    q.add_argument("-c", "--context", action="append", choices=CONTEXTS,
                   help="context(s) to write, repeatable (default: all)")
    q.set_defaults(fn=cmd_queries)

    d = sp.add_parser("demux", epilog=INPUT_HELP, formatter_class=raw, help="annotate -> inspect -> filter -> trim -> split reads into samples")
    common(d); ctx(d, "bare")
    d.add_argument("-r", "--reads", nargs="+", required=True, help="FASTQ file(s)/globs (.fastq or .fastq.gz)")
    d.add_argument("-o", "--outdir", default="demux_out")
    d.add_argument("-c", "--context", choices=CONTEXTS, default="bare",
                   help="sequence 5' of the utag: nbd, lsk, mab (end of the MAB TOPO site), or bare for reads already trimmed "
                        "of ONT adapters (default: bare)")
    d.add_argument("--single-end", action="store_true",
                   help="one utag at the read start (fwd or rev) names the sample: unique dual indexes, e.g. MAB; "
                        "the far read end is not used")
    d.add_argument("-t", "--threads", type=int, default=10)
    d.add_argument("--barbell", default="barbell")
    d.add_argument("--flank-max-errors", type=int, help="relax if too many reads are missed (default: barbell decides)")
    d.add_argument("--min-score-diff", type=float, help="min score difference between the best and second-best tag (default: barbell, 0.1); raise for fewer wrong-well calls")
    d.add_argument("--max-offset", type=int, default=250, help="max distance of utag from read end (default 250)")
    d.add_argument("--keep-merged", action="store_true", help="keep the decompressed/merged input FASTQ (deleted by default)")
    d.add_argument("--flip", dest="no_flip", action="store_false", default=True,
                   help="reverse-complement reads into forward-primer -> reverse-primer orientation")
    d.add_argument("--gzip", action="store_true")
    d.add_argument("--skip-trim", action="store_true",
                   help="write untrimmed reads to samples/ (ONT adapter/barcode, utags and primers kept; barbell trim --skip-trim)")
    d.set_defaults(fn=cmd_demux)

    t2 = sp.add_parser("twopass", epilog=INPUT_HELP, formatter_class=raw, help="barbell kit (ONT barcode) first, then utag demux inside each ONT bin")
    common(t2)
    t2.add_argument("-r", "--reads", nargs="+", required=True)
    t2.add_argument("-k", "--kit", required=True, help="barbell kit name, e.g. SQK-NBD114-96, SQK-NBD114-24, SQK-MAB114-24")
    t2.add_argument("-c", "--context", choices=CONTEXTS,
                    help="pass 2 context 5' of the utag (default: bare, pass 1 trimmed the ONT barcode)")
    t2.add_argument("--single-end", action="store_true",
                    help="pass 2 uses only the utag at the read start, fwd or rev (default for SQK-MAB kits)")
    t2.add_argument("-o", "--outdir", default="twopass_out")
    t2.add_argument("-t", "--threads", type=int, default=10)
    t2.add_argument("--barbell", default="barbell")
    t2.add_argument("--no-maximize", action="store_true", help="pass 1 without --maximize")
    t2.add_argument("--flank-max-errors", type=int, help="pass 2 (utag + primer flank) max errors (default: barbell decides)")
    t2.add_argument("--kit-flank-max-errors", type=int,
                    help="pass 1 (barbell kit) flank max errors (default: barbell decides, often 3). Raise to 4-5 if many "
                         "reads end up in no ONT bin: noisy read starts (on real MAB reads 3 -> 53%%, 4 -> 71%%, 5 -> 86%% binned)")
    t2.add_argument("--min-score-diff", type=float, help="min score difference between the best and second-best tag (default: barbell, 0.1); raise for fewer wrong-well calls")
    t2.add_argument("--max-offset", type=int, default=250)
    t2.add_argument("--flip", dest="no_flip", action="store_false", default=True,
                    help="reverse-complement reads into forward-primer -> reverse-primer orientation")
    t2.add_argument("--gzip", action="store_true")
    t2.add_argument("--skip-trim", action="store_true",
                    help="pass 2 does not trim: utags and primers stay on the reads in samples/ (pass 1, barbell kit, "
                         "always removes the ONT adapter/barcode)")
    t2.set_defaults(fn=cmd_twopass)

    cx = sp.add_parser("context", epilog=INPUT_HELP, formatter_class=raw, help="find the real kit context: align the sequence upstream of the utags to the known adapters")
    common(cx)
    cx.add_argument("-r", "--reads", nargs="+", required=True)
    cx.add_argument("-n", "--n-reads", type=int, default=20000, help="reads to scan (default 20000)")
    cx.add_argument("--upstream", type=int, default=50, help="nt upstream of the utag to report (default 50)")
    cx.add_argument("--n-align", type=int, default=2000, help="hits per read end aligned to the known adapters (default 2000)")
    cx.add_argument("--min-support", type=float, default=60, help="%% of hits needed per base for the consensus prefix (default 60)")
    cx.set_defaults(fn=cmd_context)

    t = sp.add_parser("tagqc", epilog=INPUT_HELP, formatter_class=raw, help="check the tag set, plate map layout and sample sheet coverage")
    common(t)
    t.set_defaults(fn=cmd_tagqc)

    a = ap.parse_args()
    # defaults that depend on the library: a MAB kit is read from one tag, so it needs the unique dual index design
    mab = a.cmd == "twopass" and a.kit.upper().startswith("SQK-MAB")
    a.mab_defaults = mab and a.plates is None
    if a.plates is None:
        a.plates = os.path.join(here, "data/mab-udi-platemap.csv" if mab else "data/utags-platemap.csv")
    if a.tags is None:
        udi = os.path.basename(a.plates) == "mab-udi-platemap.csv"       # its tag IDs refer to the 10 nt set
        a.tags = os.path.join(here, "data/utags-10bp-sequence.csv" if udi else "data/utags-sequence.csv")
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

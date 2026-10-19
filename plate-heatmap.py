#!/usr/bin/env python3
"""
plate-heatmap.py - plate heatmap of read counts per well, from the summary.tsv written by utag.py.

    plate-heatmap.py demux_out/summary.tsv -o plates.html
    plate-heatmap.py twopass_out/summary.tsv --scale linear --per-plate
    plate-heatmap.py demux_out/summary.tsv --palette viridis                # also cividis, magma, inferno, YlGnBu, Blues
    plate-heatmap.py old_summary.tsv --map plate_map.csv            # summaries without plate/well columns

Writes ONE self-contained HTML file (no internet, no dependencies) that shows every plate as a well grid coloured
by read count. In the page you mark negative controls (wells that should stay empty: an unused plate, NTCs) and wells
to exclude (not loaded, known dropouts); the sidebar shows sensitivity / specificity of calling a well at a read
cutoff, the cutoff for 100% specificity, AUC, reads in negative wells and which utags collect them. All other wells
are positives. The selection stays in the browser; results download as CSV.

Input (tab or comma separated; extra columns are ignored):
    sample        e.g. Plate1-A1 (plate and well are read from names like Plate1-A1, P2_B12, A1)
    reads         read count
    plate, well   where the sample sits (written by utag.py, also for custom sample names)
    ont_barcode   optional (twopass): every ONT bin gets its own plates
    utag_fwd, utag_rev   optional, shown in the tooltip and as well labels (utags)

Without plate/well columns, samples are placed by their name (Plate1-A1, P2_B12, A1) or with --map, a CSV/TSV with
columns `sample,well` and optionally `plate` (well like A1 or A01). Anything still unplaced is shown in an
"Unplaced samples" panel, where it can be selected like any other well.

Python >= 3.8, standard library only.
"""
import argparse
import csv
import json
import os
import re
import statistics
import sys

WELL_RE = re.compile(r"^(?:(?P<plate>.+?)[-_ :/.])?(?P<row>[A-Pa-p])(?P<col>\d{1,2})$")
ALIASES = {
    "sample": ("sample", "plate_well", "plate-well", "name", "well_id"),
    "reads": ("reads", "read_count", "count", "n_reads"),
    "ont": ("ont_barcode", "ont-barcode", "ont", "bin"),
    "fwd": ("utag_fwd", "utag-fwd", "fwd"),
    "rev": ("utag_rev", "utag-rev", "rev"),
}
USE_COLOR = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def paint(text, code):
    return f"\033[{code}m{text}\033[0m" if USE_COLOR else str(text)


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]


def read_table(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        head = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in head else ","
        rows = [{(k or "").strip().lower(): (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter=delim)]
    if not rows:
        sys.exit(f"{path}: empty")
    return rows


def pick(row, field):
    for k in ALIASES[field]:
        if k in row:
            return row[k]
    return None


def parse_well(text):
    m = re.fullmatch(r"([A-Pa-p])\s*0*(\d{1,2})", str(text).strip())
    if not m or int(m.group(2)) < 1 or int(m.group(2)) > 24:
        return None
    return ord(m.group(1).upper()) - 65, int(m.group(2)) - 1


def load_map(path):
    out = {}
    for r in read_table(path):
        s, w = pick(r, "sample"), r.get("well") or r.get("position")
        pos = parse_well(w) if w else None
        if not s or pos is None:
            sys.exit(f"{path}: need columns sample,well (e.g. A1) and optionally plate; bad row: {r}")
        out[s] = (r.get("plate", ""), pos)
    return out


def build(rows, sample_map):
    wells, unplaced_names = [], 0
    for r in rows:
        sample = pick(r, "sample")
        reads = pick(r, "reads")
        if sample is None or reads is None:
            sys.exit("summary needs a `sample` and a `reads` column (found: %s)" % ", ".join(rows[0]))
        try:
            reads = int(float(reads))
        except ValueError:
            sys.exit(f"bad read count for {sample}: {reads!r}")
        ont = pick(r, "ont") or ""
        plate, pos = None, None
        if sample in sample_map:
            plate, pos = sample_map[sample]
        elif r.get("well") and parse_well(r["well"]):        # summary.tsv from utag.py has plate,well
            plate, pos = r.get("plate", ""), parse_well(r["well"])
        else:
            m = WELL_RE.match(sample)
            if m and int(m.group("col")) <= 24 and 1 <= int(m.group("col")):
                plate, pos = m.group("plate") or "", (ord(m.group("row").upper()) - 65, int(m.group("col")) - 1)
        if pos is None:
            unplaced_names += 1
        wells.append(dict(id=f"{ont}|{sample}", sample=sample, ont=ont, plate=plate, pos=pos, reads=reads,
                          fwd=pick(r, "fwd") or "", rev=pick(r, "rev") or ""))
    ids = [w["id"] for w in wells]
    if len(set(ids)) != len(ids):
        sys.exit("duplicate sample names in the summary (within one ONT bin)")

    # panels: one per (ont bin, plate); unplaced samples of a bin share a list panel
    keys = sorted({(w["ont"], w["plate"] if w["pos"] is not None else None) for w in wells},
                  key=lambda k: (natural_key(k[0]), k[1] is None, natural_key(k[1] or "")))
    panels = []
    for i, (ont, plate) in enumerate(keys):
        members = [w for w in wells if w["ont"] == ont and (w["plate"] if w["pos"] is not None else None) == plate]
        if plate is None:
            title = " \u00b7 ".join(x for x in (ont, "Unplaced samples") if x)
            panels.append(dict(i=i, kind="list", title=title, rows=0, cols=0))
        else:
            mr = max(w["pos"][0] for w in members)
            mc = max(w["pos"][1] for w in members)
            rows_n = 8 if mr < 8 else 16
            cols_n = 12 if mc < 12 else 24
            title = " \u00b7 ".join(x for x in (ont, plate) if x) or "Plate"
            panels.append(dict(i=i, kind="plate", title=title, rows=rows_n, cols=cols_n))
        for w in members:
            w["panel"] = i
    out = []
    for w in wells:
        d = dict(id=w["id"], sample=w["sample"], ont=w["ont"], reads=w["reads"], fwd=w["fwd"], rev=w["rev"],
                 panel=w["panel"])
        if w["pos"] is not None:
            d.update(row=w["pos"][0], col=w["pos"][1], well=f"{chr(65 + w['pos'][0])}{w['pos'][1] + 1}")
        out.append(d)
    return out, panels, unplaced_names


def main():
    ap = argparse.ArgumentParser(description="Plate heatmap (HTML) of read counts per well from utag summary.tsv")
    ap.add_argument("summary", help="summary.tsv from utag.py demux / twopass (tab or comma separated)")
    ap.add_argument("-o", "--output", help="HTML file to write (default: plates.html next to the summary)")
    ap.add_argument("--title", help="page title (default: the summary's folder name)")
    ap.add_argument("--map", help="CSV/TSV with columns sample,well[,plate]; only needed for summaries without plate/well columns")
    ap.add_argument("--scale", choices=["log", "sqrt", "linear"], default="log", help="initial colour scale (default log)")
    ap.add_argument("--palette", choices=["YlOrRd", "viridis", "cividis", "magma", "inferno", "YlGnBu", "Blues"],
                    default="YlOrRd", help="initial colour palette (default YlOrRd; cividis is safest for colour-blind viewers)")
    ap.add_argument("--per-plate", action="store_true", help="start with the colour range set per plate instead of globally")
    ap.add_argument("--cutoff", type=int, default=10, help="initial read cutoff for calling a well (default 10)")
    a = ap.parse_args()

    rows = read_table(a.summary)
    wells, panels, unplaced = build(rows, load_map(a.map) if a.map else {})
    out = a.output or os.path.join(os.path.dirname(os.path.abspath(a.summary)), "plates.html")
    title = a.title or (os.path.basename(os.path.dirname(os.path.abspath(a.summary))) or "plates") + " - reads per well"
    data = dict(title=title, wells=wells, panels=panels, scale=a.scale, palette=a.palette, perPlate=a.per_plate, cutoff=a.cutoff)
    payload = json.dumps(data, separators=(",", ":")).replace("<", "\\u003c")
    html = TEMPLATE.replace("__TITLE__", title.replace("<", "&lt;")).replace("__DATA__", payload)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)

    total = sum(w["reads"] for w in wells)
    with_reads = sum(1 for w in wells if w["reads"] > 0)
    plates = sum(1 for p in panels if p["kind"] == "plate")
    lab = lambda s: paint(f"{s:<20}", "36")
    print(f"{paint('plate heatmap', '1;36')}  {paint(a.summary, '2')}")
    print(f"  {lab('plates')}{plates}" + (f"   (+ {unplaced} unplaced samples)" if unplaced else ""))
    print(f"  {lab('wells with reads')}{with_reads:,} / {len(wells):,}")
    print(f"  {lab('total reads')}{total:,}")
    if unplaced:
        print(f"  {paint('!', '1;33')} {unplaced} sample names are not plate-well names: they have no plate/well column; pass --map sample_wells.csv to place them")
    print(f"  {lab('written')}{out}")


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root{color-scheme:light;--bg:#f6f7f9;--panel:#ffffff;--panel2:#f1f3f6;--ink:#1c2230;--dim:#6b7383;--line:#dde1e7;--N:#d03b3b;--X:#6b7383;--P:#2a78d6;--accent:#1c2230;--zero:#eceef2}
*{box-sizing:border-box}
html,body{height:100%;margin:0}
body{background:var(--bg);color:var(--ink);font:13px/1.4 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;user-select:none;-webkit-user-select:none}
#app{display:grid;grid-template-columns:minmax(0,1fr) 350px;grid-template-rows:auto minmax(0,1fr);grid-template-areas:"head side" "main side";height:100vh}
header{grid-area:head;padding:14px 20px 10px;border-bottom:1px solid var(--line)}
header h1{margin:0 0 2px;font-size:17px;font-weight:650}
.sub{color:var(--dim);font-size:12px}
.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin-top:10px}
.controls label{display:flex;gap:6px;align-items:center;color:var(--dim)}
select,input,button{font:inherit;color:var(--ink);background:var(--panel2);border:1px solid var(--line);border-radius:6px;padding:4px 8px}
button{cursor:pointer}button:hover{border-color:#a9b0bc}
.legend{display:flex;align-items:center;gap:8px;margin-left:auto;color:var(--dim)}
.legend .bar{width:220px;height:12px;border-radius:6px;border:1px solid var(--line)}
.legend .zero{width:12px;height:12px;border-radius:2px;background:var(--zero);border:1px solid var(--line)}
main{grid-area:main;overflow:auto;padding:16px 20px 40px}
.panels{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:16px;align-items:start}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:12px 14px 14px;width:100%;justify-self:start;min-width:0}
.panel.p96{max-width:800px}.panel.p384{grid-column:1/-1;max-width:1400px}.panel.plist{grid-column:1/-1}
.phead{display:flex;flex-wrap:wrap;gap:6px 14px;align-items:baseline;margin-bottom:8px}
.phead b{font-size:14px}.pstats{color:var(--dim);font-size:12px}.pbtns{margin-left:auto;display:flex;gap:6px}
.pbtns button{padding:2px 8px;font-size:12px;color:var(--dim)}
.grid{display:grid;grid-template-columns:1.7em repeat(var(--cols),minmax(0,1fr));gap:3px;container-type:inline-size}
.hc,.hr{display:flex;align-items:center;justify-content:center;color:var(--dim);font-size:11px;cursor:pointer;border-radius:4px;min-height:1.4em}
.hc:hover,.hr:hover{background:var(--panel2);color:var(--ink)}
.well{position:relative;aspect-ratio:1;border-radius:2px;display:flex;align-items:center;justify-content:center;text-align:center;
  font-size:clamp(6.5px,calc(100cqw / (var(--cols) + 1.7) * .27),15px);font-weight:800;cursor:pointer;line-height:1;overflow:visible;
  font-variant-numeric:tabular-nums;transition:box-shadow .08s}
.well.empty{background:transparent;border:1px dashed var(--line);cursor:default}
.well:not(.empty):hover{box-shadow:0 0 0 2px #1c223088}
.well.selN{box-shadow:0 0 0 2px #fff,0 0 0 4px var(--N)}
.well.selX{opacity:.35;background-image:repeating-linear-gradient(45deg,transparent 0 3px,#fff 3px 5px)!important}
.well[data-g]::after{content:attr(data-g);position:absolute;top:-5px;right:-5px;min-width:13px;height:13px;border-radius:7px;font-size:9px;
  line-height:13px;text-align:center;color:#fff;font-weight:800;z-index:1}
.well.selN::after{background:var(--N)}.well.selX::after{background:var(--X)}
.tiles{display:flex;flex-wrap:wrap;gap:6px}
.well.tile{aspect-ratio:auto;border-radius:9px;padding:6px 10px;flex-direction:column;gap:2px;font-size:12px;min-width:92px}
.tile .tn{font-weight:500;opacity:.85;font-size:11px;word-break:break-all}.tile .tv{font-weight:700}
.well .ut{font-size:.74em;font-weight:500;line-height:1.1}
aside{grid-area:side;border-left:1px solid var(--line);background:var(--panel);overflow:auto;padding:14px 16px 30px}
aside h2{margin:14px 0 8px;font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--dim);font-weight:600}
aside h2:first-child{margin-top:0}
.sh{display:flex;align-items:center;gap:8px}.sh button{margin-left:auto;padding:2px 9px;font-size:12px;color:var(--dim);text-transform:none;letter-spacing:0;font-weight:400}
#app.wide aside{position:fixed;top:0;right:0;bottom:0;width:50vw;min-width:350px;z-index:5;box-shadow:-10px 0 28px #1c223026;padding:16px 24px 30px}
.modes{display:flex;gap:6px}
.modes button{flex:1;padding:7px 4px;border-radius:8px;display:flex;align-items:center;justify-content:center;gap:6px}
.modes button.on{border-color:var(--c);box-shadow:inset 0 0 0 1px var(--c);background:#1c22300d}
.modes .chip{width:10px;height:10px;border-radius:50%;background:var(--c)}
.modes kbd{color:var(--dim);font-size:11px}
.seg{display:flex;align-items:center;gap:6px;margin-top:8px;color:var(--dim);font-size:12px}
.seg button{padding:3px 10px;font-size:12px}.seg button.on{border-color:var(--accent);color:var(--ink);background:#1c22300d}
.gcard{border:1px solid var(--line);border-left:3px solid var(--c);border-radius:8px;background:var(--panel2);padding:6px 10px;margin-bottom:6px}
.gh{display:flex;gap:6px;align-items:center}
.gh button{padding:1px 7px;font-size:11px;color:var(--dim)}
.gh input{flex:1;min-width:0;background:transparent;border-color:transparent;font-weight:600;font-size:12.5px;padding:1px 3px}
.gh input:focus{border-color:var(--line);background:var(--bg);user-select:text;-webkit-user-select:text}
.big{font-size:15px;font-weight:700;font-variant-numeric:tabular-nums;color:var(--c);white-space:nowrap}
.big small{font-size:11px;font-weight:400;color:var(--dim)}
.gcard .s{color:var(--dim);font-size:12px}
dl{display:grid;grid-template-columns:1fr auto;gap:0 10px;margin:3px 0 0;font-size:11.5px}
.gcard dl{grid-template-columns:auto 1fr auto 1fr}
dt{color:var(--dim)}dd{margin:0;text-align:right;font-variant-numeric:tabular-nums}
.cmp{border:1px solid var(--line);border-radius:8px;padding:4px 10px;margin-bottom:6px}.cmp dl{margin:0}
.row{display:flex;gap:6px;align-items:center;flex-wrap:wrap;margin:6px 0}
.row input[type=number]{width:96px}
.btns{display:flex;gap:6px;flex-wrap:wrap}
.list{max-height:240px;overflow:auto;border:1px solid var(--line);border-radius:8px}
.item{display:flex;gap:8px;align-items:center;padding:3px 8px;border-bottom:1px solid var(--line);cursor:pointer;font-size:12px}
.item:last-child{border:0}.item:hover{background:var(--panel2)}
.item i{width:8px;height:8px;border-radius:50%;flex:none}.item span{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.item em{font-style:normal;color:var(--dim);font-variant-numeric:tabular-nums}
.hint{color:var(--dim);font-size:11.5px;margin-top:12px;line-height:1.55}
.res{border:1px solid var(--line);border-radius:8px;padding:6px 10px;margin-bottom:6px}
.res dl{grid-template-columns:1fr auto}.res dd b{font-size:13px}
.note{color:var(--dim);font-size:11.5px;line-height:1.5;margin:4px 0 6px}
.mini{width:100%;border-collapse:collapse;font-size:11.5px;font-variant-numeric:tabular-nums}
.mini td,.mini th{padding:2px 4px;border-bottom:1px solid var(--line);text-align:right}.mini th{color:var(--dim);font-weight:600}
.mini td:first-child,.mini th:first-child{text-align:left}.mini button{padding:0 6px;font-size:11px}
.strip text{fill:var(--dim);font-size:10px}
#tip{position:fixed;z-index:9;pointer-events:none;background:#fffffff5;border:1px solid var(--line);box-shadow:0 4px 14px #1c223022;color:var(--ink);border-radius:8px;padding:7px 10px;font-size:12px;display:none;max-width:280px}
#tip b{font-size:13px}#tip .k{color:var(--dim)}
@media (max-width:900px){#app{grid-template-columns:1fr;grid-template-areas:"head" "main" "side";height:auto}aside{border-left:0;border-top:1px solid var(--line)}
  #app.wide aside{position:static;width:auto;box-shadow:none}.sh button{display:none}}
</style>
</head>
<body>
<div id="app">
  <header>
    <h1 id="title"></h1>
    <div class="sub" id="overall"></div>
    <div class="controls">
      <label>colour scale <select id="scale"><option value="log">log</option><option value="sqrt">sqrt</option><option value="linear">linear</option></select></label>
      <label>palette <select id="palette"></select></label>
      <label>range <select id="range"><option value="global">all plates</option><option value="plate">per plate</option></select></label>
      <label>well labels <select id="labels"><option value="compact">compact</option><option value="exact">exact</option><option value="utags">utags</option><option value="none">none</option></select></label>
      <div class="legend"><span class="zero" title="0 reads"></span><span>0</span><span id="lo"></span><div class="bar" id="bar"></div><span id="hi"></span></div>
    </div>
  </header>
  <main><div class="panels" id="plates"></div></main>
  <aside id="side"></aside>
</div>
<div id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
(() => {
const D = JSON.parse(document.getElementById('data').textContent);
const wells = D.wells, panels = D.panels;
const byId = new Map(wells.map(w => [w.id, w]));
const KEY = 'plate-heatmap:neg:' + D.title + ':' + wells.length;
// groups: N = negative control (should have no reads), X = excluded (counts nowhere); every other well is a positive
const state = {mode: 'N', drag: 'box', scale: D.scale || 'log', palette: D.palette || 'YlOrRd', range: D.perPlate ? 'plate' : 'global', labels: 'compact', groups: {}, cutoff: D.cutoff || 10, wide: false};
try {
  const s = JSON.parse(localStorage.getItem(KEY) || 'null');
  if (s) { for (const [id, g] of Object.entries(s.groups || {})) if (byId.has(id) && (g === 'N' || g === 'X')) state.groups[id] = g; if (s.cutoff) state.cutoff = s.cutoff; }
} catch (e) {}
const save = () => { try { localStorage.setItem(KEY, JSON.stringify({groups: state.groups, cutoff: state.cutoff})); } catch (e) {} };

const $ = s => document.querySelector(s);
const esc = s => String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt = n => Math.round(n).toLocaleString('en-US');
const compact = n => {  // 950, 1.2k, 12k, 120k, 1.2M, 12M
  if (n < 1000) return String(n);
  for (const [d, u] of [[1e3, 'k'], [1e6, 'M'], [1e9, 'G']]) {
    const x = n / d, r = x < 10 ? Math.round(x * 10) / 10 : Math.round(x);
    if (r < 1000 || u === 'G') return String(r) + u;
  }
};
const sum = a => a.reduce((x, y) => x + y, 0);
const median = a => { if (!a.length) return 0; const s = a.slice().sort((x, y) => x - y), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const GC = {N: 'var(--N)', X: 'var(--X)'};
const gname = g => g === 'N' ? 'negative' : g === 'X' ? 'excluded' : 'positive';

// ---------- colour: sequential palettes only (read counts are magnitudes), low -> high
// matplotlib (perceptually uniform; cividis also reads under colour-blindness) and ColorBrewer ramps.
// Light-starting ramps begin one step in, so low counts stay distinct from the grey of 0-read wells.
const PALETTES = {
  YlOrRd:  ['#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026'],
  viridis: ['#440154', '#472d7b', '#3b528b', '#2c728e', '#21918c', '#28ae80', '#5ec962', '#addc30', '#fde725'],
  cividis: ['#00224e', '#123570', '#3b496c', '#575d6d', '#707173', '#8a8779', '#a69d75', '#c4b56c', '#fee838'],
  magma:   ['#000004', '#1c1044', '#4f127b', '#812581', '#b5367a', '#e55964', '#fb8761', '#fec287', '#fcfdbf'],
  inferno: ['#000004', '#1f0c48', '#550f6d', '#88226a', '#ba3655', '#e35933', '#f98e09', '#f9cb35', '#fcffa4'],
  YlGnBu:  ['#edf8b1', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#253494', '#081d58'],
  Blues:   ['#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'],
};
const hex2rgb = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
const PAL_RGB = Object.fromEntries(Object.entries(PALETTES).map(([k, v]) => [k, v.map(hex2rgb)]));
if (!PAL_RGB[state.palette]) state.palette = 'YlOrRd';
function cmap(t) {
  const P = PAL_RGB[state.palette], n = P.length - 1;
  t = Math.max(0, Math.min(1, t)) * n;
  const i = Math.min(n - 1, Math.floor(t)), f = t - i;
  return P[i].map((v, k) => Math.round(v + (P[i + 1][k] - v) * f));
}
const rgb = c => `rgb(${c[0]},${c[1]},${c[2]})`;
const TF = {log: v => Math.log10(v + 1), sqrt: v => Math.sqrt(v), linear: v => v};
const norm = (v, lo, hi) => { const f = TF[state.scale] || TF.log, a = f(lo), b = f(hi); return b > a ? (f(v) - a) / (b - a) : 1; };
const rng = a => a.length ? {lo: Math.min(...a), hi: Math.max(...a)} : {lo: 0, hi: 1};
const globalR = rng(wells.filter(w => w.reads > 0).map(w => w.reads));
const panelR = panels.map(p => rng(wells.filter(w => w.panel === p.i && w.reads > 0).map(w => w.reads)));
const rangeFor = w => state.range === 'plate' ? panelR[w.panel] : globalR;
function colorFor(w) {
  if (w.reads <= 0) return {bg: 'var(--zero)', fg: '#9aa1ad'};
  const r = rangeFor(w), c = cmap(norm(w.reads, r.lo, r.hi));
  return {bg: rgb(c), fg: (0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]) > 140 ? '#111' : '#fff'};
}
const label = w => state.labels === 'none' ? '' : state.labels === 'exact' ? fmt(w.reads) : compact(w.reads);
const utags = (w, sep) => [w.fwd, w.rev].filter(Boolean).map(esc).join(sep);

// ---------- build plates
function panelHTML(p) {
  const ws = wells.filter(w => w.panel === p.i);
  const tot = sum(ws.map(w => w.reads)), pos = ws.filter(w => w.reads > 0).length;
  const stats = `${fmt(tot)} reads \u00b7 ${pos}/${ws.length} wells with reads \u00b7 median ${fmt(median(ws.map(w => w.reads)))}`;
  let h = `<section class="panel ${p.kind === 'list' ? 'plist' : p.cols > 12 ? 'p384' : 'p96'}" data-p="${p.i}"><div class="phead"><b>${esc(p.title)}</b><span class="pstats">${stats}</span>` +
    `<span class="pbtns"><button data-act="selall">mark all</button><button data-act="clr">clear</button></span></div>`;
  if (p.kind === 'list') {
    h += '<div class="tiles">' + ws.map(w => `<div class="well tile" data-id="${esc(w.id)}"></div>`).join('') + '</div>';
  } else {
    const cells = new Map(ws.map(w => [w.row + '_' + w.col, w]));
    h += `<div class="grid" style="--cols:${p.cols}"><div></div>`;
    for (let c = 0; c < p.cols; c++) h += `<div class="hc" data-col="${c}">${c + 1}</div>`;
    for (let r = 0; r < p.rows; r++) {
      h += `<div class="hr" data-row="${r}">${String.fromCharCode(65 + r)}</div>`;
      for (let c = 0; c < p.cols; c++) { const w = cells.get(r + '_' + c); h += w ? `<div class="well" data-id="${esc(w.id)}"></div>` : '<div class="well empty"></div>'; }
    }
    h += '</div>';
  }
  return h + '</section>';
}
$('#plates').innerHTML = panels.map(panelHTML).join('');

function decorate(el) {
  const w = byId.get(el.dataset.id), c = colorFor(w), g = state.groups[w.id];
  el.style.background = c.bg; el.style.color = c.fg;
  if (el.classList.contains('tile')) el.innerHTML = `<span class="tn">${esc(w.sample)}</span><span class="tv">${state.labels === 'utags' ? utags(w, '/') : state.labels === 'exact' || state.labels === 'none' ? fmt(w.reads) : compact(w.reads)}</span>`;
  else if (state.labels === 'utags') el.innerHTML = `<span class="ut">${utags(w, '<br>')}</span>`;
  else el.textContent = label(w);
  el.classList.toggle('selN', g === 'N'); el.classList.toggle('selX', g === 'X');
  if (g) el.dataset.g = g === 'N' ? 'N' : '\u00d7'; else delete el.dataset.g;
}
const decorateAll = () => document.querySelectorAll('.well[data-id]').forEach(decorate);
const elMap = new Map(); document.querySelectorAll('.well[data-id]').forEach(e => elMap.set(e.dataset.id, e));
const elOf = id => elMap.get(id) || null;
const panelWells = panels.map(p => wells.filter(w => w.panel === p.i));

// ---------- header
$('#title').textContent = D.title;
const allReads = sum(wells.map(w => w.reads));
$('#overall').textContent = `${fmt(allReads)} reads \u00b7 ${wells.filter(w => w.reads > 0).length}/${wells.length} wells with reads \u00b7 ${panels.length} panel${panels.length > 1 ? 's' : ''}`;
$('#palette').innerHTML = Object.keys(PALETTES).map(k => `<option>${k}</option>`).join('');
$('#scale').value = state.scale; $('#range').value = state.range; $('#palette').value = state.palette;
function legend() {
  const stops = Array.from({length: 17}, (_, i) => rgb(cmap(i / 16)) + ' ' + (i / 16 * 100).toFixed(0) + '%').join(',');
  $('#bar').style.background = `linear-gradient(90deg,${stops})`;
  $('#lo').textContent = state.range === 'plate' ? 'min' : fmt(globalR.lo);
  $('#hi').textContent = state.range === 'plate' ? 'max (per plate)' : fmt(globalR.hi);
  $('#bar').title = {log: 'log10 scale', sqrt: 'square-root scale', linear: 'linear scale'}[state.scale];
}
for (const id of ['scale', 'palette', 'range', 'labels']) $('#' + id).addEventListener('change', e => { state[id] = e.target.value; legend(); decorateAll(); });

// ---------- selection
let paintTo = null, painting = false, rafSide = 0, drag = null;
function assign(w, g) {
  if (g) state.groups[w.id] = g; else delete state.groups[w.id];
  const el = elOf(w.id); if (el) decorate(el);
}
const modeGroup = () => state.mode === 'E' ? null : state.mode;
function applyEl(el) {                      // freehand: paint single wells
  const w = byId.get(el.dataset.id), cur = state.groups[w.id] || null;
  if (cur === paintTo) return;
  assign(w, paintTo); scheduleSide();
}
function boxApply(tw) {                     // box: everything between the start well and `tw`, inside the start panel
  const p = panels[drag.panel], ws = panelWells[drag.panel];
  let inside;
  if (p.kind === 'plate') {
    const r0 = Math.min(drag.start.row, tw.row), r1 = Math.max(drag.start.row, tw.row), c0 = Math.min(drag.start.col, tw.col), c1 = Math.max(drag.start.col, tw.col);
    inside = w => w.row >= r0 && w.row <= r1 && w.col >= c0 && w.col <= c1;
  } else {                                  // unplaced list: range in list order
    const i0 = ws.indexOf(drag.start), i1 = ws.indexOf(tw), lo = Math.min(i0, i1), hi = Math.max(i0, i1);
    const pos = new Map(ws.map((w, i) => [w.id, i])); inside = w => pos.get(w.id) >= lo && pos.get(w.id) <= hi;
  }
  for (const w of ws) {
    const want = inside(w) ? paintTo : drag.before.get(w.id), cur = state.groups[w.id] || null;
    if (cur !== want) assign(w, want);
  }
  drag.last = tw.id; scheduleSide();
}
function scheduleSide() { if (!rafSide) rafSide = requestAnimationFrame(() => { rafSide = 0; renderSide(); }); }
const plates = $('#plates');
plates.addEventListener('pointerdown', e => {
  if (e.button !== 0) return;
  const el = e.target.closest('.well[data-id]');
  if (el) {
    e.preventDefault();
    const w = byId.get(el.dataset.id), cur = state.groups[w.id] || null, m = modeGroup();
    paintTo = cur === m ? null : m; painting = true;
    drag = {panel: w.panel, start: w, last: w.id, before: new Map(panelWells[w.panel].map(x => [x.id, state.groups[x.id] || null]))};
    if (state.drag === 'box') boxApply(w); else applyEl(el);
    return;
  }
  const hdr = e.target.closest('.hc,.hr');
  if (hdr) {
    const p = +hdr.closest('.panel').dataset.p, line = hdr.classList.contains('hc') ? ['col', +hdr.dataset.col] : ['row', +hdr.dataset.row];
    lineSelect(wells.filter(w => w.panel === p && w[line[0]] === line[1]));
  }
});
window.addEventListener('pointermove', e => {
  if (!painting) return;
  const el = document.elementFromPoint(e.clientX, e.clientY);
  const t = el && el.closest && el.closest('.well[data-id]'); if (!t) return;
  const w = byId.get(t.dataset.id);
  if (state.drag === 'box') { if (w.panel === drag.panel && w.id !== drag.last) boxApply(w); }
  else applyEl(t);
});
window.addEventListener('pointerup', () => { if (painting) { painting = false; drag = null; save(); renderSide(); } });
function lineSelect(ws) {
  const m = modeGroup(), all = m && ws.every(w => state.groups[w.id] === m);
  ws.forEach(w => assign(w, all ? null : m)); save(); renderSide();
}
plates.addEventListener('click', e => {
  const b = e.target.closest('button[data-act]'); if (!b) return;
  const p = +b.closest('.panel').dataset.p, ws = wells.filter(w => w.panel === p);
  ws.forEach(w => assign(w, b.dataset.act === 'selall' ? modeGroup() : null)); save(); renderSide();
});

// tooltip
const tip = $('#tip');
function tipHTML(w) {
  const p = panels[w.panel], pt = sum(wells.filter(x => x.panel === w.panel).map(x => x.reads)), g = state.groups[w.id];
  const rows = [['reads', `<b>${fmt(w.reads)}</b>` + (pt ? ` <span class="k">(${(100 * w.reads / pt).toFixed(2)}% of panel)</span>` : '')]];
  if (w.fwd || w.rev) rows.push(['utags', `${esc(w.fwd)} / ${esc(w.rev)}`]);
  if (w.ont) rows.push(['ONT bin', esc(w.ont)]);
  rows.push(['marked', g ? `<span style="color:${GC[g]}">${gname(g)}</span>` : 'positive']);
  return `<b>${esc(w.sample)}</b>` + (w.well ? ` <span class="k">${esc(p.title)} ${w.well}</span>` : '') +
    '<div>' + rows.map(r => `<div><span class="k">${r[0]}</span> ${r[1]}</div>`).join('') + '</div>';
}
plates.addEventListener('pointerover', e => {
  const el = e.target.closest('.well[data-id]'); if (!el) { tip.style.display = 'none'; return; }
  tip.innerHTML = tipHTML(byId.get(el.dataset.id)); tip.style.display = 'block';
});
plates.addEventListener('pointermove', e => {
  if (tip.style.display === 'none') return;
  const x = Math.min(e.clientX + 14, innerWidth - tip.offsetWidth - 8), y = Math.min(e.clientY + 14, innerHeight - tip.offsetHeight - 8);
  tip.style.left = x + 'px'; tip.style.top = y + 'px';
});
plates.addEventListener('pointerleave', () => { tip.style.display = 'none'; });

// ---------- sidebar: negative controls
function auc(pos, neg) {                    // Mann-Whitney: P(positive well > negative well), ties count half
  const all = pos.map(v => [v, 1]).concat(neg.map(v => [v, 0])).sort((x, y) => x[0] - y[0]);
  let rs = 0;
  for (let i = 0; i < all.length;) { let j = i; while (j < all.length && all[j][0] === all[i][0]) j++;
    const rank = (i + 1 + j) / 2; for (let k = i; k < j; k++) if (all[k][1]) rs += rank; i = j; }
  return (rs - pos.length * (pos.length + 1) / 2) / (pos.length * neg.length);
}
const callAt = (pos, neg, t) => { const tp = pos.filter(v => v >= t).length, fp = neg.filter(v => v >= t).length;
  return {t, tp, fn: pos.length - tp, fp, tn: neg.length - fp, sens: tp / pos.length, spec: 1 - fp / neg.length}; };
const isPos = w => !state.groups[w.id];
function analysis() {
  const P = wells.filter(isPos), N = wells.filter(w => state.groups[w.id] === 'N'), X = wells.filter(w => state.groups[w.id] === 'X');
  const r = {P, N, X};
  if (!P.length || !N.length) return r;
  const pos = P.map(w => w.reads), neg = N.map(w => w.reads), maxNeg = Math.max(...neg);
  const cands = [...new Set([1, maxNeg + 1, ...pos, ...neg].filter(v => v >= 1))].sort((x, y) => x - y);
  let best = null; for (const t of cands) { const c = callAt(pos, neg, t); if (!best || c.sens + c.spec > best.sens + best.spec) best = c; }
  Object.assign(r, {pos, neg, maxNeg, negReads: sum(neg), auc: auc(pos, neg), clean: callAt(pos, neg, maxNeg + 1), best, at: callAt(pos, neg, state.cutoff)});
  // what the negatives can detect: utags no positive uses catch misreads of that tag only (specificity = lower bound)
  const uf = new Set(P.map(w => w.fwd).filter(Boolean)), ur = new Set(P.map(w => w.rev).filter(Boolean));
  if (uf.size) {
    const fF = N.filter(w => !uf.has(w.fwd)), fR = N.filter(w => !ur.has(w.rev));
    const kind = fF.length === N.length && !fR.length ? 'fwd' : fR.length === N.length && !fF.length ? 'rev' : !fF.length && !fR.length ? 'blank' : 'mixed';
    const k = new Set((kind === 'fwd' ? fF.map(w => w.fwd) : fR.map(w => w.rev))).size, m = kind === 'fwd' ? uf.size : ur.size;
    r.kind = {kind, k, m, factor: kind === 'fwd' || kind === 'rev' ? (m - 1 + k) / k : null};
  }
  return r;
}
function strip(r) {                         // reads per well, positives vs negatives, log scale, cutoff line
  const W = Math.max(318, Math.round($('#side').clientWidth - (state.wide ? 48 : 32))), lw = 58, pw = W - lw - 8, H = state.wide ? 96 : 78, mx = Math.max(10, ...wells.map(w => w.reads));
  const x = v => lw + pw * Math.log10(v + 1) / Math.log10(mx + 1);
  let h = `<svg class="strip" viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="reads per well, positives and negatives">`;
  for (const t of [0, 1, 10, 100, 1e3, 1e4, 1e5, 1e6].filter(t => t <= mx)) h += `<line x1="${x(t)}" x2="${x(t)}" y1="0" y2="${H - 14}" stroke="var(--line)"/><text x="${x(t)}" y="${H - 3}" text-anchor="middle">${compact(t)}</text>`;
  h += `<text x="${lw - 6}" y="20" text-anchor="end">positive</text><text x="${lw - 6}" y="48" text-anchor="end">negative</text>`;
  const dot = (ws, y0, c, rr, op) => ws.map((w, j) => `<circle cx="${x(w.reads).toFixed(1)}" cy="${(y0 + ((j * 2654435761) % 1000) / 1000 * 16).toFixed(1)}" r="${rr}" fill="${c}" fill-opacity="${op}"><title>${esc(w.sample)}: ${fmt(w.reads)}</title></circle>`).join('');
  h += dot(r.P, 8, 'var(--P)', 2.2, .45) + dot(r.N, 36, 'var(--N)', 3, .85);
  h += `<line x1="${x(state.cutoff - .5)}" x2="${x(state.cutoff - .5)}" y1="0" y2="${H - 14}" stroke="var(--ink)" stroke-width="1.5" stroke-dasharray="3 3"/>`;
  return h + '</svg>';
}
function renderSide() {
  const r = analysis(), pp = v => (100 * v).toFixed(1) + '%', cut = state.cutoff;
  const btn = (m, c, lab, k) => `<button data-mode="${m}" class="${state.mode === m ? 'on' : ''}" style="--c:${c}"><span class="chip"></span>${lab} <kbd>${k}</kbd></button>`;
  let h = `<h2 class="sh">Mark wells as<button id="wide" title="${state.wide ? 'back to the narrow sidebar (Esc)' : 'expand the sidebar over half the page'}">${state.wide ? '\u2192 collapse' : '\u2190 expand'}</button></h2><div class="modes">${btn('N', 'var(--N)', 'negative', 'n')}${btn('X', 'var(--X)', 'exclude', 'x')}${btn('E', 'var(--P)', 'positive', 'p')}</div>` +
    `<div class="seg">drag <button data-drag="box" class="${state.drag === 'box' ? 'on' : ''}" title="drag a rectangle">box</button><button data-drag="free" class="${state.drag === 'free' ? 'on' : ''}" title="paint any wells you pass over">freehand</button></div>` +
    `<div class="row"><input id="selQ" placeholder="Plate4, *-H12, NB02 · Plate4 …" style="flex:1;min-width:0"><button id="selGo">mark</button></div>` +
    `<div class="row"><select id="thOp"><option value="le">≤</option><option value="ge">≥</option></select><input id="thN" type="number" min="0" value="0"><button id="thGo">mark by reads</button></div>` +
    `<div class="note">${r.N.length} negative · ${r.X.length} excluded · ${r.P.length} positive wells. Negative = should have no reads (unused plate, NTC, empty well); excluded = not loaded or known dropout.</div>` +
    `<h2>Negative controls</h2><div class="row"><label style="color:var(--dim)">call a well with ≥ <input id="cut" type="number" min="1" value="${cut}" style="width:70px"> reads</label></div>`;
  if (!r.auc && r.auc !== 0) h += `<div class="note">Mark at least one negative well${r.P.length ? '' : ' and keep some positives'} to see sensitivity and specificity.</div>`;
  else {
    const kt = !r.kind ? '' : r.kind.kind === 'fwd' || r.kind.kind === 'rev'
      ? `Negatives use ${r.kind.k} ${r.kind.kind} utag(s) no positive uses: they catch ${r.kind.kind}-tag misreads only, so specificity is a lower bound. Estimated all ${r.kind.kind}-tag misassignment (× ${r.kind.factor.toFixed(2)}): <b>${(100 * r.negReads * r.kind.factor / allReads).toFixed(3)}%</b>.`
      : r.kind.kind === 'blank' ? 'Negatives use utags that positives use too (in-design blanks): they catch errors on either tag and contamination.' : 'Negatives mix unused-tag wells and in-design blanks.';
    h += `<div class="res"><dl><dt>sensitivity at ≥ ${fmt(cut)}</dt><dd><b>${pp(r.at.sens)}</b></dd><dt>specificity at ≥ ${fmt(cut)}</dt><dd><b>${pp(r.at.spec)}</b></dd>` +
      `<dt>TP / FN / TN / FP</dt><dd>${r.at.tp} / ${r.at.fn} / ${r.at.tn} / ${r.at.fp}</dd></dl></div>` +
      `<div class="res"><dl><dt>cutoff for 100% specificity</dt><dd><b>≥ ${fmt(r.clean.t)}</b></dd><dt>… sensitivity there</dt><dd>${pp(r.clean.sens)}</dd>` +
      `<dt>best cutoff (Youden)</dt><dd>≥ ${fmt(r.best.t)} · ${pp(r.best.sens)} / ${pp(r.best.spec)}</dd><dt>AUC</dt><dd>${r.auc.toFixed(4)}</dd>` +
      `<dt>reads in negative wells</dt><dd>${fmt(r.negReads)} (${fmt(1e6 * r.negReads / allReads)} ppm)</dd><dt>max in one negative well</dt><dd>${fmt(r.maxNeg)}</dd></dl></div>` +
      strip(r) + (kt ? `<div class="note">${kt}</div>` : '');
    const acc = new Map();                  // reads in negative wells per utag
    for (const w of r.N) for (const [side, t] of [['fwd', w.fwd], ['rev', w.rev]]) if (t) acc.set(side + ' ' + t, (acc.get(side + ' ' + t) || 0) + w.reads);
    const top = [...acc.entries()].filter(e => e[1] > 0).sort((x, y) => y[1] - x[1]).slice(0, state.wide ? 24 : 8);
    if (top.length) h += `<h2>Reads in negatives per utag</h2><table class="mini"><tr><th>utag</th><th>reads</th><th>share</th></tr>` +
      top.map(([t, n]) => `<tr><td>${esc(t)}</td><td>${fmt(n)}</td><td>${(100 * n / r.negReads).toFixed(0)}%</td></tr>`).join('') + '</table>';
    const fn = r.P.filter(w => w.reads < cut).sort((x, y) => x.reads - y.reads);
    if (fn.length) h += `<h2>Positives below the cutoff (${fn.length})</h2><div class="note">failed samples or demux misses; exclude the ones you know were not loaded or failed. <button id="fnAll">exclude all</button></div>` +
      `<div class="list"><table class="mini">` + fn.slice(0, 200).map(w => `<tr><td>${esc(w.sample)}${w.ont ? ' <small style="color:var(--dim)">' + esc(w.ont) + '</small>' : ''}</td><td>${fmt(w.reads)}</td><td><button data-x="${esc(w.id)}">exclude</button></td></tr>`).join('') + '</table></div>';
  }
  h += `<div class="btns" style="margin-top:12px"><button id="clrAll">reset all to positive</button><button id="exp">download results CSV</button></div>` +
    `<div class="hint">click a well to mark or unmark it · drag a box (or freehand) · click a row or column header for the whole line, “mark all” for a plate · starting on a well already in the active mode unmarks instead · keys: <b>n</b> negative, <b>x</b> exclude, <b>p</b> positive</div>`;
  const act = document.activeElement && document.activeElement.id, val = act === 'selQ' ? document.activeElement.value : null;
  $('#side').innerHTML = h;
  if (act === 'cut') { const i = $('#cut'); i.focus(); i.setSelectionRange(i.value.length, i.value.length); }
  if (act === 'selQ') { const i = $('#selQ'); i.value = val; i.focus(); }
  window.__ph.analysis = r;
}
const side = $('#side');
const glob = g => new RegExp('^' + g.replace(/[.+^${}()|[\]\\]/g, '\\$&').replace(/\*/g, '.*').replace(/\?/g, '.') + '$', 'i');
function markSelector() {                   // sample names, ONT/sample or panel titles; * and ? wildcards
  const res = ($('#selQ').value || '').split(/[\s,;]+/).filter(Boolean).map(glob); if (!res.length) return;
  wells.filter(w => res.some(re => re.test(w.sample) || re.test(panels[w.panel].title) || (w.ont && re.test(w.ont + '/' + w.sample))))
    .forEach(w => assign(w, modeGroup())); save(); renderSide();
}
side.addEventListener('click', e => {
  const t = e.target;
  if (t.closest('[data-drag]')) { state.drag = t.closest('[data-drag]').dataset.drag; renderSide(); return; }
  if (t.closest('[data-mode]')) { state.mode = t.closest('[data-mode]').dataset.mode; renderSide(); return; }
  if (t.closest('[data-x]')) { assign(byId.get(t.closest('[data-x]').dataset.x), 'X'); save(); renderSide(); return; }
  if (t.id === 'fnAll') { wells.filter(w => isPos(w) && w.reads < state.cutoff).forEach(w => assign(w, 'X')); save(); renderSide(); return; }
  if (t.id === 'selGo') { markSelector(); return; }
  if (t.id === 'wide') { state.wide = !state.wide; $('#app').classList.toggle('wide', state.wide); renderSide(); return; }
  if (t.id === 'clrAll') { wells.forEach(w => assign(w, null)); save(); renderSide(); return; }
  if (t.id === 'thGo') {
    const n = +$('#thN').value, ge = $('#thOp').value === 'ge', m = modeGroup();
    wells.filter(w => ge ? w.reads >= n : w.reads <= n).forEach(w => assign(w, m)); save(); renderSide(); return;
  }
  if (t.id === 'exp') {
    const q = v => /[",\n]/.test(String(v)) ? '"' + String(v).replace(/"/g, '""') + '"' : v, r = analysis(), lines = [];
    if (r.at) {
      const pp = v => (100 * v).toFixed(2);
      lines.push('metric,value', `positive_wells,${r.P.length}`, `negative_wells,${r.N.length}`, `excluded_wells,${r.X.length}`,
        `cutoff,${state.cutoff}`, `sensitivity_pct,${pp(r.at.sens)}`, `specificity_pct,${pp(r.at.spec)}`, `tp,${r.at.tp}`, `fn,${r.at.fn}`, `tn,${r.at.tn}`, `fp,${r.at.fp}`,
        `cutoff_100pct_specificity,${r.clean.t}`, `sensitivity_at_that_cutoff_pct,${pp(r.clean.sens)}`, `youden_cutoff,${r.best.t}`,
        `youden_sensitivity_pct,${pp(r.best.sens)}`, `youden_specificity_pct,${pp(r.best.spec)}`, `auc,${r.auc.toFixed(4)}`,
        `reads_in_negative_wells,${r.negReads}`, `negative_reads_ppm_of_all,${Math.round(1e6 * r.negReads / allReads)}`, `max_reads_in_negative_well,${r.maxNeg}`, '');
    }
    lines.push('state,sample,ont_barcode,panel,well,reads,utag_fwd,utag_rev');
    for (const w of wells) lines.push([gname(state.groups[w.id]), w.sample, w.ont, panels[w.panel].title, w.well || '', w.reads, w.fwd, w.rev].map(q).join(','));
    const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], {type: 'text/csv'})); a.download = 'negative_controls.csv'; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }
});
side.addEventListener('input', e => { if (e.target.id === 'cut') { state.cutoff = Math.max(1, +e.target.value || 1); save(); renderSide(); } });
side.addEventListener('keydown', e => { if (e.target.id === 'selQ' && e.key === 'Enter') markSelector(); });
window.addEventListener('keydown', e => {
  if (e.key === 'Escape' && state.wide) { state.wide = false; $('#app').classList.remove('wide'); renderSide(); return; }
  if (/^(INPUT|SELECT|TEXTAREA)$/.test(document.activeElement.tagName) || e.ctrlKey || e.metaKey || e.altKey) return;
  const m = {'n': 'N', '1': 'N', 'x': 'X', '2': 'X', 'p': 'E', '0': 'E', 'e': 'E'}[e.key.toLowerCase()]; if (m) { state.mode = m; renderSide(); }
});

window.__ph = {state, wells, panels};   // for debugging in the console
legend(); decorateAll(); renderSide();
})();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()

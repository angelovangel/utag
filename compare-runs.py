#!/usr/bin/env python3
"""
compare-runs.py - compare utag runs side by side, from the utag.log files written by utag.py.

    compare-runs.py demux_out/utag.log run2/utag.log
    compare-runs.py demux_out run2 twopass_out -o nbd_vs_lsk           # folders: <folder>/utag.log
    compare-runs.py runs/ --labels default fme2 lsk --baseline fme2    # a folder holding several runs

Every argument is a utag.log or a folder. A folder with a utag.log is one run, otherwise every utag.log below it
is a run (the per-bin logs of twopass, barbell/pass2_utag/<bin>/utag.log, are skipped).

Read from each log: start time, command, context/kit, barbell step times, warnings, the run summary (demux) or
two-pass summary, the read-pattern breakdown ("found"), the top patterns and, for twopass, the per-ONT-bin table.
twopass logs have no pattern tables: when the per-bin logs are still on disk, their patterns are summed.
When a summary.tsv sits next to a log, the reads per sample are compared too (twopass samples as <bin>/<sample>,
or by sample name with --ignore-bins).

demux and twopass name two metrics differently; both are reported under the twopass names:
    reads_assigned_to_wells -> reads_assigned_to_samples, wells_with_reads -> samples_with_reads,
    wells_expected -> samples_expected

Writes, into the output folder (-o, default compare_out/)
    run_comparison.csv    one row per metric, one column per run, then <run>_vs_<baseline> difference columns
    samples.csv           reads per sample and run (only when a summary.tsv was found)
    run_comparison.html   ONE self-contained page (no internet): read fate and pattern charts, metrics table with
                          a selectable baseline, top patterns, ONT bins, per-sample scatter of any two runs, and
                          negative controls: mark empty wells (e.g. an unused plate) on the plate grids and read
                          sensitivity / specificity per run (the selection stays in the browser; results download as CSV)

Python >= 3.8, standard library only.
"""
import argparse
import csv
import glob
import html
import json
import os
import re
import sys
import time

ALIAS = {"reads_assigned_to_wells": "reads_assigned_to_samples", "wells_with_reads": "samples_with_reads",
         "wells_expected": "samples_expected"}
FOUND = {"utag at both ends": "tag_both_ends", "utag at one end, primer-flank": "tag_plus_flank_other_end",
         "utag at ONE end only": "tag_one_end_only", "primer-flank only": "flank_only"}
FOUND_LABEL = {"tag_both_ends": "utag at both ends", "tag_plus_flank_other_end": "utag + primer-flank",
               "tag_one_end_only": "utag at one end only", "flank_only": "primer-flank only"}
# read counts that are also shown as % of total reads (same rule as utag.py's run_summary.csv)
PCT = lambda k: k.startswith(("reads", "pattern_reads", "unexpected_reads"))

STEP_RE = re.compile(r"^── (.+?) ─+\s*$")
NUM = r"(\d[\d,]*)"
SUMMARY_RE = re.compile(rf"^ {{4}}([a-z][a-z0-9_]+)\s+{NUM}(?:\s+[\d.]+%)?\s*$")
PATTERN_RE = re.compile(rf"^ {{4}}(\d+)\s+{NUM}\s+([\d.]+)%\s+(\S.*?)\s*$")
FOUND_RE = re.compile(rf"^ {{4}}(.+?)\s+{NUM}\s+([\d.]+)%")
BIN_RE = re.compile(rf"^ {{4}}(\S+)\s+{NUM}\s+{NUM}\s+{NUM}\s+([\d.]+)%\s+\S*\s+(\d+)/(\d+)\s*$")
BIN_INFO_RE = re.compile(rf"{NUM} reads in bin, {NUM} assigned")
TIME_RE = re.compile(r"^\s+✓ ([\d.]+)s")
UNIQ_RE = re.compile(r"pattern \(top \d+ of (\d+)\)")


def num(s):
    return int(s.replace(",", ""))


def parse_log(path):
    r = {"log": path, "started": "", "cwd": "", "command": "", "subcommand": "", "context": "", "kit": "",
         "status": "incomplete", "barbell_seconds": 0.0, "steps": {}, "warnings": [], "errors": [],
         "summary": {}, "found": {}, "patterns": {}, "unique_patterns": None, "bins": {}, "info": []}
    step, pass2_bin = "", None
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    for ln in lines:
        if ln.startswith("# ") and not r["started"]:
            m = re.match(r"# (\S+ \S+)\s+cwd: (.*)", ln)
            if m:
                r["started"], r["cwd"] = m.groups()
                continue
        if ln.startswith("# ") and not r["command"]:
            r["command"] = ln[2:].strip()
            toks = r["command"].split()
            r["subcommand"] = next((t for t in toks if t in ("demux", "twopass", "check", "inspect")), "")
            for opt, key in (("-c", "context"), ("--context", "context"), ("-k", "kit"), ("--kit", "kit")):
                if opt in toks and toks.index(opt) + 1 < len(toks):
                    r[key] = toks[toks.index(opt) + 1]
            continue
        m = STEP_RE.match(ln)
        if m:
            step = m.group(1)
            pm = re.match(r"\[pass 2\] (\S+):", step)
            pass2_bin = pm.group(1) if pm else None
            qm = re.match(r"\[1/6\] query FASTAs \((\S+)", step)
            if qm and not r["context"]:
                r["context"] = qm.group(1)
            continue
        m = TIME_RE.match(ln)
        if m:
            sec = float(m.group(1))
            r["barbell_seconds"] += sec
            r["steps"][step] = r["steps"].get(step, 0) + sec
            continue
        s = ln.strip()
        if s.startswith("! "):
            r["warnings"].append(s[2:])
            continue
        if s.startswith(("ERROR", "Traceback")) or s.startswith("error:"):
            r["errors"].append(s)
            continue
        if s.startswith("• "):
            r["info"].append(s[2:])
            m = BIN_INFO_RE.search(s)
            if m and pass2_bin:
                r["bins"].setdefault(pass2_bin, {}).update(reads=num(m.group(1)), assigned=num(m.group(2)))
            continue
        if step.startswith("[3/6] inspect"):
            m = UNIQ_RE.search(ln)
            if m:
                r["unique_patterns"] = int(m.group(1))
                continue
            m = PATTERN_RE.match(ln)
            if m:
                r["patterns"][m.group(4)] = r["patterns"].get(m.group(4), 0) + num(m.group(2))
                continue
            m = FOUND_RE.match(ln)
            if m:
                k = next((v for p, v in FOUND.items() if m.group(1).startswith(p)), None)
                if k:
                    r["found"][k] = num(m.group(2))
                continue
        if step in ("run summary", "two-pass summary"):
            r["status"] = "ok"
            m = BIN_RE.match(ln)
            if m:
                b = m.group(1)
                r["bins"][b] = {"reads": num(m.group(2)), "with_hit": num(m.group(3)), "assigned": num(m.group(4)),
                                "samples_with_reads": int(m.group(6)), "samples_expected": int(m.group(7))}
                continue
            m = SUMMARY_RE.match(ln)
            if m:
                r["summary"][ALIAS.get(m.group(1), m.group(1))] = num(m.group(2))
    if r["errors"] and r["status"] != "ok":
        r["status"] = "failed"
    if r["subcommand"] == "twopass" and not r["patterns"]:     # sum the per-bin demux logs, if still there
        sub = sorted(glob.glob(os.path.join(os.path.dirname(path), "barbell", "pass2_utag", "*", "utag.log")))
        for p in sub:
            b = parse_log(p)
            for k, v in b["found"].items():
                r["found"][k] = r["found"].get(k, 0) + v
            for k, v in b["patterns"].items():
                r["patterns"][k] = r["patterns"].get(k, 0) + v
            r["barbell_seconds"] += b["barbell_seconds"]
        if sub:
            r["info"].append(f"patterns summed over {len(sub)} per-bin logs (top patterns of each bin only)")
    return r


def read_samples(log, subcommand, ignore_bins=False):
    """({sample: reads}, {sample: [plate, well, utag_fwd, utag_rev]}) from the summary.tsv next to the log"""
    p = os.path.join(os.path.dirname(log), "summary.tsv")
    if not os.path.exists(p):
        return None, {}
    out, meta = {}, {}
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                n = int(row.get("reads") or 0)
            except ValueError:
                continue
            key, plate = row.get("sample", ""), row.get("plate", "")
            if subcommand == "twopass" and row.get("ont_barcode") and not ignore_bins:
                key, plate = f"{row['ont_barcode']}/{key}", f"{row['ont_barcode']}/{plate}"
            out[key] = out.get(key, 0) + n
            meta.setdefault(key, [plate, row.get("well", ""), row.get("utag_fwd", ""), row.get("utag_rev", "")])
    return out, meta


def find_logs(args):
    logs = []
    for a in args:
        if os.path.isfile(a):
            logs.append(a)
        elif os.path.isdir(a):
            if os.path.isfile(os.path.join(a, "utag.log")):
                logs.append(os.path.join(a, "utag.log"))
            else:
                found = sorted(p for p in glob.glob(os.path.join(a, "**", "utag.log"), recursive=True)
                               if f"{os.sep}pass2_utag{os.sep}" not in p)
                if not found:
                    sys.exit(f"no utag.log in {a}")
                logs += found
        else:
            sys.exit(f"not found: {a}")
    seen, uniq = set(), []
    for p in logs:
        if os.path.realpath(p) not in seen:
            seen.add(os.path.realpath(p))
            uniq.append(p)
    return uniq


def make_labels(logs):
    """the run folder name; parent folders are added until the labels are unique"""
    parts = [os.path.normpath(os.path.dirname(os.path.abspath(p))).split(os.sep) for p in logs]
    for depth in range(1, max(len(x) for x in parts) + 1):
        labels = ["/".join(x[-depth:]) or "run" for x in parts]
        if len(set(labels)) == len(labels):
            return labels
    return [f"run{i + 1}" for i in range(len(logs))]


def rows_for_csv(runs, top):
    """[(section, field, {label: value})], same order as shown in the HTML"""
    labels = [r["label"] for r in runs]
    rows = []
    for k in ("log", "started", "command", "subcommand", "context", "kit", "status"):
        rows.append(("run", k, {r["label"]: r[k] for r in runs}))
    rows.append(("run", "barbell_seconds", {r["label"]: round(r["barbell_seconds"], 1) for r in runs}))
    rows.append(("run", "warnings", {r["label"]: " | ".join(r["warnings"]) for r in runs}))
    keys = []
    for r in runs:
        keys += [k for k in r["summary"] if k not in keys]
    for k in keys:
        rows.append(("summary", k, {r["label"]: r["summary"].get(k, "") for r in runs}))
        if PCT(k):
            rows.append(("summary", f"{k}_pct", {r["label"]: pct(r["summary"].get(k), r["summary"].get("total_reads"))
                                                 for r in runs}))
    for k in FOUND_LABEL:
        if any(k in r["found"] for r in runs):
            rows.append(("found", k, {r["label"]: r["found"].get(k, "") for r in runs}))
            rows.append(("found", f"{k}_pct", {r["label"]: pct(r["found"].get(k), sum(r["found"].values()))
                                               for r in runs}))
    for p in top_patterns(runs, top):
        rows.append(("pattern", p, {r["label"]: r["patterns"].get(p, "") for r in runs}))
    bins = sorted({b for r in runs for b in r["bins"]})
    for b in bins:
        for f in ("reads", "with_hit", "assigned", "samples_with_reads", "samples_expected"):
            vals = {r["label"]: r["bins"].get(b, {}).get(f, "") for r in runs}
            if any(v != "" for v in vals.values()):
                rows.append(("ont_bin", f"{b}:{f}", vals))
        rows.append(("ont_bin", f"{b}:assigned_pct", {r["label"]: pct(r["bins"].get(b, {}).get("assigned"),
                                                                       r["bins"].get(b, {}).get("reads")) for r in runs}))
    return labels, rows


def pct(a, b):
    return round(100 * a / b, 2) if isinstance(a, (int, float)) and b else ""


def top_patterns(runs, top):
    """patterns in the top `top` of any run, ordered by their reads over all runs"""
    keep = {p for r in runs for p, _ in sorted(r["patterns"].items(), key=lambda x: -x[1])[:top]}
    return sorted(keep, key=lambda p: -sum(r["patterns"].get(p, 0) for r in runs))


def write_csv(path, labels, rows, base):
    others = [l for l in labels if l != base]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["section", "field"] + labels + [f"{l}_vs_{base}" for l in others])
        for sec, field, vals in rows:
            diffs = []
            for l in others:
                a, b = vals.get(l), vals.get(base)
                ok = all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in (a, b))
                diffs.append(round(a - b, 2) if ok and (sec != "run" or field == "barbell_seconds") else "")
            w.writerow([sec, field] + [vals.get(l, "") for l in labels] + diffs)


def write_samples_csv(path, runs):
    have = [r for r in runs if r["samples"] is not None]
    keys = []
    seen = set()
    for r in have:
        for k in r["samples"]:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample"] + [r["label"] for r in have])
        for k in keys:
            w.writerow([k] + [r["samples"].get(k, "") for r in have])
    return len(keys)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logs", nargs="+", help="utag.log files or folders holding them")
    ap.add_argument("-o", "--outdir", default="compare_out",
                    help="output folder for run_comparison.csv, samples.csv and run_comparison.html (default compare_out)")
    ap.add_argument("--labels", nargs="+", help="one name per run (default: the run folder names)")
    ap.add_argument("--baseline", help="run label the differences are taken against (default: the first run)")
    ap.add_argument("--title", default="utag run comparison", help="page title")
    ap.add_argument("--theme", choices=["light", "dark", "auto"], default="light",
                    help="page colours; auto follows the viewer's system setting (default light)")
    ap.add_argument("--top", type=int, default=8, help="patterns from the top N of each run (default 8)")
    ap.add_argument("--ignore-bins", action="store_true",
                    help="twopass: match samples by name only, reads summed over ONT bins (to compare with demux runs)")
    ap.add_argument("--no-samples", action="store_true", help="do not read the summary.tsv next to each log")
    a = ap.parse_args()

    logs = find_logs(a.logs)
    if a.labels and len(a.labels) != len(logs):
        sys.exit(f"--labels: {len(a.labels)} names for {len(logs)} runs ({', '.join(logs)})")
    labels = a.labels or make_labels(logs)
    if len(set(labels)) != len(labels):
        sys.exit("--labels must be unique")
    runs = []
    for p, l in zip(logs, labels):
        r = parse_log(p)
        r["label"] = l
        r["samples"], r["meta"] = (None, {}) if a.no_samples else read_samples(p, r["subcommand"], a.ignore_bins)
        runs.append(r)
    base = a.baseline or labels[0]
    if base not in labels:
        sys.exit(f"--baseline {base}: not one of {', '.join(labels)}")

    os.makedirs(a.outdir, exist_ok=True)
    metrics_csv, samples_csv, page = (os.path.join(a.outdir, f) for f in
                                      ("run_comparison.csv", "samples.csv", "run_comparison.html"))
    _, rows = rows_for_csv(runs, a.top)
    write_csv(metrics_csv, labels, rows, base)
    written = [metrics_csv]
    n_samples = 0
    if any(r["samples"] is not None for r in runs):
        n_samples = write_samples_csv(samples_csv, runs)
        written.append(samples_csv)
    elif os.path.exists(samples_csv):
        os.remove(samples_csv)        # from an earlier comparison in this folder: it would not match these runs
    with open(page, "w", encoding="utf-8") as fh:
        fh.write(render_html(runs, base, a.title, a.top, a.theme))
    written.append(page)

    w = max(len(l) for l in labels)
    for r in runs:
        s = r["summary"]
        tot, asg = s.get("total_reads"), s.get("reads_assigned_to_samples")
        frac = f"{100 * asg / tot:5.1f}% assigned" if tot and asg is not None else "no summary"
        print(f"  {r['label']:<{w}}  {r['status']:<10}  {r['subcommand'] or '?':<8} "
              f"{(f'{tot:,}' if tot else '-'):>10} reads  {frac}"
              + (f"  ({r['errors'][-1][:60]})" if r["status"] == "failed" else ""))
    if n_samples:
        print(f"  {n_samples:,} samples compared")
    for p in written:
        print(f"  -> {p}")


def render_html(runs, base, title, top, theme="light"):
    meta = {}
    for r in runs:
        for k, v in r["meta"].items():
            meta.setdefault(k, v)
    data = {
        "base": base, "meta": meta, "generated": time.strftime("%Y-%m-%d %H:%M"),
        "foundLabels": FOUND_LABEL, "patterns": top_patterns(runs, top),
        "runs": [{**{k: r[k] for k in ("label", "log", "started", "cwd", "command", "subcommand", "context", "kit",
                                       "status", "steps", "warnings", "errors", "summary", "found", "patterns",
                                       "unique_patterns", "bins", "info", "samples")},
                  "barbell_seconds": round(r["barbell_seconds"], 1)} for r in runs],
    }
    js = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    attr = "" if theme == "auto" else f' data-theme="{theme}"'
    return HTML.replace("__THEME__", attr).replace("__TITLE__", html.escape(title)).replace("__DATA__", js)


HTML = r"""<!doctype html>
<html lang="en"__THEME__><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{color-scheme:light;--bg:#f9f9f7;--panel:#fcfcfb;--panel2:#f1f0ec;--ink:#0b0b0b;--dim:#52514e;--muted:#898781;
  --grid:#e1e0d9;--axis:#c3c2b7;--ring:rgba(11,11,11,.10);
  --s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--s4:#eda100;--good:#006300;--bad:#d03b3b;
  --q0:#f1f0ec;--q1:#cde2fb;--q2:#9ec5f4;--q3:#6da7ec;--q4:#3987e5;--q5:#256abf}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])){color-scheme:dark;--bg:#0d0d0d;--panel:#1a1a19;
  --panel2:#252523;--ink:#f0efec;--dim:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;--ring:rgba(255,255,255,.10);
  --s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--good:#0ca30c;--bad:#e66767;
  --q0:#252523;--q1:#0d366b;--q2:#104281;--q3:#184f95;--q4:#1c5cab;--q5:#256abf}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#0d0d0d;--panel:#1a1a19;--panel2:#252523;--ink:#f0efec;--dim:#c3c2b7;
  --muted:#898781;--grid:#2c2c2a;--axis:#383835;--ring:rgba(255,255,255,.10);
  --s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--good:#0ca30c;--bad:#e66767;
  --q0:#252523;--q1:#0d366b;--q2:#104281;--q3:#184f95;--q4:#1c5cab;--q5:#256abf}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{padding:18px 24px 12px;border-bottom:1px solid var(--grid)}
h1{margin:0;font-size:19px;font-weight:650}
.sub{color:var(--dim);font-size:12px;margin-top:2px}
main{max-width:1280px;margin:0 auto;padding:16px 24px 48px;display:grid;gap:16px}
.panel{background:var(--panel);border:1px solid var(--ring);border-radius:12px;padding:14px 16px 16px;min-width:0}
h2{margin:0 0 10px;font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--dim);font-weight:600;display:flex;gap:12px;align-items:center;flex-wrap:wrap}
h2 .ctl{margin-left:auto;text-transform:none;letter-spacing:0;font-weight:400;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
select,input,button{font:inherit;color:var(--ink);background:var(--panel2);border:1px solid var(--ring);border-radius:6px;padding:3px 8px}
label{color:var(--dim);display:inline-flex;gap:5px;align-items:center}
.badge{display:inline-block;font-size:11px;padding:0 6px;border-radius:9px;border:1px solid currentColor;margin-left:6px;vertical-align:1px}
.ok{color:var(--good)}.failed{color:var(--bad)}.incomplete{color:var(--muted)}.base{color:var(--s1)}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
th,td{padding:4px 8px;border-bottom:1px solid var(--grid);text-align:right;white-space:nowrap}
th{color:var(--dim);font-weight:600;position:sticky;top:0;background:var(--panel)}
th:first-child,td:first-child{text-align:left}
td.pat{font-family:ui-monospace,Menlo,monospace;font-size:11.5px}
tr.sec td{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.05em;padding-top:12px;border-bottom:none}
.d{font-size:11px;margin-left:6px}.up{color:var(--good)}.down{color:var(--bad)}.neu{color:var(--muted)}
th.sortable{cursor:pointer}th.sortable:hover{color:var(--ink)}
.legend{display:flex;gap:4px 12px;flex-wrap:wrap;color:var(--dim);font-size:11.5px;margin:0 0 8px}
.legend span{white-space:nowrap}
.two.charts{grid-template-columns:minmax(0,1fr) minmax(0,1fr)}
.legend i{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:4px;vertical-align:-1px}
svg text{fill:var(--dim);font-size:11px}
.tip{position:fixed;pointer-events:none;background:var(--panel);color:var(--ink);border:1px solid var(--ring);border-radius:8px;
  padding:6px 9px;font-size:12px;box-shadow:0 4px 16px rgba(0,0,0,.18);display:none;z-index:9;max-width:420px}
.tip .t{color:var(--dim)}
.two{display:grid;grid-template-columns:minmax(0,1.25fr) minmax(0,1fr);gap:16px}
@media (max-width:900px){.two,.two.charts{grid-template-columns:1fr}header,main{padding-left:16px;padding-right:16px}}
.mono{font-family:ui-monospace,Menlo,monospace;font-size:11.5px;white-space:pre-wrap;word-break:break-all;color:var(--dim)}
.warn{color:var(--dim)}.warn::before{content:"! ";color:var(--s4);font-weight:700}
.err::before{content:"✕ ";color:var(--bad);font-weight:700}
.stat{color:var(--dim);font-size:12px;margin-top:6px}
.cell{border-radius:4px;padding:1px 6px}
details summary{cursor:pointer;color:var(--dim)}
.hidden{display:none}
.negbar{display:flex;gap:8px 14px;flex-wrap:wrap;align-items:center;margin-bottom:10px}
.modes{display:inline-flex;gap:4px}.modes button.on{border-color:var(--ink);box-shadow:inset 0 0 0 1px var(--ink)}
.modes i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px;vertical-align:0}
.plates{display:flex;flex-wrap:wrap;gap:14px;margin:6px 0 4px}
.pl{min-width:0}.pl .pt{font-weight:600;font-size:12px;cursor:pointer;margin-bottom:3px;display:flex;gap:8px;align-items:baseline}
.pl .pt small{color:var(--muted);font-weight:400}.pl .pt:hover{color:var(--s1)}
.pg{display:grid;gap:2px;user-select:none;-webkit-user-select:none}
.pg .h{font-size:9px;color:var(--muted);text-align:center;cursor:pointer;line-height:var(--cs)}.pg .h:hover{color:var(--ink)}
.w{width:var(--cs);height:var(--cs);border-radius:3px;background:var(--q0);cursor:pointer;position:relative}
.w.none{background:transparent;border:1px dashed var(--grid);cursor:default}
.w.neg{box-shadow:0 0 0 2px var(--panel),0 0 0 3.5px var(--bad)}
.w.neg::after{content:"";position:absolute;inset:30%;border-radius:50%;background:var(--panel)}
.w.exc{opacity:.3;background-image:repeating-linear-gradient(45deg,transparent 0 2px,var(--panel) 2px 3.5px)}
.w.pre{outline:2px solid var(--ink);outline-offset:0}
.kpi td.b{font-weight:600}
.note{color:var(--muted);font-size:11.5px;margin-top:6px;line-height:1.5}
.lowlist button{padding:0 6px;font-size:11px}
</style></head><body>
<header><h1>__TITLE__</h1><div class="sub" id="sub"></div></header>
<main>
  <div class="two charts">
    <section class="panel"><h2>Read fate <span class="ctl" style="font-size:12px;color:var(--muted)">% of total reads</span></h2>
      <div class="legend" id="fateLegend"></div><div id="fate"></div></section>
    <section class="panel"><h2>Tags found per read <span class="ctl" style="font-size:12px;color:var(--muted)">% of reads with a hit</span></h2>
      <div class="legend" id="foundLegend"></div><div id="found"></div></section>
  </div>
  <section class="panel"><h2>Metrics<span class="ctl">
      <label>baseline <select id="baseSel"></select></label>
      <label><input type="checkbox" id="diffOnly"> differences only</label>
      <label><input type="checkbox" id="showPct" checked> show % rows</label></span></h2>
    <div class="scroll"><table id="metrics"></table></div></section>
  <section class="panel" id="patPanel"><h2>Top read patterns<span class="ctl">
      <label><select id="patMode"><option value="share">share of reads with a hit</option><option value="reads">reads</option></select></label></span></h2>
    <div class="scroll"><table id="pats"></table></div></section>
  <section class="panel hidden" id="binPanel"><h2>ONT bins (twopass)<span class="ctl">
      <label><select id="binMode"><option value="pct">assigned %</option><option value="assigned">assigned reads</option><option value="reads">reads in bin</option><option value="samples">samples with reads</option></select></label></span></h2>
    <div class="scroll"><table id="bins"></table></div></section>
  <section class="panel hidden" id="smpPanel"><h2>Reads per sample<span class="ctl">
      <label>x <select id="sx"></select></label><label>y <select id="sy"></select></label>
      <label><select id="scale"><option value="log">log scale</option><option value="linear">linear scale</option></select></label>
      <label><input type="checkbox" id="hideZero"> hide samples with 0 in both</label></span></h2>
    <div class="two"><div><div id="scatter"></div><div class="stat" id="sstat"></div></div>
      <div><div style="display:flex;gap:8px;margin-bottom:6px"><input id="q" placeholder="filter samples…" style="flex:1">
        <select id="sortBy"><option value="diff">sort: largest difference</option><option value="fc">sort: largest fold change</option><option value="x">sort: x reads</option><option value="y">sort: y reads</option><option value="name">sort: name</option></select></div>
        <div class="scroll" style="max-height:420px;overflow-y:auto"><table id="stab"></table></div></div></div></section>
  <section class="panel hidden" id="negPanel"><h2>Negative controls<span class="ctl">
      <label>colour wells by <select id="ncol"></select></label></span></h2>
    <div class="negbar">
      <span class="modes" id="nmodes">
        <button data-m="neg" class="on"><i style="background:var(--bad)"></i>negative</button>
        <button data-m="exc"><i style="background:var(--muted)"></i>exclude</button>
        <button data-m="pos"><i style="background:var(--s1)"></i>positive</button></span>
      <input id="nsel" placeholder="Plate4, *-H12, NB02/Plate4-* …" style="width:220px"><button id="nselGo">mark</button>
      <button id="nclear">reset all to positive</button>
      <label>call a well with ≥ <input id="ncut" type="number" min="1" value="10" style="width:64px"> reads</label>
      <button id="ndl" style="margin-left:auto">download results CSV</button>
    </div>
    <div class="note" style="margin-top:0">Click or drag over wells, or click a plate name, row or column, to mark them with the active mode (click again to undo).
      Negative = should have no reads (empty wells, an unused plate, NTCs). Excluded wells count nowhere (not loaded, known dropouts). All other wells are positives.</div>
    <div class="plates" id="nplates"></div>
    <div id="nres"></div>
  </section>
  <section class="panel"><h2>Commands, timings and warnings</h2><div id="cmds"></div></section>
</main>
<div class="tip" id="tip"></div>
<script>
const D = __DATA__;
const R = D.runs, $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const fmt = v => v === undefined || v === null || v === "" ? "–" : (typeof v === "number" ? (Number.isInteger(v) ? v.toLocaleString("en-US") : v.toLocaleString("en-US", {maximumFractionDigits: 1})) : v);
const pc = (a, b) => (b ? 100 * a / b : null);
const pf = v => v === null || v === undefined ? "–" : v.toFixed(1) + "%";
const tip = $("tip");
function showTip(e, h) { tip.innerHTML = h; tip.style.display = "block"; const w = tip.offsetWidth, hh = tip.offsetHeight;
  let x = e.clientX + 14, y = e.clientY + 14; if (x + w > innerWidth - 8) x = e.clientX - w - 14; if (y + hh > innerHeight - 8) y = e.clientY - hh - 14;
  tip.style.left = x + "px"; tip.style.top = y + "px"; }
const hideTip = () => tip.style.display = "none";
let base = D.base;
const assigned = r => r.summary.reads_assigned_to_samples;

$("sub").textContent = `${R.length} runs · generated ${D.generated} · baseline ${base}`;

// ---- horizontal 100% stacked bars, one row per run
function stacked(el, legendEl, segs, rowsOf, note) {
  const colors = ["var(--s1)", "var(--s2)", "var(--s3)", "var(--s4)"];
  legendEl.innerHTML = segs.map((s, i) => `<span><i style="background:${colors[i]}"></i>${esc(s)}</span>`).join("");
  const lw = Math.min(170, 12 + 7 * Math.max(...R.map(r => r.label.length))), W = 640, rh = 26, top = 4, H = top + R.length * rh + 22;
  const pw = W - lw - 22;
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img">`;
  for (const t of [0, 25, 50, 75, 100]) { const x = lw + pw * t / 100;
    svg += `<line x1="${x}" x2="${x}" y1="${top}" y2="${H - 18}" stroke="var(--grid)"/><text x="${x}" y="${H - 5}" text-anchor="middle">${t}%</text>`; }
  R.forEach((r, i) => { const y = top + i * rh + 5, vals = rowsOf(r), tot = vals.reduce((a, b) => a + b, 0);
    svg += `<text x="${lw - 8}" y="${y + 12}" text-anchor="end" style="fill:var(--ink)">${esc(r.label.length > 24 ? "…" + r.label.slice(-23) : r.label)}</text>`;
    if (!tot) { svg += `<text x="${lw + 4}" y="${y + 12}">${note}</text>`; return; }
    let x = lw;
    vals.forEach((v, j) => { const w = pw * v / tot; if (w <= 0) return;
      const h = `<b>${esc(r.label)}</b><br><span class="t">${esc(segs[j])}</span><br>${fmt(v)} reads · ${pf(100 * v / tot)}`;
      svg += `<rect x="${x}" y="${y}" width="${Math.max(0, w - 2)}" height="16" rx="3" fill="${colors[j]}" data-h="${esc(h)}"/>`;
      if (w > 38) svg += `<text x="${x + w / 2 - 1}" y="${y + 12}" text-anchor="middle" style="fill:#fff;font-weight:600;pointer-events:none">${Math.round(100 * v / tot)}%</text>`;
      x += w; }); });
  el.innerHTML = svg + "</svg>";
  el.querySelectorAll("rect[data-h]").forEach(n => { n.addEventListener("mousemove", e => showTip(e, n.dataset.h)); n.addEventListener("mouseleave", hideTip); });
}
stacked($("fate"), $("fateLegend"), ["assigned to a sample", "tag or flank hit, not assigned", "no hit"], r => { const s = r.summary, t = s.total_reads || 0;
  const a = assigned(r) || 0, h = s.reads_with_ge1_tag_or_flank_hit ?? 0; return t ? [a, Math.max(0, h - a), Math.max(0, t - h)] : []; }, "no run summary in log");
const fk = Object.keys(D.foundLabels);
stacked($("found"), $("foundLegend"), fk.map(k => D.foundLabels[k]), r => fk.map(k => r.found[k] || 0), "no pattern table in log");

// ---- metrics table
const GOOD_UP = new Set(["reads_with_ge1_tag_or_flank_hit", "reads_assigned_to_samples", "pattern_reads_tag_both_ends", "samples_with_reads", "reads_in_ont_bins", "tag_both_ends"]);
const GOOD_DOWN = new Set(["pattern_reads_tag_one_end_only", "unexpected_reads", "unexpected_reads_decoy_wells", "unexpected_reads_other", "unexpected_or_unparsed_files", "tag_one_end_only", "flank_only", "barbell_seconds"]);
const isPctRow = k => /^(reads|pattern_reads|unexpected_reads)/.test(k);
function metricRows() {
  const rows = [], keys = [];
  R.forEach(r => Object.keys(r.summary).forEach(k => keys.includes(k) || keys.push(k)));
  rows.push(["sec", "run summary"]);
  for (const k of keys) { rows.push([k, k, R.map(r => r.summary[k]), false]);
    if (isPctRow(k) && k !== "total_reads") rows.push([k, k + " %", R.map(r => r.summary[k] === undefined ? undefined : pc(r.summary[k], r.summary.total_reads)), true]); }
  if (R.some(r => Object.keys(r.found).length)) { rows.push(["sec", "tags found per read"]);
    for (const k of fk) { rows.push([k, D.foundLabels[k], R.map(r => r.found[k]), false]);
      rows.push([k, D.foundLabels[k] + " %", R.map(r => r.found[k] === undefined ? undefined : pc(r.found[k], Object.values(r.found).reduce((a, b) => a + b, 0))), true]); } }
  rows.push(["sec", "run"]);
  rows.push(["unique_patterns", "unique read patterns", R.map(r => r.unique_patterns ?? undefined), false]);
  rows.push(["barbell_seconds", "barbell time (s)", R.map(r => r.barbell_seconds), false]);
  rows.push(["warnings", "warnings", R.map(r => r.warnings.length), false]);
  return rows;
}
function drawMetrics() {
  const bi = R.findIndex(r => r.label === base), diffOnly = $("diffOnly").checked, showPct = $("showPct").checked;
  let h = `<tr><th>metric</th>${R.map(r => `<th>${esc(r.label)}${r.label === base ? ' <span class="base">●</span>' : ""}</th>`).join("")}</tr>`;
  let pendingSec = null;
  for (const row of metricRows()) {
    if (row[0] === "sec") { pendingSec = row[1]; continue; }
    const [key, name, vals, isPct] = row;
    if (isPct && !showPct) continue;
    const defined = vals.filter(v => v !== undefined && v !== null);
    if (diffOnly && defined.every(v => Math.abs(v - defined[0]) < (isPct ? 0.05 : 1e-9))) continue;
    if (pendingSec) { h += `<tr class="sec"><td colspan="${R.length + 1}">${esc(pendingSec)}</td></tr>`; pendingSec = null; }
    const b = vals[bi];
    h += `<tr><td>${esc(name)}</td>` + vals.map((v, i) => { let cell = v === undefined || v === null ? "–" : (isPct ? pf(v) : fmt(v));
      if (i !== bi && v !== undefined && v !== null && b !== undefined && b !== null) { const d = v - b;
        if (Math.abs(d) >= (isPct ? 0.05 : 1e-9)) { const sameTotal = R[i].summary.total_reads === R[bi].summary.total_reads;
          const dir = !isPct && !sameTotal ? "neu" : GOOD_UP.has(key) ? (d > 0 ? "up" : "down") : GOOD_DOWN.has(key) ? (d < 0 ? "up" : "down") : "neu";
          const rel = !isPct && b ? ` ${d > 0 ? "+" : ""}${(100 * d / b).toFixed(1)}%` : "";
          cell += `<span class="d ${dir}">${d > 0 ? "▲" : "▼"} ${d > 0 ? "+" : ""}${isPct ? d.toFixed(1) + " pt" : fmt(Math.round(d * 10) / 10) + rel}</span>`; } }
      return `<td>${cell}</td>`; }).join("") + "</tr>"; }
  $("metrics").innerHTML = h;
  $("sub").textContent = `${R.length} runs · generated ${D.generated} · baseline ${base} · ▲▼ = difference to the baseline (green better, red worse, gray neutral; read counts are gray when the total reads differ)`;
}
$("baseSel").innerHTML = R.map(r => `<option ${r.label === base ? "selected" : ""}>${esc(r.label)}</option>`).join("");
$("baseSel").onchange = e => { base = e.target.value; drawMetrics(); };
$("diffOnly").onchange = drawMetrics; $("showPct").onchange = drawMetrics;
drawMetrics();

// ---- patterns
const qv = ["--q0", "--q1", "--q2", "--q3", "--q4", "--q5"];
function heat(f) { const i = f <= 0 ? 0 : Math.min(5, 1 + Math.floor(f * 5 / 0.4)); return `background:var(${qv[i]});${i >= 4 ? "color:#fff" : ""}`; }
function drawPats() {
  if (!D.patterns.length) { $("patPanel").classList.add("hidden"); return; }
  const mode = $("patMode").value;
  const tot = R.map(r => Object.values(r.found).reduce((a, b) => a + b, 0) || Object.values(r.patterns).reduce((a, b) => a + b, 0) || 1);
  let h = `<tr><th>pattern</th>${R.map(r => `<th>${esc(r.label)}</th>`).join("")}</tr>`;
  for (const p of D.patterns) { h += `<tr><td class="pat">${esc(p)}</td>` + R.map((r, i) => { const n = r.patterns[p];
      if (n === undefined) return `<td><span class="neu">–</span></td>`; const f = n / tot[i];
      return `<td><span class="cell" style="${heat(f)}" title="${fmt(n)} reads">${mode === "share" ? pf(100 * f) : fmt(n)}</span></td>`; }).join("") + "</tr>"; }
  $("pats").innerHTML = h + `<tr><td colspan="${R.length + 1}" class="neu" style="text-align:left;white-space:normal;font-size:11px">Only the top patterns are in a log; – = not in that run's top list.</td></tr>`;
}
$("patMode").onchange = drawPats; drawPats();

// ---- ONT bins
const bins = [...new Set(R.flatMap(r => Object.keys(r.bins)))].sort();
function drawBins() {
  const mode = $("binMode").value, br = R.filter(r => Object.keys(r.bins).length);
  let h = `<tr><th>ONT bin</th>${br.map(r => `<th>${esc(r.label)}</th>`).join("")}</tr>`;
  for (const b of bins) { h += `<tr><td>${esc(b)}</td>` + br.map(r => { const x = r.bins[b]; if (!x) return `<td class="neu">–</td>`;
      if (mode === "pct") { const f = x.reads ? x.assigned / x.reads : 0; return `<td><span class="cell" style="${heat(f * 0.4)}">${x.reads ? pf(100 * f) : "–"}</span></td>`; }
      if (mode === "samples") return `<td>${x.samples_with_reads === undefined ? "–" : fmt(x.samples_with_reads) + " / " + fmt(x.samples_expected)}</td>`;
      return `<td>${fmt(x[mode])}</td>`; }).join("") + "</tr>"; }
  $("bins").innerHTML = h;
}
if (bins.length) { $("binPanel").classList.remove("hidden"); $("binMode").onchange = drawBins; drawBins(); }

// ---- samples
const SR = R.filter(r => r.samples);
if (SR.length) {
  $("smpPanel").classList.remove("hidden");
  const opts = SR.map((r, i) => `<option value="${i}">${esc(r.label)}</option>`).join("");
  $("sx").innerHTML = opts; $("sy").innerHTML = opts;
  const bi = Math.max(0, SR.findIndex(r => r.label === D.base));
  $("sx").value = bi; $("sy").value = SR.length > 1 ? (bi === 0 ? 1 : 0) : 0;
  const keys = [...new Set(SR.flatMap(r => Object.keys(r.samples)))];
  let pts = [];
  function pearson(xs, ys) { const n = xs.length; if (n < 3) return null; const mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
    let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { const a = xs[i] - mx, b = ys[i] - my; sxy += a * b; sxx += a * a; syy += b * b; }
    return sxx && syy ? sxy / Math.sqrt(sxx * syy) : null; }
  function drawSamples() {
    const X = SR[+$("sx").value], Y = SR[+$("sy").value], log = $("scale").value === "log", hz = $("hideZero").checked;
    pts = keys.map(k => ({k, x: X.samples[k], y: Y.samples[k]})).filter(p => p.x !== undefined || p.y !== undefined)
      .map(p => ({...p, x: p.x ?? 0, y: p.y ?? 0, miss: p.x === undefined ? "x" : p.y === undefined ? "y" : ""}))
      .filter(p => !(hz && p.x === 0 && p.y === 0));
    const W = 560, H = 440, m = {l: 58, r: 14, t: 12, b: 42}, pw = W - m.l - m.r, ph = H - m.t - m.b;
    const mx = Math.max(1, ...pts.map(p => Math.max(p.x, p.y)));
    const tf = v => log ? Math.log10(v + 1) : v, tmax = tf(mx) * 1.04;
    const sx = v => m.l + pw * tf(v) / tmax, sy = v => m.t + ph - ph * tf(v) / tmax;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="reads per sample, ${esc(X.label)} vs ${esc(Y.label)}">`;
    const ticks = log ? [0, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7].filter(t => t <= mx * 1.05) : [0, .25, .5, .75, 1].map(f => Math.round(f * mx));
    for (const t of ticks) { svg += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(t)}" y2="${sy(t)}" stroke="var(--grid)"/><text x="${m.l - 6}" y="${sy(t) + 4}" text-anchor="end">${fmt(t)}</text>`;
      svg += `<line x1="${sx(t)}" x2="${sx(t)}" y1="${m.t}" y2="${m.t + ph}" stroke="var(--grid)"/><text x="${sx(t)}" y="${m.t + ph + 15}" text-anchor="middle">${fmt(t)}</text>`; }
    svg += `<line x1="${sx(0)}" y1="${sy(0)}" x2="${sx(mx * 1.04)}" y2="${sy(mx * 1.04)}" stroke="var(--axis)" stroke-dasharray="4 4"/>`;
    svg += `<text x="${m.l + pw / 2}" y="${H - 6}" text-anchor="middle" style="fill:var(--ink)">${esc(X.label)} (reads)</text>`;
    svg += `<text transform="translate(14 ${m.t + ph / 2}) rotate(-90)" text-anchor="middle" style="fill:var(--ink)">${esc(Y.label)} (reads)</text>`;
    for (const p of pts) svg += `<circle cx="${sx(p.x).toFixed(1)}" cy="${sy(p.y).toFixed(1)}" r="4" fill="var(--s1)" fill-opacity="${p.miss ? .25 : .7}" stroke="var(--panel)" stroke-width="1" data-k="${esc(p.k)}"/>`;
    $("scatter").innerHTML = svg + "</svg>";
    $("scatter").querySelectorAll("circle").forEach(n => { const p = pts.find(q => q.k === n.dataset.k);
      n.addEventListener("mousemove", e => showTip(e, `<b>${esc(p.k)}</b><br><span class="t">${esc(X.label)}</span> ${fmt(p.x)}<br><span class="t">${esc(Y.label)}</span> ${fmt(p.y)}${p.miss ? `<br><span class="t">not in ${esc((p.miss === "x" ? X : Y).label)}</span>` : ""}`));
      n.addEventListener("mouseleave", hideTip); });
    const both = pts.filter(p => !p.miss), r = pearson(both.map(p => Math.log10(p.x + 1)), both.map(p => Math.log10(p.y + 1)));
    const tx = both.reduce((a, p) => a + p.x, 0), ty = both.reduce((a, p) => a + p.y, 0);
    $("sstat").innerHTML = `${fmt(pts.length)} samples · ${fmt(both.filter(p => p.x > 0).length)} with reads in ${esc(X.label)}, ${fmt(both.filter(p => p.y > 0).length)} in ${esc(Y.label)}` +
      ` · Pearson r (log reads) ${r === null ? "–" : r.toFixed(3)} · total ${fmt(tx)} vs ${fmt(ty)}` + (pts.some(p => p.miss) ? " · faint dots: sample missing in one run" : "");
    drawSTable(X, Y);
  }
  function drawSTable(X, Y) {
    const q = $("q").value.toLowerCase(), by = $("sortBy").value;
    const fc = p => Math.log2((p.y + 1) / (p.x + 1));
    let rows = pts.filter(p => p.k.toLowerCase().includes(q));
    rows.sort(by === "diff" ? (a, b) => Math.abs(b.y - b.x) - Math.abs(a.y - a.x) : by === "fc" ? (a, b) => Math.abs(fc(b)) - Math.abs(fc(a)) : by === "x" ? (a, b) => b.x - a.x : by === "y" ? (a, b) => b.y - a.y : (a, b) => a.k.localeCompare(b.k, undefined, {numeric: true}));
    const shown = rows.slice(0, 500);
    $("stab").innerHTML = `<tr><th>sample</th><th>${esc(X.label)}</th><th>${esc(Y.label)}</th><th>log2 fold change</th></tr>` +
      shown.map(p => { const f = fc(p); return `<tr><td>${esc(p.k)}</td><td>${fmt(p.x)}</td><td>${fmt(p.y)}</td><td class="${Math.abs(f) < .01 ? "neu" : ""}">${f > 0 ? "+" : ""}${f.toFixed(2)}</td></tr>`; }).join("") +
      (rows.length > 500 ? `<tr><td colspan="4" class="neu">${fmt(rows.length - 500)} more – use the filter</td></tr>` : "");
  }
  ["sx", "sy", "scale", "hideZero"].forEach(id => $(id).onchange = drawSamples);
  $("q").oninput = () => drawSTable(SR[+$("sx").value], SR[+$("sy").value]);
  $("sortBy").onchange = $("q").oninput;
  drawSamples();
}

// ---- negative controls: wells marked in the page, analysed per run (selection kept in this browser only)
if (SR.length) (() => {
  $("negPanel").classList.remove("hidden");
  const M = D.meta, keys = [...new Set(SR.flatMap(r => Object.keys(r.samples)))];
  const SKEY = "compare-runs:neg:" + document.title + ":" + keys.length;
  let st = {}; try { st = JSON.parse(localStorage.getItem(SKEY) || "{}") || {}; } catch (e) {}
  const save = () => { try { localStorage.setItem(SKEY, JSON.stringify(st)); } catch (e) {} };
  let mode = "neg";
  const stateOf = k => st[k] || "pos";
  // plates from the summary.tsv plate/well columns; samples without a well go into an "other" panel
  const WRE = /^([A-Pa-p])0?(\d{1,2})$/, panels = new Map();
  for (const k of keys) { const m = M[k] || [], pw = WRE.exec(m[1] || ""), pn = pw && m[0] ? m[0] : "other";
    if (!panels.has(pn)) panels.set(pn, []);
    panels.get(pn).push({k, r: pw ? pw[1].toUpperCase().charCodeAt(0) - 65 : -1, c: pw ? +pw[2] - 1 : -1}); }
  const RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"].map(h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16)));
  function ramp(t) { t = Math.max(0, Math.min(1, t)) * (RAMP.length - 1); const i = Math.min(RAMP.length - 2, Math.floor(t)), f = t - i;
    return "rgb(" + RAMP[i].map((v, j) => Math.round(v + (RAMP[i + 1][j] - v) * f)).join(",") + ")"; }
  $("ncol").innerHTML = SR.map((r, i) => `<option value="${i}">${esc(r.label)}</option>`).join("");
  $("ncol").value = Math.max(0, SR.findIndex(r => r.label === D.base));

  function drawPlates() {
    let h = "";
    [...panels.entries()].sort((a, b) => a[0] === "other" ? 1 : b[0] === "other" ? -1 : a[0].localeCompare(b[0], undefined, {numeric: true})).forEach(([pn, ws], pi) => {
      const nr = pn === "other" ? 0 : Math.max(8, ...ws.map(w => w.r + 1)), nc = pn === "other" ? 0 : Math.max(12, ...ws.map(w => w.c + 1));
      const cs = nc > 12 ? 11 : 16, nN = ws.filter(w => stateOf(w.k) === "neg").length, nE = ws.filter(w => stateOf(w.k) === "exc").length;
      h += `<div class="pl" data-p="${pi}"><div class="pt" data-act="plate">${esc(pn)}<small>${nN ? nN + " neg" : ""}${nN && nE ? " · " : ""}${nE ? nE + " excl" : ""}</small></div>`;
      if (pn === "other") { h += `<div class="pg" style="--cs:${cs}px;grid-template-columns:repeat(24,${cs}px)">` + ws.map(w => `<div class="w" data-k="${esc(w.k)}"></div>`).join("") + "</div></div>"; return; }
      const at = new Map(ws.map(w => [w.r + "_" + w.c, w]));
      h += `<div class="pg" style="--cs:${cs}px;grid-template-columns:12px repeat(${nc},${cs}px)"><div></div>`;
      for (let c = 0; c < nc; c++) h += `<div class="h" data-act="col" data-c="${c}">${c + 1}</div>`;
      for (let r = 0; r < nr; r++) { h += `<div class="h" data-act="row" data-r="${r}">${String.fromCharCode(65 + r)}</div>`;
        for (let c = 0; c < nc; c++) { const w = at.get(r + "_" + c); h += w ? `<div class="w" data-k="${esc(w.k)}" data-r="${r}" data-c="${c}"></div>` : `<div class="w none"></div>`; } }
      h += "</div></div>";
    });
    $("nplates").innerHTML = h; paint();
  }
  function paint() {
    const run = SR[+$("ncol").value], mx = Math.max(1, ...Object.values(run.samples));
    $("nplates").querySelectorAll(".w[data-k]").forEach(el => { const k = el.dataset.k, v = run.samples[k], s = stateOf(k);
      el.style.background = v > 0 ? ramp(Math.log10(v + 1) / Math.log10(mx + 1)) : "";
      el.classList.toggle("neg", s === "neg"); el.classList.toggle("exc", s === "exc"); });
    $("nplates").querySelectorAll(".pl").forEach(pl => { const ks = [...pl.querySelectorAll(".w[data-k]")].map(e => e.dataset.k);
      const nN = ks.filter(k => stateOf(k) === "neg").length, nE = ks.filter(k => stateOf(k) === "exc").length;
      pl.querySelector(".pt small").textContent = (nN ? nN + " neg" : "") + (nN && nE ? " · " : "") + (nE ? nE + " excl" : ""); });
  }
  function apply(ks, toggle) {             // all already in the active mode + toggle -> back to positive
    const m = toggle && ks.length && ks.every(k => stateOf(k) === mode) ? "pos" : mode;
    for (const k of ks) { if (m === "pos") delete st[k]; else st[k] = m; }
    save(); paint(); analyse();
  }
  // mouse: click a well, drag a box inside one plate, or click a plate name / row / column header
  let drag = null;
  const box = (pl, a, b) => [...pl.querySelectorAll(".w[data-k]")].filter(e => e.dataset.r !== undefined &&
    +e.dataset.r >= Math.min(a.r, b.r) && +e.dataset.r <= Math.max(a.r, b.r) && +e.dataset.c >= Math.min(a.c, b.c) && +e.dataset.c <= Math.max(a.c, b.c));
  $("nplates").addEventListener("mousedown", e => { const w = e.target.closest(".w[data-k]"); if (!w) return; e.preventDefault();
    drag = {pl: w.closest(".pl"), a: {r: +w.dataset.r, c: +w.dataset.c}, b: {r: +w.dataset.r, c: +w.dataset.c}, el: w}; });
  $("nplates").addEventListener("mouseover", e => { const w = e.target.closest(".w[data-k]");
    if (w) { const k = w.dataset.k, m = M[k] || [];
      showTip(e, `<b>${esc(k)}</b> <span class="t">${esc(stateOf(k) === "pos" ? "positive" : stateOf(k) === "neg" ? "negative" : "excluded")}</span>` +
        (m[2] || m[3] ? `<br><span class="t">utags</span> ${esc(m[2])} / ${esc(m[3])}` : "") +
        SR.map(r => `<br><span class="t">${esc(r.label)}</span> ${fmt(r.samples[k])}`).join("")); } else hideTip();
    if (!drag || !w || w.closest(".pl") !== drag.pl || w.dataset.r === undefined) return;
    drag.b = {r: +w.dataset.r, c: +w.dataset.c};
    drag.pl.querySelectorAll(".w.pre").forEach(x => x.classList.remove("pre")); box(drag.pl, drag.a, drag.b).forEach(x => x.classList.add("pre")); });
  $("nplates").addEventListener("mouseleave", hideTip);
  document.addEventListener("mouseup", () => { if (!drag) return; const d = drag; drag = null;
    d.pl.querySelectorAll(".w.pre").forEach(x => x.classList.remove("pre"));
    const ks = d.el.dataset.r === undefined ? [d.el.dataset.k] : box(d.pl, d.a, d.b).map(x => x.dataset.k);
    apply(ks, true); });
  $("nplates").addEventListener("click", e => { const h = e.target.closest("[data-act]"); if (!h) return; const pl = h.closest(".pl");
    const ws = [...pl.querySelectorAll(".w[data-k]")].filter(x => h.dataset.act === "plate" || (h.dataset.act === "row" ? x.dataset.r === h.dataset.r : x.dataset.c === h.dataset.c));
    apply(ws.map(x => x.dataset.k), true); });
  $("nmodes").addEventListener("click", e => { const b = e.target.closest("button"); if (!b) return; mode = b.dataset.m;
    $("nmodes").querySelectorAll("button").forEach(x => x.classList.toggle("on", x === b)); });
  // selectors: sample names or plate names, * and ? wildcards, case-insensitive
  const glob = g => new RegExp("^" + g.replace(/[.+^${}()|[\]\\]/g, "\\$&").replace(/\*/g, ".*").replace(/\?/g, ".") + "$", "i");
  $("nselGo").onclick = () => { const res = $("nsel").value.split(/[\s,;]+/).filter(Boolean).map(glob);
    const ks = keys.filter(k => res.some(re => re.test(k) || re.test((M[k] || [])[0] || ""))); apply(ks, false); };
  $("nsel").addEventListener("keydown", e => { if (e.key === "Enter") $("nselGo").onclick(); });
  $("nclear").onclick = () => { st = {}; save(); paint(); analyse(); };
  $("ncol").onchange = paint;
  $("ncut").oninput = () => analyse();

  // ---- statistics
  function auc(pos, neg) {                // Mann-Whitney: P(positive well > negative well), ties count half
    const all = pos.map(v => [v, 1]).concat(neg.map(v => [v, 0])).sort((a, b) => a[0] - b[0]);
    let rs = 0;
    for (let i = 0; i < all.length;) { let j = i; while (j < all.length && all[j][0] === all[i][0]) j++;
      const rank = (i + 1 + j) / 2; for (let k = i; k < j; k++) if (all[k][1]) rs += rank; i = j; }
    return (rs - pos.length * (pos.length + 1) / 2) / (pos.length * neg.length);
  }
  const callAt = (pos, neg, t) => { const tp = pos.filter(v => v >= t).length, fp = neg.filter(v => v >= t).length;
    return {t, tp, fn: pos.length - tp, fp, tn: neg.length - fp, sens: tp / pos.length, spec: 1 - fp / neg.length}; };
  function stats(run, cut) {
    const ks = keys.filter(k => run.samples[k] !== undefined);
    const pos = ks.filter(k => stateOf(k) === "pos").map(k => run.samples[k]), neg = ks.filter(k => stateOf(k) === "neg").map(k => run.samples[k]);
    if (!pos.length || !neg.length) return null;
    const total = ks.reduce((a, k) => a + run.samples[k], 0), negReads = neg.reduce((a, b) => a + b, 0), maxNeg = Math.max(...neg);
    const cands = [...new Set([1, maxNeg + 1, ...pos, ...neg].filter(v => v >= 1))].sort((a, b) => a - b);
    let best = null; for (const t of cands) { const c = callAt(pos, neg, t); if (!best || c.sens + c.spec > best.sens + best.spec) best = c; }
    return {pos, neg, total, negReads, maxNeg, auc: auc(pos, neg), clean: callAt(pos, neg, maxNeg + 1), best, at: callAt(pos, neg, cut)};
  }
  // what the negatives can detect, from their utags: tags not used by any positive well catch misreads of that tag only
  function controlKind() {
    const pk = keys.filter(k => stateOf(k) === "pos"), nk = keys.filter(k => stateOf(k) === "neg");
    const uf = new Set(pk.map(k => (M[k] || [])[2]).filter(Boolean)), ur = new Set(pk.map(k => (M[k] || [])[3]).filter(Boolean));
    if (!nk.length || !uf.size) return null;
    const freeF = nk.filter(k => !uf.has((M[k] || [])[2])), freeR = nk.filter(k => !ur.has((M[k] || [])[3]));
    const kind = freeF.length === nk.length && !freeR.length ? "fwd" : freeR.length === nk.length && !freeF.length ? "rev"
      : !freeF.length && !freeR.length ? "blank" : "mixed";
    const kF = new Set(freeF.map(k => M[k][2])).size, kR = new Set(freeR.map(k => M[k][3])).size;
    return {kind, factor: kind === "fwd" ? (uf.size - 1 + kF) / kF : kind === "rev" ? (ur.size - 1 + kR) / kR : null,
      k: kind === "fwd" ? kF : kR, m: kind === "fwd" ? uf.size : ur.size};
  }
  let lastRows = [];
  function analyse() {
    const cut = Math.max(1, +$("ncut").value || 1), S = SR.map(r => stats(r, cut)), ck = controlKind();
    const nN = keys.filter(k => stateOf(k) === "neg").length, nE = keys.filter(k => stateOf(k) === "exc").length;
    if (!nN) { $("nres").innerHTML = `<div class="note">Mark at least one negative well above to see sensitivity and specificity.</div>`; lastRows = []; return; }
    const pp = v => v === null || v === undefined ? "–" : (100 * v).toFixed(1) + "%";
    const rows = [
      ["positive wells", s => fmt(s.pos.length)], ["negative wells", s => fmt(s.neg.length)],
      ["AUC (positive well > negative well)", s => s.auc.toFixed(4), true],
      ["reads in negative wells", s => fmt(s.negReads)],
      ["… per million assigned reads", s => fmt(Math.round(1e6 * s.negReads / s.total))],
      ["max reads in one negative well", s => fmt(s.maxNeg), true],
      [`at ≥ ${cut} reads: sensitivity`, s => pp(s.at.sens), true], [`at ≥ ${cut} reads: specificity`, s => pp(s.at.spec), true],
      [`at ≥ ${cut} reads: TP / FN / TN / FP`, s => `${s.at.tp} / ${s.at.fn} / ${s.at.tn} / ${s.at.fp}`],
      ["cutoff for 100% specificity", s => "≥ " + fmt(s.clean.t)], ["… sensitivity there", s => pp(s.clean.sens), true],
      ["best cutoff (Youden)", s => "≥ " + fmt(s.best.t)], ["… sensitivity / specificity", s => pp(s.best.sens) + " / " + pp(s.best.spec)],
    ];
    if (ck && ck.factor) rows.push([`est. all ${ck.kind}-tag misassignment (× ${ck.factor.toFixed(2)})`, s => (100 * s.negReads * ck.factor / s.total).toFixed(3) + "%"]);
    lastRows = [["metric", ...SR.map(r => r.label)], ...rows.map(([n, f]) => [n, ...S.map(s => s ? f(s).replace(/^≥ /, ">=") : "")])];
    let h = `<div class="scroll"><table class="kpi"><tr><th>metric</th>${SR.map(r => `<th>${esc(r.label)}</th>`).join("")}</tr>` +
      rows.map(([n, f, b]) => `<tr><td>${esc(n)}</td>${S.map(s => `<td class="${b ? "b" : ""}">${s ? f(s) : '<span class="neu">–</span>'}</td>`).join("")}</tr>`).join("") + "</table></div>";
    const kindTxt = !ck ? "" : ck.kind === "fwd" || ck.kind === "rev"
      ? `Negatives use ${ck.k} ${ck.kind} utag(s) that no positive well uses: they catch misreads of the ${ck.kind} tag only, so specificity here is a lower bound (${ck.kind === "fwd" ? "rev" : "fwd"}-tag errors land in positive wells unseen). Estimate: a misread of one of the ${ck.m} used tags lands on an unused one in ${ck.k}/${ck.m - 1 + ck.k} of cases.`
      : ck.kind === "blank" ? "Negatives use utags that positive wells use too (in-design blanks): they catch errors on either tag and physical contamination."
      : "Negatives mix unused-tag wells and in-design blanks; no misassignment estimate.";
    h += `<div class="note">${nN} negative, ${nE} excluded, ${keys.length - nN - nE} positive wells. ${kindTxt} Wells with 0 reads are not called at any cutoff.</div>`;
    h += `<div class="two" style="margin-top:12px"><div><div class="legend"><span><i style="background:var(--s1)"></i>positive wells</span><span><i style="background:var(--s2)"></i>negative wells</span><span style="color:var(--muted)">│ cutoff</span></div><div id="nstrip"></div></div><div id="ntags"></div></div>`;
    h += `<div id="nlow" style="margin-top:12px"></div>`;
    $("nres").innerHTML = h;
    strip(S, cut); tags(); low(cut);
  }
  function strip(S, cut) {                // reads per well, one row per run, log scale
    const lw = Math.min(150, 12 + 7 * Math.max(...SR.map(r => r.label.length))), W = 620, rh = 44, H = SR.length * rh + 30, pw = W - lw - 16;
    const mx = Math.max(10, ...SR.flatMap(r => Object.values(r.samples))), x = v => lw + pw * Math.log10(v + 1) / Math.log10(mx + 1);
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="reads per well, positives and negatives">`;
    for (const t of [0, 1, 10, 100, 1e3, 1e4, 1e5].filter(t => t <= mx)) svg += `<line x1="${x(t)}" x2="${x(t)}" y1="0" y2="${H - 18}" stroke="var(--grid)"/><text x="${x(t)}" y="${H - 5}" text-anchor="middle">${fmt(t)}</text>`;
    let pts = [];
    SR.forEach((r, i) => { const y0 = i * rh + 6;
      svg += `<text x="${lw - 8}" y="${y0 + 18}" text-anchor="end" style="fill:var(--ink)">${esc(r.label.length > 20 ? "…" + r.label.slice(-19) : r.label)}</text>`;
      keys.forEach((k, j) => { const s = stateOf(k), v = r.samples[k]; if (s === "exc" || v === undefined) return;
        const jit = ((j * 2654435761) % 1000) / 1000, y = y0 + (s === "neg" ? 22 + jit * 10 : 2 + jit * 16);
        pts.push(`<circle cx="${x(v).toFixed(1)}" cy="${y.toFixed(1)}" r="${s === "neg" ? 3.5 : 2.5}" fill="var(${s === "neg" ? "--s2" : "--s1"})" fill-opacity="${s === "neg" ? .9 : .5}" data-h="${esc(`<b>${esc(k)}</b><br><span class="t">${esc(r.label)}</span> ${fmt(v)} reads · ${s === "neg" ? "negative" : "positive"}`)}"/>`); }); });
    svg += pts.join("") + `<line x1="${x(cut - .5)}" x2="${x(cut - .5)}" y1="0" y2="${H - 18}" stroke="var(--ink)" stroke-width="1.5" stroke-dasharray="3 3"/>`;
    $("nstrip").innerHTML = svg + "</svg>";
    $("nstrip").querySelectorAll("circle").forEach(n => { n.addEventListener("mousemove", e => showTip(e, n.dataset.h)); n.addEventListener("mouseleave", hideTip); });
  }
  function tags() {                       // which utags collect the reads in negative wells
    const nk = keys.filter(k => stateOf(k) === "neg"), acc = new Map();
    for (const k of nk) { const m = M[k] || []; for (const [side, t] of [["fwd", m[2]], ["rev", m[3]]]) { if (!t) continue;
      const id = side + " " + t; if (!acc.has(id)) acc.set(id, SR.map(() => 0)); SR.forEach((r, i) => acc.get(id)[i] += r.samples[k] || 0); } }
    const rows = [...acc.entries()].sort((a, b) => b[1].reduce((x, y) => x + y, 0) - a[1].reduce((x, y) => x + y, 0)).slice(0, 12);
    $("ntags").innerHTML = rows.length ? `<div class="note" style="margin:0 0 4px">reads in negative wells per utag (top ${rows.length})</div><div class="scroll" style="max-height:300px;overflow-y:auto"><table><tr><th>utag</th>${SR.map(r => `<th>${esc(r.label)}</th>`).join("")}</tr>` +
      rows.map(([id, v]) => `<tr><td>${esc(id)}</td>${v.map(n => `<td>${fmt(n)}</td>`).join("")}</tr>`).join("") + "</table></div>" : "";
  }
  function low(cut) {                     // positives below the cutoff in every run: likely real dropouts
    const ks = keys.filter(k => stateOf(k) === "pos" && SR.every(r => r.samples[k] === undefined || r.samples[k] < cut));
    if (!ks.length) { $("nlow").innerHTML = ""; return; }
    $("nlow").innerHTML = `<div class="note" style="margin:0 0 4px">${ks.length} positive well(s) below ${cut} reads in every run: more likely failed samples than demux misses.
      <button id="nlowAll">exclude all</button></div><div class="scroll lowlist" style="max-height:220px;overflow-y:auto"><table><tr><th>sample</th>${SR.map(r => `<th>${esc(r.label)}</th>`).join("")}<th></th></tr>` +
      ks.slice(0, 200).map(k => `<tr><td>${esc(k)}</td>${SR.map(r => `<td>${fmt(r.samples[k])}</td>`).join("")}<td><button data-x="${esc(k)}">exclude</button></td></tr>`).join("") + "</table></div>";
    $("nlowAll").onclick = () => { ks.forEach(k => st[k] = "exc"); save(); paint(); analyse(); };
    $("nlow").querySelectorAll("button[data-x]").forEach(b => b.onclick = () => { st[b.dataset.x] = "exc"; save(); paint(); analyse(); });
  }
  $("ndl").onclick = () => {
    const q = v => /[",\n]/.test(v) ? '"' + String(v).replace(/"/g, '""') + '"' : v;
    const lines = lastRows.map(r => r.map(q).join(","));
    lines.push("", "sample,state");
    for (const k of keys) if (stateOf(k) !== "pos") lines.push(q(k) + "," + (stateOf(k) === "neg" ? "negative" : "excluded"));
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], {type: "text/csv"}));
    a.download = "negative_controls.csv"; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  drawPlates(); analyse();
})();

// ---- commands
$("cmds").innerHTML = R.map(r => `<details ${r.status !== "ok" ? "open" : ""} style="margin-bottom:8px"><summary><b style="color:var(--ink)">${esc(r.label)}</b><span class="badge ${r.status}">${r.status}</span> · ${esc(r.log)}</summary>
  <div style="padding:6px 0 4px 14px"><div class="mono">$ ${esc(r.command)}</div><div class="k" style="color:var(--muted);font-size:11px">cwd ${esc(r.cwd)} · started ${esc(r.started)}</div>
  ${Object.keys(r.steps).length ? `<div class="stat">${Object.entries(r.steps).map(([s, t]) => `${esc(s)} <b>${t.toFixed(1)} s</b>`).join(" · ")}</div>` : ""}
  ${r.info.map(s => `<div class="stat">• ${esc(s)}</div>`).join("")}
  ${r.warnings.map(s => `<div class="warn">${esc(s)}</div>`).join("")}
  ${r.errors.map(s => `<div class="err mono">${esc(s)}</div>`).join("")}</div></details>`).join("");
</script></body></html>
"""

if __name__ == "__main__":
    main()

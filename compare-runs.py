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
    run_comparison.html   ONE self-contained page (no internet): grouped metrics table, top patterns, ONT bins, per-sample scatter of any two runs, and
                          negative controls: mark empty wells (e.g. an unused plate) on the plate grids and read
                          sensitivity / specificity per run (the selection stays in the browser; results download as CSV)

With one run the page is a run summary (no baseline or differences; reads per sample ranked): utag.py demux
writes it as run_summary.html through render_html().

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
BIN1_RE = re.compile(rf"^ {{4}}(\S+)\s+{NUM}\s+([\d.]+)%\s+[█░]+\s*$")     # twopass --pass1-only: bin, reads, %, bar
BIN_INFO_RE = re.compile(rf"{NUM} reads in bin, {NUM} assigned")
TIME_RE = re.compile(r"^\s+✓ ([\d.]+)s")
UNIQ_RE = re.compile(r"pattern \(top \d+ of (\d+)\)")


def num(s):
    return int(s.replace(",", ""))


def parse_log(path):
    r = {"log": path, "started": "", "cwd": "", "command": "", "subcommand": "", "context": "", "kit": "",
         "status": "incomplete", "barbell_seconds": 0.0, "steps": {}, "warnings": [], "errors": [],
         "summary": {}, "found": {}, "patterns": {}, "unique_patterns": None, "bins": {}, "info": [], "pass2": [],
         "pass1": None, "pattern_reads": None}
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
            if step.startswith("[pass 1] "):     # twopass: how the reads were split by ONT barcode
                r["pass1"] = {"method": step[9:], "lines": []}
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
        if step.startswith("[pass 1] ") and s.startswith(("• ", "! ")):
            r["pass1"]["lines"].append(s)
        if s.startswith("! "):
            r["warnings"].append(s[2:])
            continue
        if s.startswith(("ERROR", "Traceback")) or s.startswith("error:"):
            r["errors"].append(s)
            continue
        if s.startswith("• "):
            r["info"].append(s[2:])
            m = re.match(rf"• {NUM} reads with a pattern", s)      # twopass --pass1-only: the share denominator
            if m:
                r["pattern_reads"] = num(m.group(1))
            m = BIN_INFO_RE.search(s)
            if m and pass2_bin:
                r["bins"].setdefault(pass2_bin, {}).update(reads=num(m.group(1)), assigned=num(m.group(2)))
            continue
        if step.startswith(("[3/6] inspect", "pass 1 inspect")):   # demux / twopass --pass1-only patterns
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
        if step == "pass 1 summary":             # twopass --pass1-only
            r["status"] = "ok"
            m = BIN1_RE.match(ln)
            if m:
                r["bins"][m.group(1)] = {"reads": num(m.group(2))}
                continue
            m = SUMMARY_RE.match(ln)
            if m:
                r["summary"][ALIAS.get(m.group(1), m.group(1))] = num(m.group(2))
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
    if r["subcommand"] == "twopass":            # the per-bin demux logs (pass 2), if still there: kept per bin,
        sub = sorted(glob.glob(os.path.join(os.path.dirname(path), "barbell", "pass2_utag", "*", "utag.log")))
        summed = not r["patterns"]               # and their patterns / tag counts summed for the whole run
        for p in sub:
            b = parse_log(p)
            b["label"] = os.path.basename(os.path.dirname(p))
            r["pass2"].append(b)
            if summed:
                for k, v in b["found"].items():
                    r["found"][k] = r["found"].get(k, 0) + v
                for k, v in b["patterns"].items():
                    r["patterns"][k] = r["patterns"].get(k, 0) + v
                r["barbell_seconds"] += b["barbell_seconds"]
        if sub and summed:
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


def harmonize_bins(runs):
    """twopass sample keys are <ONT bin>/<sample>: name a bin the same in every run when only the naming differs
    (barbell NB01 = Dorado barcode01, matched by number), using the name of the first run that has the bin"""
    num = lambda b: int(re.search(r"(\d+)$", b).group(1)) if re.search(r"(\d+)$", b) else b
    first = {}
    for r in runs:
        for k in r["samples"] or {}:
            if "/" in k:
                b = k.split("/", 1)[0]
                first.setdefault(num(b), b)
    for r in runs:
        if not r["samples"]:
            continue
        ren = lambda k: first.get(num(k.split("/", 1)[0]), k.split("/", 1)[0]) + "/" + k.split("/", 1)[1] if "/" in k else k
        r["samples"] = {ren(k): v for k, v in r["samples"].items()}
        r["meta"] = {ren(k): [ren(m[0])] + m[1:] if m and "/" in m[0] else m for k, m in r["meta"].items()}


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
    ap.add_argument("--baseline", help="run label the CSV differences are taken against, and the x run of the sample "
                                       "scatter (default: the first run)")
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
    if len(runs) > 1:
        harmonize_bins(runs)
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
        "foundLabels": FOUND_LABEL, "patterns": top_patterns(runs + [b for r in runs for b in r["pass2"]], top),
        "runs": [{**{k: r[k] for k in ("label", "log", "started", "cwd", "command", "subcommand", "context", "kit",
                                       "status", "steps", "warnings", "errors", "summary", "found", "patterns",
                                       "unique_patterns", "bins", "info", "samples", "pass1", "pattern_reads")},
                  "barbell_seconds": round(r["barbell_seconds"], 1),
                  "pass2": [{**{k: b[k] for k in ("label", "status", "steps", "warnings", "errors", "summary", "found",
                                                  "patterns", "unique_patterns", "info")},
                             "barbell_seconds": round(b["barbell_seconds"], 1)} for b in r["pass2"]]} for r in runs],
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
  --q0:#f1f0ec;--q1:#cde2fb;--q2:#9ec5f4;--q3:#6da7ec;--q4:#3987e5;--q5:#256abf;--hbg:#e7eef8;--hink:#1c3d66}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])){color-scheme:dark;--bg:#0d0d0d;--panel:#1a1a19;
  --panel2:#252523;--ink:#f0efec;--dim:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;--ring:rgba(255,255,255,.10);
  --s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--good:#0ca30c;--bad:#e66767;
  --q0:#252523;--q1:#0d366b;--q2:#104281;--q3:#184f95;--q4:#1c5cab;--q5:#256abf;--hbg:#1c2a3d;--hink:#b7cdea}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#0d0d0d;--panel:#1a1a19;--panel2:#252523;--ink:#f0efec;--dim:#c3c2b7;
  --muted:#898781;--grid:#2c2c2a;--axis:#383835;--ring:rgba(255,255,255,.10);
  --s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--good:#0ca30c;--bad:#e66767;
  --q0:#252523;--q1:#0d366b;--q2:#104281;--q3:#184f95;--q4:#1c5cab;--q5:#256abf;--hbg:#1c2a3d;--hink:#b7cdea}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{padding:18px 24px 12px;border-bottom:1px solid var(--grid)}
h1{margin:0;font-size:19px;font-weight:650}
.sub{color:var(--dim);font-size:12px;margin-top:2px}
main{max-width:1280px;margin:0 auto;padding:16px 24px 48px;display:grid;gap:16px}
.panel{background:var(--panel);border:1px solid var(--ring);border-radius:12px;padding:14px 16px 16px;min-width:0}
h2{margin:-14px -16px 12px;padding:8px 16px;min-height:38px;background:var(--hbg);border-bottom:1px solid var(--ring);border-radius:12px 12px 0 0;
  font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--hink);font-weight:650;display:flex;gap:12px;align-items:center;flex-wrap:wrap}
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
td.pat{font-family:ui-monospace,Menlo,monospace;font-size:11.5px;text-align:left}th.pat{text-align:left}
tr.sec td{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.05em;padding-top:12px;border-bottom:none}
.d{font-size:11px;margin-left:6px}.up{color:var(--good)}.down{color:var(--bad)}.neu{color:var(--muted)}
th.sortable{cursor:pointer}th.sortable:hover{color:var(--ink)}
.legend{display:flex;gap:4px 12px;flex-wrap:wrap;color:var(--dim);font-size:11.5px;margin:0 0 8px}
.legend span{white-space:nowrap}
.legend i{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:4px;vertical-align:-1px}
svg text{fill:var(--dim);font-size:11px}
.tip{position:fixed;pointer-events:none;background:var(--panel);color:var(--ink);border:1px solid var(--ring);border-radius:8px;
  padding:6px 9px;font-size:12px;box-shadow:0 4px 16px rgba(0,0,0,.18);display:none;z-index:9;max-width:420px}
.tip .t{color:var(--dim)}
.two{display:grid;grid-template-columns:minmax(0,1.25fr) minmax(0,1fr);gap:16px}
@media (max-width:900px){.two{grid-template-columns:1fr}header,main{padding-left:16px;padding-right:16px}}
.mono{font-family:ui-monospace,Menlo,monospace;font-size:11.5px;white-space:pre-wrap;word-break:break-all;color:var(--dim)}
.warn{color:var(--dim)}.warn::before{content:"! ";color:var(--s4);font-weight:700}
.err::before{content:"✕ ";color:var(--bad);font-weight:700}
.stat{color:var(--dim);font-size:12px;margin-top:6px}
.cell{border-radius:4px;padding:1px 6px}
details summary{cursor:pointer;color:var(--dim)}
.hidden{display:none}body.single .multi{display:none}
.negbar{display:flex;gap:8px 14px;flex-wrap:wrap;align-items:center;margin-bottom:10px}
.modes{display:inline-flex;gap:4px}.modes button.on{border-color:var(--ink);box-shadow:inset 0 0 0 1px var(--ink)}
.modes i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px;vertical-align:0}
.plates{display:flex;flex-wrap:wrap;gap:14px;margin:6px 0 4px}
#nplates{display:block}.bgrp+.bgrp{margin-top:10px;padding-top:8px;border-top:1px solid var(--grid)}
.bt{font-weight:650;font-size:12px;color:var(--hink);letter-spacing:.04em}
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
.kpi td.b{font-weight:600}table.kpi{width:auto}.kpi th:first-child,.kpi td:first-child{min-width:280px}
.kpi th:not(:first-child),.kpi td:not(:first-child){min-width:110px}   /* numbers next to the metric, as in the metrics table */
table.narrow{width:auto}.narrow th:first-child,.narrow td:first-child{min-width:160px}.narrow th:not(:first-child),.narrow td:not(:first-child){min-width:90px}
.p1m{font-size:14px;font-weight:600}.p1m .tag{display:inline-block;font-size:11px;font-weight:600;padding:1px 8px;border-radius:9px;margin-right:8px;
  background:var(--hbg);color:var(--hink);vertical-align:1px}.p1l{color:var(--dim);margin-top:3px}.p1l.pw::before{content:"! ";color:var(--s4);font-weight:700}
#metrics table.term{width:auto}   /* laid out as the terminal run summary */
.term th:first-child,.term td:first-child{padding-right:24px}.term th.run{text-align:center;color:var(--ink)}
.term td.of{color:var(--dim)}.term td.gap,.term th.gap{padding-right:18px}
.note{color:var(--muted);font-size:11.5px;margin-top:6px;line-height:1.5}
</style></head><body>
<header><h1>__TITLE__</h1><div class="sub" id="sub"></div></header>
<main>
  <section class="panel hidden" id="p1Panel"><h2>Pass 1: ONT barcodes</h2><div id="p1"></div></section>
  <section class="panel"><h2>Metrics<span class="ctl">
      <label class="multi"><input type="checkbox" id="diffOnly"> differences only</label>
      <label><input type="checkbox" id="showPct" checked> show % of total</label></span></h2>
    <div id="metrics"></div></section>
  <section class="panel" id="patPanel"><h2>Top read patterns<span class="ctl">
      <label><select id="patMode"><option value="share">share of reads with a hit</option><option value="reads">reads</option></select></label></span></h2>
    <div class="scroll"><table id="pats"></table></div></section>
  <section class="panel hidden" id="binPanel"><h2>ONT bins (twopass)</h2>
    <div class="scroll"><table id="bins"></table></div></section>
  <section class="panel hidden" id="smpPanel"><h2>Reads per sample<span class="ctl">
      <label class="multi">x <select id="sx"></select></label><label class="multi">y <select id="sy"></select></label>
      <label><select id="scale"><option value="log">log scale</option><option value="linear">linear scale</option></select></label>
      <label class="multi"><input type="checkbox" id="hideZero"> hide samples with 0 in both</label></span></h2>
    <div class="two"><div><div id="scatter"></div><div class="stat" id="sstat"></div></div>
      <div><div style="display:flex;gap:8px;margin-bottom:6px"><input id="q" placeholder="filter samples…" style="flex:1">
        <select id="sortBy"><option value="diff">sort: largest difference</option><option value="fc">sort: largest fold change</option><option value="x">sort: x reads</option><option value="y">sort: y reads</option><option value="name">sort: name</option></select></div>
        <div class="scroll" style="max-height:420px;overflow-y:auto"><table id="stab"></table></div></div></div></section>
  <section class="panel hidden" id="negPanel"><h2>Negative controls<span class="ctl">
      <label class="multi">colour wells by <select id="ncol"></select></label></span></h2>
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
const SINGLE = R.length === 1;              // one run: a run summary, nothing to compare against
if (SINGLE) document.body.classList.add("single");
const subLine = () => SINGLE
  ? `${R[0].subcommand || "run"}${/\s--dorado\b/.test(R[0].command) ? " --dorado" : ""}${/\s--pass1-only\b/.test(R[0].command) ? " --pass1-only" : ""}${R[0].context ? " -c " + R[0].context : ""}${R[0].kit ? " -k " + R[0].kit : ""} · started ${R[0].started} · generated ${D.generated}`
  : `${R.length} runs · generated ${D.generated}`;
$("sub").textContent = subLine();
// columns of the run views: the runs, or for one twopass run "all bins" + each ONT bin (pass 2)
// (bins with < 0.1% of the reads, at least 10, are counted in "all bins" only)
const P2 = SINGLE && R[0].pass2 ? R[0].pass2 : [], minBin = Math.max(10, 0.001 * ((R[0] && R[0].summary.total_reads) || 0));
const P2small = P2.filter(b => (b.summary.total_reads || 0) < minBin);
const C = P2.length ? [{...R[0], label: "all bins"}, ...P2.filter(b => !P2small.includes(b))] : R;
if (P2small.length && R[0].samples)          // their samples (all wells of the plate map, ~0 reads) would swamp the sample views
  for (const k of Object.keys(R[0].samples)) if (P2small.some(b => k.startsWith(b.label + "/"))) delete R[0].samples[k];

const fk = Object.keys(D.foundLabels);

// ---- twopass pass 1: Dorado / MinKNOW bins or barbell kit
const p1kind = r => !r.pass1 ? "" : /dorado/i.test(r.pass1.method) ? "Dorado / MinKNOW bins" : "barbell kit";
if (R.some(r => r.pass1)) {
  $("p1Panel").classList.remove("hidden");
  $("p1").innerHTML = R.filter(r => r.pass1).map(r => `<div style="margin-bottom:10px">` +
    `<div class="p1m">${SINGLE ? "" : `<span style="color:var(--dim);font-weight:400">${esc(r.label)}</span> · `}` +
    `<span class="tag">${esc(p1kind(r))}</span>${esc(r.pass1.method)}</div>` +
    (/dorado/i.test(r.pass1.method)
      ? `<div class="p1l">The reads were already split by ONT barcode when they were basecalled (Dorado / MinKNOW barcodeNN folders or read headers); utag used those bins as they are.</div>`
      : /\s--pass1-only\b/.test(r.command)
      ? `<div class="p1l">utag split the reads by ONT barcode with barbell kit (and trimmed the barcode) into ont_bins/; --pass1-only: no utag demux.</div>`
      : `<div class="p1l">utag split the reads by ONT barcode with barbell kit (and trimmed the barcode) before searching the utags in each bin.</div>`) +
    r.pass1.lines.map(l => `<div class="p1l ${l.startsWith("! ") ? "pw" : ""}">${esc(l.replace(/^[•!] /, ""))}</div>`).join("") + `</div>`).join("");
  const kinds = [...new Set(R.map(p1kind).filter(Boolean))];
  $("binPanel").querySelector("h2").firstChild.textContent = `ONT bins (${kinds.join(", ")})`;
}

// ---- metrics table
// background by share (0..1), the blue steps of the top read patterns table (--q0..--q5); metrics % and ONT bins tables
const shade = f => { const i = f <= 0 ? 0 : Math.min(5, 1 + Math.floor(f * 5)); return `background:var(--q${i});${i >= 4 ? "color:#fff" : ""}`; };
const isPctRow = k => /^(reads|pattern_reads|unexpected_reads)/.test(k);
// the run summary as utag prints it: metric, value, % of total reads (for the read counts)
function drawMetrics() {
  const diffOnly = $("diffOnly").checked, showPct = $("showPct").checked, keys = [];
  C.forEach(r => Object.keys(r.summary).forEach(k => keys.includes(k) || keys.push(k)));
  const ofTotal = k => k === "total_reads" || isPctRow(k);
  const many = C.length > 1, span = showPct ? 2 : 1;
  let h = many ? `<tr><th></th>${C.map(r => `<th class="run" colspan="${span}">${esc(r.label)}</th>`).join("")}</tr>` : "";
  h += `<tr class="rule"><th>metric</th>${C.map(() => `<th class="${showPct ? "" : "gap"}">value</th>` + (showPct ? `<th class="gap">of total</th>` : "")).join("")}</tr>`;
  for (const k of keys) {
    const vals = C.map(r => r.summary[k]), defined = vals.filter(v => v !== undefined && v !== null);
    if (diffOnly && defined.every(v => v === defined[0])) continue;
    h += `<tr><td>${esc(k)}</td>` + C.map(r => { const v = r.summary[k], t = r.summary.total_reads;
      return `<td class="${showPct ? "" : "gap"}">${v === undefined ? "–" : fmt(v)}</td>` +
        (showPct ? `<td class="of gap">${v !== undefined && ofTotal(k) && t ? `<span class="cell" style="${shade(v / t)}">${pf(100 * v / t)}</span>` : ""}</td>` : ""); }).join("") + "</tr>";
  }
  let out = `<div class="scroll"><table class="term">${h}</table></div>`;
  if (P2small.length) out += `<div class="note">${P2small.length} small ONT bin(s) not shown as columns (< ${fmt(Math.round(minBin))} reads; counted in all bins): ` +
    `${P2small.map(b => esc(b.label) + " " + fmt(b.summary.total_reads || 0)).join(", ")}</div>`;
  $("metrics").innerHTML = out;
  $("sub").textContent = subLine();
}
$("diffOnly").onchange = drawMetrics; $("showPct").onchange = drawMetrics;
drawMetrics();

// ---- patterns
const qv = ["--q0", "--q1", "--q2", "--q3", "--q4", "--q5"];
function heat(f) { const i = f <= 0 ? 0 : Math.min(5, 1 + Math.floor(f * 5 / 0.4)); return `background:var(${qv[i]});${i >= 4 ? "color:#fff" : ""}`; }
// a pattern as a read: grey line, a block at each end where something was found (as in utag context's report),
// and in the middle for @prev_left hits (a second barcode inside the read: chimeric / concatenated reads);
// F = utag-fwd (blue), R = utag-rev (orange); filled = utag, outlined = primer-flank only; the point = strand (fw right, rc left)
function patRead(p) {
  const W = 130, H = 16, bw = 24, y = 2, h = 12;
  let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(p)}"><rect x="0" y="${H / 2 - 1.5}" width="${W}" height="3" rx="1.5" fill="var(--axis)"/>`;
  const parts = p.split("__").map(x => /^(F|R)(tag|flank)\[(fw|rc)\b.*@(left|right|prev_left)\b/.exec(x)).filter(Boolean);
  const mids = parts.filter(m => m[4] === "prev_left").length;
  let j = 0;
  for (const m of parts) {
    const [, side, kind, strand, end] = m, col = side === "F" ? "var(--s1)" : "var(--s2)";
    const x0 = end === "left" ? 0 : end === "right" ? W - bw : (W - bw) / 2 + (j++ - (mids - 1) / 2) * (bw + 2), x1 = x0 + bw, t = 6;   // arrow: body + point towards the strand
    const pts = strand === "fw" ? [[x0, y], [x1 - t, y], [x1, y + h / 2], [x1 - t, y + h], [x0, y + h]]
                                : [[x0 + t, y], [x1, y], [x1, y + h], [x0 + t, y + h], [x0, y + h / 2]];
    svg += `<polygon points="${pts.map(q => q.join(",")).join(" ")}" fill="${kind === "tag" ? col : "var(--panel)"}" stroke="${col}" stroke-width="1.5"/>` +
      `<text x="${(x0 + x1) / 2 + (strand === "fw" ? -2 : 2)}" y="${y + 9.5}" text-anchor="middle" style="font-size:9px;font-weight:700;fill:${kind === "tag" ? "#fff" : col}">${side}</text>`;
  }
  return svg + "</svg>";
}
function drawPats() {
  if (!D.patterns.length) { $("patPanel").classList.add("hidden"); return; }
  const mode = $("patMode").value;
  const tot = C.map(r => Object.values(r.found).reduce((a, b) => a + b, 0) || r.pattern_reads || Object.values(r.patterns).reduce((a, b) => a + b, 0) || 1);
  let h = `<tr><th>read</th><th class="pat">pattern</th>${C.map(r => `<th>${esc(r.label)}</th>`).join("")}</tr>`;
  for (const p of D.patterns) { h += `<tr><td style="line-height:0">${patRead(p)}</td><td class="pat">${esc(p)}</td>` + C.map((r, i) => { const n = r.patterns[p];
      if (n === undefined) return `<td><span class="neu">–</span></td>`; const f = n / tot[i];
      return `<td><span class="cell" style="${heat(f)}" title="${fmt(n)} reads">${mode === "share" ? pf(100 * f) : fmt(n)}</span></td>`; }).join("") + "</tr>"; }
  const kit = R.every(r => /\s--pass1-only\b/.test(r.command));   // barbell kit patterns: F / R = the ONT barcode
  const what = kit ? ["ONT barcode with its front flank", "ONT barcode with its rear flank", "barcode", "flank"] : ["utag-fwd", "utag-rev", "utag", "primer-flank"];
  $("pats").innerHTML = h + `<tr><td colspan="${C.length + 2}" class="neu" style="text-align:left;white-space:normal;font-size:11px">` +
    `${kit ? "barbell kit patterns (pass 1). " : ""}Only the top patterns are in a log; – = not in that run's top list. ` +
    `Read: <b style="color:var(--s1)">F</b> ${what[0]}, <b style="color:var(--s2)">R</b> ${what[1]} at that read end; filled = ${what[2]}, ` +
    `outlined = ${what[3]} only; the point shows the strand (fw → , rc ←).</td></tr>`;
}
$("patMode").onchange = drawPats; drawPats();

// ---- ONT bins
// bins matched across runs by their number (NB01 = barcode01, as in the sample sheet); row label = all names used
const binKey = b => { const m = /(\d+)$/.exec(b); return m ? "#" + +m[1] : b; };
const binRows = new Map();
R.forEach(r => Object.keys(r.bins).forEach(b => { const k = binKey(b); if (!binRows.has(k)) binRows.set(k, new Set()); binRows.get(k).add(b); }));
const bins = [...binRows.keys()].sort((a, b) => a.localeCompare(b, undefined, {numeric: true}));
const binOfRun = (r, k) => { const n = Object.keys(r.bins).find(b => binKey(b) === k); return n ? r.bins[n] : undefined; };
function drawBins() {                     // one row per bin; reads, hits, assigned, assigned %, samples (per run)
  const br = R.filter(r => Object.keys(r.bins).length), multi = br.length > 1;
  const p1only = br.every(r => Object.values(r.bins).every(x => x.assigned === undefined));   // twopass --pass1-only
  const cols = p1only ? ["reads in bin", "of total reads"] : ["reads in bin", "with hit", "assigned", "assigned %", "samples with reads"];
  let h = multi ? `<tr><th></th>${br.map(r => `<th colspan="${cols.length}" style="text-align:center;border-bottom:1px solid var(--axis)">${esc(r.label)}</th>`).join("")}</tr>` : "";
  h += `<tr><th>ONT bin</th>${br.map(() => cols.map(c => `<th>${c}</th>`).join("")).join("")}</tr>`;
  for (const b of bins) { h += `<tr><td>${esc([...binRows.get(b)].join(" · "))}</td>` + br.map(r => { const x = binOfRun(r, b);
      if (!x) return cols.map(() => `<td class="neu">–</td>`).join("");
      if (p1only) { const t = r.summary.total_reads, f = t ? x.reads / t : null;
        return `<td>${fmt(x.reads)}</td><td>${f === null ? "–" : `<span class="cell" style="${shade(f)}">${pf(100 * f)}</span>`}</td>`; }
      const f = x.reads && x.assigned !== undefined ? x.assigned / x.reads : null, sf = x.samples_expected ? x.samples_with_reads / x.samples_expected : null;
      return `<td>${fmt(x.reads)}</td><td>${fmt(x.with_hit)}</td><td>${fmt(x.assigned)}</td>` +
        `<td>${f === null ? "–" : `<span class="cell" style="${shade(f)}">${pf(100 * f)}</span>`}</td>` +
        `<td>${x.samples_with_reads === undefined ? "–" : `<span class="cell" style="${sf === null ? "" : shade(sf)}">${fmt(x.samples_with_reads)} / ${fmt(x.samples_expected)}</span>`}</td>`; }).join("") + "</tr>"; }
  $("bins").innerHTML = h;
}
if (bins.length) { $("binPanel").classList.remove("hidden"); drawBins(); }

// ---- samples
const SR = R.filter(r => r.samples);
if (SR.length === 1) {                      // one run: samples ranked by reads, and how even they are
  $("smpPanel").classList.remove("hidden");
  const X = SR[0], ALL = Object.entries(X.samples).map(([k, v]) => ({k, v}));
  let all = ALL, tot = 0;
  const setBin = b => { all = b ? ALL.filter(p => p.k.startsWith(b + "/")) : ALL; tot = all.reduce((a, p) => a + p.v, 0); };
  const sbins = X.subcommand === "twopass" ? [...new Set(ALL.map(p => p.k.split("/")[0]).filter((b, i, a) => ALL[i].k.includes("/")))].sort() : [];
  if (sbins.length > 1) {                   // twopass: one ONT bin or all of them
    $("smpPanel").querySelector("h2 .ctl").insertAdjacentHTML("afterbegin",
      `<label>ONT bin <select id="sbin"><option value="">all bins</option>${sbins.map(b => `<option>${esc(b)}</option>`).join("")}</select></label>`);
    $("sbin").onchange = () => { setBin($("sbin").value); drawRank(); };
  }
  setBin("");
  // sort shared by the table (click a header) and the plot: reads, sample name, or plate and well (row-major)
  const srt = {key: "reads", dir: -1}, WRE = /^([A-Pa-p])0?(\d{1,2})$/;
  const wellKey = p => { const mm = D.meta[p.k] || [], w = WRE.exec(mm[1] || "");
    return [mm[0] || "~", w ? w[1].toUpperCase().charCodeAt(0) : 999, w ? +w[2] : 999]; };
  const cmp = (a, b) => { let d;
    if (srt.key === "reads") d = a.v - b.v;
    else if (srt.key === "name") d = a.k.localeCompare(b.k, undefined, {numeric: true});
    else { const x = wellKey(a), y = wellKey(b); d = x[0].localeCompare(y[0], undefined, {numeric: true}) || x[1] - y[1] || x[2] - y[2]; }
    return srt.dir * d || a.k.localeCompare(b.k, undefined, {numeric: true}); };
  const shownRows = () => { const fq = $("q").value.toLowerCase(); return all.filter(p => p.k.toLowerCase().includes(fq)).sort(cmp); };
  const q = (a, f) => { const s = [...a].sort((x, y) => x - y), i = (s.length - 1) * f, lo = Math.floor(i); return s.length ? s[lo] + (s[Math.min(s.length - 1, lo + 1)] - s[lo]) * (i - lo) : 0; };
  let med = 0;
  function drawRank() {                     // statistics over the ONT bin; plot and table follow the filter and sort
    const v = all.map(p => p.v).sort((a, b) => b - a), n = v.length;
    med = q(v, .5);
    const mean = tot / Math.max(1, n), sd = Math.sqrt(v.reduce((a, x) => a + (x - mean) ** 2, 0) / Math.max(1, n)), with0 = v.filter(x => !x).length;
    const within2 = med ? v.filter(x => x >= med / 2 && x <= med * 2).length : 0, p10 = q(v, .1), p90 = q(v, .9);
    $("sstat").innerHTML = `${fmt(n)} samples · ${fmt(n - with0)} with reads${with0 ? ` · <span style="color:var(--bad)">${fmt(with0)} with 0</span>` : ""} · ` +
      `median ${fmt(Math.round(med))} · mean ${fmt(Math.round(mean))} · min ${fmt(v[n - 1] ?? 0)} · max ${fmt(v[0] ?? 0)} · CV ${mean ? (sd / mean).toFixed(2) : "–"} · ` +
      `90th / 10th percentile ${p10 ? (p90 / p10).toFixed(1) + "×" : "–"} · within 2× of the median ${pf(pc(within2, n))}`;
    drawSorted();
  }
  function drawSorted() { const rows = shownRows(); drawPlot(rows); drawRankTable(rows); }
  function drawPlot(rows) {
    const log = $("scale").value === "log", n = rows.length;
    const W = 560, H = 360, m = {l: 58, r: 14, t: 12, b: 42}, pw = W - m.l - m.r, ph = H - m.t - m.b, mx = Math.max(1, ...all.map(p => p.v));
    const tf = x => log ? Math.log10(x + 1) : x, tmax = tf(mx) * 1.04, sx = i => m.l + pw * (i + .5) / Math.max(1, n), sy = x => m.t + ph - ph * tf(x) / tmax;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="reads per sample">`;
    const ticks = log ? [0, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7].filter(t => t <= mx * 1.05) : [0, .25, .5, .75, 1].map(f => Math.round(f * mx));
    for (const t of ticks) svg += `<line x1="${m.l}" x2="${W - m.r}" y1="${sy(t)}" y2="${sy(t)}" stroke="var(--grid)"/><text x="${m.l - 6}" y="${sy(t) + 4}" text-anchor="end">${fmt(t)}</text>`;
    for (const [f, lab] of [[1, "median"], [2, "2× median"], [.5, "½× median"]]) { if (!med) break; const y = sy(med * f);
      if (y < m.t || y > m.t + ph) continue;
      svg += `<line x1="${m.l}" x2="${W - m.r}" y1="${y}" y2="${y}" stroke="var(--axis)" stroke-dasharray="${lab === "median" ? "" : "4 4"}"/><text x="${W - m.r - 4}" y="${y - 4}" text-anchor="end">${lab}</text>`; }
    const order = srt.key === "reads" ? (srt.dir < 0 ? "most reads first" : "fewest reads first") : srt.key === "name" ? "by name" : "by plate and well";
    svg += `<text x="${m.l + pw / 2}" y="${H - 6}" text-anchor="middle" style="fill:var(--ink)">samples, ${order} (${fmt(n)})</text>`;
    svg += `<text transform="translate(14 ${m.t + ph / 2}) rotate(-90)" text-anchor="middle" style="fill:var(--ink)">reads</text>`;
    const bw = Math.max(1, pw / Math.max(1, n) - (n > 150 ? 0 : 1));
    rows.forEach((p, i) => { const y = sy(p.v);
      svg += `<rect x="${(sx(i) - bw / 2).toFixed(1)}" y="${y.toFixed(1)}" width="${bw.toFixed(2)}" height="${Math.max(0, m.t + ph - y).toFixed(1)}" fill="${p.v ? "var(--s1)" : "var(--bad)"}" data-k="${esc(p.k)}"/>`;
      if (!p.v) svg += `<rect x="${(sx(i) - bw / 2).toFixed(1)}" y="${m.t + ph - 3}" width="${bw.toFixed(2)}" height="3" fill="var(--bad)"/>`; });
    $("scatter").innerHTML = svg + "</svg>";
    $("scatter").querySelectorAll("rect[data-k]").forEach(el => { const p = all.find(x => x.k === el.dataset.k), mm = D.meta[p.k] || [];
      el.addEventListener("mousemove", e => showTip(e, `<b>${esc(p.k)}</b><br>${fmt(p.v)} reads · ${pf(pc(p.v, tot))} of assigned` +
        (med ? ` · ${(p.v / med).toFixed(2)}× median` : "") + (mm[1] ? `<br><span class="t">${esc(mm[0])} ${esc(mm[1])} · utags ${esc(mm[2])} / ${esc(mm[3])}</span>` : "")));
      el.addEventListener("mouseleave", hideTip); });
  }
  function drawRankTable(rows) {
    const shown = rows.slice(0, 500), arrow = k => srt.key === k ? (srt.dir < 0 ? " ▼" : " ▲") : "";
    const th = (k, lab) => `<th class="sortable" data-s="${k}">${lab}${arrow(k)}</th>`;
    $("stab").innerHTML = `<tr>${th("name", "sample")}${th("well", "plate · well")}${th("reads", "reads")}${th("reads", "% of assigned")}${th("reads", "× median")}</tr>` +
      shown.map(p => { const mm = D.meta[p.k] || [], r = med ? p.v / med : null;
        return `<tr><td>${esc(p.k)}</td><td>${esc([mm[0], mm[1]].filter(Boolean).join(" · "))}</td><td>${fmt(p.v)}</td><td>${pf(pc(p.v, tot))}</td>` +
          `<td class="${r === null ? "neu" : r < .5 || r > 2 ? "down" : ""}">${r === null ? "–" : r.toFixed(2)}</td></tr>`; }).join("") +
      (rows.length > 500 ? `<tr><td colspan="5" class="neu">${fmt(rows.length - 500)} more – use the filter</td></tr>` : "");
  }
  $("stab").addEventListener("click", e => { const h = e.target.closest("th[data-s]"); if (!h) return; const k = h.dataset.s;
    srt.dir = srt.key === k ? -srt.dir : k === "reads" ? -1 : 1; srt.key = k; drawSorted(); });
  $("sortBy").classList.add("multi");
  $("scale").onchange = drawSorted;
  $("q").oninput = drawSorted;
  drawRank();
} else if (SR.length) {
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
  const M = D.meta, KEYS = [...new Set(SR.flatMap(r => Object.keys(r.samples)))];
  let keys = KEYS;                          // the wells analysed: all, or one ONT bin of a twopass run
  const SKEY = "compare-runs:neg:" + document.title + ":" + KEYS.length;
  let st = {}; try { st = JSON.parse(localStorage.getItem(SKEY) || "{}") || {}; } catch (e) {}
  const save = () => { try { localStorage.setItem(SKEY, JSON.stringify(st)); } catch (e) {} };
  let mode = "neg";
  const stateOf = k => st[k] || "pos";
  // plates from the summary.tsv plate/well columns; samples without a well go into an "other" panel
  const WRE = /^([A-Pa-p])0?(\d{1,2})$/, panels = new Map();
  function buildPanels() { panels.clear();
    for (const k of keys) { const m = M[k] || [], pw = WRE.exec(m[1] || ""), pn = pw && m[0] ? m[0] : "other";
      if (!panels.has(pn)) panels.set(pn, []);
      panels.get(pn).push({k, r: pw ? pw[1].toUpperCase().charCodeAt(0) - 65 : -1, c: pw ? +pw[2] - 1 : -1}); } }
  // twopass: each ONT bin is its own set of wells, so negatives are analysed per bin (marks are kept for all bins)
  const nbins = SINGLE && SR[0].subcommand === "twopass" ? [...new Set(KEYS.filter(k => k.includes("/")).map(k => k.split("/")[0]))].sort() : [];
  if (nbins.length > 1) {
    $("negPanel").querySelector("h2 .ctl").insertAdjacentHTML("afterbegin",
      `<label>ONT bin <select id="nbin">${nbins.map(b => `<option>${esc(b)}</option>`).join("")}<option value="">all bins</option></select></label>`);
    keys = KEYS.filter(k => k.startsWith(nbins[0] + "/"));
    $("nbin").onchange = () => { const b = $("nbin").value; keys = b ? KEYS.filter(k => k.startsWith(b + "/")) : KEYS;
      buildPanels(); drawPlates(); analyse(); };
  }
  buildPanels();
  const RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"].map(h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16)));
  function ramp(t) { t = Math.max(0, Math.min(1, t)) * (RAMP.length - 1); const i = Math.min(RAMP.length - 2, Math.floor(t)), f = t - i;
    return "rgb(" + RAMP[i].map((v, j) => Math.round(v + (RAMP[i + 1][j] - v) * f)).join(",") + ")"; }
  $("ncol").innerHTML = SR.map((r, i) => `<option value="${i}">${esc(r.label)}</option>`).join("");
  $("ncol").value = Math.max(0, SR.findIndex(r => r.label === D.base));

  function drawPlates() {
    let h = "", grp = null;                 // twopass plates (<bin>/<plate>): one row of plates per ONT bin
    const binOf = pn => pn.includes("/") ? pn.slice(0, pn.indexOf("/")) : "";
    [...panels.entries()].sort((a, b) => a[0] === "other" ? 1 : b[0] === "other" ? -1 : a[0].localeCompare(b[0], undefined, {numeric: true})).forEach(([pn, ws], pi) => {
      const bin = pn === "other" ? "" : binOf(pn);
      if (bin !== grp) { h += (grp !== null ? "</div></div>" : "") + `<div class="bgrp">${bin ? `<div class="bt">${esc(bin)}</div>` : ""}<div class="plates">`; grp = bin; }
      const nr = pn === "other" ? 0 : Math.max(8, ...ws.map(w => w.r + 1)), nc = pn === "other" ? 0 : Math.max(12, ...ws.map(w => w.c + 1));
      const cs = nc > 12 ? 11 : 16, nN = ws.filter(w => stateOf(w.k) === "neg").length, nE = ws.filter(w => stateOf(w.k) === "exc").length;
      h += `<div class="pl" data-p="${pi}"><div class="pt" data-act="plate" title="${esc(pn)}">${esc(bin ? pn.slice(bin.length + 1) : pn)}<small>${nN ? nN + " neg" : ""}${nN && nE ? " · " : ""}${nE ? nE + " excl" : ""}</small></div>`;
      if (pn === "other") { h += `<div class="pg" style="--cs:${cs}px;grid-template-columns:repeat(24,${cs}px)">` + ws.map(w => `<div class="w" data-k="${esc(w.k)}"></div>`).join("") + "</div></div>"; return; }
      const at = new Map(ws.map(w => [w.r + "_" + w.c, w]));
      h += `<div class="pg" style="--cs:${cs}px;grid-template-columns:12px repeat(${nc},${cs}px)"><div></div>`;
      for (let c = 0; c < nc; c++) h += `<div class="h" data-act="col" data-c="${c}">${c + 1}</div>`;
      for (let r = 0; r < nr; r++) { h += `<div class="h" data-act="row" data-r="${r}">${String.fromCharCode(65 + r)}</div>`;
        for (let c = 0; c < nc; c++) { const w = at.get(r + "_" + c); h += w ? `<div class="w" data-k="${esc(w.k)}" data-r="${r}" data-c="${c}"></div>` : `<div class="w none"></div>`; } }
      h += "</div></div>";
    });
    $("nplates").innerHTML = h + (grp !== null ? "</div></div>" : ""); paint();
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
    $("nres").innerHTML = h;
    strip(S, cut); tags();
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
    $("ntags").innerHTML = rows.length ? `<div class="note" style="margin:0 0 4px">reads in negative wells per utag (top ${rows.length})</div><div class="scroll" style="max-height:300px;overflow-y:auto"><table class="narrow"><tr><th>utag</th>${SR.map(r => `<th>${esc(r.label)}</th>`).join("")}</tr>` +
      rows.map(([id, v]) => `<tr><td>${esc(id)}</td>${v.map(n => `<td>${fmt(n)}</td>`).join("")}</tr>`).join("") + "</table></div>" : "";
  }
  $("ndl").onclick = () => {
    const q = v => /[",\n]/.test(v) ? '"' + String(v).replace(/"/g, '""') + '"' : v;
    const lines = lastRows.map(r => r.map(q).join(","));
    lines.push("", "sample,state");
    for (const k of keys) if (stateOf(k) !== "pos") lines.push(q(k) + "," + (stateOf(k) === "neg" ? "negative" : "excluded"));
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], {type: "text/csv"}));
    a.download = `negative_controls${$("nbin") && $("nbin").value ? "_" + $("nbin").value : ""}.csv`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  drawPlates(); analyse();
})();

// ---- commands
$("cmds").innerHTML = R.map(r => `<details ${r.status !== "ok" ? "open" : ""} style="margin-bottom:8px"><summary><b style="color:var(--ink)">${esc(r.label)}</b><span class="badge ${r.status}">${r.status}</span> · ${esc(r.log)}</summary>
  <div style="padding:6px 0 4px 14px"><div class="mono">$ ${esc(r.command)}</div><div class="k" style="color:var(--muted);font-size:11px">cwd ${esc(r.cwd)} · started ${esc(r.started)}</div>
  ${Object.keys(r.steps).length ? `<div class="stat">${Object.entries(r.steps).map(([s, t]) => `${esc(s)} <b>${t.toFixed(1)} s</b>`).join(" · ")}</div>` : ""}
  ${r.info.map(s => `<div class="stat">• ${esc(s)}</div>`).join("")}
  ${r.warnings.map(s => `<div class="warn">${esc(s)}</div>`).join("")}
  ${r.errors.map(s => `<div class="err mono">${esc(s)}</div>`).join("")}
  ${(r.pass2 || []).map(b => `<div class="stat"><b style="color:var(--ink)">pass 2 · ${esc(b.label)}</b> <span class="badge ${b.status}">${b.status}</span>` +
      ` · barbell ${b.barbell_seconds.toFixed(1)} s</div>` + b.warnings.map(s => `<div class="warn">${esc(s)}</div>`).join("") +
      b.errors.map(s => `<div class="err mono">${esc(s)}</div>`).join("")).join("")}</div></details>`).join("");
</script></body></html>
"""

if __name__ == "__main__":
    main()

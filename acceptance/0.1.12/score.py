"""Scores the prose-auto 0.1.12 held-out benchmark once (PREREG.md §5, §6, §9). Frozen with its SHA-256 in FREEZE.md
before any held-out raw file is downloaded or opened.

    python score.py DIR       reads DIR/{reference,results,decisions,checks}.tsv, writes DIR/score.json, DIR/score.md

Inputs (tab-separated, one header line; an empty cell is missing):
  reference.tsv  run, class (bottom-up | top-down | no-proteome), organism_verified (yes | no), label, protease, cys,
                 prec_window_ppm, frag_window, frag_unit (the adjudicated reference annotation, §3)
  results.tsv    run, arm (auto | reference | sdrf | defaults | best_psm | ablation names), status (ok | stopped |
                 missing), psm, pep, pg, psm_fdp, pep_fdp (counts at 1 % FDR; FDP in %, combined entrapment estimate)
  decisions.tsv  run, arm, setting (label | protease | cys | prec_window_ppm | frag_window), value, status (decided |
                 default | unknown | review)
  checks.tsv     run, kind (correct | related_close | related_distant | unrelated | incomplete | wrong_label |
                 wrong_protease | natural), check (pre | post), verdict (ok | warning | database | settings |
                 not assessed)
Only numpy and the standard library; the bootstrap uses numpy.random.default_rng(20261120).
"""
import csv
import json
import math
import os
import sys

import numpy as np

SEED, B = 20261120, 10000
MIN_PSM = 100
P1_MARGIN, P3_MARGIN, FDP_MAX = 0.98, -0.05, 1.5
# §9: per kind and check, the one-sided 95 % Wilson lower bound the sensitivity (specificity for "correct") must reach
CHECK_BOUND = {("pre", "correct"): 0.75, ("pre", "unrelated"): 0.65, ("pre", "related_distant"): 0.50,
               ("post", "correct"): 0.75, ("post", "unrelated"): 0.85, ("post", "wrong_label"): 0.70,
               ("post", "wrong_protease"): 0.70}   # power.py: the largest bound with power >= 0.80 at the development
#                                                     share; wrong settings capped at 0.70 (development 24/25, 25/25)
REPORTED_ONLY = ("related_close", "incomplete", "natural")


def read(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [{k: (v if v != "" else None) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


# ---- statistics -----------------------------------------------------------------------------------------------------

def wilson(k, n, z=1.6449):
    """one-sided 95 % Wilson bounds (lower, upper) for k of n"""
    if n == 0:
        return (None, None)
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def _binom_cdf(k, n, p):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))


def clopper_lower(k, n, alpha=0.05):
    """one-sided (1 - alpha) Clopper-Pearson lower bound for k of n"""
    if k == 0:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):          # P(X >= k | p) = alpha
        mid = (lo + hi) / 2
        if 1 - _binom_cdf(k - 1, n, mid) < alpha:
            lo = mid
        else:
            hi = mid
    return lo


def boot(stat, data, q, rng=None):
    """percentile bootstrap over runs: the q-quantile of stat over B resamples of the rows of data"""
    rng = rng or np.random.default_rng(SEED)
    data = np.asarray(data)
    idx = rng.integers(0, len(data), size=(B, len(data)))
    return float(np.quantile(stat(data[idx]), q))


def median_rows(x):
    return np.median(x, axis=-1)


# ---- the run set and the ratios -------------------------------------------------------------------------------------

def eligible(ref):
    """runs in the primary analysis: bottom-up with a verified organism (§3, §5); the others are counted"""
    ok = {r["run"] for r in ref if r.get("class") == "bottom-up" and r.get("organism_verified") == "yes"}
    out = {"eligible": len(ok), "top-down": sum(r.get("class") == "top-down" for r in ref),
           "no-proteome": sum(r.get("class") == "no-proteome" for r in ref),
           "organism not verified": sum(r.get("class") == "bottom-up" and r.get("organism_verified") != "yes"
                                        for r in ref)}
    return ok, out


def _num(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def arm_table(results):
    return {(r["run"], r["arm"]): r for r in results}


def ratios(tab, runs, num, den, field="psm"):
    """paired per-run ratios num/den (§5): a stopped arm has 0; both below MIN_PSM: left out; den below MIN_PSM and
    num above it: +inf (num better); a missing arm: left out. Returns (ratios, counts)."""
    out, cnt = [], {"pairs": 0, "missing": 0, "both below": 0, "denominator below": 0}
    for run in sorted(runs):
        a, b = tab.get((run, num)), tab.get((run, den))
        if not a or not b or "missing" in (a["status"], b["status"]):
            cnt["missing"] += 1
            continue
        x = 0.0 if a["status"] == "stopped" else (_num(a, field) or 0.0)
        y = 0.0 if b["status"] == "stopped" else (_num(b, field) or 0.0)
        if x < MIN_PSM and y < MIN_PSM:
            cnt["both below"] += 1
            continue
        cnt["pairs"] += 1
        if y < MIN_PSM:
            cnt["denominator below"] += 1
            out.append(math.inf)
        else:
            out.append(x / y)
    return out, cnt


def ratio_endpoint(tab, runs, num, den, field="psm"):
    r, cnt = ratios(tab, runs, num, den, field)
    if not r:
        return {**cnt, "median": None, "ci95": [None, None]}
    lr = np.log(np.clip(np.array(r), 1e-9, 1e9))   # +inf as 1e9: the median is rank-based
    lo = boot(median_rows, lr, 0.025)
    hi = boot(median_rows, lr, 0.975)
    m = float(np.exp(np.median(lr)))
    return {**cnt, "median": math.inf if m >= 1e8 else m, "ci95": [float(np.exp(lo)), float(np.exp(hi))],
            "gain": int(sum(x > 1 for x in r)), "loss": int(sum(x < 1 for x in r)), "equal": int(sum(x == 1 for x in r))}


def compliant(row):
    f, p = _num(row, "psm_fdp"), _num(row, "pep_fdp")
    return f is not None and p is not None and f <= FDP_MAX and p <= FDP_MAX


def p3(tab, runs):
    """P3: paired difference of the share of runs at FDP <= 1.5 % (PSM and peptide level), auto minus reference, over
    the runs with >= MIN_PSM in both; one-sided 95 % lower bound (5th bootstrap percentile) >= P3_MARGIN. The absolute
    auto share with its one-sided Clopper-Pearson bound is reported."""
    pairs = []
    for run in sorted(runs):
        a, b = tab.get((run, "auto")), tab.get((run, "reference"))
        if a and b and a["status"] == b["status"] == "ok" and (_num(a, "psm") or 0) >= MIN_PSM \
                and (_num(b, "psm") or 0) >= MIN_PSM:
            pairs.append((compliant(a), compliant(b)))
    if not pairs:
        return {"runs": 0}
    d = np.array([[x, y] for x, y in pairs], dtype=float)
    diff = lambda s: s[..., 0].mean(axis=-1) - s[..., 1].mean(axis=-1)
    k = int(d[:, 0].sum())
    return {"runs": len(pairs), "auto_share": k / len(pairs), "reference_share": float(d[:, 1].mean()),
            "difference": float(d[:, 0].mean() - d[:, 1].mean()), "difference_lower95": boot(diff, d, 0.05),
            "auto_share_lower95": clopper_lower(k, len(pairs))}


# ---- secondary endpoints --------------------------------------------------------------------------------------------

def accuracy(ref, dec, runs, arm="auto"):
    """S1-S5: per setting, exact agreement with the reference among the runs where both are stated, Wilson bounds;
    S4 windows within 1.5x of the measured window; S5 abstentions"""
    truth = {r["run"]: r for r in ref}
    out = {}
    for s in ("label", "protease", "cys", "prec_window_ppm", "frag_window"):
        k = n = abst = 0
        for d in dec:
            if d["arm"] != arm or d["setting"] != s or d["run"] not in runs:
                continue
            t = truth[d["run"]].get(s)
            if d.get("status") in ("default", "unknown", "review"):
                abst += 1
            if t is None or d.get("value") is None:
                continue
            n += 1
            if s in ("prec_window_ppm", "frag_window"):
                v, w = float(d["value"]), float(t)
                k += w / 1.5 <= v <= w * 1.5
            else:
                k += d["value"].strip().lower() == t.strip().lower()
        lo, hi = wilson(k, n)
        out[s] = {"correct": k, "n": n, "share": k / n if n else None, "wilson95": [lo, hi], "abstained": abst}
    return out


def checks(rows):
    """S6 and the §9 check rules: per check and kind, flagged (sensitivity) or ok (specificity, kind correct) among the
    assessed scenarios, one-sided 95 % Wilson bounds, against CHECK_BOUND"""
    out = {}
    for (chk, kind) in sorted({(r["check"], r["kind"]) for r in rows}):
        rs = [r for r in rows if r["check"] == chk and r["kind"] == kind]
        ass = [r for r in rs if r["verdict"] != "not assessed"]
        hit = sum((r["verdict"] == "ok") if kind == "correct" else (r["verdict"] != "ok") for r in ass)
        exact = None
        if kind in ("wrong_label", "wrong_protease"):
            exact = sum(r["verdict"] == "settings" for r in ass)
        elif kind in ("related_close", "related_distant", "unrelated", "incomplete") and chk == "post":
            exact = sum(r["verdict"] == "database" for r in ass)
        lo, hi = wilson(hit, len(ass))
        bound = CHECK_BOUND.get((chk, kind))
        out[f"{chk}:{kind}"] = {"scenarios": len(rs), "assessed": len(ass), "hit": hit,
                                "share": hit / len(ass) if ass else None, "wilson95": [lo, hi],
                                "exact_verdict": exact, "bound": bound,
                                "pass": None if bound is None or lo is None else lo >= bound}
    return out


def cys_rule(tab, dec, never="auto_never", always="auto_always"):
    """§9: FASTag -alkylation_out instead of the competition stays the default if, on the runs whose Cys decision
    differs between the never and always arms, the summed PSM ratio never/always >= 0.99 and no run is below 0.95"""
    v = {(d["run"], d["arm"]): d["value"] for d in dec if d["setting"] == "cys"}
    runs = sorted({r for r, a in v if a == never and v.get((r, always)) not in (None, v[(r, never)])})
    pairs = [(_num(tab[(r, never)], "psm") or 0, _num(tab[(r, always)], "psm") or 0) for r in runs
             if (r, never) in tab and (r, always) in tab]
    if not pairs:
        return {"runs": 0, "pass": None}
    s = sum(a for a, _ in pairs) / max(1, sum(b for _, b in pairs))
    worst = min((a / b for a, b in pairs if b), default=1.0)
    return {"runs": len(pairs), "summed_ratio": s, "worst": worst, "pass": s >= 0.99 and worst >= 0.95}


def score(d):
    ref, res, dec, chk = (read(os.path.join(d, f"{n}.tsv")) for n in ("reference", "results", "decisions", "checks"))
    runs, counts = eligible(ref)
    tab = arm_table(res)
    out = {"runs": counts, "P1": ratio_endpoint(tab, runs, "auto", "reference"),
           "P2": ratio_endpoint(tab, runs, "auto", "defaults"), "P3": p3(tab, runs)}
    out["P1"]["pass"] = out["P1"]["ci95"][0] is not None and out["P1"]["ci95"][0] >= P1_MARGIN
    out["P2"]["pass"] = out["P2"]["ci95"][0] is not None and out["P2"]["ci95"][0] > 1
    out["P3"]["pass"] = out["P3"].get("difference_lower95") is not None and out["P3"]["difference_lower95"] >= P3_MARGIN
    out["pass"] = out["P1"]["pass"] and out["P2"]["pass"] and out["P3"]["pass"]
    out["reported"] = {f"{f} {n}/{m}": ratio_endpoint(tab, runs, n, m, f)
                       for f in ("psm", "pep", "pg") for n, m in (("auto", "reference"), ("auto", "defaults"),
                                                                ("auto", "sdrf"), ("auto", "best_psm"))
                       if f != "psm" or m in ("sdrf", "best_psm")}
    out["S1-S5"] = accuracy(ref, dec, runs)
    out["S6"] = checks(chk)
    out["S8 top-down"] = {"top-down runs": counts["top-down"],
                          "stopped": sum(1 for r in ref if r.get("class") == "top-down"
                                         and (tab.get((r["run"], "auto")) or {}).get("status") == "stopped"),
                          "bottom-up stopped": sum(1 for r in runs if (tab.get((r, "auto")) or {}).get("status")
                                                   == "stopped")}
    out["rule cys"] = cys_rule(tab, dec)
    return out


def markdown(s):
    f = lambda x: "n/a" if x is None else f"{x:.3f}" if isinstance(x, float) else str(x)
    lines = [f"# Held-out score (PREREG 0.1.12)", "", f"Runs: {s['runs']}", "", "| Endpoint | n | estimate | bound | pass |",
             "|---|---|---|---|---|"]
    for k in ("P1", "P2"):
        e = s[k]
        lines.append(f"| {k} median PSM ratio | {e.get('pairs')} | {f(e['median'])} | {f(e['ci95'][0])} | {e['pass']} |")
    e = s["P3"]
    lines.append(f"| P3 FDP share auto - reference | {e.get('runs')} | {f(e.get('difference'))} | "
                 f"{f(e.get('difference_lower95'))} | {e['pass']} |")
    lines += ["", f"**Study passes: {s['pass']}**", "", "| Check:kind | assessed | share | Wilson lower | bound | pass |",
              "|---|---|---|---|---|---|"]
    for k, e in s["S6"].items():
        lines.append(f"| {k} | {e['assessed']} | {f(e['share'])} | {f(e['wilson95'][0])} | {f(e['bound'])} | {e['pass']} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    d = sys.argv[1]
    s = score(d)
    json.dump(s, open(os.path.join(d, "score.json"), "w"), indent=1, default=str)
    open(os.path.join(d, "score.md"), "w").write(markdown(s))
    print(markdown(s))

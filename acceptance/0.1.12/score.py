"""Scores the prose-auto 0.1.12 held-out benchmark once (PREREG.md §5, §6, §9). Frozen with its SHA-256 in FREEZE.md
before any held-out raw file is downloaded or opened.

    python score.py DIR [HELDOUT]   reads DIR/{reference,results,decisions,checks,scenarios}.tsv (and timing.tsv if
                                    present) and HELDOUT (default: heldout.tsv next to this file); writes
                                    DIR/score.json and DIR/score.md

Every input is validated before anything is computed: required columns, enumerations, numbers, duplicate keys, runs
outside the draw. A blank cell is missing and is never a result. Any violation stops the scoring with the reason.

Inputs (tab-separated, one header line; written by extract.py from the arm outputs, except reference.tsv):
  reference.tsv  run, class (bottom-up | top-down | no-proteome), organism_verified (yes | no), label, protease, cys,
                 reporter_ions (yes | no), prec_window, prec_unit, frag_window, frag_unit. run = accession/file of the
                 draw (heldout.tsv); a drawn run without a row is missing (download or conversion failed).
  results.tsv    run, arm, status (ok | stopped | missing), psm, pep, pg, psm_fdp, pep_fdp. ok rows need all five
                 numbers; stopped and missing rows none. FDP in %, the combined entrapment estimate.
  decisions.tsv  run, arm, setting (label | protease | cys | prec_window | frag_window | ms2_analyzer_qc |
                 ms2_analyzer_header), value, unit (ppm | Da, windows only), status (auto | review | default | given |
                 unknown)
  checks.tsv     run, kind, check (pre | post), verdict (ok | warning | database | settings | not assessed)
  scenarios.tsv  the frozen manifest after adjudication (scenarios.py): which (run, kind) scenarios were built
  timing.tsv     run, version, input_type, step, seconds (S7, descriptive)
Only numpy and the standard library; every bootstrap uses numpy.random.default_rng(20261120) and B resamples.
"""
import csv
import json
import math
import os
import sys

import numpy as np

SEED, B, ALPHA = 20261120, 10000, 0.05
MIN_PSM = 100            # an arm counts as identifying from 100 PSMs at 1 % FDR
MIN_ELIGIBLE = 120       # fewer eligible runs: P1-P3 get no verdict (§2.2)
P1_MARGIN, P3_MARGIN, FDP_MAX = 0.98, -0.05, 1.5
P3_ABS = 0.75            # P3 absolute gate: one-sided 95 % Clopper-Pearson lower bound of auto's share (power.py)
MZ_REF = 1000.0          # windows in Da and ppm are compared at m/z 1,000
# §9 check gates: (check, kind) -> the bound H0 "share <= bound" is tested against (exact one-sided binomial, Holm
# over all gates at ALPHA). Each bound is the largest (step 0.05, from 0.10) with power >= 0.80 at ALPHA/9 in power.py,
# given at least 10 expected assessed scenarios;
# wrong_protease is capped at wrong_label's bound (development 25/25 sits at the ceiling). A kind without such a bound
# is experimental: reported with its interval, never a gate.
CHECK_BOUND = {("pre", "correct"): 0.70, ("pre", "unrelated"): 0.50, ("pre", "related_close"): 0.45,
               ("pre", "related_distant"): 0.45, ("post", "correct"): 0.75, ("post", "unrelated"): 0.75,
               ("post", "related_distant"): 0.20, ("post", "wrong_label"): 0.65, ("post", "wrong_protease"): 0.65}
EXPERIMENTAL = {("post", "related_close"), ("pre", "incomplete"), ("post", "incomplete")}

CLASSES = {"bottom-up", "top-down", "no-proteome"}
STATUS = {"ok", "stopped", "missing"}
DEC_STATUS = {"auto", "review", "default", "given", "unknown"}
ABSTAIN = {"default", "unknown", "review"}
SETTINGS = {"label", "protease", "cys", "prec_window", "frag_window", "ms2_analyzer_qc", "ms2_analyzer_header"}
UNITS = {"ppm", "Da"}
CHECKS = {"pre", "post"}
KINDS = {"correct", "related_close", "related_distant", "unrelated", "incomplete", "wrong_label", "wrong_protease",
         "natural"}
VERDICTS = {"ok", "warning", "database", "settings", "not assessed"}
ISOBARIC = ("tmt", "itraq")


class InputError(ValueError):
    pass


# ---- reading and validation -----------------------------------------------------------------------------------------

def read(path, cols, required=True):
    """rows of a TSV as dicts, values stripped, blank -> None; the columns must be present"""
    if not os.path.exists(path):
        if required:
            raise InputError(f"{path}: missing")
        return []
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        miss = [c for c in cols if c not in (r.fieldnames or [])]
        if miss:
            raise InputError(f"{path}: columns {miss} missing")
        return [{k: (v.strip() if isinstance(v, str) and v.strip() != "" else None) for k, v in row.items()}
                for row in r]


def _enum(rows, col, allowed, where, nullable=False):
    for i, r in enumerate(rows, 2):
        v = r.get(col)
        if (v is None and not nullable) or (v is not None and v not in allowed):
            raise InputError(f"{where} line {i}: {col} = {v!r}, expected one of {sorted(allowed)}")


def _unique(rows, key, where):
    seen = set()
    for i, r in enumerate(rows, 2):
        k = tuple(r.get(c) for c in key)
        if k in seen:
            raise InputError(f"{where} line {i}: duplicate {dict(zip(key, k))}")
        seen.add(k)


def _number(v, where):
    try:
        x = float(v)
    except (TypeError, ValueError):
        raise InputError(f"{where}: {v!r} is not a number")
    if not math.isfinite(x) or x < 0:
        raise InputError(f"{where}: {v!r} is not a finite non-negative number")
    return x


def load(d, heldout):
    drawn = [f"{r['accession']}/{r['file']}" for r in read(heldout, ["accession", "file"])]
    if len(set(drawn)) != len(drawn):
        raise InputError(f"{heldout}: duplicate runs")
    known = set(drawn)
    ref = read(f"{d}/reference.tsv", ["run", "class", "organism_verified", "label", "protease", "cys", "reporter_ions",
                                      "prec_window", "prec_unit", "frag_window", "frag_unit"])
    res = read(f"{d}/results.tsv", ["run", "arm", "status", "psm", "pep", "pg", "psm_fdp", "pep_fdp"])
    dec = read(f"{d}/decisions.tsv", ["run", "arm", "setting", "value", "unit", "status"])
    chk = read(f"{d}/checks.tsv", ["run", "kind", "check", "verdict"])
    scn = read(f"{d}/scenarios.tsv", ["accession", "file", "kind", "databases"])
    tim = read(f"{d}/timing.tsv", ["run", "version", "input_type", "step", "seconds"], required=False)
    _enum(ref, "class", CLASSES, "reference.tsv")
    _enum(ref, "organism_verified", {"yes", "no"}, "reference.tsv")
    _enum(ref, "reporter_ions", {"yes", "no"}, "reference.tsv", nullable=True)
    _enum(ref, "prec_unit", UNITS, "reference.tsv", nullable=True)
    _enum(ref, "frag_unit", UNITS, "reference.tsv", nullable=True)
    _unique(ref, ["run"], "reference.tsv")
    for r in ref:
        if r["run"] not in known:
            raise InputError(f"reference.tsv: run {r['run']} is not in the draw")
        for w, u in (("prec_window", "prec_unit"), ("frag_window", "frag_unit")):
            if (r[w] is None) != (r[u] is None):
                raise InputError(f"reference.tsv {r['run']}: {w} and {u} go together")
            if r[w] is not None:
                r[w] = _number(r[w], f"reference.tsv {r['run']} {w}")
    _enum(res, "status", STATUS, "results.tsv")
    _unique(res, ["run", "arm"], "results.tsv")
    for r in res:
        where = f"results.tsv {r['run']} {r['arm']}"
        if r["arm"] is None:
            raise InputError(f"{where}: arm missing")
        nums = ["psm", "pep", "pg", "psm_fdp", "pep_fdp"]
        if r["status"] == "ok":
            for c in nums:
                r[c] = _number(r[c], f"{where} {c}")
        elif any(r[c] is not None for c in nums):
            raise InputError(f"{where}: a {r['status']} arm has no counts")
    _enum(dec, "setting", SETTINGS, "decisions.tsv")
    _enum(dec, "status", DEC_STATUS, "decisions.tsv")
    _enum(dec, "unit", UNITS, "decisions.tsv", nullable=True)
    _unique(dec, ["run", "arm", "setting"], "decisions.tsv")
    for r in dec:
        if r["setting"] in ("prec_window", "frag_window") and r["value"] is not None:
            if r["unit"] is None:
                raise InputError(f"decisions.tsv {r['run']} {r['setting']}: unit missing")
            r["value"] = _number(r["value"], f"decisions.tsv {r['run']} {r['setting']}")
    _enum(chk, "check", CHECKS, "checks.tsv")
    _enum(chk, "kind", KINDS, "checks.tsv")
    _enum(chk, "verdict", VERDICTS, "checks.tsv")
    _unique(chk, ["run", "kind", "check"], "checks.tsv")
    for rows, where in ((res, "results.tsv"), (dec, "decisions.tsv"), (chk, "checks.tsv"), (tim, "timing.tsv")):
        for r in rows:
            if r["run"] not in known and not (where == "checks.tsv" and r["kind"] != "correct"):
                raise InputError(f"{where}: run {r['run']} is not in the draw")
    return drawn, ref, res, dec, chk, scn, tim


# ---- statistics (one implementation for score.py and power.py) ------------------------------------------------------

def median(x, axis=-1):
    """the sample median of extended reals in [0, inf]: the middle order statistic, or for an even count the mean of
    the two middle ones (inf if either is inf; no -inf exists, so no inf - inf)"""
    s = np.sort(np.asarray(x, dtype=float), axis=axis)
    n = s.shape[axis]
    lo, hi = np.take(s, (n - 1) // 2, axis=axis), np.take(s, n // 2, axis=axis)
    return (lo + hi) / 2


def quantile(x, q):
    """the q-quantile of bootstrap statistics: the order statistic at ceil(q * B) (numpy method inverted_cdf; no
    interpolation, so inf and 0 stay exact)"""
    return float(np.quantile(np.asarray(x, dtype=float), q, method="inverted_cdf"))


def boot(stat, data, rng=None, b=None):
    """B bootstrap replicates of stat over the rows (runs) of data"""
    rng = rng if rng is not None else np.random.default_rng(SEED)
    data = np.asarray(data)
    idx = rng.integers(0, len(data), size=(b or B, len(data)))
    return stat(data[idx])


def wilson(k, n, z=1.6449):
    """one-sided 95 % Wilson bounds (lower, upper) for k of n"""
    if n == 0:
        return (None, None)
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)"""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def clopper_lower(k, n, alpha=ALPHA):
    """one-sided (1 - alpha) Clopper-Pearson lower bound for k of n"""
    if n == 0:
        return None
    if k == 0:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):          # binom_sf(k, n, p) = alpha
        mid = (lo + hi) / 2
        if binom_sf(k, n, mid) < alpha:
            lo = mid
        else:
            hi = mid
    return lo


def holm(pvalues):
    """Holm's step-down at ALPHA: {key: rejected}"""
    order = sorted(pvalues.items(), key=lambda kv: kv[1])
    out, m = {}, len(order)
    for i, (k, p) in enumerate(order):
        if p > ALPHA / (m - i):
            out.update({k2: False for k2, _ in order[i:]})
            break
        out[k] = True
    return out


# ---- runs, ratios, P1-P3 --------------------------------------------------------------------------------------------

def accounting(drawn, ref):
    """every drawn run in one class (§2.2); eligible = bottom-up with a verified organism"""
    by = {r["run"]: r for r in ref}
    cls = {}
    for run in drawn:
        r = by.get(run)
        cls[run] = "missing" if r is None else "organism not verified" if (
            r["class"] == "bottom-up" and r["organism_verified"] != "yes") else r["class"]
    counts = {k: sum(v == k for v in cls.values()) for k in
              ("bottom-up", "top-down", "no-proteome", "organism not verified", "missing")}
    return {r for r, c in cls.items() if c == "bottom-up"}, {"drawn": len(drawn), **counts}, cls


def _count(row, field):
    return 0.0 if row["status"] == "stopped" else row[field]


def ratios(tab, runs, num, den, field="psm"):
    """paired per-run ratios num/den of field (§5). Eligibility is by PSMs: an arm missing -> left out; both arms below
    MIN_PSM PSMs -> left out; den below MIN_PSM PSMs and num at or above -> inf. A stopped arm has 0. For peptides and
    protein groups of an eligible pair: den 0 and num 0 -> left out, den 0 -> inf. Returns (ratios, counts)."""
    out, cnt = [], {"pairs": 0, "missing": 0, "both below": 0, "denominator below": 0, "both zero": 0}
    for run in sorted(runs):
        a, b = tab.get((run, num)), tab.get((run, den))
        if not a or not b or "missing" in (a["status"], b["status"]):
            cnt["missing"] += 1
            continue
        pa, pb = _count(a, "psm"), _count(b, "psm")
        if pa < MIN_PSM and pb < MIN_PSM:
            cnt["both below"] += 1
            continue
        if pb < MIN_PSM:
            cnt["pairs"] += 1
            cnt["denominator below"] += 1
            out.append(math.inf)
            continue
        x, y = _count(a, field), _count(b, field)
        if y == 0:
            if x == 0:
                cnt["both zero"] += 1
                continue
            cnt["pairs"] += 1
            out.append(math.inf)
            continue
        cnt["pairs"] += 1
        out.append(x / y)
    return out, cnt


def ratio_endpoint(tab, runs, num, den, field="psm", rng=None):
    r, cnt = ratios(tab, runs, num, den, field)
    if not r:
        return {**cnt, "median": None, "ci95": [None, None]}
    meds = boot(median, np.array(r), rng)
    return {**cnt, "median": float(median(r)), "ci95": [quantile(meds, 0.025), quantile(meds, 0.975)],
            "gain": int(sum(x > 1 for x in r)), "loss": int(sum(x < 1 for x in r)), "equal": int(sum(x == 1 for x in r))}


def compliant(row):
    """FDP <= FDP_MAX on PSM and peptide level, with >= MIN_PSM PSMs; a stopped arm is not compliant"""
    return row["status"] == "ok" and row["psm"] >= MIN_PSM and row["psm_fdp"] <= FDP_MAX and row["pep_fdp"] <= FDP_MAX


def p3_pairs(tab, runs):
    """(auto compliant, reference compliant) over the eligible runs whose reference arm identifies (>= MIN_PSM); auto
    stopped or below MIN_PSM counts as noncompliant; a missing arm is left out and counted"""
    pairs, missing = [], 0
    for run in sorted(runs):
        a, b = tab.get((run, "auto")), tab.get((run, "reference"))
        if not a or not b or "missing" in (a["status"], b["status"]):
            missing += 1
            continue
        if b["status"] == "ok" and b["psm"] >= MIN_PSM:
            pairs.append((compliant(a), compliant(b)))
    return pairs, missing


def p3(tab, runs, rng=None):
    """P3 (§5), two gates, both required: relative, the 5th bootstrap percentile of auto share minus reference share
    >= P3_MARGIN; absolute, the one-sided 95 % Clopper-Pearson lower bound of auto's share >= P3_ABS"""
    pairs, missing = p3_pairs(tab, runs)
    if not pairs:
        return {"runs": 0, "missing": missing}
    d = np.array(pairs, dtype=float)
    diffs = boot(lambda s: s[..., 0].mean(axis=-1) - s[..., 1].mean(axis=-1), d, rng)
    k = int(d[:, 0].sum())
    return {"runs": len(pairs), "missing": missing, "auto_share": k / len(pairs),
            "reference_share": float(d[:, 1].mean()), "difference": float(d[:, 0].mean() - d[:, 1].mean()),
            "difference_lower95": quantile(diffs, 0.05), "auto_share_lower95": clopper_lower(k, len(pairs)),
            "reference_share_lower95": clopper_lower(int(d[:, 1].sum()), len(pairs))}


# ---- secondary endpoints and §9 rules --------------------------------------------------------------------------------

def _ppm(v, unit):
    return v if unit == "ppm" else v / MZ_REF * 1e6


def _norm(v):
    return str(v).strip().lower()


def accuracy(ref, dec, runs, arm="auto"):
    """S1-S5 over the eligible runs with a reference value: exact agreement (windows: within 1.5x after conversion to
    ppm at m/z MZ_REF); a run without a decision row counts as wrong ("absent"); abstentions are counted"""
    truth = {r["run"]: r for r in ref}
    d = {(x["run"], x["setting"]): x for x in dec if x["arm"] == arm}
    out = {}
    for s, unit in (("label", None), ("protease", None), ("cys", None), ("prec_window", "prec_unit"),
                    ("frag_window", "frag_unit")):
        k = n = absent = abst = 0
        for run in sorted(runs):
            t = truth[run].get(s)
            if t is None:
                continue
            n += 1
            x = d.get((run, s))
            if x is None or x["value"] is None:
                absent += 1
                continue
            abst += x["status"] in ABSTAIN
            if unit:
                w, v = _ppm(t, truth[run][unit]), _ppm(x["value"], x["unit"])
                k += w / 1.5 <= v <= w * 1.5
            else:
                k += _norm(x["value"]) == _norm(t)
        out[s] = {"correct": k, "n": n, "share": k / n if n else None, "wilson95": list(wilson(k, n)),
                  "absent": absent, "abstained": abst}
    return out


def built_scenarios(scn):
    return {(f"{r['accession']}/{r['file']}", r["kind"]) for r in scn
            if r["kind"] not in ("excluded", "unresolved") and r["databases"] not in (None, "none")}


def checks(chk, scn, runs):
    """S6 and the §9 check gates. Specificity (kind correct): the auto arm's own checks on the eligible runs.
    Sensitivity: the built scenarios of the frozen manifest (wrong_label and wrong_protease: the built correct
    scenarios). A scenario without a check row is missing; not assessed is left out; both are counted. Each gate:
    exact one-sided binomial test of H0 share <= bound, Holm over all gates at ALPHA."""
    built = built_scenarios(scn)
    rows = {(r["run"], r["kind"], r["check"]): r["verdict"] for r in chk}
    out, pv = {}, {}
    for chk_ in sorted(CHECKS):
        for kind in sorted(KINDS):
            if kind == "natural":
                want = sorted({r for r, k, c in rows if k == "natural" and c == chk_})
            elif kind == "correct":
                want = sorted(runs)
            elif kind in ("wrong_label", "wrong_protease"):   # the single-organism runs with a built correct scenario
                multi = {r for r, k in built if k == "incomplete"}
                want = sorted(r for r, k in built if k == "correct" and r not in multi) if chk_ == "post" else []
            else:
                want = sorted(r for r, k in built if k == kind)
            if not want:
                continue
            v = [rows.get((r, kind, chk_)) for r in want]
            ass = [x for x in v if x not in (None, "not assessed")]
            hit = sum((x == "ok") if kind == "correct" else (x != "ok") for x in ass)
            exact = sum(x == "settings" for x in ass) if kind in ("wrong_label", "wrong_protease") else \
                sum(x == "database" for x in ass) if chk_ == "post" and kind != "correct" else None
            e = {"scenarios": len(want), "missing": v.count(None), "not_assessed": v.count("not assessed"),
                 "assessed": len(ass), "hit": hit, "share": hit / len(ass) if ass else None,
                 "wilson95": list(wilson(hit, len(ass))), "exact_verdict": exact,
                 "bound": CHECK_BOUND.get((chk_, kind)),
                 "status": "gate" if (chk_, kind) in CHECK_BOUND else
                 "experimental" if (chk_, kind) in EXPERIMENTAL else "reported"}
            if e["bound"] is not None:
                e["p_value"] = binom_sf(hit, len(ass), e["bound"]) if ass else 1.0
                pv[f"{chk_}:{kind}"] = e["p_value"]
            out[f"{chk_}:{kind}"] = e
    for k, rej in holm(pv).items():
        out[k]["pass"] = rej
    return out


def _cys_rule(tab, dec, a1, a2, runs):
    """on the eligible runs whose Cys decision differs between arms a1 and a2: summed PSM ratio a1/a2 >= 0.99 and no
    run below 0.95. A run with either arm missing is left out and counted; a stopped arm has 0 PSMs."""
    v = {(d["run"], d["arm"]): d["value"] for d in dec if d["setting"] == "cys"}
    diff = sorted(r for r in runs if (r, a1) in v and (r, a2) in v and v[(r, a1)] != v[(r, a2)])
    pairs, missing = [], 0
    for r in diff:
        x, y = tab.get((r, a1)), tab.get((r, a2))
        if not x or not y or "missing" in (x["status"], y["status"]):
            missing += 1
            continue
        pairs.append((_count(x, "psm"), _count(y, "psm")))
    if not pairs:
        return {"runs": 0, "missing": missing, "pass": None}
    s = sum(a for a, _ in pairs) / sum(b for _, b in pairs) if sum(b for _, b in pairs) else math.inf
    worst = min((a / b for a, b in pairs if b), default=math.inf)
    return {"runs": len(pairs), "missing": missing, "summed_ratio": s, "worst": worst,
            "pass": s >= 0.99 and worst >= 0.95}


def rules(tab, ref, dec, runs, acc):
    """§9 detector rules other than the check gates"""
    truth = {r["run"]: r for r in ref}
    out = {"cys: -alkylation_out (auto_never) vs competition (auto_always)":
           _cys_rule(tab, dec, "auto_never", "auto_always", runs),
           "mini-search fallback (auto) vs always (auto_always)": _cys_rule(tab, dec, "auto", "auto_always", runs)}
    lab = [r for r in runs if truth[r].get("label") is not None]
    base = sum(_norm(truth[r]["label"]) in ("none", "label free") for r in lab)
    d = {(x["run"], x["setting"], x["arm"]): x for x in dec}
    bad = [r for r in runs if truth[r].get("reporter_ions") == "no" and (r, "label", "auto") in d
           and d[(r, "label", "auto")]["value"] and _norm(d[(r, "label", "auto")]["value"]).startswith(ISOBARIC)]
    out["label: S1 >= always-label-free baseline, no isobaric call without reporter ions"] = {
        "auto": acc["label"]["correct"], "baseline": base, "n": len(lab), "isobaric_without_reporters": len(bad),
        "pass": (acc["label"]["correct"] >= base and not bad) if lab else None}
    prot = [r for r in runs if truth[r].get("protease") is not None]
    stated = sum((r, "protease", "sdrf") in d and d[(r, "protease", "sdrf")]["value"] is not None
                 and _norm(d[(r, "protease", "sdrf")]["value"]) == _norm(truth[r]["protease"]) for r in prot)
    out["protease: S2 >= stated-protease (SDRF arm) baseline"] = {
        "auto": acc["protease"]["correct"], "baseline": stated, "n": len(prot),
        "pass": acc["protease"]["correct"] >= stated if prot else None}
    qc = [(d[(r, "ms2_analyzer_qc", "auto")]["value"], d[(r, "ms2_analyzer_header", "auto")]["value"]) for r in runs
          if (r, "ms2_analyzer_qc", "auto") in d and (r, "ms2_analyzer_header", "auto") in d
          and d[(r, "ms2_analyzer_qc", "auto")]["value"] and d[(r, "ms2_analyzer_header", "auto")]["value"]]
    agree = sum(_norm(a) == _norm(b) for a, b in qc)
    out["FASTag QC: analyser class agrees with the header in >= 98 %"] = {
        "agree": agree, "n": len(qc), "pass": agree / len(qc) >= 0.98 if qc else None}
    return out


def timing(tim):
    """S7 (descriptive): median seconds per input type, version and step"""
    g = {}
    for r in tim:
        g.setdefault((r["input_type"], r["version"], r["step"]), []).append(_number(r["seconds"], "timing.tsv"))
    return {" | ".join(k): {"runs": len(v), "median_s": float(np.median(v))} for k, v in sorted(g.items())} \
        if g else "timing.tsv not provided"


def score(d, heldout=None):
    heldout = heldout or os.path.join(os.path.dirname(os.path.abspath(__file__)), "heldout.tsv")
    drawn, ref, res, dec, chk, scn, tim = load(d, heldout)
    runs, counts, cls = accounting(drawn, ref)
    tab = {(r["run"], r["arm"]): r for r in res}
    out = {"runs": counts, "class": cls, "P1": ratio_endpoint(tab, runs, "auto", "reference"),
           "P2": ratio_endpoint(tab, runs, "auto", "defaults"), "P3": p3(tab, runs)}
    enough = len(runs) >= MIN_ELIGIBLE
    out["P1"]["pass"] = (out["P1"]["ci95"][0] is not None and out["P1"]["ci95"][0] >= P1_MARGIN) if enough else None
    out["P2"]["pass"] = (out["P2"]["ci95"][0] is not None and out["P2"]["ci95"][0] > 1) if enough else None
    e = out["P3"]
    e["relative_pass"] = e.get("difference_lower95") is not None and e["difference_lower95"] >= P3_MARGIN
    e["absolute_pass"] = e.get("auto_share_lower95") is not None and e["auto_share_lower95"] >= P3_ABS
    e["pass"] = (e["relative_pass"] and e["absolute_pass"]) if enough else None
    out["pass"] = (out["P1"]["pass"] and out["P2"]["pass"] and out["P3"]["pass"]) if enough else None
    out["verdict"] = ("pass" if out["pass"] else "fail") if enough else \
        f"no verdict: {len(runs)} eligible runs, {MIN_ELIGIBLE} needed"
    arms = sorted({a for _, a in tab} - {"auto"})
    out["reported"] = {f"{f} auto/{m}": ratio_endpoint(tab, runs, "auto", m, f)
                       for m in arms for f in ("psm", "pep", "pg") if not (f == "psm" and m in ("reference", "defaults"))}
    out["S1-S5"] = acc = accuracy(ref, dec, runs)
    out["S6"] = checks(chk, scn, runs)
    out["S7"] = timing(tim)
    out["S8 top-down"] = {"top-down runs": counts["top-down"],
                          "stopped": sum(1 for r, c in cls.items() if c == "top-down"
                                         and (tab.get((r, "auto")) or {}).get("status") == "stopped"),
                          "bottom-up stopped": sum(1 for r in runs
                                                   if (tab.get((r, "auto")) or {}).get("status") == "stopped")}
    out["rules"] = rules(tab, ref, dec, runs, acc)
    return out


def markdown(s):
    f = lambda x: "n/a" if x is None else ("inf" if x == math.inf else f"{x:.3f}") if isinstance(x, float) else str(x)
    lines = ["# Held-out score (PREREG 0.1.12)", "", f"Runs: {s['runs']}", "", f"**Verdict: {s['verdict']}**", "",
             "| Endpoint | n | estimate | bound | pass |", "|---|---|---|---|---|"]
    for k in ("P1", "P2"):
        e = s[k]
        lines.append(f"| {k} median PSM ratio | {e.get('pairs')} | {f(e['median'])} | {f(e['ci95'][0])} | {e['pass']} |")
    e = s["P3"]
    lines.append(f"| P3 relative (auto - reference) | {e.get('runs')} | {f(e.get('difference'))} | "
                 f"{f(e.get('difference_lower95'))} | {e.get('relative_pass')} |")
    lines.append(f"| P3 absolute (auto share) | {e.get('runs')} | {f(e.get('auto_share'))} | "
                 f"{f(e.get('auto_share_lower95'))} | {e.get('absolute_pass')} |")
    lines += ["", "| Check:kind | status | assessed | share | Wilson lower | bound | Holm pass |",
              "|---|---|---|---|---|---|---|"]
    for k, e in s["S6"].items():
        lines.append(f"| {k} | {e['status']} | {e['assessed']} | {f(e['share'])} | {f(e['wilson95'][0])} | "
                     f"{f(e['bound'])} | {e.get('pass')} |")
    lines += ["", "| §9 rule | pass |", "|---|---|"] + [f"| {k} | {v.get('pass')} |" for k, v in s["rules"].items()]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    d = sys.argv[1]
    s = score(d, sys.argv[2] if len(sys.argv) > 2 else None)
    json.dump(s, open(os.path.join(d, "score.json"), "w"), indent=1, default=str)
    open(os.path.join(d, "score.md"), "w").write(markdown(s))
    print(markdown(s))

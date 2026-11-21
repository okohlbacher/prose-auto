"""Scores the prose-auto 0.1.12 held-out benchmark once (PREREG.md §5, §6, §9). Frozen with its SHA-256 in FREEZE.md
before any held-out raw file is downloaded or opened.

    python score.py DIR [HELDOUT] [--pilot]   reads DIR/{reference,results,decisions,checks,scenarios}.tsv (and timing.tsv if
                                    present) and HELDOUT (default: heldout.tsv next to this file); writes
                                    DIR/score.json and DIR/score.md

Every input is validated before anything is computed: required columns, enumerations, numbers, duplicate keys, runs
outside the draw. A blank cell is missing and is never a result. Any violation stops the scoring with the reason.

Inputs (tab-separated, one header line; written by extract.py from the arm outputs, except reference.tsv):
  reference.tsv  run, class (bottom-up | top-down | no-proteome), organism_verified (yes | no), label, protease, cys,
                 reporter_ions (yes | no), prec_window, prec_unit, frag_window, frag_unit. run = accession/file of the
                 draw (heldout.tsv); a drawn run without a row is missing (download or conversion failed).
  results.tsv    run, arm, status (ok | stopped | missing), psm, pep, pg, psm_fdp, pep_fdp. ok rows need all five
                 numbers; stopped and missing rows none. FDP in %, the combined entrapment estimate. A bottom-up run
                 whose auto or reference arm has no row or is missing did not execute: it is missing, before the
                 120-run minimum (§2.2).
  decisions.tsv  run, arm, setting (label | protease | cys | prec_window | frag_window | ms2_analyzer_qc |
                 ms2_analyzer_header), value, unit (ppm | Da, windows only), status (auto | review | default | given |
                 unknown)
  checks.tsv     run, source (auto | scenario | natural), arm, kind, check (pre | post), verdict (ok | warning |
                 database | settings | not assessed), seconds. source auto: the auto arm's own checks on a held-out run
                 (kind correct; specificity); scenario: a built scenario of the manifest (any kind, its run need not be
                 in the draw); natural: the check on the SDRF or defaults arm's result (arm names it); only the arms of
                 the natural cohort count (score.natural_cohort), others are counted and ignored. seconds: the check's
                 run time, may be blank.
  scenario_counts.tsv  run, kind, status (ok | missing), psm: every built scenario's search, whatever its checks gave
                 (S6 losses); written by extract.py, optional (losses are then not reported)
  scenarios.tsv  the frozen manifest after adjudication (scenarios.py): which (run, kind) scenarios were built
  timing.tsv     run, version, input_type, step, seconds (S7, descriptive)
Only numpy and the standard library; every bootstrap uses numpy.random.default_rng(20261120) and B resamples, the seed
restarted for each endpoint. score.json writes +inf as the string "inf" (strict JSON).
"""
import csv
import json
import math
from fractions import Fraction
import os
import sys

import numpy as np

SEED, B, ALPHA = 20261120, 10000, 0.05
MIN_PSM = 100            # an arm counts as identifying from 100 PSMs at 1 % FDR
MIN_INFORMATIVE = 5      # a Cys deployment rule needs this many informative runs, else no verdict
MIN_ELIGIBLE = 120       # fewer eligible runs: P1-P3 get no verdict (§2.2)
P1_MARGIN, P3_MARGIN, FDP_MAX = 0.98, -0.05, 1.5
P3_ABS = 0.75            # P3 absolute gate: one-sided 95 % Clopper-Pearson lower bound of auto's share (power.py)
MZ_REF = 1000.0          # windows in Da and ppm are compared at m/z 1,000
# §9 check gates: (check, kind) -> the bound H0 "share <= bound" is tested against (exact one-sided binomial, Holm
# over all gates at ALPHA). Each bound is the largest (step 0.05, from 0.10) with power >= 0.80 at ALPHA/9 in power.py,
# given at least 10 expected assessed scenarios (specificity at the 5th percentile of the eligible runs under high
# attrition; related_distant at the manifest's rank mix);
# wrong_protease is capped at wrong_label's bound (development 25/25 sits at the ceiling). A kind without such a bound
# is experimental: reported with its interval, never a gate.
CHECK_BOUND = {("pre", "correct"): 0.65, ("pre", "unrelated"): 0.50, ("pre", "related_close"): 0.45,
               ("pre", "related_distant"): 0.50, ("post", "correct"): 0.65, ("post", "unrelated"): 0.75,
               ("post", "related_distant"): 0.25, ("post", "wrong_label"): 0.65, ("post", "wrong_protease"): 0.65}
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
SOURCES = {"auto", "scenario", "natural"}
GATES = sorted(CHECK_BOUND)
# S1: prose-auto's label vocabulary (prose_auto/schema.py), the granularity S1 scores; the reference uses the same
# names. Every "no label" spelling is one category; TMT10/11 and TMTpro 16/18 are one category each, as prose-auto
# reports them.
LABEL_CANON = {"label free": "none", "none": "none", "no label": "none", "no isobaric label": "none",
               "no silac": "none", "tmt2plex": "TMT2plex", "tmt6plex": "TMT6plex", "tmt10/11plex": "TMT10/11plex",
               "tmt10plex": "TMT10/11plex", "tmt11plex": "TMT10/11plex", "tmtpro (16/18plex)": "TMTpro (16/18plex)",
               "tmtpro16plex": "TMTpro (16/18plex)", "tmtpro18plex": "TMTpro (16/18plex)",
               "tmt0 / tmt (126 only)": "TMT0", "itraq4plex": "iTRAQ4plex", "itraq8plex": "iTRAQ8plex",
               "silac (lys8/arg10)": "SILAC (Lys8/Arg10)", "silac (arg6)": "SILAC (Arg6)",
               "silac lys4 or dimethyl (d4)": "SILAC Lys4 or dimethyl (d4)", "silac (lys4/arg6)": "SILAC (Lys4/Arg6)",
               "silac triplex (lys4/arg6, lys8/arg10)": "SILAC triplex (Lys4/Arg6, Lys8/Arg10)",
               "15n metabolic labelling": "15N", "heavy water (2h2o) labelling": "2H2O"}
LABEL_CANON.update({v.lower(): v for v in set(LABEL_CANON.values())})   # the canonical names map to themselves
ISOBARIC = ("tmt", "itraq")
QC_NO_CALL = "no call"   # extract.py: FASTag QC ran but its 90 % rule made no MS2 analyser call


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


def validate_checks(chk, where):
    """the check rows' enumerations and keys (shared with power.py): a blank verdict is never a result"""
    _enum(chk, "source", SOURCES, where)
    _enum(chk, "check", CHECKS, where)
    _enum(chk, "kind", KINDS, where)
    _enum(chk, "verdict", VERDICTS, where)
    _unique(chk, ["run", "source", "arm", "kind", "check"], where)


def load(d, heldout):
    drawn = [f"{r['accession']}/{r['file']}" for r in read(heldout, ["accession", "file"])]
    if len(set(drawn)) != len(drawn):
        raise InputError(f"{heldout}: duplicate runs")
    known = set(drawn)
    ref = read(f"{d}/reference.tsv", ["run", "class", "organism_verified", "label", "protease", "cys", "reporter_ions",
                                      "prec_window", "prec_unit", "frag_window", "frag_unit"])
    res = read(f"{d}/results.tsv", ["run", "arm", "status", "psm", "pep", "pg", "psm_fdp", "pep_fdp"])
    dec = read(f"{d}/decisions.tsv", ["run", "arm", "setting", "value", "unit", "status"])
    chk = read(f"{d}/checks.tsv", ["run", "source", "arm", "kind", "check", "verdict", "seconds"])
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
    validate_checks(chk, "checks.tsv")
    in_manifest = {f"{r['accession']}/{r['file']}" for r in scn}
    for r in chk:
        for c in ("seconds",):
            if r[c] is not None:
                r[c] = _number(r[c], f"checks.tsv {r['run']} {c}")
        if r["source"] == "scenario":
            if r["run"] not in in_manifest:
                raise InputError(f"checks.tsv: scenario run {r['run']} is not in the manifest")
        elif r["run"] not in known:
            raise InputError(f"checks.tsv: run {r['run']} is not in the draw")
        if r["source"] == "auto" and r["kind"] != "correct":
            raise InputError(f"checks.tsv {r['run']}: the auto arm's own checks are kind correct")
        if r["source"] == "natural" and r["arm"] is None:
            raise InputError(f"checks.tsv {r['run']}: a natural-failure check names its arm")
    for rows, where in ((res, "results.tsv"), (dec, "decisions.tsv"), (tim, "timing.tsv")):
        for r in rows:
            if r["run"] not in known:
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

PRIMARY_ARMS = ("auto", "reference")   # a run executed when both have a result that is not missing


def accounting(drawn, ref, tab=None):
    """every drawn run in one class (§2.2); eligible = bottom-up with a verified organism that executed: its auto and
    reference arms each have a result that is not missing (tab given). A run without a reference row, or a verified
    bottom-up run that did not execute, is missing; "not executed" counts the latter. The 120-run minimum is applied
    to the eligible runs, so failed runs never count toward it."""
    by = {r["run"]: r for r in ref}
    cls, not_exec = {}, 0
    for run in drawn:
        r = by.get(run)
        c = "missing" if r is None else "organism not verified" if (
            r["class"] == "bottom-up" and r["organism_verified"] != "yes") else r["class"]
        if c == "bottom-up" and tab is not None and any(
                (tab.get((run, a)) or {}).get("status", "missing") == "missing" for a in PRIMARY_ARMS):
            c, not_exec = "missing", not_exec + 1
        cls[run] = c
    counts = {k: sum(v == k for v in cls.values()) for k in
              ("bottom-up", "top-down", "no-proteome", "organism not verified", "missing")}
    return {r for r, c in cls.items() if c == "bottom-up"}, {"drawn": len(drawn), **counts,
                                                              "not executed": not_exec}, cls


def eligible_before_arms(d, heldout):
    """the conditional top-up rule's count (§2.1), before any arm runs: the drawn runs of heldout that are bottom-up
    with a verified organism (reference.tsv; class no-proteome = no reference database) and a complete conversion
    (DIR/conversion.tsv: run, status ok | failed; failed covers download, conversion and the MS2-count check). It
    reads nothing of any arm or endpoint. Every drawn run must have a reference row or a failed conversion, else the
    annotation is incomplete and the rule may not run. Returns (eligible count, the accounting counts)."""
    drawn = [f"{r['accession']}/{r['file']}" for r in read(heldout, ["accession", "file"])]
    ref = read(f"{d}/reference.tsv", ["run", "class", "organism_verified"])
    conv = read(f"{d}/conversion.tsv", ["run", "status"])
    _enum(ref, "class", CLASSES, "reference.tsv")
    _enum(ref, "organism_verified", {"yes", "no"}, "reference.tsv")
    _enum(conv, "status", {"ok", "failed"}, "conversion.tsv")
    _unique(ref, ["run"], "reference.tsv")
    _unique(conv, ["run"], "conversion.tsv")
    known = set(drawn)
    for r in ref + conv:
        if r["run"] not in known:
            raise InputError(f"run {r['run']} is not in the draw")
    failed = {r["run"] for r in conv if r["status"] == "failed"}
    ok = {r["run"] for r in conv if r["status"] == "ok"}
    has_ref = {r["run"] for r in ref}
    open_ = [r for r in drawn if r not in failed and not (r in ok and r in has_ref)]
    if open_:
        raise InputError(f"{len(open_)} drawn runs are neither annotated and converted nor failed: {open_[:3]}")
    runs, counts, _ = accounting(drawn, [r for r in ref if r["run"] not in failed])
    return len(runs), counts


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


def p3(tab, runs, rng=None, margin=None):
    """P3 (§5), two gates, both required: relative, the 5th bootstrap percentile of auto share minus reference share
    >= P3_MARGIN; absolute, the one-sided 95 % Clopper-Pearson lower bound of auto's share >= P3_ABS"""
    margin = P3_MARGIN_FRAC if margin is None else margin
    pairs, missing = p3_pairs(tab, runs)
    if not pairs:
        return {"runs": 0, "missing": missing}
    n = len(pairs)
    k, kr = sum(a for a, _ in pairs), sum(b for _, b in pairs)
    lo = p3_relative_lower(pairs, rng)
    return {"runs": n, "missing": missing, "auto_share": k / n, "reference_share": kr / n, "difference": (k - kr) / n,
            "difference_lower95": lo / n, "relative_pass_exact": lo >= margin * n,
            "auto_share_lower95": clopper_lower(k, n), "reference_share_lower95": clopper_lower(kr, n)}


P3_MARGIN_FRAC = Fraction(-1, 20)   # P3_MARGIN as an exact fraction: the boundary -0.05 itself passes


def p3_relative_lower(pairs, rng=None, b=None):
    """the one-sided 95 % lower bound of the paired difference in compliant runs, as a signed count: per run +1 (auto
    compliant, reference not), -1 (the reverse) or 0; the 5th percentile (order statistic) of the bootstrap sums.
    Integers throughout, so the -0.05 boundary is exact (score.py and power.py)"""
    sgn = np.array([int(a) - int(c) for a, c in pairs], dtype=np.int64)
    sums = boot(lambda s: s.sum(axis=-1), sgn, rng, b)
    return int(np.quantile(sums, 0.05, method="inverted_cdf"))


# ---- secondary endpoints and §9 rules --------------------------------------------------------------------------------

def _ppm(v, unit):
    return v if unit == "ppm" else v / MZ_REF * 1e6


def _norm(v):
    return str(v).strip().lower()


def canon_label(v):
    """the S1 category of a label name (LABEL_CANON); an unknown name stays itself"""
    return LABEL_CANON.get(_norm(v), _norm(v))


def label_family(v):
    """the S1 family line: isobaric (TMT, TMTpro, iTRAQ), SILAC, none, else the category (15N, 2H2O)"""
    c = canon_label(v)
    return "isobaric" if c.lower().startswith(ISOBARIC) else "SILAC" if c.lower().startswith("silac") else c


def accuracy(ref, dec, runs, arm="auto"):
    """S1-S5 over the eligible runs with a reference value: exact agreement (windows: within 1.5x after conversion to
    ppm at m/z MZ_REF); a run without a decision row counts as wrong ("absent"); abstentions are counted"""
    truth = {r["run"]: r for r in ref}
    d = {(x["run"], x["setting"]): x for x in dec if x["arm"] == arm}
    out = {}
    for s, unit in (("label", None), ("protease", None), ("cys", None), ("prec_window", "prec_unit"),
                    ("frag_window", "frag_unit")):
        k = n = absent = abst = fam = 0
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
            elif s == "label":
                k += canon_label(x["value"]) == canon_label(t)
                fam += label_family(x["value"]) == label_family(t)
            else:
                k += _norm(x["value"]) == _norm(t)
        out[s] = {"correct": k, "n": n, "share": k / n if n else None, "wilson95": list(wilson(k, n)),
                  "absent": absent, "abstained": abst}
        if s == "label":
            out[s]["family_correct"] = fam
            out[s]["family_share"] = fam / n if n else None
    return out


def scenario_counts(d, scn, pilot=False):
    """scenario_counts.tsv: required when scenarios were built (the pilot predates it), unique (run, kind), and every
    built scenario covered"""
    built = built_scenarios(scn)
    rows_ = read(f"{d}/scenario_counts.tsv", ["run", "kind", "status", "psm"], required=bool(built) and not pilot)
    _unique(rows_, ["run", "kind"], "scenario_counts.tsv")
    _enum(rows_, "status", {"ok", "missing"}, "scenario_counts.tsv")
    if rows_ or not pilot:
        miss = sorted(built - {(r["run"], r["kind"]) for r in rows_})
        if miss:
            raise InputError(f"scenario_counts.tsv: built scenarios without a row: {miss[:5]}")
    return rows_


def built_scenarios(scn):
    return {(f"{r['accession']}/{r['file']}", r["kind"]) for r in scn
            if r["kind"] not in ("excluded", "unresolved") and r["databases"] not in (None, "none")}


NATURAL_ARMS = ("sdrf", "defaults")


def natural_cohort(tab, runs):
    """the natural failures (§6 S6): (run, arm) of the eligible runs where arm (SDRF or defaults) has a result with
    strictly fewer than half the reference arm's PSMs (a stopped arm has 0); derived from the search counts, never
    from which checks arrived"""
    out = set()
    for run in runs:
        ref = tab.get((run, "reference"))
        if not ref or ref["status"] != "ok":
            continue
        for arm in NATURAL_ARMS:
            a = tab.get((run, arm))
            if a and a["status"] != "missing" and _count(a, "psm") < 0.5 * ref["psm"]:
                out.add((run, arm))
    return out


def checks(chk, scn, runs, tab=None, counts=None):
    """S6 and the §9 check gates. Specificity (kind correct): the auto arm's own checks (source auto) on the eligible
    runs. Sensitivity: the scenario checks (source scenario) of the built scenarios of the frozen manifest
    (wrong_label and wrong_protease: the single-organism runs with a built correct scenario); the correct-scenario
    controls are reported per check against the manifest. Natural failures: the cohort of natural_cohort (tab given);
    natural checks outside it are counted and ignored. A scenario without a check result is missing; not assessed is
    left out; both are counted. Each of the nine gates: exact one-sided binomial test of H0 share <= bound, Holm over
    all nine at ALPHA; a gate without an assessed scenario keeps p = 1 and does not pass. Also per kind: the median
    PSM ratio scenario / correct scenario (the loss, from scenario_counts.tsv: every built scenario's search, whatever
    its checks gave; a scenario or control without an ok count is counted, a control with 0 PSMs gives no ratio) and
    the median seconds of each check."""
    built = built_scenarios(scn)
    multi = {r for r, k in built if k == "incomplete"}
    rows = {(r["run"], r["source"], r["arm"], r["kind"], r["check"]): r for r in chk}
    cnt = {(r["run"], r["kind"]): r for r in counts or []}
    ok_psm = lambda run, kind: (lambda r: _number(r["psm"], "scenario_counts.tsv") if r and r["status"] == "ok"
                                and r["psm"] is not None else None)(cnt.get((run, kind)))
    cohort = natural_cohort(tab, runs) if tab is not None else \
        {(r, a) for r, s_, a, k, c in rows if s_ == "natural"}
    out, pv = {}, {f"{c}:{k}": 1.0 for c, k in GATES}
    for chk_ in sorted(CHECKS):
        for kind in sorted(KINDS):
            if kind == "natural":
                keys = [(r, "natural", a, "natural", chk_) for r, a in sorted(cohort)] if chk_ == "post" else []
            elif kind == "correct":
                keys = [(r, "auto", None, "correct", chk_) for r in sorted(runs)]
            elif kind in ("wrong_label", "wrong_protease"):
                keys = [(r, "scenario", None, kind, chk_) for r in sorted(r for r, k in built if k == "correct"
                                                                             and r not in multi)] if chk_ == "post" else []
            else:
                keys = [(r, "scenario", None, kind, chk_) for r in sorted(r for r, k in built if k == kind)]
            gated = (chk_, kind) in CHECK_BOUND
            if not keys and not gated:
                continue
            got = [rows.get(k_) for k_ in keys]
            v = [g["verdict"] if g else None for g in got]
            ass = [x for x in v if x not in (None, "not assessed")]
            hit = sum((x == "ok") if kind == "correct" else (x != "ok") for x in ass)
            exact = sum(x == "settings" for x in ass) if kind in ("wrong_label", "wrong_protease") else \
                sum(x == "database" for x in ass) if chk_ == "post" and kind not in ("correct", "natural") else None
            secs = [g["seconds"] for g in got if g and g["seconds"] is not None]
            e = {"scenarios": len(keys), "missing": v.count(None), "not_assessed": v.count("not assessed"),
                 "assessed": len(ass), "hit": hit, "share": hit / len(ass) if ass else None,
                 "wilson95": list(wilson(hit, len(ass))), "exact_verdict": exact,
                 "seconds_median": float(np.median(secs)) if secs else None,
                 "bound": CHECK_BOUND.get((chk_, kind)),
                 "status": "gate" if gated else "experimental" if (chk_, kind) in EXPERIMENTAL else "reported"}
            if kind not in ("correct", "natural") and chk_ == "post":
                loss, no_count, no_ctl = [], 0, 0
                for r, *_ in keys:
                    x, c = ok_psm(r, kind), ok_psm(r, "correct")
                    if x is None:
                        no_count += 1
                    elif not c:
                        no_ctl += 1
                    else:
                        loss.append(x / c)
                e.update({"psm_ratio_to_correct_median": float(median(loss)) if loss else None,
                          "psm_ratio_n": len(loss), "without_count": no_count,
                          "without_control_count": no_ctl})
            if gated:
                e["p_value"] = binom_sf(hit, len(ass), e["bound"]) if ass else 1.0
                pv[f"{chk_}:{kind}"] = e["p_value"]
            out[f"{chk_}:{kind}"] = e
    ctl = sorted(r for r, k in built if k == "correct")
    for chk_ in sorted(CHECKS):
        g = [rows.get((r, "scenario", None, "correct", chk_)) for r in ctl]
        out[f"scenario correct (controls) {chk_}"] = {
            "scenarios": len(ctl), "checks": sum(x is not None for x in g),
            "ok": sum(x is not None and x["verdict"] == "ok" for x in g),
            "not_assessed": sum(x is not None and x["verdict"] == "not assessed" for x in g)}
    out["natural checks outside the cohort"] = {"checks": sum(1 for r, s_, a, k, c in rows if s_ == "natural"
                                                              and (r, a) not in cohort)}
    for k, rej in holm(pv).items():
        out[k]["pass"] = rej
    return out


def _cys_rule(tab, dec, a1, a2, runs):
    """on the eligible runs whose Cys decision differs between arms a1 and a2: summed PSM ratio a1/a2 >= 0.99 and no
    run below 0.95. A run with either arm missing is left out and counted; a stopped arm has 0 PSMs; a run with 0 in
    both arms is uninformative, left out and counted (it is never evidence); x/0 with x > 0 is +inf. Fewer than
    MIN_INFORMATIVE informative runs: no verdict."""
    v = {(d["run"], d["arm"]): d["value"] for d in dec if d["setting"] == "cys"}
    diff = sorted(r for r in runs if (r, a1) in v and (r, a2) in v and v[(r, a1)] != v[(r, a2)])
    pairs, missing, zero = [], 0, 0
    for r in diff:
        x, y = tab.get((r, a1)), tab.get((r, a2))
        if not x or not y or "missing" in (x["status"], y["status"]):
            missing += 1
            continue
        a, b = _count(x, "psm"), _count(y, "psm")
        if a == 0 and b == 0:
            zero += 1
            continue
        pairs.append((a, b))
    out = {"runs": len(pairs), "differing": len(diff), "missing": missing, "uninformative": zero}
    if len(pairs) < MIN_INFORMATIVE:
        return {**out, "pass": None}
    sa, sb = sum(a for a, _ in pairs), sum(b for _, b in pairs)
    s = sa / sb if sb else math.inf
    worst = min(a / b if b else math.inf for a, b in pairs)
    return {**out, "summed_ratio": s, "worst": worst, "pass": s >= 0.99 and worst >= 0.95}


def rules(tab, ref, dec, runs, acc):
    """§9 detector rules other than the check gates"""
    truth = {r["run"]: r for r in ref}
    out = {"cys: -alkylation_out (auto_never) vs competition (auto_always)":
           _cys_rule(tab, dec, "auto_never", "auto_always", runs),
           "mini-search fallback (auto) vs always (auto_always)": _cys_rule(tab, dec, "auto", "auto_always", runs)}
    lab = [r for r in runs if truth[r].get("label") is not None]
    base = sum(canon_label(truth[r]["label"]) == "none" for r in lab)
    d = {(x["run"], x["setting"], x["arm"]): x for x in dec}
    bad = [r for r in runs if truth[r].get("reporter_ions") == "no" and (r, "label", "auto") in d
           and d[(r, "label", "auto")]["value"] and _norm(d[(r, "label", "auto")]["value"]).startswith(ISOBARIC)]
    out["label: S1 >= always-label-free baseline, no isobaric call without reporter ions"] = {
        "auto": acc["label"]["correct"], "baseline": base, "n": len(lab), "isobaric_without_reporters": len(bad),
        "pass": (acc["label"]["correct"] >= base and not bad) if lab else None}
    ok = lambda r, arm: (r, "protease", arm) in d and d[(r, "protease", arm)]["value"] is not None
    with_ref = [r for r in runs if truth[r].get("protease") is not None]
    prot = [r for r in with_ref if ok(r, "sdrf")]   # the same runs for both: a reference and an SDRF protease
    hit = lambda r, arm: ok(r, arm) and _norm(d[(r, "protease", arm)]["value"]) == _norm(truth[r]["protease"])
    a_, s_ = sum(hit(r, "auto") for r in prot), sum(hit(r, "sdrf") for r in prot)
    out["protease: S2 >= stated-protease (SDRF arm) baseline"] = {
        "auto": a_, "baseline": s_, "n": len(prot), "without_sdrf_protease": len(with_ref) - len(prot),
        "pass": a_ >= s_ if prot else None}
    # the denominator: every eligible run with header evidence; FASTag QC's own call (none when its 90 % rule made
    # no call, or no QC ran) disagrees unless it equals the header
    head = [r for r in runs if ((d.get((r, "ms2_analyzer_header", "auto")) or {}).get("value"))]
    qcv = lambda r: (d.get((r, "ms2_analyzer_qc", "auto")) or {}).get("value")
    called = [r for r in head if qcv(r) not in (None, QC_NO_CALL)]
    agree = sum(_norm(qcv(r)) == _norm(d[(r, "ms2_analyzer_header", "auto")]["value"]) for r in called)
    out["FASTag QC: analyser class agrees with the header in >= 98 %"] = {
        "agree": agree, "n": len(head), "called": len(called), "no_call": len(head) - len(called),
        "without_header": len(runs) - len(head), "pass": agree / len(head) >= 0.98 if head else None}
    return out


def timing(tim):
    """S7 (descriptive): median seconds per input type, version and step"""
    g = {}
    for r in tim:
        g.setdefault((r["input_type"], r["version"], r["step"]), []).append(_number(r["seconds"], "timing.tsv"))
    return {" | ".join(k): {"runs": len(v), "median_s": float(np.median(v))} for k, v in sorted(g.items())} \
        if g else "timing.tsv not provided"


def primary(tab, runs, pilot=False, p3_abs=None, p3_margin=None):
    """P1, P2, P3 and the verdict (§5, §9), the one routine score.py and power.py share: every bootstrap restarts
    the frozen seed with B resamples; fewer than MIN_ELIGIBLE eligible runs (or the pilot): no verdict. p3_abs and
    p3_margin override the frozen P3 gates (power.py's design grid only)."""
    p3_abs = P3_ABS if p3_abs is None else p3_abs
    p3_margin = P3_MARGIN_FRAC if p3_margin is None else p3_margin
    out = {"P1": ratio_endpoint(tab, runs, "auto", "reference"), "P2": ratio_endpoint(tab, runs, "auto", "defaults"),
           "P3": p3(tab, runs, margin=p3_margin)}
    enough = len(runs) >= MIN_ELIGIBLE and not pilot   # the pilot (PREREG §0) never gets a verdict
    out["P1"]["pass"] = (out["P1"]["ci95"][0] is not None and out["P1"]["ci95"][0] >= P1_MARGIN) if enough else None
    out["P2"]["pass"] = (out["P2"]["ci95"][0] is not None and out["P2"]["ci95"][0] > 1) if enough else None
    e = out["P3"]
    e["relative_pass"] = bool(e.get("relative_pass_exact"))
    e["absolute_pass"] = e.get("auto_share_lower95") is not None and e["auto_share_lower95"] >= p3_abs
    e["pass"] = (e["relative_pass"] and e["absolute_pass"]) if enough else None
    out["pass"] = (out["P1"]["pass"] and out["P2"]["pass"] and out["P3"]["pass"]) if enough else None
    out["verdict"] = ("pass" if out["pass"] else "fail") if enough else \
        "pilot: no verdict" if pilot else f"no verdict: {len(runs)} eligible runs, {MIN_ELIGIBLE} needed"
    if not enough:
        e["relative_pass"] = e["absolute_pass"] = None
    return out


def score(d, heldout=None, pilot=False):
    heldout = heldout or os.path.join(os.path.dirname(os.path.abspath(__file__)), "heldout.tsv")
    drawn, ref, res, dec, chk, scn, tim = load(d, heldout)
    cnt_rows = scenario_counts(d, scn, pilot)
    tab = {(r["run"], r["arm"]): r for r in res}
    runs, counts, cls = accounting(drawn, ref, tab)
    out = {"runs": counts, "class": cls, **primary(tab, runs, pilot)}
    arms = sorted({a for _, a in tab} - {"auto"})
    out["reported"] = {f"{f} auto/{m}": ratio_endpoint(tab, runs, "auto", m, f)
                       for m in arms for f in ("psm", "pep", "pg") if not (f == "psm" and m in ("reference", "defaults"))}
    out["S1-S5"] = acc = accuracy(ref, dec, runs)
    out["S6"] = checks(chk, scn, runs, tab, cnt_rows)
    out["S7"] = timing(tim)
    out["S8 top-down"] = {"top-down runs": counts["top-down"],
                          "stopped": sum(1 for r, c in cls.items() if c == "top-down"
                                         and (tab.get((r, "auto")) or {}).get("status") == "stopped"),
                          "bottom-up stopped": sum(1 for r in runs
                                                   if (tab.get((r, "auto")) or {}).get("status") == "stopped")}
    out["rules"] = rules(tab, ref, dec, runs, acc)
    return out


def jsonable(o):
    """score.json in strict JSON: +inf as the string "inf", None as null, sets and tuples as lists, any other
    non-JSON value as its string; NaN never occurs and stops the dump (allow_nan=False)"""
    if isinstance(o, float) and math.isinf(o):
        return "inf" if o > 0 else "-inf"
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [jsonable(v) for v in o]
    if o is None or isinstance(o, (bool, int, float, str)):
        return o
    return str(o)


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
        if not isinstance(e, dict) or "status" not in e:
            continue
        lines.append(f"| {k} | {e['status']} | {e['assessed']} | {f(e['share'])} | {f(e['wilson95'][0])} | "
                     f"{f(e['bound'])} | {e.get('pass')} |")
    lines += ["", "| §9 rule | pass |", "|---|---|"] + [f"| {k} | {v.get('pass')} |" for k, v in s["rules"].items()]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    args = [x for x in sys.argv[1:] if x != "--pilot"]
    d = args[0]
    s = score(d, args[1] if len(args) > 1 else None, pilot="--pilot" in sys.argv)
    json.dump(jsonable(s), open(os.path.join(d, "score.json"), "w"), indent=1, allow_nan=False)
    open(os.path.join(d, "score.md"), "w").write(markdown(s))
    print(markdown(s))

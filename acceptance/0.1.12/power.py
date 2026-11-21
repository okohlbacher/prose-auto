"""Power of the PREREG 0.1.12 criteria, through score.py's own decision routine (score.primary).

    python power.py PRIMARY CHECKS [sims] [--grid] [--overwrite] [--workers=N]
        -> power.txt (or power_grid.txt) next to this file, and stdout

PRIMARY: the pilot's scoring inputs (PREREG §0: drawn.tsv, reference.tsv, results.tsv, decisions.tsv, checks.tsv,
scenarios.tsv). CHECKS: dev_checks.tsv (run, source, arm, kind, check, verdict; the 0.1.11 check calibration, validated
with score.validate_checks) and related_ranks.tsv (check, rank, assessed, hit; the related_distant calibration by the
rank of the common ancestor, v0112/related_strata.txt). The frozen copies are power_inputs/ next to this file.
An existing power.txt or power_grid.txt is never overwritten without --overwrite (the files are frozen, §10).

Primary endpoints. The held-out design is the manifest: heldout.tsv's runs, of which those tagged top-down
(experiment field) and those refdb_draw.tsv marks no-proteome are fixed out of P1-P3. Every other run is drawn from
the pilot's drawn runs as a whole: its class (eligible, organism not verified, missing) and all its arm results
together, so attrition and effects come from the same observed runs. The pilot's metadata strata were often wrong
(§0) and most strata have 0 or 1 pilot runs, so no stratification is used: stratum-specific power is unknown. Each
simulated set goes through score.accounting and score.primary unchanged (the 120-run minimum, eligibility, B
resamples with the frozen seed restarted per endpoint, the exact P3 boundary), so a simulated pass is a scorer pass.
Scenarios (every row reports eligible runs, mean and 5th-95th percentile, P(no verdict), P1, P2, P3 and joint power):
  conditional   the pilot as the truth;
  predictive    each simulation first resamples the 19 pilot runs (the pilot's own uncertainty, both attrition and
                effects);
  attrition     the eligible share at the pilot's two-sided 95 % Clopper-Pearson lower bound (an extra loss on top
                of the pilot's), conditional and predictive;
  shifts        auto PSMs x0.98 / x0.97, and 2 or 5 % of the compliant auto runs made noncompliant;
  projects      leave one pilot project out (range and the most influential project), and every pilot FDP within
                0.05 points of the 1.5 % limit moved across it, one at a time.
--grid: joint power for 187 (the first draw), 250, 300 (the adopted draw) and 400 drawn runs (the fixed top-down and
no-proteome shares scaled to the size), P3 bounds X in {0.70, 0.75, 0.80} and margins {-0.05, -0.10}, at the pilot
effect and two worse ones; SIMS // 4 simulations per cell (at least 200), stated in the output.

Checks. Per gate the expected assessed count is the expected built count times the development share assessed;
for specificity (kind correct) the built count is the 5th percentile of the eligible runs in the predictive
attrition scenario, so the bounds hold when attrition is high. The development specificity rows are the 0.1.11
correct-database scenario checks (source scenario; development has no auto-arm checks); the held-out gate uses the
auto arm's own checks. related_distant uses the rank mix of the held-out manifest (common_rank) and the rank-specific
development shares. Each bound is the largest (step 0.05, from 0.10) with exact single-gate power >= 0.80 at ALPHA/9
(Bonferroni), given at least 10 expected assessed scenarios. Holm rejects whatever Bonferroni rejects, so these exact
powers are lower bounds of each gate's power under Holm, whatever the dependence between the checks."""
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as SC   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
P = C = None            # set by main() from the command line (never at import: tests import this module)
GRID, WORKERS, SIMS = False, 1, 1000
FDP_NEAR = 0.05          # pilot FDPs within this many points of the limit are moved across it (sensitivity)


def pilot_units(p):
    """the pilot's drawn runs as units: (run, class after score.accounting, {arm: result row})"""
    drawn, ref, res, *_ = SC.load(p, os.path.join(p, "drawn.tsv"))
    tab = {(r["run"], r["arm"]): r for r in res}
    _, counts, cls = SC.accounting(drawn, ref, tab)
    return [(r, cls[r], {a: row for (rr, a), row in tab.items() if rr == r}) for r in drawn], counts


def design(n_total):
    """(runs resampled from the pilot, fixed top-down, fixed no-proteome) for a held-out set of n_total drawn runs:
    the exact manifests at 187 (the first draw) and 300, scaled shares otherwise"""
    held = SC.read(os.path.join(HERE, "heldout.tsv"), ["accession", "file", "experiment"])
    nop = {(r["accession"], r["file"]) for r in SC.read(os.path.join(HERE, "refdb_draw.tsv"),
                                                          ["accession", "file", "class"]) if r["class"] == "no-proteome"}
    rows = held[:187] if n_total == 187 else held
    td = sum("top-down" in (r["experiment"] or "").lower() for r in rows)
    np_ = sum((r["accession"], r["file"]) in nop and "top-down" not in (r["experiment"] or "").lower() for r in rows)
    if n_total not in (187, len(held)):
        td, np_ = round(td * n_total / len(rows)), round(np_ * n_total / len(rows))
    return n_total - td - np_, td, np_


def simulate(units, m, rng, extra_loss=0.0, shift=1.0, flip=0.0, x=None, margin=None):
    """one simulated held-out set of m resampled runs through score.accounting and score.primary"""
    tab, drawn, ref = {}, [], []
    for i in range(m):
        run, cls, arms = units[rng.integers(len(units))]
        key = f"sim{i}"
        drawn.append(key)
        if cls == "missing" or (cls == "bottom-up" and extra_loss and rng.random() < extra_loss):
            continue                                    # no reference row: missing
        ref.append({"run": key, "class": "bottom-up" if cls in ("bottom-up", "organism not verified") else cls,
                    "organism_verified": "yes" if cls == "bottom-up" else "no"})
        for a, row in arms.items():
            row = dict(row)
            if a == "auto" and row["status"] == "ok":
                if shift != 1.0:
                    row["psm"] = row["psm"] * shift
                if flip and SC.compliant(row) and rng.random() < flip:
                    row["psm_fdp"] = SC.FDP_MAX + 1
            tab[(key, a)] = row
    runs, _, _ = SC.accounting(drawn, ref, tab)
    out = SC.primary(tab, runs, p3_abs=x, p3_margin=margin)
    return len(runs), out["verdict"].startswith("no verdict"), bool(out["P1"]["pass"]), \
        bool(out["P2"]["pass"]), bool(out["P3"]["pass"]), bool(out["pass"])


def task(args):
    """one scenario: (label, units, m, sims, seed, predictive, extra_loss, shift, flip, x, margin) -> summary"""
    label, units, m, sims, seed, predictive, extra_loss, shift, flip, x, margin = args
    rng = np.random.default_rng([SC.SEED, seed])
    el, nov, p1, p2, p3, joint = [], 0, 0, 0, 0, 0
    for _ in range(sims):
        u = [units[i] for i in rng.integers(len(units), size=len(units))] if predictive else units
        n, nv, a, b, c, d = simulate(u, m, rng, extra_loss, shift, flip, x, margin)
        el.append(n)
        nov, p1, p2, p3, joint = nov + nv, p1 + a, p2 + b, p3 + c, joint + d
    q = np.percentile(el, [5, 95])
    return {"label": label, "sims": sims, "eligible": float(np.mean(el)), "eligible5": float(q[0]),
            "eligible95": float(q[1]), "no_verdict": nov / sims, "P1": p1 / sims, "P2": p2 / sims, "P3": p3 / sims,
            "joint": joint / sims}


def run_tasks(ts):
    if WORKERS > 1:
        with ProcessPoolExecutor(WORKERS) as ex:
            return list(ex.map(task, ts))
    return [task(t) for t in ts]


def line(e):
    return (f"{e['label']}: eligible {e['eligible']:.0f} (5-95 % {e['eligible5']:.0f}-{e['eligible95']:.0f}), "
            f"P(no verdict) {e['no_verdict']:.2f}; power P1 {e['P1']:.2f}, P2 {e['P2']:.2f}, P3 {e['P3']:.2f}, "
            f"joint {e['joint']:.2f} ({e['sims']} simulations)")


def attrition_loss(counts):
    """the extra loss that brings the eligible share to the two-sided 95 % Clopper-Pearson lower bound"""
    k, n = counts["bottom-up"], counts["drawn"] - counts["top-down"] - counts["no-proteome"]
    lo = SC.clopper_lower(k, n, alpha=0.025)
    return 1 - lo / (k / n), lo, k, n


def exact_power(n, p, bound, m):
    """P(the exact one-sided binomial test of share <= bound rejects at ALPHA/m) when the true share is p"""
    kc = next((x for x in range(n + 1) if SC.binom_sf(x, n, bound) <= SC.ALPHA / m), None)
    return 0.0 if kc is None else SC.binom_sf(kc, n, p)


def check_section(n_correct):
    dchk = SC.read(os.path.join(C, "dev_checks.tsv"), ["run", "source", "arm", "kind", "check", "verdict"])
    SC.validate_checks(dchk, "dev_checks.tsv")
    ranks = SC.read(os.path.join(C, "related_ranks.tsv"), ["check", "rank", "assessed", "hit"])
    mani = SC.read(os.path.join(HERE, "scenarios.tsv"), ["accession", "file", "kind", "databases", "common_rank"])
    built = SC.built_scenarios(mani)
    multi = {r for r, k_ in built if k_ == "incomplete"}
    n_built = {k_: sum(1 for _, kk in built if kk == k_) for k_ in ("related_close", "related_distant", "unrelated",
                                                                    "incomplete")}
    n_built["wrong_label"] = n_built["wrong_protease"] = sum(1 for r, k_ in built if k_ == "correct" and r not in multi)
    n_built["correct"] = n_correct
    mix = {}   # the rank of the common ancestor of each built related_distant scenario (its first species)
    for r in mani:
        if r["kind"] == "related_distant" and r["databases"] not in (None, "none"):
            rk = (r["common_rank"] or "").split(";")[0].strip()
            mix[rk] = mix.get(rk, 0) + 1
    m = len(SC.GATES)
    srcs = sorted({r["source"] for r in dchk if r["kind"] == "correct"})
    out = [f"checks: the expected specificity count is {n_correct} eligible runs (the 5th percentile of the predictive "
           f"attrition scenario); development specificity from source {', '.join(srcs)} (the correct-database "
           f"scenario checks of the 0.1.11 calibration); related_distant rank mix of the manifest: "
           + ", ".join(f"{k} {v}" for k, v in sorted(mix.items()))]
    gates = sorted((c, k_) for c in ("pre", "post") for k_ in ("correct", "unrelated", "related_close",
                   "related_distant", "incomplete", "wrong_label", "wrong_protease")
                   if not (c == "pre" and k_ in ("wrong_label", "wrong_protease")))
    for c, k_ in gates:
        v = [r["verdict"] for r in dchk if r["check"] == c and r["kind"] == k_ and r["source"] == "scenario"]
        ass = [x for x in v if x != "not assessed"]
        if not ass:
            out.append(f"{c}:{k_}: no assessed development rows; power unknown")
            continue
        hit = sum((x == "ok") if k_ == "correct" else (x != "ok") for x in ass)
        p, sa, note, worst = hit / len(ass), len(ass) / len(v), "", None
        if k_ == "related_distant":
            by = {r["rank"]: (int(r["hit"]), int(r["assessed"])) for r in ranks if r["check"] == c}
            known = {rk: n_ for rk, n_ in mix.items() if rk in by and by[rk][1]}
            if known:
                p_rank = sum(n_ * by[rk][0] / by[rk][1] for rk, n_ in known.items()) / sum(known.values())
                worst = min(by[rk][0] / by[rk][1] for rk in known)
                note = ("; by rank " + ", ".join(f"{rk} {by[rk][0]}/{by[rk][1]}" for rk in sorted(known))
                        + f", rank-weighted {p_rank:.3f}, worst rank {worst:.3f}")
                p = p_rank
        n = int(round(n_built[k_] * sa))
        best = max((b / 100 for b in range(10, 100, 5) if exact_power(n, p, b / 100, m) >= 0.8), default=None) \
            if n >= 10 else None
        frozen = SC.CHECK_BOUND.get((c, k_), "experimental" if (c, k_) in SC.EXPERIMENTAL else "reported")
        fp = ""
        if isinstance(frozen, float):
            fp = f", exact power at it {exact_power(n, p, frozen, m):.3f} (a lower bound under Holm)"
            if worst is not None:
                fp += f", {exact_power(n, worst, frozen, m):.3f} at the worst rank"
        out.append(f"{c}:{k_}: development {p:.3f} ({hit}/{len(ass)} pooled{note}), held-out n={n}; largest bound "
                   f"with exact power >= 0.80 at alpha/{m}: {best}; frozen: {frozen}{fp}")
    return out


def main(argv):
    global P, C, GRID, WORKERS, SIMS
    args = [a for a in argv if not a.startswith("--")]
    P, C = (args + [None, None])[:2]
    SIMS = int(args[2]) if len(args) > 2 else 1000
    GRID = "--grid" in argv
    WORKERS = int(next((a.split("=")[1] for a in argv if a.startswith("--workers=")), "1"))
    target = os.path.join(HERE, "power_grid.txt" if GRID else "power.txt")
    if os.path.exists(target) and "--overwrite" not in argv:
        sys.exit(f"{target} exists (frozen); pass --overwrite to regenerate it")
    units, counts = pilot_units(P)
    loss, lo, k, n = attrition_loss(counts)
    cmd = "python power.py " + " ".join(a for a in argv if a != "--overwrite")
    out = [f"command: {cmd}",
           f"pilot: {counts}; eligible share {k}/{n}, two-sided 95 % lower bound {lo:.3f} (extra loss {loss:.3f})"]
    m300, td, np_ = design(300)
    out.append(f"held-out design: 300 drawn, {td} tagged top-down and {np_} no-proteome fixed, {m300} drawn from the "
               f"{len(units)} pilot runs")
    if GRID:
        sims = max(200, SIMS // 4)
        out.append(f"design grid: joint power (P1 / P2 / P3 in brackets), {sims} simulations per cell, conditional:")
        cells = []
        for nd in (187, 250, 300, 400):
            mm = design(nd)[0]
            for x in (0.70, 0.75, 0.80):
                for mg in (SC.P3_MARGIN_FRAC, Fraction(-1, 10)):
                    for j, (sh, fl) in enumerate(((1.0, 0.0), (0.98, 0.02), (0.97, 0.05))):
                        cells.append((f"{nd}|{x}|{float(mg)}|{j}", units, mm, sims, len(cells), False, 0.0, sh, fl,
                                      x, mg))
        res = {r["label"]: r for r in run_tasks(cells)}
        f = lambda r: f"{r['joint']:.2f} ({r['P1']:.2f}/{r['P2']:.2f}/{r['P3']:.2f})"
        for nd in (187, 250, 300, 400):
            for x in (0.70, 0.75, 0.80):
                for mg in (-0.05, -0.1):
                    c = [res[f"{nd}|{x}|{mg}|{j}"] for j in range(3)]
                    out.append(f"drawn {nd}, X {x:.2f}, margin {mg:+.2f}: pilot effect {f(c[0])}; x0.98 and 2 % fewer "
                               f"compliant {f(c[1])}; x0.97 and 5 % fewer {f(c[2])}; P(no verdict) "
                               f"{c[0]['no_verdict']:.2f}")
    else:
        base = lambda label, seed, **kw: (label, kw.get("units", units), m300, kw.get("sims", SIMS), seed,
                                          kw.get("pred", False), kw.get("loss", 0.0), kw.get("shift", 1.0),
                                          kw.get("flip", 0.0), None, None)
        ts = [base("conditional (pilot as truth)", 1), base("predictive (pilot resampled)", 2, pred=True),
              base("attrition at the lower bound, conditional", 3, loss=loss),
              base("attrition at the lower bound, predictive", 4, pred=True, loss=loss),
              base("auto PSMs x0.98", 5, shift=0.98), base("auto PSMs x0.97", 6, shift=0.97),
              base("2 % of compliant auto runs noncompliant", 7, flip=0.02),
              base("5 % of compliant auto runs noncompliant", 8, flip=0.05)]
        loo = [base(f"without {u[0].split('/')[0]}", 100 + i, units=[v for v in units if v is not u],
                    sims=max(200, SIMS // 5)) for i, u in enumerate(units)]
        near = []
        for i, (run, cls, arms) in enumerate(units):
            for a in ("auto", "reference"):
                row = arms.get(a)
                if cls == "bottom-up" and row and row["status"] == "ok" and row["psm"] >= SC.MIN_PSM:
                    for fld in ("psm_fdp", "pep_fdp"):
                        if abs(row[fld] - SC.FDP_MAX) <= FDP_NEAR:
                            moved = dict(row)
                            moved[fld] = SC.FDP_MAX if row[fld] > SC.FDP_MAX else SC.FDP_MAX + 0.001
                            u2 = list(units)
                            u2[i] = (run, cls, {**arms, a: moved})
                            near.append(base(f"{run.split('/')[0]} {a} {fld} {row[fld]:.4f} -> {moved[fld]:.4f}",
                                             200 + len(near), units=u2))
        res = run_tasks(ts + loo + near)
        out += [line(e) for e in res[:len(ts)]]
        lr = res[len(ts):len(ts) + len(loo)]
        lo_e, hi_e = min(lr, key=lambda e: e["joint"]), max(lr, key=lambda e: e["joint"])
        out.append(f"leave one pilot project out (conditional, {lr[0]['sims']} simulations each): joint power "
                   f"{lo_e['joint']:.2f}-{hi_e['joint']:.2f}; lowest {lo_e['label']}, highest {hi_e['label']}")
        out += [f"  {line(e)}" for e in sorted(lr, key=lambda e: e["joint"])[:3]]
        out += (["pilot FDPs within 0.05 points of 1.5 % moved across it (conditional):"]
                + [f"  {line(e)}" for e in res[len(ts) + len(loo):]]) if near else \
            ["no pilot FDP within 0.05 points of 1.5 %"]
        out += check_section(int(math.floor(res[3]["eligible5"])))
    txt = "\n".join(out) + "\n"
    open(target, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])

"""Power of the PREREG 0.1.12 criteria, with score.py's own functions.

    python power.py PRIMARY CHECKS [sims] [--grid]   -> power.txt next to this file (or power_grid.txt), and stdout

PRIMARY: scoring inputs (drawn.tsv, reference.tsv, results.tsv, decisions.tsv, checks.tsv, scenarios.tsv) and
pool.tsv (run, pool) of the runs whose paired outcomes stand in for the held-out ones: the pilot (PREREG §0), or the
development ID corpus (v0112/dev_inputs.py: prose-auto 0.1.3 against the optimal settings).
CHECKS: dev_checks.tsv (run, source, arm, kind, check, verdict), the 0.1.11 check calibration, validated with
score.validate_checks (a blank verdict stops the script).

Primary endpoints. Each simulation draws the held-out design by quota: a run of a heldout.tsv stratum is drawn from the
PRIMARY runs of the matching pool (stratum_pool); a pool without runs draws from all runs (said in the output). The
runs tagged top-down are removed and each other run is lost with probability MISSING: the PRIMARY data's share of
drawn runs that were missing, not verified or without a proteome, at least 0.10. score.py's rules then apply
unchanged: the 120-run minimum, PSM eligibility, the ratio conventions, P3's signed-count relative gate and its
Clopper-Pearson absolute gate. Shifts show the margins: auto PSMs x0.98 / x0.97, and 2 or 5 % of the compliant auto
runs made noncompliant. P2 needs a defaults arm on every PRIMARY run, and at least 10 of them (the pilot), else its
power is unknown (the development corpus: 16 data sets of their own, apart from the 228 runs).
--grid: joint power (P1, P2, P3 both gates) for draws of 187 to 400 runs and P3 bounds X in {0.70, 0.75, 0.80} and
margins in {-0.05, -0.10}, at the PRIMARY effect and at the worse ones.

Checks. Per gate, the held-out scenario count (the manifest) times the development share assessed, and the development
share hit. Each bound is the largest (step 0.05, from 0.10) with exact single-gate power >= 0.80 at ALPHA/9 (an exact
binomial sum, no simulation), given at least 10 expected assessed scenarios. The joint power under Holm at the frozen
bounds is simulated (its Monte Carlo standard error is printed)."""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as SC   # noqa: E402

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
P, C = ARGS[0], ARGS[1]
SIMS = int(ARGS[2]) if len(ARGS) > 2 else 1000
GRID = "--grid" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
B_SIM = 2000
rng = np.random.default_rng(SC.SEED)
out = []

drawn, ref, res, dec, chk, scn, tim = SC.load(P, os.path.join(P, "drawn.tsv"))
tab = {(r["run"], r["arm"]): r for r in res}
pool_of = {r["run"]: r["pool"] for r in SC.read(os.path.join(P, "pool.tsv"), ["run", "pool"])}
eligible, counts, _ = SC.accounting(drawn, ref)
runs_ok = [r for r in pool_of if r in eligible and (r, "auto") in tab and (r, "reference") in tab
           and "missing" not in (tab[(r, "auto")]["status"], tab[(r, "reference")]["status"])]
# the share of drawn runs lost before the endpoints (missing, organism not verified, no proteome): from the PRIMARY
# data, at least 10 %; the 5 runs tagged top-down are removed apart
MISSING = max(0.10, 1 - len(eligible) / max(1, counts["drawn"] - counts["top-down"]))
pools = {}
for r in runs_ok:
    pools.setdefault(pool_of[r], []).append(r)
with_defaults = [r for r in runs_ok if (r, "defaults") in tab and tab[(r, "defaults")]["status"] != "missing"]


def stratum_pool(s):
    if s.startswith("protease:"):
        return "nontryptic"
    return {"acquisition:Orbitrap DIA": "orbitrap DIA", "instrument:Astral": "astral",
            "instrument:timsTOF DDA": "tof DDA", "instrument:timsTOF DIA": "tof DIA", "instrument:TOF": "tof DDA",
            "instrument:ion trap MS2": "iontrap"}.get(s, "orbitrap DDA")


held = SC.read(os.path.join(HERE, "heldout.tsv"), ["stratum", "experiment"])
held = [h for h in held if "top-down" not in (h["experiment"] or "").lower()]
want = [stratum_pool(h["stratum"]) for h in held]
empty = sorted({p for p in want if not pools.get(p)})
out.append("primary pools: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(pools.items()))
           + (f"; no runs in {', '.join(empty)}: drawn from all runs" if empty else ""))
out.append(f"lost before the endpoints: {MISSING:.3f} of the drawn runs ({counts})")
out.append(f"held-out draw: {len(held)} runs not tagged top-down, by pool: " +
           ", ".join(f"{k} {want.count(k)}" for k in sorted(set(want))))


def simulate(n_drawn, shift=1.0, flip=0.0):
    """one simulated held-out set of n_drawn bottom-up runs (the quota mix of heldout.tsv, scaled)"""
    t, runs = {}, set()
    order = [want[i % len(want)] for i in range(n_drawn)]
    for i, p in enumerate(order):
        if rng.random() < MISSING:
            continue
        src_pool = pools.get(p) or runs_ok
        src = src_pool[rng.integers(len(src_pool))]
        a, b = dict(tab[(src, "auto")]), dict(tab[(src, "reference")])
        if a["status"] == "ok" and shift != 1.0:
            a["psm"] = a["psm"] * shift
        if flip and SC.compliant(a) and rng.random() < flip:
            a["psm_fdp"] = SC.FDP_MAX + 1
        t[(i, "auto")], t[(i, "reference")] = a, b
        if (src, "defaults") in tab:
            t[(i, "defaults")] = dict(tab[(src, "defaults")])
        runs.add(i)
    return t, runs


def primary(n_drawn=None, shift=1.0, flip=0.0, x=None, margin=None, sims=SIMS):
    n_drawn = len(want) if n_drawn is None else n_drawn
    x = SC.P3_ABS if x is None else x
    margin = SC.P3_MARGIN_FRAC if margin is None else margin
    p1 = p2 = p3 = joint = n_el = pairs1 = pairs3 = 0
    p2_known = len(with_defaults) >= 10 and len(with_defaults) == len(runs_ok)
    for _ in range(sims):
        t, runs = simulate(n_drawn, shift, flip)
        n_el += len(runs)
        if len(runs) < SC.MIN_ELIGIBLE:
            continue
        r, cnt = SC.ratios(t, runs, "auto", "reference")
        pairs1 += cnt["pairs"]
        ok1 = SC.quantile(SC.boot(SC.median, np.array(r), rng, B_SIM), 0.025) >= SC.P1_MARGIN
        ok2 = False
        if p2_known:
            r2, _ = SC.ratios(t, runs, "auto", "defaults")
            ok2 = bool(r2) and SC.quantile(SC.boot(SC.median, np.array(r2), rng, B_SIM), 0.025) > 1
        pr, _ = SC.p3_pairs(t, runs)
        pairs3 += len(pr)
        n = len(pr)
        ok3 = n > 0 and SC.p3_relative_lower(pr, rng, B_SIM) >= margin * n and \
            SC.clopper_lower(sum(a for a, _ in pr), n) >= x
        p1, p2, p3 = p1 + ok1, p2 + ok2, p3 + ok3
        joint += ok1 and ok2 and ok3
    return {"P1": p1 / sims, "P2": p2 / sims if p2_known else None, "P3": p3 / sims,
            "joint": joint / sims if p2_known else None, "eligible": n_el / sims, "P1 pairs": pairs1 / sims,
            "P3 pairs": pairs3 / sims}


def fmt(v):
    return "unknown" if v is None else f"{v:.2f}"


# the PRIMARY effects
pr, _ = SC.p3_pairs(tab, set(runs_ok))
k, kr, n = sum(a for a, _ in pr), sum(b for _, b in pr), len(pr)
r1, _ = SC.ratios(tab, set(runs_ok), "auto", "reference")
r2, _ = SC.ratios(tab, set(with_defaults), "auto", "defaults")
out.append(f"primary data: {len(runs_ok)} runs; median PSM ratio auto/reference {SC.median(r1):.3f} (n={len(r1)}), "
           f"auto/defaults " + (f"{SC.median(r2):.3f} (n={len(r2)})" if r2 else "n/a")
           + f"; FDP-compliant auto {k}/{n} = {k / n:.3f}, reference {kr}/{n} = {kr / n:.3f}; one-sided 95 % lower "
           f"bounds auto {SC.clopper_lower(k, n):.3f}, reference {SC.clopper_lower(kr, n):.3f}")
if not GRID:
    for shift, flip in ((1.0, 0.0), (0.98, 0.0), (0.97, 0.0), (1.0, 0.02), (1.0, 0.05)):
        e = primary(shift=shift, flip=flip)
        out.append(f"auto PSMs x{shift}, {flip:.0%} of compliant auto runs made noncompliant: eligible "
                   f"{e['eligible']:.0f}, P1 pairs {e['P1 pairs']:.0f}, P3 pairs {e['P3 pairs']:.0f}; power P1 "
                   f"{fmt(e['P1'])}, P2 {fmt(e['P2'])}, P3 both gates {fmt(e['P3'])}, joint {fmt(e['joint'])}")
else:
    out.append("design grid (joint power of P1, P2 and P3 with both gates; P1 / P2 / P3 in brackets):")
    for nd in (187, 250, 300, 400):
        for x in (0.70, 0.75, 0.80):
            for m_ in (SC.P3_MARGIN_FRAC, SC.Fraction(-1, 10)):
                cells = []
                for shift, flip in ((1.0, 0.0), (0.98, 0.02), (0.97, 0.05)):
                    e = primary(nd, shift, flip, x, m_, sims=max(200, SIMS // 4))
                    cells.append(f"{fmt(e['joint'])} ({fmt(e['P1'])}/{fmt(e['P2'])}/{fmt(e['P3'])})")
                out.append(f"drawn {nd}, X {x:.2f}, margin {float(m_):+.2f}: pilot effect {cells[0]}; "
                           f"x0.98 and 2 % fewer compliant {cells[1]}; x0.97 and 5 % fewer {cells[2]}")
    txt = "\n".join(out) + "\n"
    open(os.path.join(HERE, "power_grid.txt"), "w").write(txt)
    print(txt)
    sys.exit(0)

# checks
dchk = SC.read(os.path.join(C, "dev_checks.tsv"), ["run", "source", "arm", "kind", "check", "verdict"])
SC.validate_checks(dchk, "dev_checks.tsv")
mani = SC.read(os.path.join(HERE, "scenarios.tsv"), ["accession", "file", "kind", "databases"])
built = SC.built_scenarios(mani)
multi = {r for r, k_ in built if k_ == "incomplete"}
n_built = {k_: sum(1 for _, kk in built if kk == k_) for k_ in ("related_close", "related_distant", "unrelated",
                                                                "incomplete")}
n_built["wrong_label"] = n_built["wrong_protease"] = sum(1 for r, k_ in built if k_ == "correct" and r not in multi)
n_built["correct"] = int(round(len(want) * (1 - MISSING)))
gates = [(c, k_) for c in ("pre", "post") for k_ in ("correct", "unrelated", "related_close", "related_distant",
                                                     "incomplete", "wrong_label", "wrong_protease")
         if not (c == "pre" and k_ in ("wrong_label", "wrong_protease"))]
rate = {}
for c, k_ in gates:
    v = [r["verdict"] for r in dchk if r["check"] == c and r["kind"] == k_]
    ass = [x for x in v if x != "not assessed"]
    if ass:
        hit = sum((x == "ok") if k_ == "correct" else (x != "ok") for x in ass)
        rate[(c, k_)] = (hit / len(ass), len(ass) / len(v), hit, len(ass))
m = len(SC.GATES)


def exact_power(n, p, bound):
    """P(the exact one-sided binomial test of share <= bound rejects at ALPHA/m) when the true share is p"""
    kc = next((x for x in range(n + 1) if SC.binom_sf(x, n, bound) <= SC.ALPHA / m), None)
    return 0.0 if kc is None else SC.binom_sf(kc, n, p)


for (c, k_), (p, sa, hit, na) in sorted(rate.items()):
    n = int(round(n_built[k_] * sa))
    best = max((b / 100 for b in range(10, 100, 5) if exact_power(n, p, b / 100) >= 0.8), default=None) \
        if n >= 10 else None
    frozen = SC.CHECK_BOUND.get((c, k_), "experimental" if (c, k_) in SC.EXPERIMENTAL else "reported")
    fp = f", exact power at it {exact_power(n, p, frozen):.3f}" if isinstance(frozen, float) else ""
    out.append(f"{c}:{k_}: development {p:.3f} ({hit}/{na}), held-out n={n}; largest bound with exact power >= 0.80 "
               f"at alpha/{m}: {best}; frozen: {frozen}{fp}")
hits = {g: 0 for g in SC.GATES}
for _ in range(SIMS):
    pv = {}
    for g in SC.GATES:
        p, sa, _, _ = rate[g]
        n = int(round(n_built[g[1]] * sa))
        pv[g] = SC.binom_sf(int(rng.binomial(n, p)), n, SC.CHECK_BOUND[g])
    for g, rej in SC.holm(pv).items():
        hits[g] += rej
se = lambda q: math.sqrt(q * (1 - q) / SIMS)
out.append(f"joint power under Holm at the frozen bounds ({SIMS} simulations, standard error <= "
           f"{max(se(hits[g] / SIMS) for g in hits):.3f}): " +
           ", ".join(f"{c}:{k_} {hits[(c, k_)] / SIMS:.2f}" for c, k_ in sorted(hits)))
txt = "\n".join(out) + "\n"
open(os.path.join(HERE, "power.txt"), "w").write(txt)
print(txt)

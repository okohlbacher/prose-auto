"""Power of the PREREG 0.1.12 criteria, simulated from the development data with score.py's own functions.

    python power.py DEVDIR [sims]   -> power.txt next to this file, and stdout

DEVDIR is written by v0112/dev_inputs.py (in-sample development data):
- the ID corpus, prose-auto 0.1.3 auto against the optimal settings (stands in for the reference arm);
- the 16-dataset arm search with the defaults arm;
- the 0.1.11 check calibration, related pairs split by the rank of their common ancestor.

Primary endpoints. Each simulation draws the 187 held-out runs by quota: a run of a stratum is drawn from the matching
development pool (protease strata: the non-tryptic sets; Orbitrap DIA, Astral, timsTOF DDA/DIA, TOF, ion trap: their
instrument pools; all others: Orbitrap DDA). The 5 runs tagged top-down are removed, and each other run is missing
(download, conversion, organism not verified) with probability 0.10. score.py's rules then apply unchanged: the
120-run minimum, PSM eligibility, the ratio conventions, P3's two gates. Effects are shifted to show the margins:
P1 ratios times 0.98 and 0.97; P3 with 2 or 5 % of the compliant auto runs made noncompliant, and with auto exactly
as compliant as the reference arm. The development auto arm is prose-auto 0.1.3: 8 of its runs fall below 100 PSMs
where the optimal settings identify, which counts as noncompliant in P3.
P2: only 16 development data sets have a defaults arm; its power is reported as unknown.

Checks. Per gate, the held-out scenario count (from the manifest) times the development share assessed; the
development share hit; all gates simulated jointly with Holm at 0.05. A kind is gated at the largest bound (step
0.05, from 0.10) whose single-gate power at alpha/9 is >= 0.80, given at least 10 expected assessed scenarios;
otherwise it is experimental. The development
data for related_distant are 9 pairs (order or class): its power is uncertain."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as SC   # noqa: E402

D = sys.argv[1]
SIMS = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
HERE = os.path.dirname(os.path.abspath(__file__))
B_SIM = 2000
MISSING = 0.10
rng = np.random.default_rng(SC.SEED)
out = []

drawn, ref, res, dec, chk, scn, tim = SC.load(D, os.path.join(D, "drawn.tsv"))
tab = {(r["run"], r["arm"]): r for r in res}
pool_of = {r["run"]: r["pool"] for r in SC.read(os.path.join(D, "dev_pool.tsv"), ["run", "pool"])}
dev_runs = [r for r in pool_of if pool_of[r] != "defaults16" and (r, "auto") in tab and (r, "reference") in tab
            and "missing" not in (tab[(r, "auto")]["status"], tab[(r, "reference")]["status"])]
pools = {}
for r in dev_runs:
    pools.setdefault(pool_of[r], []).append(r)


def stratum_pool(s):
    if s.startswith("protease:"):
        return "nontryptic"
    return {"acquisition:Orbitrap DIA": "orbitrap DIA", "instrument:Astral": "astral",
            "instrument:timsTOF DDA": "tof DDA", "instrument:timsTOF DIA": "tof DIA", "instrument:TOF": "tof DDA",
            "instrument:ion trap MS2": "iontrap"}.get(s, "orbitrap DDA")


held = SC.read(os.path.join(HERE, "heldout.tsv"), ["stratum", "experiment"])
held = [h for h in held if "top-down" not in (h["experiment"] or "").lower()]
want = [stratum_pool(h["stratum"]) for h in held]
out.append("development pools: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(pools.items())))
out.append(f"held-out draw: {len(held)} runs not tagged top-down, by pool: " +
           ", ".join(f"{k} {want.count(k)}" for k in sorted(set(want))))


def simulate(shift=1.0, flip=0.0, as_reference=False):
    """one simulated held-out set: {(run, arm): row}, eligible runs; as_reference: auto's FDP and status are the
    reference arm's (auto exactly as compliant as the reference)"""
    t, runs = {}, set()
    for i, p in enumerate(want):
        if rng.random() < MISSING:
            continue
        src = pools[p][rng.integers(len(pools[p]))]
        a, b = dict(tab[(src, "auto")]), dict(tab[(src, "reference")])
        if as_reference:
            a = dict(b)
        if a["status"] == "ok" and shift != 1.0:
            a["psm"] = a["psm"] * shift
        if flip and SC.compliant(a) and rng.random() < flip:
            a["psm_fdp"] = SC.FDP_MAX + 1
        t[(i, "auto")], t[(i, "reference")] = a, b
        runs.add(i)
    return t, runs


def primary(shift=1.0, flip=0.0, as_reference=False):
    p1 = p3r = p3a = n_el = pairs1 = pairs3 = 0
    for _ in range(SIMS):
        t, runs = simulate(shift, flip, as_reference)
        n_el += len(runs)
        if len(runs) < SC.MIN_ELIGIBLE:
            continue
        r, cnt = SC.ratios(t, runs, "auto", "reference")
        pairs1 += cnt["pairs"]
        p1 += SC.quantile(SC.boot(SC.median, np.array(r), rng, B_SIM), 0.025) >= SC.P1_MARGIN
        pr, _ = SC.p3_pairs(t, runs)
        pairs3 += len(pr)
        dd = np.array(pr, dtype=float)
        diffs = SC.boot(lambda s: s[..., 0].mean(axis=-1) - s[..., 1].mean(axis=-1), dd, rng, B_SIM)
        rel = SC.quantile(diffs, 0.05) >= SC.P3_MARGIN
        p3r += rel
        p3a += rel and SC.clopper_lower(int(dd[:, 0].sum()), len(dd)) >= SC.P3_ABS
    return p1 / SIMS, p3r / SIMS, p3a / SIMS, n_el / SIMS, pairs1 / SIMS, pairs3 / SIMS


for shift, flip in ((1.0, 0.0), (0.98, 0.0), (0.97, 0.0), (1.0, 0.02), (1.0, 0.05)):
    p1, p3r, p3a, ne, n1, n3 = primary(shift, flip)
    out.append(f"auto PSMs x{shift}, {flip:.0%} of compliant auto runs made noncompliant: eligible {ne:.0f}, "
               f"P1 pairs {n1:.0f}, P3 pairs {n3:.0f}; power P1 {p1:.2f}, P3 relative {p3r:.2f}, "
               f"P3 both gates {p3a:.2f}")
p1, p3r, p3a, ne, n1, n3 = primary(as_reference=True)
out.append(f"auto exactly as compliant as the reference arm (its counts and FDP): P3 relative {p3r:.2f}, P3 both gates "
           f"{p3a:.2f}")
out.append("P2: power unknown (16 development data sets with a defaults arm)")

pr, _ = SC.p3_pairs(tab, set(dev_runs))
dd = np.array(pr, dtype=float)
out.append(f"P3 development: auto share {dd[:, 0].mean():.3f}, reference share {dd[:, 1].mean():.3f} (n={len(dd)}); "
           f"one-sided 95 % lower bounds auto {SC.clopper_lower(int(dd[:, 0].sum()), len(dd)):.3f}, reference "
           f"{SC.clopper_lower(int(dd[:, 1].sum()), len(dd)):.3f}; P3_ABS = {SC.P3_ABS}")

# checks
dchk = SC.read(os.path.join(D, "dev_checks.tsv"), ["run", "kind", "check", "verdict"])
mani = SC.read(os.path.join(HERE, "scenarios.tsv"), ["accession", "file", "kind", "databases"])
built = SC.built_scenarios(mani)
multi = {r for r, k in built if k == "incomplete"}
n_built = {k: sum(1 for _, kk in built if kk == k) for k in ("related_close", "related_distant", "unrelated",
                                                             "incomplete")}
n_built["wrong_label"] = n_built["wrong_protease"] = sum(1 for r, k in built if k == "correct" and r not in multi)
n_built["correct"] = int(round(len(want) * (1 - MISSING)))
gates = [(c, k) for c in ("pre", "post") for k in ("correct", "unrelated", "related_close", "related_distant",
                                                   "incomplete", "wrong_label", "wrong_protease")
         if not (c == "pre" and k in ("wrong_label", "wrong_protease"))]
rate = {}
for c, k in gates:
    v = [r["verdict"] for r in dchk if r["check"] == c and r["kind"] == k]
    ass = [x for x in v if x != "not assessed"]
    if ass:
        hit = sum((x == "ok") if k == "correct" else (x != "ok") for x in ass)
        rate[(c, k)] = (hit / len(ass), len(ass) / len(v), hit, len(ass))
m = len(SC.CHECK_BOUND)
for (c, k), (p, sa, hit, na) in sorted(rate.items()):
    n = int(round(n_built[k] * sa))
    ks = rng.binomial(n, p, SIMS)
    single = lambda bnd: np.mean([SC.binom_sf(int(x), n, bnd) <= SC.ALPHA / m for x in ks])
    best = max((b / 100 for b in range(10, 100, 5) if single(b / 100) >= 0.8), default=None) if n >= 10 else None
    frozen = SC.CHECK_BOUND.get((c, k), "experimental" if (c, k) in SC.EXPERIMENTAL else "reported")
    out.append(f"{c}:{k}: development {p:.3f} ({hit}/{na}), held-out n={n}; largest bound with power >= 0.80 at "
               f"alpha/{m}: {best}; frozen: {frozen}")
hits = {g: 0 for g in SC.CHECK_BOUND}
for _ in range(SIMS):
    pv = {}
    for g, bnd in SC.CHECK_BOUND.items():
        p, sa, _, _ = rate[g]
        n = int(round(n_built[g[1]] * sa))
        pv[g] = SC.binom_sf(int(rng.binomial(n, p)), n, bnd)
    for g, rej in SC.holm(pv).items():
        hits[g] += rej
out.append("joint power under Holm at the frozen bounds: " +
           ", ".join(f"{c}:{k} {hits[(c, k)] / SIMS:.2f}" for c, k in sorted(hits)))
txt = "\n".join(out) + "\n"
open(os.path.join(HERE, "power.txt"), "w").write(txt)
print(txt)

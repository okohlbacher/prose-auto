"""Power of the PREREG 0.1.12 criteria, simulated from the development data (score.py's inputs built by
v0112/dev_inputs.py: the ID corpus, 0.1.3 auto against optimal settings; the 16-dataset defaults search; the 0.1.11
check calibration). usage: power.py DEVDIR [sims] -> power.txt next to this file and on stdout.

Every simulation draws the held-out sizes from the development distribution with replacement and applies score.py's
rule (same bootstrap, B reduced to 2000 for speed). Expected held-out sizes: 187 runs minus 5 tagged top-down minus
10 % missing or unverified = 164 pairs; per wrong-database kind the 40 runs minus the excluded ones (specificity:
all 164), times the development share assessed. The check bounds in score.py are the largest (step 0.05) with power
>= 0.80 at the development share."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as SC   # noqa: E402

D = sys.argv[1]
SIMS = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
N_PAIRS = 164
SC.B = 2000
rng = np.random.default_rng(SC.SEED)
out = []

ref, res, chk = (SC.read(os.path.join(D, f"{n}.tsv")) for n in ("reference", "results", "checks"))
runs, _ = SC.eligible(ref)
tab = SC.arm_table(res)


def ratio_power(num, den, margin, strict, shifts=(1.0,)):
    r, _ = SC.ratios(tab, runs, num, den)
    lr = np.log(np.clip(np.array(r), 1e-9, 1e9))
    for s in shifts:
        hits = 0
        for _ in range(SIMS):
            x = rng.choice(lr, N_PAIRS) + np.log(s)
            lo = SC.boot(SC.median_rows, x, 0.025, rng)
            hits += (lo > np.log(margin)) if strict else (lo >= np.log(margin))
        out.append(f"{num}/{den}: development median {np.exp(np.median(lr)):.3f} (n={len(lr)}), true ratios x{s}: "
                   f"power {hits / SIMS:.2f} for lower 95 % bound {'>' if strict else '>='} {margin}")


ratio_power("auto", "reference", SC.P1_MARGIN, False, shifts=(1.0, 0.98, 0.97, 0.965))
ratio_power("auto", "defaults", 1.0, True)

pairs = []
for run in sorted(runs):
    a, b = tab.get((run, "auto")), tab.get((run, "reference"))
    if a and b and (SC._num(a, "psm") or 0) >= SC.MIN_PSM and (SC._num(b, "psm") or 0) >= SC.MIN_PSM:
        pairs.append((SC.compliant(a), SC.compliant(b)))
d = np.array(pairs, dtype=float)
diff = lambda s: s[..., 0].mean(axis=-1) - s[..., 1].mean(axis=-1)
hits = abs95 = 0
for _ in range(SIMS):
    x = d[rng.integers(0, len(d), N_PAIRS)]
    hits += SC.boot(diff, x, 0.05, rng) >= SC.P3_MARGIN
    abs95 += x[:, 0].mean() >= 0.95
out.append(f"P3: development FDP <= {SC.FDP_MAX} % share auto {d[:, 0].mean():.3f}, reference {d[:, 1].mean():.3f} "
           f"(n={len(d)}); power {hits / SIMS:.2f} for lower 95 % bound of the difference >= {SC.P3_MARGIN}; "
           f"the 0.1.11 rule (auto share >= 0.95) passes in {abs95 / SIMS:.2f}")

# the checks: per kind, binomial at the development share, held-out scenario counts
N_KIND = {"correct": N_PAIRS,   # specificity: every held-out run (the auto arm's own checks)
          "unrelated": 38, "related_distant": 38, "related_close": 38, "incomplete": 8,
          "wrong_label": 29, "wrong_protease": 29}
for (chk_, kind), bound in sorted(SC.CHECK_BOUND.items()):
    e = SC.checks(chk)[f"{chk_}:{kind}"]
    p, share_ass = e["share"], e["assessed"] / e["scenarios"]
    n = int(round(N_KIND[kind] * share_ass))
    k = rng.binomial(n, p, SIMS)
    lows = np.array([SC.wilson(int(x), n)[0] for x in k])
    pw = np.mean(lows >= bound)
    best = max((b / 100 for b in range(30, 100, 5) if np.mean(lows >= b / 100) >= 0.8), default=None)
    out.append(f"{chk_}:{kind}: development {p:.3f} ({e['hit']}/{e['assessed']}), held-out n={n}: power {pw:.2f} "
               f"for the one-sided Wilson lower bound >= {bound}; largest bound (step 0.05) with power >= 0.80: {best}")
txt = "\n".join(out) + "\n"
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "power.txt"), "w").write(txt)
print(txt)

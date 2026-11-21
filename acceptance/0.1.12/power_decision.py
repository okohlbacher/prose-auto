"""Decision table for the author (not frozen, not in SHA256SUMS): PREDICTIVE power (the 19 pilot runs resampled in
each simulation, power.py's predictive scenario) for 300 / 400 / 500 drawn runs x P3 relative margin -0.05 / -0.10,
X = 0.75, through score.primary unchanged.

    python power_decision.py ACC [workers]   ACC = acceptance/0.1.12 (power.py, pilot/score, heldout.tsv, refdb_draw.tsv)
        -> power_decision.txt next to this script

Per cell: joint, P1, P2, P3 power, P(no verdict), eligible 5-95 %; the variants with 2 % and 5 % of the compliant
auto runs made noncompliant; the leave-one-pilot-project-out minimum (predictive, 200 simulations each). 400 and 500
scale the manifest's top-down and no-proteome shares (power.design, as the grid). Also the predictive figure of the
adopted design (300, -0.05) at the pilot's lower-bound eligibility. 1,000 simulations at 300, 500 otherwise."""
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction

ACC = os.path.abspath(sys.argv[1])
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 1
sys.path.insert(0, ACC)
import power as PW   # noqa: E402

MARGINS = (Fraction(-1, 20), Fraction(-1, 10))


def main():
    units, counts = PW.pilot_units(os.path.join(ACC, "pilot", "score"))
    loss, lo, k, n = PW.attrition_loss(counts)
    tasks, seed = [], 1000
    for nd in (300, 400, 500):
        m = PW.design(nd)[0]
        sims = 1000 if nd == 300 else 500
        for mg in MARGINS:
            for var, fl in (("base", 0.0), ("2%", 0.02), ("5%", 0.05)):
                tasks.append((f"{nd}|{float(mg)}|{var}", units, m, sims, seed, True, 0.0, 1.0, fl, 0.75, mg))
                seed += 1
            for i, u in enumerate(units):
                tasks.append((f"{nd}|{float(mg)}|loo|{u[0].split('/')[0]}", [v for v in units if v is not u], m, 200,
                              seed, True, 0.0, 1.0, 0.0, 0.75, mg))
                seed += 1
    tasks.append(("300|-0.05|attrition", units, PW.design(300)[0], 1000, seed, True, loss, 1.0, 0.0, 0.75, MARGINS[0]))
    with ProcessPoolExecutor(WORKERS) as ex:
        res = {r["label"]: r for r in ex.map(PW.task, tasks)}
    f = lambda e: (f"joint {e['joint']:.2f} (P1 {e['P1']:.2f} / P2 {e['P2']:.2f} / P3 {e['P3']:.2f}), "
                   f"P(no verdict) {e['no_verdict']:.2f}, eligible {e['eligible5']:.0f}-{e['eligible95']:.0f}")
    out = [f"command: python power_decision.py {' '.join(sys.argv[1:])}",
           f"pilot: {counts}; lower-bound eligible share {lo:.3f} (extra loss {loss:.3f}); X 0.75; predictive model"]
    for nd in (300, 400, 500):
        for mg in MARGINS:
            key = f"{nd}|{float(mg)}"
            loo = [v for kk, v in res.items() if kk.startswith(key + "|loo|")]
            worst = min(loo, key=lambda e: e["joint"])
            out.append(f"drawn {nd}, margin {float(mg):+.2f} ({res[key + '|base']['sims']} simulations): "
                       f"{f(res[key + '|base'])}; 2 % fewer compliant: joint {res[key + '|2%']['joint']:.2f}; "
                       f"5 % fewer: joint {res[key + '|5%']['joint']:.2f}; leave one project out: minimum joint "
                       f"{worst['joint']:.2f} ({worst['label'].split('|')[-1]}, 200 simulations each)")
    out.append(f"adopted design (300, -0.05) at the lower-bound eligibility, predictive: {f(res['300|-0.05|attrition'])} "
               f"({res['300|-0.05|attrition']['sims']} simulations)")
    txt = "\n".join(out) + "\n"
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "power_decision.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()

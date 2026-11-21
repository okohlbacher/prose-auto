"""Top up the held-out draw of PREREG 0.1.12 from 187 to 300 runs (§0, §2.1), seed 20261122, with the strata,
queries, predicates and file rule of draw_sets.py. Every PXD in exclude.txt, heldout.tsv (the 187 runs), wrongdb.tsv,
candidates.tsv, pilot/pilot.tsv and pilot/pilot_candidates.tsv is left out. Only project-level fields and file names
and sizes are read.

Quotas: the 187-run quotas scaled to 300 (largest remainder; ties by stratum name), minus the runs already drawn. A
stratum whose pool runs out passes its unfilled places on, one at a time, to the strata in draw order that still
have candidates (logged). One run per project; a project counts for the first stratum that takes it.

python3 draw_topup.py  ->  topup.tsv, topup_candidates.tsv, draw_topup.log, and the 113 rows appended to heldout.tsv
(only when heldout.tsv is still the 187-run draw: its SHA-256 below; otherwise the script stops)."""
import collections
import hashlib
import math
import random
import re
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import draw_sets as D   # noqa: E402

SEED = 20261122
N_TOTAL = 300
HELDOUT_187 = "c9f398249ada92aa5f8f17a9c264e4d693a2cf276fab94b63743ce30993012b6"
COLS = ["stratum", "accession", "file", "mb", "instruments", "organisms", "experiment", "quantification"]


def quotas(have):
    tot = sum(have.values())
    ex = {k: Fraction(v * N_TOTAL, tot) for k, v in have.items()}
    q = {k: math.floor(v) for k, v in ex.items()}
    for k in sorted(ex, key=lambda k: (-(ex[k] - q[k]), k))[:N_TOTAL - sum(q.values())]:
        q[k] += 1
    return {k: q[k] - have[k] for k in have}


def main():
    held = HERE / "heldout.tsv"
    if hashlib.sha256(held.read_bytes()).hexdigest() != HELDOUT_187:
        sys.exit("heldout.tsv is not the 187-run draw; restore it before drawing the top-up")
    rng = random.Random(SEED)
    lines = []
    log = lambda s: (lines.append(s), print(s, flush=True))
    excl = set()
    for f in ("exclude.txt", "heldout.tsv", "wrongdb.tsv", "candidates.tsv", "pilot/pilot.tsv",
              "pilot/pilot_candidates.tsv"):
        excl |= set(re.findall(r"PXD\d{6}", (HERE / f).read_text()))
    have = collections.Counter(l.split("\t")[0] for l in held.read_text().splitlines()[1:])
    want = quotas(have)
    log(f"seed {SEED}; {len(excl)} excluded accessions; top-up {sum(want.values())} runs to {N_TOTAL}")
    seen, chosen, cand, pools, got = set(excl), [], [], {}, collections.Counter()

    def take(name, p):
        run = D.pick_run(p["accession"], rng)
        cand.append((name, p["accession"], run[0] if run else "-"))
        if not run:
            return False
        seen.add(p["accession"])
        chosen.append(dict(stratum=name, accession=p["accession"], file=run[0], mb=run[1],
                           **{k: "; ".join(map(str, p.get(f) or [])) for k, f in (
                               ("instruments", "instruments"), ("organisms", "organisms"),
                               ("experiment", "experimentTypes"), ("quantification", "quantificationMethods"))}))
        got[name] += 1
        return True

    def fill(name, n):
        pool = pools[name]
        while n and pool:
            p = pool.pop(0)
            if p["accession"] not in seen and take(name, p):
                n -= 1
        return n

    short = 0
    for name, _, query, pred in D.STRATA:
        found = [p for q in ((query,) if isinstance(query, str) else query) for p in D.search(q, log)]
        pool = [p for p in found if p.get("accession") not in seen
                and str(p.get("publicationDate", "")) >= D.MIN_DATE and pred(p)]
        pool = sorted({p["accession"]: p for p in pool}.values(), key=lambda p: p["accession"])
        rng.shuffle(pool)
        pools[name] = pool
        n0 = len(pool)
        left = fill(name, want[name])
        short += left
        log(f"{name}: {got[name]} of {want[name]} (pool {n0})")
    order = [s[0] for s in D.STRATA]
    while short:
        moved = False
        for name in order:
            if short and pools[name] and fill(name, 1) == 0:
                short -= 1
                moved = True
                log(f"unfilled place passed to {name}")
        if not moved:
            log(f"{short} places left unfilled: every pool is empty")
            break
    rows = "".join("\t".join(str(r[k]) for k in COLS) + "\n" for r in chosen)
    (HERE / "topup.tsv").write_text("\t".join(COLS) + "\n" + rows)
    (HERE / "topup_candidates.tsv").write_text("stratum\taccession\tfile\n" + "".join("\t".join(r) + "\n" for r in cand))
    with open(held, "a") as fh:
        fh.write(rows)
    log(f"top-up {len(chosen)} runs, {sum(r['mb'] for r in chosen)} MB; heldout.tsv now "
        f"{len(held.read_text().splitlines()) - 1} runs, sha256 {hashlib.sha256(held.read_bytes()).hexdigest()}")
    (HERE / "draw_topup.log").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    assert sum(quotas({"a": 10, "b": 7}).values()) == N_TOTAL - 17
    main()

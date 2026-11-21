"""Top up the held-out draw of PREREG 0.1.12 from 187 to 300 runs (§0, §2.1), seed 20261122, with the strata,
queries, predicates and file rule of draw_sets.py. Every PXD in exclude.txt, heldout.tsv (the 187 runs), wrongdb.tsv,
candidates.tsv, pilot/pilot.tsv and pilot/pilot_candidates.tsv is left out. Only project-level fields and file names
and sizes are read.

Quotas: the 187-run quotas scaled to 300 (largest remainder; ties by stratum name), minus the runs already drawn. A
stratum whose pool runs out passes its unfilled places on, one at a time, to the strata in draw order that still
have candidates (logged). One run per project; a project counts for the first stratum that takes it.

python3 draw_topup.py  ->  topup.tsv, topup_candidates.tsv, draw_topup.log, and the 113 rows appended to heldout.tsv
(only when heldout.tsv is still the 187-run draw: its SHA-256 below; otherwise the script stops).

python3 draw_topup.py --replay  ->  topup_replay.tsv, draw_topup_replay.log: the same draw against today's PRIDE, with
the first 187 rows of heldout.tsv as the 187-run draw; heldout.tsv is not touched. The draw reads PRIDE's live search
results and file lists, which were not archived, so a replay can only show that today's PRIDE gives the same runs;
the frozen record is the append itself (topup.tsv = rows 188-300 of heldout.tsv).

python3 draw_topup.py --batch K DIR  ->  batchK.tsv, batchK_candidates.tsv, draw_batchK.log, rows appended to
heldout.tsv: the conditional top-up of PREREG §2.1 (amendment 2026-11-22b). DIR holds reference.tsv and conversion.tsv
of every drawn run, after the blind annotation and conversion and before any arm. The batch is drawn only when
score.eligible_before_arms counts fewer than MIN_ELIGIBLE (120) eligible runs, heldout.tsv holds exactly the runs of
the batches before (300 + 50 (K - 1)), and the cap (500 drawn) is not reached. Batch K draws min(50, 500 - drawn) runs
with seed BATCH_SEEDS[K], fixed in advance. Its quotas are the 187-run quotas scaled to the new total, minus the runs
each stratum already has (at least 0); unfilled places pass on as above. Its exclusions are those of the top-up plus
every drawn run and every candidate of the top-up and of earlier batches. An empty batch (every pool exhausted) is
logged, and the study then gets no verdict, as at the cap."""
import glob
import os
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
BATCH, CAP = 50, 500
BATCH_SEEDS = {1: 20261201, 2: 20261202, 3: 20261203, 4: 20261204}   # fixed in PREREG §2.1 before any annotation
HELDOUT_187 = "c9f398249ada92aa5f8f17a9c264e4d693a2cf276fab94b63743ce30993012b6"
COLS = ["stratum", "accession", "file", "mb", "instruments", "organisms", "experiment", "quantification"]


def quotas(have, total=N_TOTAL, current=None):
    """the 187-run quotas (have) scaled to total by largest remainder (ties by name), minus the runs each stratum
    already has (current, default have), at least 0"""
    tot = sum(have.values())
    ex = {k: Fraction(v * total, tot) for k, v in have.items()}
    q = {k: math.floor(v) for k, v in ex.items()}
    for k in sorted(ex, key=lambda k: (-(ex[k] - q[k]), k))[:total - sum(q.values())]:
        q[k] += 1
    cur = have if current is None else current
    return {k: max(0, q[k] - cur.get(k, 0)) for k in have}


def draw(seed, want, excl, log, n_target):
    """the seeded draw of draw_sets.STRATA by quota (want), passing unfilled places on; returns (chosen, candidates)"""
    rng = random.Random(seed)
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

    short = n_target - sum(want.values())   # places the clipped quotas leave open pass on too
    for name, _, query, pred in D.STRATA:
        found = [p for q in ((query,) if isinstance(query, str) else query) for p in D.search(q, log)]
        pool = [p for p in found if p.get("accession") not in seen
                and str(p.get("publicationDate", "")) >= D.MIN_DATE and pred(p)]
        pool = sorted({p["accession"]: p for p in pool}.values(), key=lambda p: p["accession"])
        rng.shuffle(pool)
        pools[name] = pool
        n0 = len(pool)
        short += fill(name, want.get(name, 0))
        log(f"{name}: {got[name]} of {want.get(name, 0)} (pool {n0})")
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
    return chosen, cand


def _rows(chosen):
    return "".join("\t".join(str(r[k]) for k in COLS) + "\n" for r in chosen)


def _excl(*extra):
    excl = set()
    for f in ("exclude.txt", "wrongdb.tsv", "candidates.tsv", "pilot/pilot.tsv", "pilot/pilot_candidates.tsv") + extra:
        excl |= set(re.findall(r"PXD\d{6}", (HERE / f).read_text()))
    return excl


def batch(k, d):
    """batch k of the conditional top-up (PREREG §2.1): only by the eligibility rule, before any arm"""
    sys.path.insert(0, str(HERE))
    import score as SC
    held = HERE / "heldout.tsv"
    h = os.path.dirname(os.path.dirname(os.path.abspath(d)))   # DIR is H/annotation/topup
    fz = os.environ.get("PROSE_AUTO_FREEZE", str(HERE / "FREEZE.md"))
    if (os.path.exists(fz) and "Step 4: the annotation sealed" in open(fz).read()) or \
            glob.glob(f"{h}/runs/*/*/*/*.auto.json"):
        sys.exit(f"batch {k}: refused: the annotation is sealed or an arm has run (the top-up comes before any arm)")
    rows = held.read_text().splitlines()[1:]
    if k not in BATCH_SEEDS or len(rows) != N_TOTAL + BATCH * (k - 1):
        sys.exit(f"batch {k}: heldout.tsv has {len(rows)} runs, expected {N_TOTAL + BATCH * (k - 1)}")
    n_el, counts = SC.eligible_before_arms(d, str(held))
    if n_el >= SC.MIN_ELIGIBLE:
        sys.exit(f"batch {k}: {n_el} eligible runs, {SC.MIN_ELIGIBLE} needed: no top-up")
    if len(rows) >= CAP:
        sys.exit(f"batch {k}: the cap of {CAP} drawn runs is reached: no verdict")
    size = min(BATCH, CAP - len(rows))
    lines = []
    log = lambda s: (lines.append(s), print(s, flush=True))
    log(f"batch {k}: {n_el} eligible of {len(rows)} drawn ({counts}); drawing {size}, seed {BATCH_SEEDS[k]}")
    first = collections.Counter(l.split("\t")[0] for l in rows[:187])
    cur = collections.Counter(l.split("\t")[0] for l in rows)
    want = quotas(first, len(rows) + size, cur)
    prev = ("topup_candidates.tsv",) + tuple(f"batch{j}_candidates.tsv" for j in range(1, k))
    excl = _excl(*prev) | set(re.findall(r"PXD\d{6}", held.read_text()))
    chosen, cand = draw(BATCH_SEEDS[k], want, excl, log, size)
    (HERE / f"batch{k}.tsv").write_text("\t".join(COLS) + "\n" + _rows(chosen))
    (HERE / f"batch{k}_candidates.tsv").write_text("stratum\taccession\tfile\n" + "".join("\t".join(r) + "\n" for r in cand))
    with open(held, "a") as fh:
        fh.write(_rows(chosen))
    if not chosen:
        log(f"batch {k}: every pool is exhausted: no verdict")
    log(f"batch {k}: {len(chosen)} runs; heldout.tsv now {len(held.read_text().splitlines()) - 1} runs, sha256 "
        f"{hashlib.sha256(held.read_bytes()).hexdigest()}")
    (HERE / f"draw_batch{k}.log").write_text("\n".join(lines) + "\n")


def main(replay=False):
    held = HERE / "heldout.tsv"
    first = "".join(held.read_text().splitlines(keepends=True)[:188])
    if hashlib.sha256((first if replay else held.read_text()).encode()).hexdigest() != HELDOUT_187:
        sys.exit("heldout.tsv is not the 187-run draw; restore it before drawing the top-up")
    lines = []
    log = lambda s: (lines.append(s), print(s, flush=True))
    excl = _excl() | set(re.findall(r"PXD\d{6}", first))
    have = collections.Counter(l.split("\t")[0] for l in first.splitlines()[1:])
    want = quotas(have)
    log(f"seed {SEED}; {len(excl)} excluded accessions; top-up {sum(want.values())} runs to {N_TOTAL}")
    chosen, cand = draw(SEED, want, excl, log, sum(want.values()))
    rows = _rows(chosen)
    if replay:
        (HERE / "topup_replay.tsv").write_text("\t".join(COLS) + "\n" + rows)
        (HERE / "draw_topup_replay.log").write_text("\n".join(lines) + "\n")
        return
    (HERE / "topup.tsv").write_text("\t".join(COLS) + "\n" + rows)
    (HERE / "topup_candidates.tsv").write_text("stratum\taccession\tfile\n" + "".join("\t".join(r) + "\n" for r in cand))
    with open(held, "a") as fh:
        fh.write(rows)
    log(f"top-up {len(chosen)} runs, {sum(r['mb'] for r in chosen)} MB; heldout.tsv now "
        f"{len(held.read_text().splitlines()) - 1} runs, sha256 {hashlib.sha256(held.read_bytes()).hexdigest()}")
    (HERE / "draw_topup.log").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    assert sum(quotas({"a": 10, "b": 7}).values()) == N_TOTAL - 17
    if "--batch" in sys.argv:
        i = sys.argv.index("--batch")
        batch(int(sys.argv[i + 1]), sys.argv[i + 2])
    else:
        main(replay="--replay" in sys.argv)

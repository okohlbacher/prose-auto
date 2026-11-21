"""Draw the pilot of PREREG 0.1.12 §0: 20 runs, one per PRIDE project, by quota, seed 20261121, with the predicates and
queries of the held-out draw (../draw_sets.py). Every PXD in ../exclude.txt, ../heldout.tsv, ../wrongdb.tsv and
../candidates.tsv (every project the held-out draw looked at) is left out, so the pilot shares no project with the
held-out sets. Only project-level fields and file names and sizes are read.

python3 draw_pilot.py   ->  pilot.tsv, pilot_candidates.tsv, draw_pilot.log (this folder)
Every future draw excludes pilot.tsv's projects as well (PREREG §0)."""
import random
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import draw_sets as D   # noqa: E402

SEED = 20261121
D.MAX_MB = 3000          # the pilot keeps files small enough to download and search in a day
TARGET = {"protease:Glu-C": 1, "protease:Asp-N": 1, "protease:chymotrypsin": 1, "protease:Lys-C": 1,
          "label:TMTpro": 1, "label:TMT": 1, "label:iTRAQ": 1, "label:SILAC": 1, "enrichment:phospho": 1,
          "instrument:Astral": 1, "instrument:timsTOF DDA": 1, "instrument:timsTOF DIA": 1, "instrument:TOF": 1,
          "instrument:ion trap MS2": 2, "acquisition:Orbitrap DIA": 1, "organism:plant": 1,
          "organism:bacterium": 1, "organism:mouse or rat": 1, "baseline:Orbitrap DDA label-free human": 1}


def main():
    rng = random.Random(SEED)
    lines = []
    log = lambda s: (lines.append(s), print(s, flush=True))
    excl = set()
    for f in ("exclude.txt", "heldout.tsv", "wrongdb.tsv", "candidates.tsv"):
        excl |= set(re.findall(r"PXD\d{6}", (HERE.parent / f).read_text()))
    log(f"seed {SEED}; {len(excl)} excluded accessions (exclude.txt, heldout.tsv, wrongdb.tsv, candidates.tsv)")
    chosen, cand, seen = [], [], set(excl)
    for name, _, query, pred in D.STRATA:
        target = TARGET.get(name, 0)
        if not target:
            continue
        found = [p for q in ((query,) if isinstance(query, str) else query) for p in D.search(q, log)]
        pool = [p for p in found if p.get("accession") not in seen
                and str(p.get("publicationDate", "")) >= D.MIN_DATE and pred(p)]
        pool = sorted({p["accession"]: p for p in pool}.values(), key=lambda p: p["accession"])
        rng.shuffle(pool)
        got = 0
        for p in pool:
            if got >= target:
                break
            run = D.pick_run(p["accession"], rng)
            cand.append((name, p["accession"], run[0] if run else "-"))
            if not run:
                continue
            seen.add(p["accession"])
            chosen.append(dict(stratum=name, accession=p["accession"], file=run[0], mb=run[1],
                               instruments="; ".join(map(str, p.get("instruments") or [])),
                               organisms="; ".join(map(str, p.get("organisms") or [])),
                               experiment="; ".join(map(str, p.get("experimentTypes") or [])),
                               quantification="; ".join(map(str, p.get("quantificationMethods") or []))))
            got += 1
        log(f"{name}: {got} of {target} (pool {len(pool)})")
    cols = ["stratum", "accession", "file", "mb", "instruments", "organisms", "experiment", "quantification"]
    (HERE / "pilot.tsv").write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[k]) for k in cols) + "\n"
                                                                   for r in chosen))
    (HERE / "pilot_candidates.tsv").write_text("stratum\taccession\tfile\n" + "".join("\t".join(r) + "\n" for r in cand))
    log(f"pilot {len(chosen)} runs, {sum(r['mb'] for r in chosen)} MB")
    (HERE / "draw_pilot.log").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

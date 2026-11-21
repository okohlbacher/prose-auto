"""The reference-database manifest of PREREG 0.1.12 §3.1: for every drawn run the organisms, taxon ids and UniProt
reference proteomes that all primary arms (auto, reference, SDRF, defaults, best_psm) search, plus cRAP and the
entrapment twins. Frozen with its SHA-256 after adjudication and before any arm runs.

    python refdb.py RUNS.tsv TAXDIR candidates.tsv [DRAW.tsv] > refdb.tsv

RUNS.tsv: accession, file, organisms (the adjudicated organisms; before adjudication the PRIDE field), optionally
class. The organism names follow scenarios.species_names; TAXDIR and candidates.tsv are resolve_uniprot.py's output
for RUNS.tsv (with --snapshot after adjudication). A run is no-proteome when no species is named, a name is not in
the taxonomy, or a species has no reference proteome of its own; it then has no database and leaves P1-P3 (§2.2).
A species named twice (two names, one taxon) is listed once. DRAW.tsv (heldout.tsv): every drawn run must have a row
in RUNS.tsv, else the script stops (adjudicated.tsv covers every drawn run, §3)."""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenarios import at_rank, load_tax, species_names   # noqa: E402


def main(runs, taxdir, candfile, draw=None):
    names, rank, parent = load_tax(taxdir)
    own = {int(r["target"]): r for r in csv.DictReader(open(candfile), delimiter="\t") if r["level"] == "self"}
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "file", "class", "organisms", "taxids", "proteomes", "note"])
    rows = list(csv.DictReader(open(runs), delimiter="\t"))
    if draw:
        have = {(r["accession"], r["file"]) for r in rows}
        miss = [f"{r['accession']}/{r['file']}" for r in csv.DictReader(open(draw), delimiter="\t")
                if (r["accession"], r["file"]) not in have]
        if miss:
            raise SystemExit(f"{len(miss)} drawn runs have no row in {runs}: {', '.join(miss[:5])}")
    for r in rows:
        sp, unres = [], []
        for n in species_names(r["organisms"]):
            t = names.get(n.lower())
            if t is None:
                unres.append(n)
            elif at_rank(t, rank, parent, "species") is not None or rank.get(t) == "species":
                st = at_rank(t, rank, parent, "species") or t
                if all(st != x for _, x in sp):   # one taxon once (two names of one species)
                    sp.append((names[f"#{t}"], st))
        cls = r.get("class") or "bottom-up"
        no_own = [n for n, t in sp if t not in own]
        why = f"not in the taxonomy: {', '.join(unres)}" if unres else "no species named" if not sp else \
            f"no reference proteome of its own: {', '.join(no_own)}" if no_own else None
        if why and cls == "bottom-up":
            cls = "no-proteome"
        w.writerow([r["accession"], r["file"], cls, "; ".join(n for n, _ in sp), "; ".join(str(t) for _, t in sp),
                    "" if why else "; ".join(own[t]["proteome"] for _, t in sp), why or ""])


if __name__ == "__main__":
    main(*sys.argv[1:5])

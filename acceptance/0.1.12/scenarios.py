"""The wrong-database scenarios of PREREG 0.1.12 §2.2, by a taxonomy rule (replaces the hand-made map of
v0111/dbcal.py). Frozen with its SHA-256 and its output in FREEZE.md before any held-out raw file is downloaded.

    python scenarios.py wrongdb.tsv TAXDIR candidates.tsv [adjudicated.tsv] > scenarios.tsv

After the reference annotation, adjudicated.tsv (accession, file, organisms, class) replaces the organisms of its runs
and excludes a run whose class is not bottom-up; the frozen scripts then rebuild the manifest (resolve_uniprot.py with
the same file first), and the rebuilt manifest is hashed before any arm runs.

TAXDIR holds names.dmp and nodes.dmp (NCBI format; resolve_uniprot.py writes the subset from UniProt); candidates.tsv
(resolve_uniprot.py) lists per species (target) its own reference proteome (level self) and the reference proteomes of
other species at the family, order and class level. Per run of wrongdb.tsv:
- organisms: the organisms field split at ";", parenthesised common names and strain suffixes dropped, the first two
  words kept; a name whose rank is above species (Fungi, Bacteria) or that names a metagenome is dropped;
- excluded (no scenarios): the adjudicated class is not bottom-up; the experiment field says top-down; no species is
  left; a name is not in the taxonomy; or a species has no reference proteome of its own (its correct control would
  be a relative: Dendrobaena veneta);
- correct: the species themselves, each with its own reference proteome;
- related_close: per species, the candidate of another genus in the same family; related_distant: of another family
  in the same order, else the same class; each records the rank of the common ancestor;
- unrelated: Arabidopsis thaliana for animals, fungi and other eukaryotes, Homo sapiens for plants, bacteria and
  archaea (the group of the first species);
- incomplete (two or more species): the first species alone.
Ties: the candidate with the most proteins, then the lowest taxid. A scenario without a candidate is "none" and is
reported, not built. Every built row names the UniProt proteome of each of its species; a built row without one stops
the script."""
import csv
import re
import sys

UNRELATED = {"Metazoa": "Arabidopsis thaliana", "Fungi": "Arabidopsis thaliana", "Eukaryota": "Arabidopsis thaliana",
             "Viridiplantae": "Homo sapiens", "Bacteria": "Homo sapiens", "Archaea": "Homo sapiens"}
GROUP_ORDER = ("Metazoa", "Fungi", "Viridiplantae", "Bacteria", "Archaea", "Eukaryota")


def load_tax(d):
    names, rank, parent = {}, {}, {}
    for line in open(f"{d}/names.dmp"):
        f = [x.strip() for x in line.split("|")]
        if f[3] == "scientific name":
            names[f[1].lower()] = int(f[0])
            names.setdefault(f"#{f[0]}", f[1])
    for line in open(f"{d}/nodes.dmp"):
        f = [x.strip() for x in line.split("|")]
        parent[int(f[0])], rank[int(f[0])] = int(f[1]), f[2]
    return names, rank, parent


def lineage(t, parent):
    out = [t]
    while t in parent and parent[t] != t:
        t = parent[t]
        out.append(t)
    return out


def at_rank(t, rank, parent, r):
    return next((x for x in lineage(t, parent) if rank.get(x) == r), None)


def species_names(field):
    out = []
    for o in (field or "").split(";"):
        o = re.sub(r"\(.*?\)", " ", o)
        w = o.split()
        if not w or "metagenome" in o.lower():
            continue
        out.append(" ".join(w[:2]) if len(w) > 1 else w[0])
    return out


def pick(t, cands, level):
    """the candidate for species t: close = the family level (another genus, same family); distant = the order level
    (another family, same order), else the class level (another order, same class). Ties: most proteins, then the
    lowest taxid. Returns (candidate, rank of the common ancestor) or (None, None)."""
    for lvl in (("family",) if level == "close" else ("order", "class")):
        rows = [c for c in cands if c["target"] == t and c["level"] == lvl]
        if rows:
            return min(rows, key=lambda c: (-c["proteins"], c["taxid"])), lvl
    return None, None


def group(t, names, parent):
    lin = {names.get(f"#{x}") for x in lineage(t, parent)}
    return next((g for g in GROUP_ORDER if g in lin), None)


_ADJ = {}


def organisms_of(r, adjudicated=None):
    """the organisms field of a wrongdb.tsv row, replaced by the adjudicated one when adjudicated.tsv lists the run"""
    if adjudicated and adjudicated not in _ADJ:
        _ADJ[adjudicated] = {(x["accession"], x["file"]): x for x in csv.DictReader(open(adjudicated), delimiter="\t")}
    a = _ADJ.get(adjudicated, {}).get((r["accession"], r["file"]))
    return a["organisms"] if a else r["organisms"]


def main(wrongdb, taxdir, candfile, adjudicated=None):
    names, rank, parent = load_tax(taxdir)
    cands = [{"target": int(r["target"]), "level": r["level"], "species": r["species"], "taxid": int(r["taxid"]),
              "proteins": int(r["proteins"]), "proteome": r["proteome"]}
             for r in csv.DictReader(open(candfile), delimiter="\t")]
    own = {c["target"]: c for c in cands if c["level"] == "self"}
    adj = {}
    if adjudicated:
        organisms_of({"accession": "", "file": "", "organisms": ""}, adjudicated)
        adj = _ADJ[adjudicated]
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "file", "kind", "databases", "taxids", "proteomes", "common_rank", "note"])

    def built(kind, species):
        """a built scenario row: every species needs its own reference proteome (else the frozen manifest is wrong)"""
        miss = [n for n, t, _ in species if t not in own]
        if miss:
            raise SystemExit(f"{acc} {kind}: no reference proteome for {', '.join(miss)}")
        w.writerow([acc, f, kind, "; ".join(n for n, _, _ in species), "; ".join(str(t) for _, t, _ in species),
                    "; ".join(own[t]["proteome"] for _, t, _ in species),
                    "; ".join(x for _, _, x in species if x), ""])

    for r in csv.DictReader(open(wrongdb), delimiter="\t"):
        acc, f = r["accession"], r["file"]
        sp, unres = [], []
        for n in species_names(organisms_of(r, adjudicated)):
            t = names.get(n.lower())
            if t is None:
                unres.append(n)
                continue
            if rank.get(t) not in ("species", "subspecies", "strain", "no rank") and \
                    at_rank(t, rank, parent, "species") is None:
                continue   # a higher taxon (Fungi, Bacteria)
            sp.append((names[f"#{t}"], at_rank(t, rank, parent, "species") or t))
        a = adj.get((acc, f))
        no_own = [n for n, t in sp if t not in own]
        why = f"class {a['class']} (adjudication)" if a and a.get("class") not in (None, "", "bottom-up") else \
            "top-down" if "top-down" in (r.get("experiment") or "").lower() else \
            f"not in the taxonomy: {', '.join(unres)}" if unres else "no species named" if not sp else \
            f"no reference proteome of its own: {', '.join(no_own)} (the correct control is undefined)" if no_own \
            else None
        if why:   # excluded runs are reported and never replaced; a scenario is never built from a partial list
            w.writerow([acc, f, "excluded", "", "", "", "", why])
            continue
        built("correct", [(n, t, "") for n, t in sp])
        for lvl in ("close", "distant"):
            got = [pick(t, cands, lvl) for _, t in sp]
            if any(c is None for c, _ in got):
                w.writerow([acc, f, f"related_{lvl}", "none", "", "", "", "no candidate"])
            else:
                for c, _ in got:   # a candidate is its own reference proteome
                    own.setdefault(c["taxid"], c)
                built(f"related_{lvl}", [(c["species"], c["taxid"], x) for c, x in got])
        u = UNRELATED.get(group(sp[0][1], names, parent))
        ut = names.get(u.lower()) if u else None
        if ut is None:
            w.writerow([acc, f, "unrelated", "none", "", "", "", "unrelated organism not resolved"])
        else:
            built("unrelated", [(u, ut, "")])
        if len(sp) > 1:
            built("incomplete", [(sp[0][0], sp[0][1], "")])


if __name__ == "__main__":
    main(*sys.argv[1:5])

"""The wrong-database scenarios of PREREG 0.1.12 §2.2, by a taxonomy rule (replaces the hand-made map of
v0111/dbcal.py). Frozen with its SHA-256 and its output in FREEZE.md before any held-out raw file is downloaded.

    python scenarios.py wrongdb.tsv TAXDIR candidates.tsv > scenarios.tsv

TAXDIR holds names.dmp and nodes.dmp (NCBI format; resolve_uniprot.py writes the subset from UniProt); candidates.tsv
(resolve_uniprot.py) lists per species (target) its own reference proteome (level self) and the reference proteomes of
other species at the family, order and class level. Per run of wrongdb.tsv:
- organisms: the organisms field split at ";", parenthesised common names and strain suffixes dropped, the first two
  words kept; a name whose rank is above species (Fungi, Bacteria) or that names a metagenome is dropped;
- excluded (no scenarios): the experiment field says top-down, no species is left, or a name is not in the taxonomy;
- correct: the species themselves (the database of the reference annotation);
- related_close: per species, the candidate of another genus in the same family; related_distant: of another family
  in the same order, else the same class; each records the rank of the common ancestor;
- unrelated: Arabidopsis thaliana for animals, fungi and other eukaryotes, Homo sapiens for plants, bacteria and
  archaea (the group of the first species);
- incomplete (two or more species): the first species alone.
Ties: the candidate with the most proteins, then the lowest taxid. A scenario without a candidate is "none" and is
reported, not built."""
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


def main(wrongdb, taxdir, candfile):
    names, rank, parent = load_tax(taxdir)
    cands = [{"target": int(r["target"]), "level": r["level"], "species": r["species"], "taxid": int(r["taxid"]),
              "proteins": int(r["proteins"]), "proteome": r["proteome"]}
             for r in csv.DictReader(open(candfile), delimiter="\t")]
    own = {c["target"]: c for c in cands if c["level"] == "self"}
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "file", "kind", "databases", "taxids", "common_rank", "note"])
    for r in csv.DictReader(open(wrongdb), delimiter="\t"):
        acc, f = r["accession"], r["file"]
        sp, unres = [], []
        for n in species_names(r["organisms"]):
            t = names.get(n.lower())
            if t is None:
                unres.append(n)
                continue
            if rank.get(t) not in ("species", "subspecies", "strain", "no rank") and \
                    at_rank(t, rank, parent, "species") is None:
                continue   # a higher taxon (Fungi, Bacteria)
            sp.append((names[f"#{t}"], at_rank(t, rank, parent, "species") or t))
        why = "top-down" if "top-down" in (r.get("experiment") or "").lower() else \
            f"not in the taxonomy: {', '.join(unres)}" if unres else "no species named" if not sp else None
        if why:   # a name the frozen taxonomy lacks excludes the run: a scenario is never built from a partial list
            w.writerow([acc, f, "excluded", "", "", "", why])
            continue
        w.writerow([acc, f, "correct", "; ".join(n for n, _ in sp), "; ".join(str(t) for _, t in sp), "",
                    "proteomes " + "; ".join(own[t]["proteome"] if t in own else "none (nearest relative, as prose-auto)"
                                             for _, t in sp)])
        for lvl in ("close", "distant"):
            got = [pick(t, cands, lvl) for _, t in sp]
            if any(c is None for c, _ in got):
                w.writerow([acc, f, f"related_{lvl}", "none", "", "", "no candidate"])
            else:
                w.writerow([acc, f, f"related_{lvl}", "; ".join(c["species"] for c, _ in got),
                            "; ".join(str(c["taxid"]) for c, _ in got), "; ".join(x for _, x in got),
                            "proteomes " + "; ".join(c["proteome"] for c, _ in got)])
        u = UNRELATED.get(group(sp[0][1], names, parent))
        w.writerow([acc, f, "unrelated", u or "none", names.get(u.lower(), "") if u else "", "", ""])
        if len(sp) > 1:
            w.writerow([acc, f, "incomplete", sp[0][0], sp[0][1], "", f"of {len(sp)} species"])


if __name__ == "__main__":
    main(*sys.argv[1:4])

"""The inputs of scenarios.py from UniProt (taxonomy and reference-proteome metadata only, no sequences):

    python3 resolve_uniprot.py wrongdb.tsv OUTDIR [adjudicated.tsv]

adjudicated.tsv (accession, file, organisms, class), written after the reference annotation, replaces the organisms
field of wrongdb.tsv for its runs (scenarios.py reads the same file). The unrelated organisms of scenarios.py
(Arabidopsis thaliana, Homo sapiens) are always resolved.

OUTDIR/names.dmp, nodes.dmp   the lineages of every organism named in wrongdb.tsv (NCBI dump format)
OUTDIR/candidates.tsv          per species: its own reference proteome (level self) and the reference proteomes of
                               other species of another genus in its family (family), of another family in its order
                               (order) and of another order in its class (class); one per species (most proteins)
OUTDIR/queries.tsv             every request: URL, UniProt release and date, SHA-256 of the response body
Standard library only."""
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scenarios import UNRELATED, organisms_of, species_names   # noqa: E402

API = "https://rest.uniprot.org"
LOG = []


def get(url):
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                body = r.read()
                LOG.append((url, r.headers.get("X-UniProt-Release", ""), r.headers.get("X-UniProt-Release-Date", ""),
                            hashlib.sha256(body).hexdigest()))
                nxt = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("Link", "") or "")
                return json.loads(body), (nxt.group(1) if nxt else None)
        except OSError:
            if k == 3:
                raise
            time.sleep(10 * (k + 1))


def taxon(name):
    q = urllib.parse.quote(f'scientific:"{name}"')
    js, _ = get(f"{API}/taxonomy/search?query={q}&fields=id,scientific_name,rank,lineage&format=json&size=25")
    hit = [h for h in js["results"] if h["scientificName"].lower() == name.lower()]
    return hit[0] if hit else None


def references(taxid):
    """the reference proteomes under taxid: {species name: (taxid, upid, proteins)}, the largest per species"""
    q = urllib.parse.quote(f"(taxonomy_id:{taxid}) AND (reference:true)")
    url, out = f"{API}/proteomes/search?query={q}&fields=upid,organism,protein_count&format=json&size=500", {}
    while url:
        js, url = get(url)
        for h in js["results"]:
            sp = " ".join(h["taxonomy"]["scientificName"].split()[:2])
            row = (h["taxonomy"]["taxonId"], h["id"], h.get("proteinCount", 0))
            if sp not in out or (row[2], -row[0]) > (out[sp][2], -out[sp][0]):
                out[sp] = row
    return out


def main(wrongdb, out, adjudicated=None):
    os.makedirs(out, exist_ok=True)
    nodes, names, cands, done = {}, {}, [], set()
    rows = list(csv.DictReader(open(wrongdb), delimiter="\t"))
    wanted = [n for r in rows for n in species_names(organisms_of(r, adjudicated))] + sorted(set(UNRELATED.values()))
    for n in wanted:
        if n.lower() in done:
            continue
        done.add(n.lower())
        t = taxon(n)
        if not t:
            print(f"unresolved: {n}", file=sys.stderr)
            continue
        chain = [x for x in t.get("lineage", [])] + [{"taxonId": t["taxonId"], "scientificName": t["scientificName"],
                                                     "rank": t["rank"]}]
        parent = 1
        for x in chain:
            nodes[x["taxonId"]] = (parent if x["taxonId"] != parent else 1, x["rank"])
            names[x["taxonId"]] = x["scientificName"]
            parent = x["taxonId"]
        if t["rank"] != "species":
            continue
        sp = t["scientificName"]
        rank_of = {x["rank"]: x["taxonId"] for x in chain}
        levels = [("self", t["taxonId"], None), ("family", rank_of.get("family"), rank_of.get("genus")),
                  ("order", rank_of.get("order"), rank_of.get("family")),
                  ("class", rank_of.get("class"), rank_of.get("order"))]
        got = {}
        for lvl, inner, outer in levels:
            if inner is None or (lvl == "class" and got.get("order")):
                continue
            pool = references(inner)
            if lvl == "self":
                got[lvl] = {k: v for k, v in pool.items() if k == sp}
            else:
                near = references(outer) if outer else {}
                got[lvl] = {k: v for k, v in pool.items() if k not in near and k != sp}
            for k, (tid, up, n_) in got[lvl].items():
                cands.append((t["taxonId"], lvl, k, tid, up, n_))
    names.setdefault(1, "root")
    nodes.setdefault(1, (1, "no rank"))
    with open(f"{out}/nodes.dmp", "w") as fh:
        fh.writelines(f"{k}\t|\t{p}\t|\t{rk}\t|\n" for k, (p, rk) in sorted(nodes.items()))
    with open(f"{out}/names.dmp", "w") as fh:
        fh.writelines(f"{k}\t|\t{v}\t|\t\t|\tscientific name\t|\n" for k, v in sorted(names.items()))
    with open(f"{out}/candidates.tsv", "w") as fh:
        fh.write("target\tlevel\tspecies\ttaxid\tproteome\tproteins\n")
        fh.writelines("\t".join(map(str, c)) + "\n" for c in cands)
    with open(f"{out}/queries.tsv", "w") as fh:
        fh.write("url\tuniprot_release\trelease_date\tsha256\n")
        fh.writelines("\t".join(x) + "\n" for x in LOG)
    print(f"{len(done)} names, {len(nodes)} taxa, {len(cands)} candidates, {len(LOG)} requests", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])

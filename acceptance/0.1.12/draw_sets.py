"""Draw the prose-auto 0.1.11 held-out and wrong-database sets (PREREG.md) from PRIDE, by quota, with a fixed seed.

Only project-level PRIDE fields (instruments, quantification methods, experiment types, organisms, keyword hits in
the project texts) and file names and sizes are read; no SDRF, no result file, no spectrum. Every PXD named in
exclude.txt is left out (the release, ID-corpus, mass-truth, mixture, non-tryptic, HLA, intact-protein and FASTag
acceptance sets). One run per project. Strata are filled in the order below, rarest first; a project counts for
the first stratum that takes it. The chosen files go into the ABI PRIDE mirror (pxd get) only after the freeze.

python3 draw_sets.py exclude.txt   ->  heldout.tsv, wrongdb.tsv, candidates.tsv, draw.log (this folder)
"""
import json, random, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

SEED = 20261120
API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
HERE = Path(__file__).resolve().parent
PAGES = 6            # 100 projects per page and query
MIN_DATE = "2014-01-01"
RAW_EXT = (".raw", ".d.zip", ".d.tar", ".d.tar.gz", ".wiff")   # Thermo, Bruker, SCIEX vendor files
MIN_MB, MAX_MB = 50, 6000

L = lambda p: " ".join(str(x) for x in p.get("instruments") or []).lower()
TXT = lambda p: " ".join(str(p.get(k) or "") for k in ("title", "projectDescription", "sampleProcessingProtocol",
                                                          "dataProcessingProtocol", "keywords")).lower()
QM = lambda p: " ".join(str(x) for x in p.get("quantificationMethods") or []).lower()
EXP = lambda p: " ".join(str(x) for x in p.get("experimentTypes") or []).lower()
ORG = lambda p: " ".join(str(x) for x in p.get("organisms") or []).lower()
DIA = lambda p: "data-independent" in EXP(p) or re.search(r"\bdia\b|swath|diapasef", TXT(p)) is not None
NOTRYP = lambda p: not re.search(r"(?<!chymo)trypsin|(?<!chymo)tryptic", TXT(p))   # chymotrypsin is no trypsin
ORBI = ("orbitrap", "q exactive", "exploris", "lumos", "eclipse", "fusion", "ascend")
ASTRAL = lambda p: "astral" in L(p)
TIMS = lambda p: "timstof" in L(p)
TOF = lambda p: any(k in L(p) for k in ("tripletof", "zenotof", "synapt", "xevo", "6600", "5600", "impact", "maxis"))
IT = lambda p: any(k in L(p) for k in ("ltq orbitrap xl", "ltq orbitrap velos", "orbitrap velos", "ltq velos",
                                        "orbitrap elite", "ltq orbitrap", "ltq")) and "astral" not in L(p)
ORBITRAP = lambda p: any(k in L(p) for k in ORBI) and not ASTRAL(p)
SINGLE = lambda p: len(p.get("organisms") or []) == 1
PLANT = ("arabidopsis", "oryza", "zea", "triticum", "glycine", "solanum", "nicotiana", "brassica", "hordeum",
         "medicago", "vitis", "populus", "chlamydomonas", "physcomitr")
FUNGUS = ("saccharomyces", "candida", "schizosaccharomyces", "aspergillus", "neurospora", "cryptococcus",
          "pichia", "komagataella", "fusarium", "yarrowia")
BACT = ("escherichia", "bacillus", "pseudomonas", "mycobacterium", "staphylococcus", "streptococcus", "salmonella",
        "listeria", "clostridi", "vibrio", "streptomyces", "bacteroides", "corynebacterium", "lactobacillus",
        "synechocystis", "klebsiella", "acinetobacter", "helicobacter", "campylobacter")
VERT = ("danio", "gallus", "xenopus", "bos taurus", "sus scrofa", "ovis", "equus", "canis", "felis", "oryctolagus",
        "macaca", "salmo", "oncorhynchus")
INVERT = ("drosophila", "caenorhabditis", "apis", "bombyx", "aedes", "anopheles", "crassostrea", "mytilus")
has = lambda p, names: any(n in ORG(p) for n in names)

# (stratum, target runs, PRIDE keyword query or queries, predicate)
STRATA = [
    ("protease:Glu-C", 10, "Glu-C", lambda p: re.search(r"glu-?c\b", TXT(p)) and NOTRYP(p)),
    ("protease:Asp-N", 10, "Asp-N", lambda p: re.search(r"asp-?n\b", TXT(p)) and NOTRYP(p)),
    ("protease:chymotrypsin", 10, "chymotrypsin", lambda p: "chymotrypsin" in TXT(p) and NOTRYP(p)),
    ("protease:Lys-C", 10, ("Lys-C", "LysC", "endoproteinase Lys-C"), lambda p: re.search(r"lys-?c\b", TXT(p)) and NOTRYP(p)),
    ("label:TMTpro SPS-MS3", 6, ("TMTpro SPS", "SPS-MS3", "MS3"), lambda p: "tmtpro" in TXT(p) and re.search(r"\bsps\b|ms3", TXT(p))),
    ("label:TMTpro", 8, "TMTpro", lambda p: "tmtpro" in TXT(p)),
    ("label:TMT", 10, "TMT", lambda p: "tmt" in QM(p) + TXT(p) and "tmtpro" not in TXT(p)),
    ("label:iTRAQ", 8, "iTRAQ", lambda p: "itraq" in QM(p) + TXT(p)),
    ("label:SILAC", 8, "SILAC", lambda p: "silac" in QM(p) + TXT(p)),
    ("enrichment:glyco", 6, "glycoproteomics", lambda p: re.search(r"glycopeptide|glycoproteom", TXT(p))),
    ("enrichment:phospho", 8, "phosphoproteomics", lambda p: re.search(r"phosphopeptide|phosphoproteom|tio2|imac", TXT(p))),
    ("instrument:Astral", 8, "Astral", ASTRAL),
    ("instrument:timsTOF DDA", 8, "timsTOF", lambda p: TIMS(p) and not DIA(p)),
    ("instrument:timsTOF DIA", 8, "diaPASEF", lambda p: TIMS(p) and DIA(p)),
    ("instrument:TOF", 6, "TripleTOF", TOF),
    ("instrument:ion trap MS2", 8, "LTQ Orbitrap", IT),
    ("acquisition:Orbitrap DIA", 10, "data-independent acquisition", lambda p: ORBITRAP(p) and DIA(p)),
    ("organism:plant", 6, "plant proteome", lambda p: SINGLE(p) and has(p, PLANT)),
    ("organism:fungus", 5, "yeast proteome", lambda p: SINGLE(p) and has(p, FUNGUS)),
    ("organism:bacterium", 6, "bacterial", lambda p: SINGLE(p) and has(p, BACT)),
    ("organism:other vertebrate", 5, "zebrafish proteome", lambda p: SINGLE(p) and has(p, VERT)),
    ("organism:invertebrate", 4, "drosophila proteome", lambda p: SINGLE(p) and has(p, INVERT)),
    ("organism:mouse or rat", 6, "mouse proteome", lambda p: SINGLE(p) and has(p, ("mus musculus", "rattus"))),
    ("baseline:Orbitrap DDA label-free human", 13, "label-free",
     lambda p: ORBITRAP(p) and not DIA(p) and SINGLE(p) and has(p, ("homo sapiens",))),
]
# wrong-database set: 30 single-organism held-out runs (seeded) and 10 runs of multi-species projects
MIX = ("mixed species", 10, "benchmark", lambda p: len({o.split("(")[0].strip() for o in
                                                                        p.get("organisms") or []}) >= 2)


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:   # PRIDE returns sporadic 5xx
            err = e
            time.sleep(10 * (i + 1))
    raise err


def search(q, log):
    out = []
    for page in range(PAGES):
        url = f"{API}/search/projects?keyword={urllib.parse.quote(q)}&pageSize=100&page={page}"
        res = get(url)
        out += res
        if len(res) < 100:
            break
    log(f"query {q!r}: {len(out)} projects")
    return out


def raw_files(acc):
    files, page = [], 0
    while True:
        res = get(f"{API}/projects/{acc}/files?pageSize=500&page={page}")
        files += res
        if len(res) < 500:
            return files
        page += 1


def pick_run(acc, rng):
    """one vendor file of a reasonable size, seeded; None if the project has none"""
    ok = []
    for f in raw_files(acc):
        n, mb = f.get("fileName", ""), int(f.get("fileSizeBytes") or 0) / 1e6
        if n.lower().endswith(RAW_EXT) and MIN_MB <= mb <= MAX_MB and not n.startswith("."):
            ok.append((n, round(mb)))
    return rng.choice(sorted(ok)) if ok else None


def main(exclude_file):
    rng = random.Random(SEED)
    log_lines = []
    log = lambda s: (log_lines.append(s), print(s, flush=True))
    excl = set(re.findall(r"PXD\d{6}", Path(exclude_file).read_text()))
    log(f"seed {SEED}; {len(excl)} excluded accessions")
    chosen, cand_rows, seen = [], [], set(excl)
    for name, target, query, pred in STRATA + [MIX]:
        found = [p for q in ((query,) if isinstance(query, str) else query) for p in search(q, log)]
        pool = [p for p in found if p.get("accession") not in seen and str(p.get("publicationDate", ""))
                >= MIN_DATE and pred(p)]
        pool = sorted({p["accession"]: p for p in pool}.values(), key=lambda p: p["accession"])
        rng.shuffle(pool)
        got = 0
        for p in pool:
            if got >= target:
                break
            run = pick_run(p["accession"], rng)
            cand_rows.append((name, p["accession"], run[0] if run else "-"))
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
    held = [c for c in chosen if c["stratum"] != MIX[0]]
    mix = [c for c in chosen if c["stratum"] == MIX[0]]
    single = [c for c in held if ";" not in c["organisms"]]
    wrong = sorted(rng.sample(single, 30), key=lambda c: c["accession"]) + mix
    cols = ["stratum", "accession", "file", "mb", "instruments", "organisms", "experiment", "quantification"]
    for fn, rows in (("heldout.tsv", held), ("wrongdb.tsv", wrong)):
        with open(HERE / fn, "w") as fh:
            fh.write("\t".join(cols) + "\n" + "".join("\t".join(str(r[k]) for k in cols) + "\n" for r in rows))
    with open(HERE / "candidates.tsv", "w") as fh:
        fh.write("stratum\taccession\tfile\n" + "".join("\t".join(r) + "\n" for r in cand_rows))
    log(f"held-out {len(held)} runs, wrong-database {len(wrong)} runs ({len(mix)} multi-species)")
    (HERE / "draw.log").write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])

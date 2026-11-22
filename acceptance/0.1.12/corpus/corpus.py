"""ProSE identifications on the whole corpus: prose-auto 0.1.3 settings (auto) against the reference settings (optimal)
and the authors' SDRF settings (sdrf). Same ProSE build, entrapment design and protocol as prose_window_sweep/sweep.py.

  corpus.py prep             plan.json: per run the three arms (cli, database), staged inputs, entrapment databases
  corpus.py run W [H]        all searches, W workers (8 threads each), H of them for unspecific-cleavage arms (CORPUS_DONE)
  corpus.py table            corpus_ids/corpus.tsv: per run x arm PSMs, peptides, protein groups, FDP

Arms (only settings change; protocol as the sweep: <= 1 variable mod per peptide, isotope error -2, Percolator):
  auto     prose_cli of the final 0.1.3 trees (*_v9f where they exist: v9d + 9 runs re-decided in v9e; else *_v9d), database = the decided organism(s)
  optimal  the 0.1.3 decisions with label, Cys, enzyme, enrichment and organism replaced by the reference (SDRF corrected
           by score_fix.tsv, evaluate.py) where it differs; HLA: unspecific cleavage; extra decided variable mods dropped;
           windows by 0.1.3's per-class rule on the measured error (= the sweep's best arms); database = true proteome(s)
  sdrf     the authors' SDRF: windows, enzyme, modifications, organism (missing fields: the auto value)
Protein groups: picked-protein FDR (Savitski 2015) over the 1 % FDR PSMs, protein score = best Percolator score,
indistinguishable proteins (same peptide set) merged; entrapment FDP on groups whose proteins are all ENTRAP_.
"""
import collections, glob, json, os, re, shutil, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed

T = "/ceph/ibmi/abi/oliver/prose-auto-20261029"
X = "/scratch/kohlbach/corpus"
OUT = f"{T}/corpus_ids"
sys.path.insert(0, T)
sys.path.insert(0, f"{T}/harness/v9e_code")
import entrap  # noqa: E402

PR = "/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr"
DEPS = "/ceph/ibmi/abi/oliver/opt/prose-check-20261028-deps/lib"
PROTO = ["-Search:modifications:variable_max_per_peptide", "1", "-Search:precursor:isotope_error_min", "-2",
         "-Search:report:top_hits", "10", "-Search:report:isotope_error_convention", "theoretical_minus_observed",
         "-Search:rescore", "true", "-no_summary"]
SETS = [("benchmark", "valset_v9f", f"{T}/valset_all/run.list"), ("benchmark", "borg_ho_v9f", f"{T}/borg_ho_v3/run.list"),
        ("nontryptic", "nontryptic_v9d", f"{T}/nontryptic/run.list"),
        ("nontryptic_dia", "nontryptic_dia_v9d", f"{T}/nontryptic_dia/run.list"),
        ("hla", "hla_v9d", f"{T}/hla_extra.list"), ("mixset", "mixset_v9f", f"{T}/mixset/run.list")]
EXCLUDE = {("PXD000900", "HeLa_ArgC_HCD_3"), ("PXD040407", None), ("PXD051338", None)}   # nontryptic exclusions
VALID = json.load(open(f"{X}/prose_valid.json")) if os.path.exists(f"{X}/prose_valid.json") else {}

# reference value (SDRF / score_fix) -> prose-auto decision vocabulary, used where the decision does not match
LABEL = {"label free": "label free", "tmt": "TMT6/10/11plex", "TMTpro": "TMTpro (16/18plex)",
         "silac": "SILAC (Lys8/Arg10)", "itraq": "iTRAQ4plex"}
CYS = {"carbamid": "Carbamidomethyl (C)", "iodoacet": "Carbamidomethyl (C)", "methylth": "Methylthio (C)",
       "propiona": "Propionamide (C)", "nethylma": "Nethylmaleimide (C)", "none": "none (free Cys)"}
ENZ = {"trypsin": "Trypsin", "lysc": "Lys-C", "argc": "Arg-C", "aspn": "Asp-N", "chymotry": "Chymotrypsin",
       "gluc": "Glu-C", "glutamyl": "Glu-C", "unspecif": "unspecific cleavage", "lysn": "Lys-N"}
ENR = {"phospho": "Phospho (STY)", "none": "none", "acetyl": "Acetyl (K)", "glygly": "GlyGly (K)", "ubiquit": "GlyGly (K)"}


def norm(s):
    return re.sub(r"[^a-z0-9+]", "", str(s).lower())


def stem(s):
    return re.sub(r"\.(raw|mzml|d|d\.zip|d\.tar\.gz|wiff)$", "", os.path.basename(s.rstrip("/")).lower())


def species(org):
    """'Homo sapiens+Escherichia coli' / 'Homo_sapiens + Saccharomyces_cerevisiae' -> list of names"""
    return [re.sub(r"\s*\(.*\)", "", x).replace("_", " ").strip() for x in re.split(r"\s*\+\s*", str(org or "")) if x.strip()]


def jobs():
    import evaluate as ev
    out = []
    for setname, tree, rl in SETS:
        orgs = {}
        for l in open(rl):
            path, org, run, acc = l.rstrip("\n").split("\t")
            orgs[(run.split(":", 1)[0], stem(run.split(":", 1)[-1]))] = (org, path)   # '<acc>[-pure]:<file>'
            orgs.setdefault((acc, stem(run.split(":", 1)[-1])), (org, path))
        for js in sorted(glob.glob(f"{T}/{tree}/out/*/auto/*.auto.json")):
            acc = js.split("/")[-3]; pxd = acc.split("-")[0]; run = os.path.basename(js)[:-len(".auto.json")]
            if (pxd, run) in EXCLUDE or (pxd, None) in EXCLUDE:
                continue
            org, _ = orgs.get((acc, stem(run)), (None, None))
            sd = f"{T}/{tree}/meta/{ev.SDRF.get(pxd, pxd + '.sdrf.tsv')}"
            truth = ev.sdrf_truth(sd, run) if os.path.exists(sd) else {}
            ref = {p: ev.FIX.get((pxd, p), truth.get(p)) for p in ("label", "cys_alkylation", "enzyme", "enrichment")}
            out.append(dict(set=setname, tree=tree, acc=acc, run=run, json=js, true_org=org, ref=ref,
                            sdrf=sd if os.path.exists(sd) else None))
    return out


def correct(setname, dec, ref):
    """decision values of the optimal arm: {param: (value, note)} for the parameters that change"""
    import evaluate as ev
    ch = {}
    v = lambda p: (dec.get(p) or {}).get("value")
    for p, table in (("label", LABEL), ("cys_alkylation", CYS), ("enzyme", ENZ), ("enrichment", ENR)):
        t = ref.get(p)
        if p == "enzyme" and setname == "hla":
            t = "unspecific cleavage"
        if t in (None, "", "none stated"):
            continue
        m = ev.match(p, ev.pred_str(v(p)), t)
        if m is True or m == "partial":
            continue
        key = ev.label_class(t) if p == "label" else norm(t)
        new = next((val for k, val in table.items() if key.startswith(norm(k)) or key == k), None)
        if p == "enzyme" and "+" in t:   # Lys-C+Trypsin and the like: the first named enzyme
            new = table.get(next((k for k in table if norm(t).startswith(k)), ""), new)
        ch[p] = (new, f"{v(p)} -> {new or 'unmapped ' + str(t)}")
    return ch


def optimal_cli(a, setname, ref):
    from prose_auto import prose_params as pp
    from prose_auto.evidence import Decision
    dec = a["decisions"]
    ch = correct(setname, dec, ref)
    D = {}
    for p, d in dec.items():
        D[p] = Decision(parameter=p, value=d.get("value"), confidence=d.get("confidence", 1.0), source=d.get("source", ""),
                        status=d.get("status", "auto"))
    for p, (val, _) in ch.items():
        if val:
            D[p].value, D[p].source, D[p].confidence = val, "reference", 1.0
    if "variable_mods" in D:
        D["variable_mods"].value = []
    for p in ("precursor_tolerance", "fragment_tolerance"):   # the per-class rule on the measured error
        if p in D and D[p].source in ("user", "repository", "text", "llm"):
            D[p].source = "reference-measured"
    return pp.to_cli(pp.to_prose(D)["params"]), {p: n for p, (val, n) in ch.items()}


MOD_ALIAS = {"iodoacetamide derivatized": "Carbamidomethyl", "deamidation": "Deamidated", "met-loss": "Met-loss",
             "lys8": "Label:13C(6)15N(2)", "arg10": "Label:13C(6)15N(4)", "lys4": "Label:2H(4)", "arg6": "Label:13C(6)"}


def sdrf_settings(path, run):
    """authors' settings of this file's SDRF row(s): windows, enzyme, fixed/variable mods, organisms"""
    import csv
    rows = list(csv.reader(open(path, encoding="utf-8", errors="replace"), delimiter="\t"))
    h = [c.lower().strip() for c in rows[0]]
    body = [r + [""] * (len(h) - len(r)) for r in rows[1:]]
    fi = [i for i, c in enumerate(h) if "comment[data file]" in c]
    sel = [r for r in body if fi and stem(r[fi[0]]) == stem(run)] or body
    cols = lambda key: [i for i, c in enumerate(h) if key in c]
    vals = lambda key: [r[i] for r in sel for i in cols(key) if r[i] and r[i].lower() not in ("not available", "not applicable")]
    first = lambda key: (vals(key) or [None])[0]
    fixed, var, notes = [], [], []
    for m in dict.fromkeys(vals("comment[modification parameters]")):
        kv = dict(x.split("=", 1) for x in m.split(";") if "=" in x)
        name = (kv.get("NT") or "").strip(); mt = (kv.get("MT") or "variable").lower(); pp_ = (kv.get("PP") or "").lower()
        name = MOD_ALIAS.get(name.lower(), name)
        ta = re.sub(r"(?i)protein|[nc]-term", "", kv.get("TA") or "")
        sites = re.sub(r"[^A-Z]", "", ta.upper())
        if "protein n-term" in pp_:
            cands = [f"{name} (Protein N-term)"] + [f"{name} (Protein N-term {s})" for s in sites]
        elif "n-term" in pp_:
            cands = [f"{name} (N-term)"] if not sites else [f"{name} (N-term {s})" for s in sites]
        elif "c-term" in pp_:
            cands = [f"{name} (C-term)"] if not sites else [f"{name} (C-term {s})" for s in sites]
        else:
            cands = [f"{name} ({s})" for s in sites] or [f"{name} (N-term)"]
        if mt.startswith("fixed") and "protein n-term" in pp_:   # ProSE refuses fixed protein-terminal mods
            mt = "variable"; notes.append(f"fixed protein N-term mod as variable: {name}")
        ok = {x.lower(): x for x in VALID.get("fixed" if mt.startswith("fixed") else "variable") or []}
        for c in cands:   # exact (any case), else the N-terminal form of the residue (Gln->pyro-Glu (Q) -> (N-term Q))
            alt = re.sub(r" \(([A-Z])\)$", r" (N-term \1)", c)
            hit = ok.get(c.lower()) or ok.get(alt.lower())
            if hit:
                (fixed if mt.startswith("fixed") else var).append(hit)
            else:
                notes.append(f"mod not in ProSE: {c}")
    enz = first("comment[cleavage agent details]")
    enz = (re.search(r"NT=([^;]+)", enz) or [None, enz])[1] if enz else None
    if enz:
        e = norm(enz)
        enz = next((val for k, val in ENZ.items() if e.startswith(k)), None) or enz
    orgs = list(dict.fromkeys(v for v in vals("characteristics[organism]")))
    return dict(prec=first("comment[precursor mass tolerance]"), frag=first("comment[fragment mass tolerance]"), enzyme=enz,
                fixed=list(dict.fromkeys(fixed)), variable=list(dict.fromkeys(var)), organisms=orgs, notes=notes)


def setopt(cli, k, vals):
    """replace option k (a list value takes all following non-option tokens) or append it"""
    cli = list(cli)
    vals = [str(x) for x in (vals if isinstance(vals, list) else [vals])]
    if k in cli:
        i = cli.index(k); j = i + 1
        while j < len(cli) and not cli[j].startswith("-Search:"):
            j += 1
        cli[i + 1:j] = vals
    else:
        cli += [k, *vals]
    return cli


def tol(v):   # as sweep.tol: (value, 'Da' | 'ppm'), None for missing or zero
    m = re.search(r"([\d.]+)\s*(ppm|da|mmu|th)?", str(v or ""), re.I)
    if not m or float(m.group(1)) == 0:
        return None
    x, u = float(m.group(1)), (m.group(2) or "ppm").lower()
    return (x / 1000, "Da") if u == "mmu" else (x, "Da" if u in ("da", "th") else "ppm")


def sdrf_cli(auto_cli, s):
    cli, notes = list(auto_cli), list(s["notes"])
    for key, opts in (("prec", ("-Search:precursor:mass_tolerance_lower", "-Search:precursor:mass_tolerance_upper",
                                "-Search:precursor:mass_tolerance_unit")),
                      ("frag", ("-Search:fragment:mass_tolerance", "-Search:fragment:mass_tolerance_unit"))):
        t = tol(s[key])
        if t:
            for o in opts[:-1]:
                cli = setopt(cli, o, round(t[0], 6))
            cli = setopt(cli, opts[-1], t[1])
        else:
            notes.append(f"{key} window missing: auto")
    if s["enzyme"]:
        e = {"Trypsin": "Trypsin", "Glu-C": "glutamyl endopeptidase"}.get(s["enzyme"], s["enzyme"])
        if e in (VALID.get("enzyme") or [e]):
            cli = setopt(cli, "-Search:enzyme", e)
        else:
            notes.append(f"enzyme not in ProSE: {e}: auto")
    else:
        notes.append("enzyme missing: auto")
    if s["fixed"] or s["variable"]:
        cli = setopt(cli, "-Search:modifications:fixed", s["fixed"] or [])
        cli = setopt(cli, "-Search:modifications:variable", s["variable"] or [])
    else:
        notes.append("mods missing: auto")
    return cli, notes


ALIAS = {"cercopithecus aethiops": "Chlorocebus sabaeus", "physcomitrella patens": "Physcomitrium patens",
         "lactobacillus plantarum": "Lactiplantibacillus plantarum"}   # synonyms of the reference proteomes


def dbfor(names):
    """entrapment database (seed 0) of the concatenated proteomes; None if one is missing"""
    sys.path.insert(0, f"{T}/org/rca")
    import rca
    hits, missing = [], []
    for n in names:
        n = ALIAS.get(" ".join(n.lower().split()[:2]), n)
        h = rca.truth_db(n) or rca.shortlist_db(n)
        (hits.append(h) if h else missing.append(n))
    if not hits:
        return None, f"no proteome for {'+'.join(missing)}"
    if missing:   # a mixture: search the species that have a proteome
        print("  missing proteome:", missing, flush=True)
    key = "+".join(sorted({re.sub(r"\W+", "_", h[0]) for h in hits}))
    ent = f"{X}/db/{key}_entrap_s0.fasta"
    if not os.path.exists(ent):
        os.makedirs(f"{X}/db", exist_ok=True)
        td, seen = ent.replace("_entrap_s0", "_td"), set()
        with open(td + ".tmp", "w") as fh:
            for h in hits:
                for hd, s in entrap.read_fasta(h[1]):
                    acc = hd.split()[0]
                    if hd.startswith("DECOY_") or acc in seen:
                        continue
                    seen.add(acc)
                    fh.write(f">{hd}\n{s}\n>DECOY_{hd}\n{s[::-1]}\n")
        os.replace(td + ".tmp", td)
        entrap.makedb(td, ent + ".tmp", 0)
        os.replace(ent + ".tmp", ent)
    return ent, key


def prep():
    plan, skipped = [], []
    for j in jobs():
        a = json.load(open(j["json"]))
        i = a["input"]
        dia = str((a["decisions"].get("acquisition") or {}).get("value") or "").startswith("DIA")
        # DIA: prose-auto's DIAspeXtract pseudo-spectra; ddaPASEF: the .d (ProSE reads it); else the given file
        src = i["pseudo_spectra"] if dia and i.get("pseudo_spectra") and os.path.exists(i["pseudo_spectra"]) else i["given"]
        if not os.path.exists(src):
            skipped.append((j["acc"], j["run"], f"input missing {src}")); continue
        uid = f"{j['set']}__{j['acc']}__{re.sub(r'[^A-Za-z0-9._-]+', '_', j['run'])}"[:180]
        ext = ".d" if src.rstrip("/").endswith(".d") else ".mzML"
        mz = f"{X}/in/{uid}{ext}"
        if not os.path.exists(mz):
            os.makedirs(f"{X}/in", exist_ok=True)
            (shutil.copytree if ext == ".d" else shutil.copy)(src, mz + ".tmp" + ext); os.replace(mz + ".tmp" + ext, mz)
        dec = a["decisions"]
        arms = {}
        db, key = dbfor(species(dec.get("organism", {}).get("value")))
        arms["auto"] = dict(cli=a["prose_cli"], db=db, dbkey=key, notes=[] if db else [key])
        cli, ch = optimal_cli(a, j["set"], j["ref"])
        db, key = dbfor(species(j["true_org"])) if j["true_org"] else (None, "no reference organism")
        arms["optimal"] = dict(cli=cli, db=db, dbkey=key, changed=ch, notes=[] if db else [key])
        if j["sdrf"]:
            s = sdrf_settings(j["sdrf"], j["run"])
            cli, notes = sdrf_cli(a["prose_cli"], s)
            db, key = dbfor([x for o in s["organisms"] for x in species(o)]) if s["organisms"] else (None, "no SDRF organism")
            arms["sdrf"] = dict(cli=cli, db=db, dbkey=key, notes=notes + ([] if db else [key]), sdrf=s)
        v = lambda p: (dec.get(p) or {}).get("value")
        lowres2 = str(v("ms2_analyzer")) in ("IT", "LIT", "PAULIONTRAP", "RADIALEJECTIONLINEARIONTRAP", "AXIALEJECTIONLINEARIONTRAP")
        ins = str(v("instrument") or "").lower()
        cls = "iontrap" if lowres2 else "astral" if "astral" in ins else \
            "tof" if ("TOF" in str(v("ms2_analyzer")).upper() or "tims" in ins or "tof" in ins) else "orbitrap"
        plan.append(dict(uid=uid, set=j["set"], acc=j["acc"], run=j["run"], cls=cls, acq=v("acquisition"), input=mz,
                         auto_org=v("organism"), true_org=j["true_org"], auto_enzyme=v("enzyme"),
                         ref_enzyme=("unspecific cleavage" if j["set"] == "hla" else j["ref"].get("enzyme")), arms=arms))
        print("prep", uid, cls, {k: (x["dbkey"], bool(x["db"])) for k, x in arms.items()}, flush=True)
    for s in skipped:
        print("SKIPPED", *s, flush=True)
    json.dump(plan, open(f"{X}/plan.json", "w"), indent=1)
    print(len(plan), "runs", sum(len(p["arms"]) for p in plan), "arms", flush=True)


def proteins(path, db):
    """protein groups at 1 % picked-protein FDR from the rank-1 PSMs at q <= 0.01 (ENTRAP_ groups counted apart)"""
    txt = open(path).read()
    prot = dict(re.findall(r'<ProteinHit id="([^"]+)" accession="([^"]+)"', txt))
    best, peps = {}, collections.defaultdict(set)
    for m in entrap.PID.finditer(txt):
        h = re.search(r'<PeptideHit score="([^"]+)" sequence="([^"]+)"[^>]*?protein_refs="([^"]*)"(.*?)</PeptideHit>', m.group(0), re.S)
        if not h or float(h.group(1)) > 0.01:
            continue
        sc = float(entrap.UP(h.group(4), "percolator_score") or "nan")
        for ref in h.group(3).split():
            acc = prot.get(ref)
            if acc:
                best[acc] = max(best.get(acc, -1e300), sc if sc == sc else -1e300)
                peps[acc].add(entrap.strip(h.group(2)))
    # picked: target vs its decoy (DECOY_<acc>); the better one competes
    pairs = {}
    for acc, sc in best.items():
        base, dec = (acc[6:], True) if acc.startswith("DECOY_") else (acc, False)
        if base not in pairs or sc > pairs[base][0]:
            pairs[base] = (sc, dec, acc)
    ranked = sorted(pairs.values(), reverse=True)
    nd = nt = 0; qs = []
    for sc, dec, acc in ranked:
        nd += dec; nt += not dec
        qs.append(nd / max(1, nt))
    for i in range(len(qs) - 2, -1, -1):
        qs[i] = min(qs[i], qs[i + 1])
    acc_t = [acc for (sc, dec, acc), q in zip(ranked, qs) if not dec and q <= 0.01]
    groups = {}
    for acc in acc_t:   # indistinguishable proteins: identical peptide set
        groups.setdefault(frozenset(peps[acc]), []).append(acc)
    ne = sum(all(a.startswith("ENTRAP_") for a in g) for g in groups.values())
    n = len(groups)
    return dict(groups=n, groups_T=n - ne, groups_E=ne, fdp=100 * 2 * ne / n if n else 0.0)   # r ~ 1 at protein level


def one(task):
    p, arm, threads = task
    a = p["arms"][arm]
    o = f"{X}/out/{p['uid']}/{arm}"
    js = f"{o}/eval.json"
    if os.path.exists(js) or not a["db"]:
        return p["uid"], arm, "skip"
    os.makedirs(o, exist_ok=True)
    cli = list(a["cli"])
    if "-Search:precursor:max_charge" in cli:            # 0.1.3 can emit > 20 (PXD062749: 56); ProSE accepts 1..20
        i = cli.index("-Search:precursor:max_charge") + 1
        cli[i] = str(min(int(cli[i]), 20))
    enz = cli[cli.index("-Search:enzyme") + 1] if "-Search:enzyme" in cli else "Trypsin"
    proto = PROTO + ([] if enz == "unspecific cleavage" else ["-Search:peptide:missed_cleavages", "2"])
    idx = f"{o}/out.idXML"
    cmd = [f"{PR}/bin/ProSE", "-in", p["input"], "-database", a["db"], "-out_idxml", idx, *cli, *proto, "-threads", str(threads)]
    tmp = f"/dev/shm/kohlbach-tmp/corpus_{os.getpid()}"
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, OPENMS_DATA_PATH=f"{PR}/share/OpenMS", LD_LIBRARY_PATH=f"{PR}/lib:{DEPS}", OPENMS_TMPDIR=tmp)
    json.dump(cmd, open(f"{o}/cmd.json", "w"))
    with open(f"{o}/log", "w") as log:
        rc = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    shutil.rmtree(tmp, ignore_errors=True)
    if rc:
        return p["uid"], arm, f"FAILED rc={rc}"
    entrap.evaluate(js + ".tmp", a["db"], [idx])
    e = json.load(open(js + ".tmp"))
    e["protein"] = proteins(idx, a["db"])
    json.dump(e, open(js + ".tmp", "w"), indent=1)
    os.replace(js + ".tmp", js)
    os.unlink(idx)
    return p["uid"], arm, "done"


def heavy(p, arm):
    cli = p["arms"][arm]["cli"]
    return "-Search:enzyme" in cli and cli[cli.index("-Search:enzyme") + 1] == "unspecific cleavage"


def run(workers, heavy_workers):
    plan = json.load(open(f"{X}/plan.json"))
    tasks = [(p, a) for p in plan for a in p["arms"] if p["arms"][a]["db"]]
    light = [(p, a, 8) for p, a in tasks if not heavy(p, a)]
    hv = [(p, a, 16) for p, a in tasks if heavy(p, a)]
    print(len(light), "searches,", len(hv), "unspecific-cleavage searches", flush=True)
    # ponytail: two pools; the unspecific ones need up to ~500 GB each, so few at a time
    with ProcessPoolExecutor(workers) as ex, ProcessPoolExecutor(heavy_workers) as eh:
        fs = [ex.submit(one, t) for t in light] + [eh.submit(one, t) for t in hv]
        for f in as_completed(fs):
            try:
                print(*f.result(), flush=True)
            except Exception as e:   # one bad search must not stop the corpus
                print("ERROR", repr(e)[:300], flush=True)
    print("CORPUS_DONE", flush=True)


def table():
    plan = json.load(open(f"{X}/plan.json"))
    os.makedirs(OUT, exist_ok=True)
    cols = ["set", "acc", "run", "class", "acq", "arm", "database", "psm", "psm_fdp", "pep", "pep_fdp", "protein_groups",
            "protein_fdp", "auto_org", "true_org", "auto_enzyme", "ref_enzyme", "changed", "notes"]
    with open(f"{OUT}/corpus.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in plan:
            for arm, a in p["arms"].items():
                js = f"{X}/out/{p['uid']}/{arm}/eval.json"
                e = json.load(open(js)) if os.path.exists(js) else None
                ps, pp_, pr = (e["psm"], e["peptide"], e["protein"]) if e else ({}, {}, {})
                row = [p["set"], p["acc"], p["run"], p["cls"], p["acq"], arm, a["dbkey"], ps.get("accepted", ""),
                       round(ps["combined"], 3) if e else "", pp_.get("accepted", ""), round(pp_["combined"], 3) if e else "",
                       pr.get("groups", ""), round(pr["fdp"], 3) if e else "", p["auto_org"], p["true_org"], p["auto_enzyme"],
                       p["ref_enzyme"], "; ".join(f"{k}: {v}" for k, v in (a.get("changed") or {}).items()),
                       "; ".join(a.get("notes") or [])]
                fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "prep":
        prep()
    elif c == "run":
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 14, int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    elif c == "table":
        table()

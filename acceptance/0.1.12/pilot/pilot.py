"""The pilot of PREREG 0.1.12 §0 on the IBMI nodes: the four primary arms on the 20 pilot runs, scored with score.py
--pilot. Development tooling: the corpus protocol of corpus_ids/corpus.py (ProSE build, entrapment seed 0, counting).

    python pilot.py sheet      per run: the methods texts, an SDRF if PRIDE has one, the raw-data indicators (MS
                               levels, analysers, reporter-ion and SILAC-pair shares): the annotator's sheet, no
                               prose-auto output
    python pilot.py refwin     the reference windows: a wide-window mini-search with the annotated settings (DDA mzML;
                               DIA and timsTOF take the instrument standard), |median| + 4 robust SD, Gaussian fit
    python pilot.py auto W     prose-auto 0.1.12 default on every run (the auto arm), W in parallel
    python pilot.py prep       the arms' ProSE commands and entrapment databases -> plan.json
    python pilot.py search W   every search, W in parallel (8 threads each)
    python pilot.py tree       the output tree of extract.py, then extract.py and score.py --pilot

Paths: P (this folder on Ceph) holds pilot.tsv, reference.tsv (the annotation) and the results; X (node-local
scratch) the inputs and searches. The arms share one database per run: the proteomes of the annotated organisms as
prose-auto resolves them for --organism (its report's search_database), with reversed decoys and the seed-0
entrapment twins."""
import csv, glob, json, os, re, shutil, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

T = "/ceph/ibmi/abi/oliver/prose-auto-20261029"
P = os.path.dirname(os.path.abspath(__file__))
X = "/scratch/kohlbach/pilot0112"
CODE = f"{T}/prose-auto"                       # the 0.1.12 release commit
sys.path.insert(0, CODE)
sys.path.insert(0, f"{T}/corpus_ids")
sys.path.insert(0, T)
PXD = ["python3", "-I", os.path.expanduser("~/.claude/skills/pride/scripts/pxd.py")]
FASTAG = "/ceph/ibmi/abi/oliver/opt/fastag-1.15.3/FASTag/FASTag"
CRAP = f"{T}/search/db/crap.fasta"
ENV = dict(os.environ, PYTHONPATH=CODE, PATH=f"{T}/env/bin:" + os.environ["PATH"],
           PROSE_AUTO_FILECONVERTER="/ceph/ibmi/abi/oliver/spextract/diaweaver-7ed1169/bin/FileConverter",
           PROSE_AUTO_TEXT_CACHE=f"{T}/text_cache_v019", PROSE_AUTO_JEV_CACHE=f"{T}/v0112/jev_cache_pilot",
           TMPDIR="/scratch/kohlbach/tmp", OPENMS_TMPDIR="/dev/shm/kohlbach-tmp",
           OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4", MKL_NUM_THREADS="4")
REPORTERS = {"TMT": [126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 130.134825, 131.138180],
             "TMTpro": [126.127726, 127.124761, 128.128116, 129.131471, 130.134825, 131.138180, 132.141535,
                        133.144890, 134.148245],
             "iTRAQ4": [114.1112, 115.1083, 116.1116, 117.1150], "iTRAQ8": [113.1078, 118.1112, 119.1150, 121.1220]}


def runs():
    return list(csv.DictReader(open(f"{P}/pilot.tsv"), delimiter="\t"))


def key(r):
    return f"{r['accession']}/{r['file']}"


def stem(f):
    return re.sub(r"\.(raw|wiff|d\.zip|d\.tar\.gz|d\.tar|d)$", "", f, flags=re.I)


def local_input(r):
    """the staged input of a run: <X>/in/<acc>/<stem>.mzML, or .d for Bruker"""
    d = f"{X}/in/{r['accession']}"
    s = stem(r["file"])
    for p in (f"{d}/{s}.mzML", f"{d}/{s}.d"):
        if os.path.exists(p):
            return p
    return None


def reference():
    return {r["run"]: r for r in csv.DictReader(open(f"{P}/reference.tsv"), delimiter="\t")}


# ---- the annotator's sheet --------------------------------------------------------------------------------------------

def indicators(path):
    """MS levels, analysers (filter strings), reporter-ion shares in 500 MS2/MS3 spectra (10 ppm), SILAC-like K/R pairs
    are not measured here (the methods decide)"""
    import pyopenms as oms
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    levels, filt, rep, seen = {}, {}, {k: 0 for k in REPORTERS}, 0
    for s in exp:
        lv = s.getMSLevel()
        levels[lv] = levels.get(lv, 0) + 1
        fs = s.getMetaValue("filter string") if s.metaValueExists("filter string") else ""
        if fs:
            filt[str(fs).split(" ")[0]] = filt.get(str(fs).split(" ")[0], 0) + 1
        if lv >= 2 and seen < 500:
            seen += 1
            mz, _ = s.get_peaks()
            for k, ions in REPORTERS.items():
                rep[k] += sum(any(abs(m - x) <= x * 1e-5 for m in mz) for x in ions) >= len(ions) - 1
    return {"ms_levels": levels, "filter_prefix": filt, "reporter_share": {k: round(v / max(1, seen), 3)
                                                                           for k, v in rep.items()}}


def sheet():
    sys.path.insert(0, os.path.join(CODE))
    from prose_auto import repository as repo
    out = []
    for r in runs():
        p = local_input(r)
        txt = {}
        try:
            js = json.load(urlopen(f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{r['accession']}"))
            txt = {k: js.get(k) for k in ("title", "sampleProcessingProtocol", "dataProcessingProtocol",
                                           "instruments", "organisms", "quantificationMethods", "keywords")}
        except Exception as e:
            txt = {"error": repr(e)}
        sdrf = sorted(glob.glob(f"/ceph/ibmi/abi/data/PXD/{r['accession']}/raw/*sdrf*"))
        ind = indicators(p) if p and p.endswith(".mzML") else {"note": "Bruker .d: analysis.tdf"}
        out.append({"run": key(r), "stratum": r["stratum"], "input": p, "pride": txt, "sdrf": sdrf, "indicators": ind})
        print("sheet", key(r), flush=True)
    json.dump(out, open(f"{P}/annotation_sheet.json", "w"), indent=1)


def urlopen(u):
    import urllib.request
    return urllib.request.urlopen(u, timeout=120)


# ---- organism verification (pilot rule, PREREG §0) ----------------------------------------------------------------

TAX = f"{T}/opt/FASTag-tax-ext/share-FASTag-taxonomy"


def verify():
    """FASTag's database-free species ranking (genus rank, the extended index) on every run: the stated organisms are
    verified when each stated genus is among the top 5 genera. A genus the index does not hold, or a run FASTag cannot
    assess (it fails on DIA and Bruker input, ranks nothing on small runs, and an isobaric label it was not told about
    garbles its tags), is left to the reference arm (>= 100 PSMs at 1 % FDR on the stated proteomes; decided in tree).
    -> verify.json"""
    from prose_auto import external as ext
    from scenarios import species_names
    index = {l.split("|")[1].strip().lower() for l in open(f"{TAX}/names.dmp") if "scientific name" in l}
    ref, out = reference(), {}
    for r in runs():
        k, p = key(r), local_input(r)
        genera = [n.split()[0] for n in species_names(ref[k]["organisms"])]
        try:
            ev = ext.fastag_species(FASTAG, p, threads=16, taxonomy_dir=TAX, rank="genus", log=lambda *x: None,
                                    subsample=0.3)
            top = [str(t["name"]).split()[0] for t in (ev[0].stats or {}).get("taxa", [])][:5]
        except Exception as e:
            top = [f"error: {e!r}"[:200]]
        status = {g: ("in top 5" if g in top else "not in the index" if g.lower() not in index else "not in top 5")
                  for g in genera}
        out[k] = {"stated_genera": genera, "fastag_top5": top, "status": status}
        print("verify", k, status, top, flush=True)
    json.dump(out, open(f"{P}/verify.json", "w"), indent=1)


# ---- reference windows -------------------------------------------------------------------------------------------------

STD = {"orbitrap": (10.0, (20.0, "ppm")), "astral": (10.0, (20.0, "ppm")), "tof": (20.0, (40.0, "ppm")),
       "iontrap": (10.0, (0.5, "Da"))}


def refwin_one(r, ref):
    import numpy as np
    from prose_auto import tagsearch as ts
    from prose_auto.hypotheses import fragment_errors
    from prose_auto.detect_tolerance import fit_gauss_uniform
    from prose_auto.spectra import load_mzml
    a = ref[key(r)]
    p = local_input(r)
    cls = a["instrument_class"]
    if not p or not p.endswith(".mzML") or a.get("acquisition", "DDA") != "DDA":
        prec, frag = STD[cls]
        return {"prec_window": prec, "prec_unit": "ppm", "frag_window": frag[0], "frag_unit": frag[1],
                "source": "instrument standard (DIA / timsTOF)"}
    hr2 = cls != "iontrap"
    mods = ts.Mods()
    alk = a["cys"]
    mods.cys = ts.CYS_MODS.get(alk, mods.cys) if alk != "none (free Cys)" else 0.0
    lab = a["label"]
    k_ = next((k for p_, k in (("TMTpro", "TMTpro"), ("TMT", "TMT"), ("iTRAQ8", "iTRAQ8plex"), ("iTRAQ", "iTRAQ4plex"))
               if lab.startswith(p_)), None)
    if k_:
        mods.k = mods.nterm = ts.LABEL_MASS[k_]
    enz = a["protease"] if a["protease"] in ts.ENZYMES else "Trypsin"
    db = json.load(open(glob.glob(f"{P}/out/{r['accession']}/{stem(r['file'])}/auto/*.auto.json")[0]))[
        "search_database"]["fasta"]
    run = load_mzml(p, max_ms2=20000)
    scans = ts.select_scans(run, 5000, "random", 0)
    peaks = [ts.prepare_peaks(s) for s in scans]
    idx = ts.index_for("ref_" + re.sub(r"\W+", "_", "+".join(map(os.path.basename, db))), db[0], mods, enz)
    pw, fw = 50.0, (50.0 if hr2 else 1.0)
    ok = ts.fdr_filter(ts.score_search(scans, idx, pw, fw, hr2, mods, peaks=peaks))[0]
    pe = np.array([x[4] for x in ok])
    fe = fragment_errors(scans, peaks, ok, mods, fw, hr2)
    if pe.size < 20 or fe.size < 100:
        prec, frag = STD[cls]
        return {"prec_window": prec, "prec_unit": "ppm", "frag_window": frag[0], "frag_unit": frag[1],
                "source": f"instrument standard ({pe.size} PSMs in the wide search)"}
    med = float(np.median(pe))
    sd = float(1.4826 * np.median(np.abs(pe - med)))
    g = fit_gauss_uniform(fe, fw)
    fwin = abs(g[0]) + 4 * g[1] if g else STD[cls][1][0]
    return {"prec_window": round(abs(med) + 4 * sd, 3), "prec_unit": "ppm", "frag_window": round(fwin, 4),
            "frag_unit": "ppm" if hr2 else "Da", "source": f"wide search, {pe.size} PSMs"}


def refwin():
    ref = reference()
    out = {}
    for r in runs():
        if ref.get(key(r), {}).get("class") != "bottom-up":
            continue
        try:
            out[key(r)] = refwin_one(r, ref)
        except Exception as e:
            out[key(r)] = {"error": repr(e)[:300]}
        print("refwin", key(r), out[key(r)], flush=True)
    json.dump(out, open(f"{P}/refwin.json", "w"), indent=1)


# ---- the auto arm ------------------------------------------------------------------------------------------------------

def auto_one(r, organisms):
    from scenarios import species_names
    o = f"{P}/out/{r['accession']}/{stem(r['file'])}/auto"
    if glob.glob(f"{o}/*.auto.json"):
        return key(r), "exists"
    os.makedirs(o, exist_ok=True)
    p = local_input(r)
    args = ["python", "-m", "prose_auto.cli", "analyze", p, "--panel", f"{T}/panel", "--fastag", FASTAG,
            "--fastag-taxonomy", f"{T}/opt/FASTag-tax-ext/share-FASTag-taxonomy",
            "--diaspextract", f"{T}/opt/bin/diaspextract", "--threads", "8", "-o", o,
            "--jev", "--sdrf", "off", "--contaminants", CRAP, "--accession", r["accession"]]
    if p.endswith(".d"):
        args += ["--fileconverter", ENV["PROSE_AUTO_FILECONVERTER"]]
    for sp in species_names(organisms):
        args += ["--organism", sp]
    t = time.time()
    with open(f"{o}/run.log", "w") as fh:
        rc = subprocess.run(["/usr/bin/time", "-v", "-o", f"{o}/time.txt"] + args, stdout=fh, stderr=subprocess.STDOUT,
                            cwd=o, env=ENV).returncode
    return key(r), f"rc={rc} {time.time() - t:.0f} s"


def auto(workers):
    sys.path.insert(0, os.path.dirname(P))
    ref = reference()
    with ProcessPoolExecutor(workers) as ex:
        fs = [ex.submit(auto_one, r, ref[key(r)]["organisms"]) for r in runs() if local_input(r)]
        for f in as_completed(fs):
            print("auto", *f.result(), flush=True)


# ---- the arms ----------------------------------------------------------------------------------------------------------

def ref_cli(s, prec, frag):
    """a ProSE command from settings (label, protease, cys, enrichment) and windows; unstated fields: ProSE defaults"""
    from prose_auto import prose_params as pp
    fixed = ([s["cys"]] if s.get("cys") and s["cys"] != "none (free Cys)" else [])
    var = ["Oxidation (M)"]
    lab = s.get("label") or "label free"
    fixed += pp.LABEL_MODS.get(lab, [])
    var += pp.SILAC_VARIABLE.get(lab, [])
    var += [m for m in pp.ENRICHMENT_MODS.get(s.get("enrichment") or "none", []) if m not in var]
    cli = []
    if prec:
        cli += ["-Search:precursor:mass_tolerance_lower", str(prec[0]), "-Search:precursor:mass_tolerance_upper",
                str(prec[0]), "-Search:precursor:mass_tolerance_unit", prec[1]]
    if frag:
        cli += ["-Search:fragment:mass_tolerance", str(frag[0]), "-Search:fragment:mass_tolerance_unit", frag[1]]
    if s.get("protease"):
        cli += ["-Search:enzyme", pp.ENZYME_NAMES.get(s["protease"], s["protease"])]
    return cli + ["-Search:modifications:fixed", *fixed, "-Search:modifications:variable", *var]


def entrap_db(fastas):
    import entrap
    k = re.sub(r"\W+", "_", "+".join(sorted(os.path.basename(f) for f in fastas)))[:150]
    ent = f"{X}/db/{k}_entrap_s0.fasta"
    if not os.path.exists(ent):
        os.makedirs(f"{X}/db", exist_ok=True)
        td, seen = ent.replace("_entrap_s0", "_td"), set()
        with open(td + ".tmp", "w") as fh:
            for f in fastas:
                for hd, sq in entrap.read_fasta(f):
                    acc = hd.split()[0]
                    if hd.startswith("DECOY_") or acc in seen:
                        continue
                    seen.add(acc)
                    fh.write(f">{hd}\n{sq}\n>DECOY_{hd}\n{sq[::-1]}\n")
        os.replace(td + ".tmp", td)
        entrap.makedb(td, ent + ".tmp", 0)
        os.replace(ent + ".tmp", ent)
    return ent


def prep():
    ref, win = reference(), json.load(open(f"{P}/refwin.json"))
    plan = []
    for r in runs():
        a = ref.get(key(r))
        js = glob.glob(f"{P}/out/{r['accession']}/{stem(r['file'])}/auto/*.auto.json")
        if not a or a["class"] != "bottom-up" or not js:
            continue
        rep = json.load(open(js[0]))
        if "no_digest" in rep:
            plan.append(dict(run=key(r), stopped=True))
            continue
        i = rep["input"]
        dia = str((rep["decisions"].get("acquisition") or {}).get("value") or "").startswith("DIA")
        src = i["pseudo_spectra"] if dia and i.get("pseudo_spectra") and os.path.exists(i["pseudo_spectra"]) \
            else local_input(r)
        db = entrap_db(rep["search_database"]["fasta"])
        w = win[key(r)]
        prec, frag = (w["prec_window"], w["prec_unit"]), (w["frag_window"], w["frag_unit"])
        sdrf = json.loads(a["sdrf_settings"]) if a.get("sdrf_settings") else {}
        arms = {"auto": rep["prose_cli"],
                "reference": ref_cli(a, prec, frag),
                "sdrf": ref_cli(sdrf, (sdrf["prec_window"], sdrf["prec_unit"]) if sdrf.get("prec_window") else None,
                                (sdrf["frag_window"], sdrf["frag_unit"]) if sdrf.get("frag_window") else None),
                "defaults": ["-Search:enzyme", "Trypsin/P", "-Search:modifications:fixed", "Carbamidomethyl (C)",
                             "-Search:modifications:variable", "Oxidation (M)"]}
        settings = {"reference": {k: a[k] for k in ("label", "protease", "cys")} | {"prec_window": prec[0],
                    "prec_unit": prec[1], "frag_window": frag[0], "frag_unit": frag[1]},
                    "sdrf": {k: sdrf.get(k) for k in ("label", "protease", "cys", "prec_window", "prec_unit",
                                                     "frag_window", "frag_unit")},
                    "defaults": {"protease": "Trypsin/P", "cys": "Carbamidomethyl (C)", "label": "label free"}}
        plan.append(dict(run=key(r), input=src, db=db, arms=arms, settings=settings))
    json.dump(plan, open(f"{P}/plan.json", "w"), indent=1)
    print(len(plan), "runs", flush=True)


def search_one(p, arm):
    import corpus as C
    import entrap
    o = f"{X}/out/{p['run']}/{arm}"
    js = f"{o}/eval.json"
    if os.path.exists(js):
        return p["run"], arm, "exists"
    os.makedirs(o, exist_ok=True)
    cli = list(p["arms"][arm])
    if "-Search:precursor:max_charge" in cli:
        i = cli.index("-Search:precursor:max_charge") + 1
        cli[i] = str(min(int(cli[i]), 20))
    enz = cli[cli.index("-Search:enzyme") + 1] if "-Search:enzyme" in cli else "Trypsin"
    proto = C.PROTO + ([] if enz == "unspecific cleavage" else ["-Search:peptide:missed_cleavages", "2"])
    idx = f"{o}/out.idXML"
    cmd = [f"{C.PR}/bin/ProSE", "-in", p["input"], "-database", p["db"], "-out_idxml", idx, *cli, *proto,
           "-threads", "8"]
    tmp = f"/dev/shm/kohlbach-tmp/pilot_{os.getpid()}"
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, OPENMS_DATA_PATH=f"{C.PR}/share/OpenMS", LD_LIBRARY_PATH=f"{C.PR}/lib:{C.DEPS}",
               OPENMS_TMPDIR=tmp)
    json.dump(cmd, open(f"{o}/cmd.json", "w"))
    t = time.time()
    with open(f"{o}/log", "w") as log:
        rc = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    shutil.rmtree(tmp, ignore_errors=True)
    if rc:
        return p["run"], arm, f"FAILED rc={rc}"
    entrap.evaluate(js + ".tmp", p["db"], [idx])
    e = json.load(open(js + ".tmp"))
    e["protein"] = C.proteins(idx, p["db"])
    e["seconds"] = round(time.time() - t, 1)
    json.dump(e, open(js + ".tmp", "w"), indent=1)
    os.replace(js + ".tmp", js)
    os.unlink(idx)   # scratch is short; the counts are in eval.json
    return p["run"], arm, "done"


def search(workers):
    plan = [p for p in json.load(open(f"{P}/plan.json")) if not p.get("stopped")]
    with ProcessPoolExecutor(workers) as ex:
        fs = [ex.submit(search_one, p, a) for p in plan for a in p["arms"]]
        for f in as_completed(fs):
            try:
                print("search", *f.result(), flush=True)
            except Exception as e:
                print("ERROR", repr(e)[:300], flush=True)
    print("PILOT_SEARCH_DONE", flush=True)


# ---- the output tree, extraction and the score -------------------------------------------------------------------------

def tree():
    out, ref = f"{P}/tree", reference()
    shutil.rmtree(out, ignore_errors=True)
    plan = {p["run"]: p for p in json.load(open(f"{P}/plan.json"))}
    for r in runs():
        k = key(r)
        base = f"{out}/runs/{r['accession']}/{r['file']}"
        src = f"{P}/out/{r['accession']}/{stem(r['file'])}/auto"
        for arm in ("auto", "reference", "sdrf", "defaults"):
            d = f"{base}/{arm}"
            os.makedirs(d, exist_ok=True)
            if arm == "auto":
                for f in glob.glob(f"{src}/*.auto.json"):
                    shutil.copy(f, d)
            elif k in plan and not plan[k].get("stopped"):
                json.dump(plan[k]["settings"][arm], open(f"{d}/settings.json", "w"))
            e = f"{X}/out/{k}/{arm}/eval.json"
            if os.path.exists(e):
                shutil.copy(e, d)
            else:   # extract.py: a report with no_digest is stopped whatever this says
                open(f"{d}/status", "w").write("missing")
    inp = f"{P}/score"
    os.makedirs(inp, exist_ok=True)
    subprocess.run([f"{T}/env/bin/python", f"{os.path.dirname(P)}/extract.py", out, inp], check=True)
    with open(f"{inp}/reference.tsv", "w") as fh:
        cols = ["run", "class", "organism_verified", "label", "protease", "cys", "reporter_ions", "prec_window",
                "prec_unit", "frag_window", "frag_unit"]
        win = json.load(open(f"{P}/refwin.json"))
        ver = json.load(open(f"{P}/verify.json"))
        fh.write("\t".join(cols) + "\n")
        reported = {key(r) for r in runs() if glob.glob(f"{P}/out/{r['accession']}/{stem(r['file'])}/auto/*.auto.json")}
        for k, a in ref.items():
            if k not in reported:   # no tool reads the input (DIAspeXtract, OpenTIMS): a conversion failure, missing
                continue
            w = win.get(k, {})
            e = f"{X}/out/{k}/reference/eval.json"
            ref_psm = json.load(open(e))["psm"]["accepted"] if os.path.exists(e) else 0
            v_ = ver.get(k, {})
            st, top = v_.get("status", {}), v_.get("fastag_top5", [])
            # FASTag cannot assess a run when it fails (DIA and Bruker input), ranks nothing (small runs), or the run
            # carries an isobaric label it was not told about: then the reference arm decides (>= 100 PSMs)
            assessable = bool(top) and not str(top[0]).startswith("error") and \
                not str(a.get("label") or "").startswith(("TMT", "iTRAQ"))
            ok = bool(st) and all(x == "in top 5" or ((x == "not in the index" or not assessable) and ref_psm >= 100)
                                  for x in st.values())
            row = {**a, "organism_verified": "yes" if ok else "no",
                   **{c: w.get(c, "") for c in ("prec_window", "prec_unit", "frag_window", "frag_unit")}}
            fh.write("\t".join(str(row.get(c, "") or "") for c in cols) + "\n")
    with open(f"{inp}/drawn.tsv", "w") as fh:
        fh.write("accession\tfile\n" + "".join(f"{r['accession']}\t{r['file']}\n" for r in runs()))
    open(f"{inp}/scenarios.tsv", "w").write("accession\tfile\tkind\tdatabases\n")
    with open(f"{inp}/pool.tsv", "w") as fh:
        fh.write("run\tpool\n" + "".join(f"{key(r)}\t{pool_of(r['stratum'])}\n" for r in runs()))
    subprocess.run([f"{T}/env/bin/python", f"{os.path.dirname(P)}/score.py", inp, f"{inp}/drawn.tsv", "--pilot"],
                   check=True)


def pool_of(stratum):
    """the power.py pool of a stratum (power.stratum_pool, repeated here so this script needs no power.py import)"""
    if stratum.startswith("protease:"):
        return "nontryptic"
    return {"acquisition:Orbitrap DIA": "orbitrap DIA", "instrument:Astral": "astral",
            "instrument:timsTOF DDA": "tof DDA", "instrument:timsTOF DIA": "tof DIA", "instrument:TOF": "tof DDA",
            "instrument:ion trap MS2": "iontrap"}.get(stratum, "orbitrap DDA")


if __name__ == "__main__":
    cmd = sys.argv[1]
    sys.path.insert(0, os.path.dirname(P))
    {"sheet": sheet, "verify": verify, "refwin": refwin, "prep": prep, "tree": tree,
     "auto": lambda: auto(int(sys.argv[2])), "search": lambda: search(int(sys.argv[2]))}[cmd]()

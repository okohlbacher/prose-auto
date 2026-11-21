"""The held-out harness of PREREG 0.1.12 §10 and §10.1, frozen with FREEZE.md (§10 step 2) before any raw file is
downloaded. Each stage is resumable: a finished unit (a FASTA, a conversion, an analysis, a search) is not redone.

    python harness.py STAGE [options]

  setup                  the run's Jev cache: a copy of jev_cache_v0111 (§1)
  fasta MANIFEST...      UniProt proteome FASTAs of every proteome id in the manifests (refdb_draw.tsv or refdb.tsv,
                         scenarios.tsv) and of the verification panel -> H/db/<UP>.fasta, H/db/requests.tsv; the
                         frozen reference panel H/panel/<Genus_species>.fasta (§3.1) from the refdb manifest's
                         species (the first manifest)
  get                    pxd get of every drawn and wrong-database file (a .wiff with its .wiff.scan) into the shared mirror (§10 step 3)
  convert [--workers W]  Thermo .raw -> mzML (pxd convert, the mirror's FileConverter recipe); .wiff -> mzML (pxd
                         convert --tool mzpeak-convert); Bruker .d.zip / .d.tar -> H/in/<acc>/<stem>.d (unpacked,
                         zip CRCs checked); then the MS2-count check (§2.2) -> H/annotation/conversion.tsv
  sheet                  the annotators' sheets (§3, blind): PRIDE project texts, SDRF paths, raw-data indicators ->
                         H/annotation/sheets/<acc>.json; no prose-auto output exists or is read
  verify --index I       §3.1 organism verification of run I (a 5,000-MS2 subset, seed 0, searched with the
                         adjudicated settings against the labelled verification database) -> H/annotation/verify/
  windows --index I      §3.3 reference windows of run I -> H/annotation/windows/
  reference              reference.tsv (score.py's schema; drawn runs) from adjudicated.tsv, verify and windows
  topup-count            the conditional top-up's count (§2.1) on the drawn runs; a batch is then drawn with
                         draw_topup.py --batch K H/annotation/topup
  analyze [--workers W]  every prose-auto analyze (§10.1 step 5a, IBMI, Jev): per eligible run auto, spectra_only,
                         texts, auto_never, auto_always, fast_0110; per built scenario its pre-search analysis
  plan                   the searches per run and scenario -> H/plan/<acc>__<file>.json, H/plan/release.txt
  task --index I         array B (BinAC2): every search and check of run I of the release list
  scenario --index I     array C (BinAC2): every scenario search and check of wrong-database run I
  tree                   the frozen output tree of extract.py (H/out) from the task outputs

Paths: T (the frozen tree), H = T/v0112/heldout (inputs, plans, outputs on Ceph; the BinAC2 copy keeps the same
paths inside the container), X = node-local scratch ($TMPDIR on BinAC2). Every stage first checks the SHA-256 of the
executables it runs against FREEZE.md (check_frozen) and stops on a mismatch. Jev: JEV_MODEL must name a pinned model
(never *latest*); the key comes from the environment, which only ssh stdin fills."""
import csv
import glob
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import time
import urllib.request
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
ACC = os.path.dirname(HERE)
T = os.environ.get("PROSE_AUTO_T", "/ceph/ibmi/abi/oliver/prose-auto-20261029")
H = os.environ.get("PROSE_AUTO_H", f"{T}/v0112/heldout")
# a dry run (development runs only) points these elsewhere; the held-out run uses the defaults
RUNS_TSV = os.environ.get("PROSE_AUTO_RUNS", f"{ACC}/heldout.tsv")
WRONGDB_TSV = os.environ.get("PROSE_AUTO_WRONGDB", f"{ACC}/wrongdb.tsv")
FREEZE = os.environ.get("PROSE_AUTO_FREEZE", f"{ACC}/FREEZE.md")
X = os.environ.get("TMPDIR_HELDOUT") or os.environ.get("TMPDIR") or "/scratch/kohlbach/heldout0112"
CODE = os.environ.get("PROSE_AUTO_CODE", f"{T}/prose-auto")
sys.path.insert(0, HERE)
sys.path.insert(0, ACC)
sys.path.insert(0, CODE)   # the frozen prose_auto first: the environment holds an older installed copy
import arms as A   # noqa: E402


PRE_FREEZE_STAGES = ("setup", "fasta")   # step 2 fetches the FASTA files, then writes FREEZE.md


def frozen_prose_auto(stage):
    """every stage (IBMI and BinAC2): the prose_auto this process imports is the frozen tree (CODE), and, once
    FREEZE.md exists, its tree hash equals FREEZE.md's line 'tree <sha256>  <CODE>/prose_auto'"""
    import prose_auto
    where = os.path.realpath(os.path.dirname(prose_auto.__file__))
    if where != os.path.realpath(f"{CODE}/prose_auto"):
        raise SystemExit(f"prose_auto imported from {where}, not the frozen {CODE}/prose_auto")
    if not os.path.exists(FREEZE):
        if stage in PRE_FREEZE_STAGES:
            return
        raise SystemExit(f"no {FREEZE}: stage {stage} runs only after the freeze")
    m = re.search(rf"^tree ([0-9a-f]{{64}})  {re.escape(CODE)}/prose_auto\b", open(FREEZE).read(), re.M)
    import freeze
    if not m or m.group(1) != freeze.tree_hash(f"{CODE}/prose_auto"):
        raise SystemExit(f"{CODE}/prose_auto differs from its FREEZE.md tree hash")

TOOLS = {"python": f"{T}/env/bin/python", "panel": f"{H}/panel", "crap": f"{T}/search/db/crap.fasta",
         "fastag": "/ceph/ibmi/abi/oliver/opt/fastag-1.15.3/FASTag/FASTag",
         "fastag_taxonomy": f"{T}/opt/FASTag-tax-ext/share-FASTag-taxonomy",
         "diaspextract": f"{T}/opt/bin/diaspextract",
         "fileconverter": "/ceph/ibmi/abi/oliver/spextract/diaweaver-7ed1169/bin/FileConverter",
         "prose": "/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr/bin/ProSE"}
TOOLS_0110 = dict(TOOLS, fastag="/ceph/ibmi/abi/oliver/opt/fastag-1.10.0/FASTag/FASTag")
CODE_0110 = f"{T}/harness/v0110c_code"
PROSE_ENV = {"OPENMS_DATA_PATH": "/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr/share/OpenMS",
             "LD_LIBRARY_PATH": "/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr/lib:"
                                "/ceph/ibmi/abi/oliver/opt/prose-check-20261028-deps/lib"}
PXD = ["python3", "-I", os.path.expanduser("~/.claude/skills/pride/scripts/pxd.py")]
MIRROR = "/ceph/ibmi/abi/data/PXD"
# the verification database (§3.1): the stated proteome(s) plus these, and cRAP
VERIFY_PANEL = {"human": "UP000005640", "mouse": "UP000000589", "yeast": "UP000002311", "ecoli": "UP000000625",
                "arabidopsis": "UP000006548"}
SUBSET, SUBSET_SEED = 5000, 0
VERIFY_MIN_PSM, VERIFY_STATED, VERIFY_COMPONENT = 100, 0.80, 0.05
MS2_TOLERANCE = 0     # the mzML's MS2 count must equal the vendor count exactly
LBL = "LBL_"          # label prefix of the verification database's headers
THREADS = int(os.environ.get("HELDOUT_THREADS", "8"))
FROZEN_TOOLS = ("fastag", "diaspextract", "fileconverter", "prose")


# ---- runs, manifests, the freeze ------------------------------------------------------------------------------------

def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def runs():
    return rows(RUNS_TSV)


def all_runs():
    """the drawn runs and the wrong-database runs outside the draw (they are downloaded, converted and annotated too)"""
    seen, out = set(), []
    for r in runs() + rows(WRONGDB_TSV):
        if key(r) not in seen:
            seen.add(key(r))
            out.append(r)
    return out


CONV = f"{H}/annotation/conversion.tsv"


def key(r):
    return f"{r['accession']}/{r['file']}"


def stem(f):
    return re.sub(r"\.(raw|wiff|d\.zip|d\.tar\.gz|d\.tar|d)$", "", f, flags=re.I)


def fname(k):
    return k.replace("/", "__")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def frozen_hashes():
    """{path: sha256} of FREEZE.md (lines '<sha256>  <path>'); a path listed twice with two digests stops (a later
    section never silently supersedes an earlier one)"""
    p = FREEZE
    if not os.path.exists(p):
        return {}
    out = {}
    for m in re.finditer(r"^([0-9a-f]{64})  (\S.*)$", open(p).read(), re.M):
        if out.setdefault(m.group(2), m.group(1)) != m.group(1):
            raise SystemExit(f"FREEZE.md lists {m.group(2)} with two digests")
    return out


def frozen_tree(path):
    """stop unless the tree hash of a directory equals FREEZE.md's 'tree <sha256>  <path>' line"""
    import freeze
    m = re.search(rf"^tree ([0-9a-f]{{64}})  {re.escape(path)}\b", open(FREEZE).read(), re.M)
    if not m or m.group(1) != freeze.tree_hash(path):
        raise SystemExit(f"{path} differs from its FREEZE.md tree hash")


def files_of(path):
    """the files of a spectra input (an mzML, or every file of a .d folder)"""
    if os.path.isdir(path):
        return sorted(os.path.join(d, f) for d, _, fs in os.walk(path) for f in fs)
    return [path]


def child_import_origin(e, code=None):
    """the prose_auto a child 'python -m prose_auto.cli' imports with environment e is the frozen tree"""
    code = code or CODE
    out = subprocess.run([TOOLS["python"], "-c", "import prose_auto, os; print(os.path.realpath(prose_auto.__file__))"],
                         env=e, cwd=X, capture_output=True, text=True).stdout.strip()
    if os.path.dirname(out) != os.path.realpath(f"{code}/prose_auto"):
        raise SystemExit(f"a child imports prose_auto from {out}, not {code}/prose_auto")


# §10.1 reruns: a failed unit (search, check, analysis, verification, conversion) is rerun once, then missing. Attempts
# are counted in <dir>/<unit>.attempts across reruns of the stage or Slurm task (a killed task counts: the attempt is
# written before the run). A stage with a unit left to rerun exits 1, so Slurm records the task failed and the rerun
# (2x memory and time) is due; a unit out of attempts is final (missing, <unit>.failed names the step).
ATTEMPTS = 2
FAILED = []      # units of this process that failed and may be rerun
PUBLISH = []     # paths under H this process wrote (heldout.slurm copies only these back)


def run_unit(d, unit, fn):
    """'ok', 'retry' or 'missing'; fn() returns True on success"""
    os.makedirs(d, exist_ok=True)
    f = f"{d}/{unit}.attempts"
    n = int(open(f).read()) if os.path.exists(f) else 0
    if n >= ATTEMPTS:
        return "missing"
    open(f, "w").write(f"{n + 1}\n")
    try:
        ok = fn()
    except Exception as e:   # noqa: BLE001 (the unit's failure, recorded)
        print(f"{unit} {d}: {e!r}"[:300], flush=True)
        ok = False
    if ok:
        return "ok"
    if n + 1 >= ATTEMPTS:
        open(f"{d}/{unit}.failed", "w").write(f"{unit}\n")
        return "missing"
    FAILED.append(f"{unit} {d}")
    return "retry"


def check_frozen(paths):
    """stop unless every path's SHA-256 equals its FREEZE.md entry (§10.1: a mismatch is a failed task)"""
    fz = frozen_hashes()
    bad = [p for p in paths if fz.get(p) != sha256(p)]
    if bad:
        raise SystemExit(f"not frozen or changed since the freeze: {', '.join(bad)}")


def no_jev():
    """a dry run on development runs may leave Jev out (HELDOUT_NO_JEV=1); never with the real FREEZE.md"""
    if os.environ.get("HELDOUT_NO_JEV") != "1":
        return False
    if FREEZE == f"{ACC}/FREEZE.md":
        raise SystemExit("HELDOUT_NO_JEV is for dry runs only (set PROSE_AUTO_FREEZE elsewhere)")
    return True


JEV_FIELDS = {"JEV_MODEL": "Jev model", "JEV_URL": "Jev endpoint", "PROSE_AUTO_JEV_VOTES": "Jev votes"}


def jev_frozen():
    """the Jev model, endpoint and vote count FREEZE.md records ('- Jev model: X' ...)"""
    txt = open(FREEZE).read() if os.path.exists(FREEZE) else ""
    out = {}
    for var, label in JEV_FIELDS.items():
        m = re.search(rf"^- {label}: (\S+)\s*$", txt, re.M)
        if not m:
            raise SystemExit(f"FREEZE.md records no {label}")
        out[var] = m.group(1)
    if "latest" in out["JEV_MODEL"]:
        raise SystemExit("FREEZE.md's Jev model is not pinned (*latest*)")
    return out


def jev_env():
    """the frozen Jev model, endpoint and votes; an environment that names others stops"""
    fz = jev_frozen()
    for var, v in fz.items():
        if os.environ.get(var) not in (None, "", v):
            raise SystemExit(f"{var}={os.environ[var]} differs from FREEZE.md ({v})")
    if not (os.environ.get("JEV_API_KEY") or os.environ.get("TYPESAFE_API_KEY")):
        raise SystemExit("no Jev key in the environment (it comes over ssh stdin only)")
    return fz


def env(extra=None):
    e = dict(os.environ, PYTHONPATH=CODE, PATH=f"{T}/env/bin:" + os.environ["PATH"],
             PROSE_AUTO_FILECONVERTER=TOOLS["fileconverter"], PROSE_AUTO_TEXT_CACHE=f"{T}/text_cache_v019",
             PROSE_AUTO_JEV_CACHE=f"{H}/jev_cache", OPENMS_TMPDIR=f"{X}/openms",
             OMP_NUM_THREADS=str(THREADS), OPENBLAS_NUM_THREADS=str(THREADS), MKL_NUM_THREADS=str(THREADS))
    os.makedirs(e["OPENMS_TMPDIR"], exist_ok=True)
    return {**e, **(extra or {})}


def freeze_append(title, paths):
    """append a step's hashes to FREEZE.md (freeze.section)"""
    import freeze
    with open(FREEZE, "a") as fh:
        fh.write(freeze.section(title, paths))


JEV_CACHE_FROZEN = f"{T}/jev_cache_v0111"


def setup():
    """the Jev cache of the run: a copy of jev_cache_v0111 (§1); new requests go only into the copy"""
    if not os.path.exists(f"{H}/jev_cache"):
        shutil.copytree(JEV_CACHE_FROZEN, f"{H}/jev_cache")
    print("Jev cache", f"{H}/jev_cache")


# ---- FASTA files and the reference panel (§3.1, §10 step 2 and 4) ---------------------------------------------------

UNIPROT = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3A{up}%29"


def proteome_ids(manifests):
    ids = set(VERIFY_PANEL.values())
    for m in manifests:
        for r in rows(m):
            ids |= {x.strip() for x in (r.get("proteomes") or "").split(";") if x.strip().startswith("UP")}
    return sorted(ids)


def panel_entries(refdb):
    """{Genus_species: UP id} of the refdb manifest's species (one FASTA per species, as prose-auto's --panel)"""
    out = {}
    for r in rows(refdb):
        for n, up in zip((r.get("organisms") or "").split(";"), (r.get("proteomes") or "").split(";")):
            n, up = n.strip(), up.strip()
            if n and up.startswith("UP"):
                out["_".join(n.split()[:2])] = up
    return out


def fasta(manifests):
    os.makedirs(f"{H}/db", exist_ok=True)
    log = f"{H}/db/requests.tsv"
    for up in proteome_ids(manifests):
        dst = f"{H}/db/{up}.fasta"
        if os.path.exists(dst):
            continue
        url = UNIPROT.format(up=up)
        with urllib.request.urlopen(url, timeout=600) as r:
            body, rel = r.read(), r.headers.get("X-UniProt-Release", "")
        if not body.startswith(b">"):
            raise SystemExit(f"{up}: no FASTA from UniProt")
        open(dst + ".part", "wb").write(body)
        os.replace(dst + ".part", dst)
        with open(log, "a") as fh:
            fh.write(f"{url}\t{rel}\t{hashlib.sha256(body).hexdigest()}\n")
        print("fasta", up, len(body), flush=True)
    os.makedirs(f"{H}/panel", exist_ok=True)
    for name, up in panel_entries(manifests[0]).items():
        p = f"{H}/panel/{name}.fasta"
        if not os.path.exists(p):
            shutil.copyfile(f"{H}/db/{up}.fasta", p)


def search_db(proteomes, crap=True):
    """target + reversed decoys + seed-0 entrapment twins of the given UniProt proteomes (and cRAP), once per set
    (pilot.entrap_db); the FASTA files are the frozen H/db copies"""
    sys.path.insert(0, T)
    import entrap
    srcs = [f"{H}/db/{u}.fasta" for u in proteomes] + ([TOOLS["crap"]] if crap else [])
    k = "+".join(sorted(proteomes)) + ("+crap" if crap else "")
    ent = f"{H}/db/search/{k}_entrap_s0.fasta"
    if os.path.exists(ent):
        return ent
    os.makedirs(os.path.dirname(ent), exist_ok=True)
    td, seen = ent.replace("_entrap_s0", "_td"), set()
    with open(td + ".part", "w") as fh:
        for f in srcs:
            for hd, sq in entrap.read_fasta(f):
                acc = hd.split()[0]
                if hd.startswith("DECOY_") or acc in seen:
                    continue
                seen.add(acc)
                fh.write(f">{hd}\n{sq}\n>DECOY_{hd}\n{sq[::-1]}\n")
    os.replace(td + ".part", td)
    entrap.makedb(td, ent + ".part", 0)
    os.replace(ent + ".part", ent)
    return ent


# ---- download, conversion, the MS2-count check (§2.2, §10 step 3) ----------------------------------------------------

def mirror_raw(r):
    return f"{MIRROR}/{r['accession']}/raw/{r['file']}"


def get():
    for r in all_runs():
        files = [r["file"]] + ([r["file"] + ".scan"] if r["file"].lower().endswith(".wiff") else [])
        rc = subprocess.run(PXD + ["get", r["accession"], "--files", *files]).returncode
        print("get", key(r), "ok" if rc == 0 else f"FAILED rc={rc}", flush=True)


def mzml_counts(path):
    """(spectra of every MS level, MS2-or-higher spectra) of an mzML, from its ms level cvParams (MS:1000511)"""
    n = n2 = 0
    with open(path, "rb") as fh:
        for line in fh:
            if b'accession="MS:1000511"' in line:
                m = re.search(rb'value="(\d+)"', line)
                if m:
                    n += 1
                    n2 += int(m.group(1)) >= 2
    return n, n2


def ms2_mzml(path):
    return mzml_counts(path)[1]


def vendor_counts(path, work):
    """(spectra of every MS level, MS2 spectra) of the vendor file itself, read by the frozen FASTag from the .raw
    (Thermo RawFileReader) or the .d (Bruker TDF), never from the converted mzML: one full pass with -qc_out
    (count_ms2) and -census (input_spectra); None for a format FASTag does not read (.wiff)"""
    if path.lower().endswith(".wiff"):
        return None
    os.makedirs(work, exist_ok=True)
    subprocess.run([TOOLS["fastag"], "-in", path, "-out", f"{work}/tags.tsv", "-qc_out", f"{work}/qc.tsv", "-census",
                    "-threads", str(THREADS)], capture_output=True, check=True)
    val = lambda f, sec, k: next((int(float(x[2])) for x in (l.rstrip("\n").split("\t") for l in open(f))
                                  if len(x) > 2 and x[0] == sec and x[1] == k and x[2]), None)
    total = val(f"{work}/tags.census_summary.tsv", "census", "input_spectra")
    ms2 = val(f"{work}/qc.tsv", "qc", "count_ms2")
    if total is None or ms2 is None:
        raise SystemExit(f"{path}: FASTag reports no spectrum counts")
    return total, ms2


def conversion_status(counts, vendor, fmt):
    """(status, note) of the conversion check (§2.2, amendment 2026-11-22c): counts = (spectra, MS2 spectra) of the
    mzML, vendor = the same two counts of the vendor file itself. ok when the mzML has MS2 spectra and both counts
    agree within MS2_TOLERANCE (a partial conversion loses spectra); .wiff has no independent reader, so it rests on
    the converter's validated index"""
    if counts is None:
        return "failed", "no converted spectra"
    total, ms2 = counts
    if ms2 == 0:
        return "failed", "no MS2 spectra"
    if vendor is None:
        return ("ok", f"MS2 {ms2}; no independent vendor count (.wiff): the converter's validated index") \
            if fmt == "wiff" else ("failed", "no vendor count")
    vt, v2 = vendor
    if abs(total - vt) > MS2_TOLERANCE or abs(ms2 - v2) > MS2_TOLERANCE:
        return "failed", f"converted {total} spectra / {ms2} MS2, vendor file {vt} / {v2}"
    return "ok", f"{total} spectra, MS2 {ms2} (vendor file {vt} / {v2})"


def unpack(src, dst):
    """a Bruker .d archive to dst (a .d directory); zip CRCs are checked; members escaping dst are refused"""
    if os.path.isdir(dst):
        return dst
    tmp = dst + ".part"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    if src.lower().endswith(".zip"):
        with zipfile.ZipFile(src) as z:
            if z.testzip():
                raise SystemExit(f"{src}: bad CRC")
            for m in z.namelist():
                if m.startswith("/") or ".." in m.split("/"):
                    raise SystemExit(f"{src}: member {m} escapes")
            z.extractall(tmp)
    else:
        with tarfile.open(src) as t:
            t.extractall(tmp, filter="data")
    found = [d for d, _, fs in os.walk(tmp) if d.endswith(".d") and "analysis.tdf" in fs]
    if len(found) != 1:
        raise SystemExit(f"{src}: {len(found)} .d folders with analysis.tdf")
    os.replace(found[0], dst)
    shutil.rmtree(tmp, ignore_errors=True)
    return dst


def convert_one(r):
    k, raw, low = key(r), mirror_raw(r), r["file"].lower()
    try:
        if not os.path.exists(raw):
            return k, "failed", "", "download missing", None, None
        if low.endswith((".d.zip", ".d.tar", ".d.tar.gz")):
            inp = unpack(raw, f"{H}/in/{r['accession']}/{stem(r['file'])}.d")
            fmt, conv = "d", None
        else:
            fmt = "wiff" if low.endswith(".wiff") else "raw"
            extra = ["--tool", "mzpeak-convert"] if fmt == "wiff" else []
            subprocess.run(PXD + ["convert", r["accession"], "--files", r["file"], "--to", "mzML", *extra],
                           check=True, capture_output=True)
            inp = subprocess.run(PXD + ["path", r["accession"], r["file"], "--kind", "mzML"], check=True,
                                 capture_output=True, text=True).stdout.strip()
            conv = mzml_counts(inp)
        ven = vendor_counts(raw if fmt != "d" else inp, f"{X}/qc/{fname(k)}")
        if fmt == "d":   # the .d is read natively by every tool: the archive's integrity and readable MS2 spectra
            st, note = ("ok", f"{ven[0]} spectra, MS2 {ven[1]} (unpacked .d)") if ven and ven[1] else \
                ("failed", "no MS2 spectra in the .d")
        else:
            st, note = conversion_status(conv, ven, fmt)
        return k, st, inp, note, conv, ven
    except Exception as e:   # one failure is that run's failure (§10.1: rerun once, then missing)
        return k, "failed", "", repr(e)[:200], None, None


PXD_TOOLS = os.path.expanduser("~/.claude/skills/pride/tools.json")
STEP3 = "Step 3: the converted inputs (mzML, unpacked .d)"


def convert(workers):
    """the conversions, each tried at most ATTEMPTS times (§10.1); once every run is final, the inputs' SHA-256 are
    appended to FREEZE.md (step 3) and checked again before every use"""
    check_frozen([TOOLS["fastag"], PXD_TOOLS])
    out = CONV
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = {x["run"]: x for x in rows(out)} if os.path.exists(out) else {}
    todo = [r for r in all_runs() if done.get(key(r), {}).get("status") != "ok"
            and int(done.get(key(r), {}).get("attempts") or 0) < ATTEMPTS]
    res = dict(done)
    with ProcessPoolExecutor(workers) as ex:
        for f in as_completed([ex.submit(convert_one, r) for r in todo]):
            k, st, inp, note, conv, ven = f.result()
            prev = done.get(k, {})
            res[k] = {"run": k, "status": st, "input": inp, "note": note, "ms2_converted": conv, "ms2_vendor": ven,
                      "attempts": int(prev.get("attempts") or 0) + 1}
            print("convert", k, st, note, flush=True)
            if st != "ok" and res[k]["attempts"] < ATTEMPTS:
                FAILED.append(f"convert {k}")
    with open(out + ".part", "w") as fh:
        cols = ["run", "status", "input", "note", "ms2_converted", "ms2_vendor", "attempts"]
        fh.write("\t".join(cols) + "\n" + "".join("\t".join(str(v.get(c, "") if v.get(c) is not None else "")
                                                             for c in cols) + "\n" for v in res.values()))
    os.replace(out + ".part", out)
    if not FAILED and STEP3 not in (open(FREEZE).read() if os.path.exists(FREEZE) else ""):
        freeze_append(STEP3, [f for v in inputs().values() for f in files_of(v)])


def conversions_final():
    """stop while a conversion may still be rerun (§10.1)"""
    open_ = [x["run"] for x in rows(CONV) if x["status"] != "ok" and int(x.get("attempts") or 0) < ATTEMPTS]
    if open_:
        raise SystemExit(f"{len(open_)} conversions are due a rerun (harness.py convert): {open_[:3]}")


def check_inputs(paths):
    """the converted inputs (and pseudo-spectra) before use: their SHA-256 as appended to FREEZE.md"""
    check_frozen([f for p in paths for f in files_of(p)])


def inputs():
    return {x["run"]: x["input"] for x in rows(CONV) if x["status"] == "ok"}


# ---- the annotators' sheets (§3, blind) ------------------------------------------------------------------------------

REPORTERS = {"TMT": [126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 130.134825, 131.138180],
             "TMTpro": [126.127726, 127.124761, 128.128116, 129.131471, 130.134825, 131.138180, 132.141535,
                        133.144890, 134.148245],
             "iTRAQ4": [114.1112, 115.1083, 116.1116, 117.1150], "iTRAQ8": [113.1078, 118.1112, 119.1150, 121.1220]}


def indicators(path):
    """the raw-data indicators of §3.2 for the annotators: MS levels, analysers (filter strings), reporter-ion shares
    in 200 MS2/MS3 spectra (10 ppm; §3.2's manual check); Bruker .d: analysis.tdf frame types and acquisition"""
    if path.rstrip("/").endswith(".d"):
        db = sqlite3.connect(f"file:{path}/analysis.tdf?mode=ro", uri=True)
        frames = dict(db.execute("select MsMsType, count(*) from Frames group by MsMsType").fetchall())
        meta = dict(db.execute("select Key, Value from GlobalMetadata where Key in "
                               "('AcquisitionSoftware','InstrumentName','MzAcqRangeLower','MzAcqRangeUpper')"))
        return {"frames_by_msms_type": frames, "metadata": meta}
    import pyopenms as oms
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    levels, filt, rep, seen = {}, {}, {k: 0 for k in REPORTERS}, 0
    for s in exp:
        lv = s.getMSLevel()
        levels[lv] = levels.get(lv, 0) + 1
        fs = str(s.getMetaValue("filter string")) if s.metaValueExists("filter string") else ""
        if fs:
            filt[fs.split(" ")[0]] = filt.get(fs.split(" ")[0], 0) + 1
        if lv >= 2 and seen < 200:
            seen += 1
            mz, _ = s.get_peaks()
            for k, ions in REPORTERS.items():
                rep[k] += sum(any(abs(m - x) <= x * 1e-5 for m in mz) for x in ions) >= len(ions) - 1
    return {"ms_levels": levels, "filter_prefix": filt,
            "reporter_share": {k: round(v / max(1, seen), 3) for k, v in rep.items()}}


def sheet():
    conversions_final()
    os.makedirs(f"{H}/annotation/sheets", exist_ok=True)
    inp = inputs()
    for r in all_runs():
        dst = f"{H}/annotation/sheets/{fname(key(r))}.json"
        if os.path.exists(dst):
            continue
        try:
            with urllib.request.urlopen(f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{r['accession']}",
                                        timeout=120) as u:
                js = json.load(u)
            txt = {k: js.get(k) for k in ("title", "projectDescription", "sampleProcessingProtocol",
                                           "dataProcessingProtocol", "instruments", "organisms",
                                           "quantificationMethods", "keywords")}
        except Exception as e:
            txt = {"error": repr(e)}
        sdrf = sorted(glob.glob(f"{MIRROR}/{r['accession']}/raw/*sdrf*"))
        # a run whose conversion failed is annotated from its texts (adjudicated.tsv needs a row for every run, §3)
        ind = indicators(inp[key(r)]) if key(r) in inp else {"conversion": "failed: no raw-data indicators"}
        json.dump({"run": key(r), "input": inp.get(key(r)), "pride": txt, "sdrf": sdrf, "indicators": ind},
                  open(dst, "w"), indent=1)
        print("sheet", key(r), flush=True)


# ---- organism verification (§3.1) and windows (§3.3) ----------------------------------------------------------------

def adjudicated():
    """adjudicated.tsv (annotate.py final): accession, file, organisms, class, label, protease, cys, enrichment,
    reporter_ions, acquisition, instrument_class, sdrf_settings (JSON), organism_disputed, disputed; refdb.tsv gives
    each run's proteomes"""
    adj = {key(r): r for r in rows(f"{H}/annotation/adjudicated.tsv")}
    ref = {key(r): r for r in rows(f"{H}/annotation/refdb.tsv")}
    return adj, ref


def verify_verdict(psms, stated):
    """§3.1 from the subset search's target PSMs at 1 % FDR, each the set of database labels its peptide falls on
    (cRAP is 'crap'): verified when >= 100 PSMs pass and >= 80 % of those not explained by cRAP fall on a stated
    proteome (a shared peptide counts as stated); a mixture needs each stated component in >= 5 % of them uniquely.
    The most frequent other label among the rest is named when the stated proteomes fail."""
    n = len(psms)
    rest = [s for s in psms if s and not s <= {"crap"}]
    on = [s for s in rest if s & set(stated)]
    share = len(on) / len(rest) if rest else 0.0
    uniq = {c: sum(1 for s in rest if s - {"crap"} == {c}) / len(rest) if rest else 0.0 for c in stated}
    other = {}
    for s in rest:
        if not s & set(stated):
            for x in s - {"crap"}:
                other[x] = other.get(x, 0) + 1
    ok = n >= VERIFY_MIN_PSM and share >= VERIFY_STATED and \
        (len(stated) < 2 or all(v >= VERIFY_COMPONENT for v in uniq.values()))
    return {"psms": n, "not_crap": len(rest), "stated_share": share, "component_unique_share": uniq,
            "verified": ok, "points_to": max(other, key=other.get) if other and not ok else None}


def subset_mzml(src, dst, n=SUBSET, seed=SUBSET_SEED):
    """every MS1 and a seeded random n of the MS2 (MS3) spectra of src; a .d is first converted with the frozen
    FileConverter (the conversion prose-auto uses for ddaPASEF)"""
    import pyopenms as oms
    if os.path.exists(dst):
        return dst
    if src.rstrip("/").endswith(".d"):
        mz = f"{X}/conv/{os.path.basename(src)}.mzML"
        if not os.path.exists(mz):
            os.makedirs(os.path.dirname(mz), exist_ok=True)
            subprocess.run([TOOLS["fileconverter"], "-in", src, "-out", mz, "-out_type", "mzML", "-no_progress"],
                           check=True, capture_output=True)
        src = mz
    exp, out = oms.MSExperiment(), oms.MSExperiment()
    oms.MzMLFile().load(src, exp)
    ms2 = [i for i, s in enumerate(exp) if s.getMSLevel() >= 2]
    keep = set(random.Random(seed).sample(ms2, min(n, len(ms2))))
    for i, s in enumerate(exp):
        if s.getMSLevel() == 1 or i in keep:
            out.addSpectrum(s)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    oms.MzMLFile().store(dst + ".part.mzML", out)
    os.replace(dst + ".part.mzML", dst)
    return dst


def labelled_db(stated_ups, dst):
    """the verification database: each proteome's headers prefixed LBL_<label>|, then reversed decoys"""
    if os.path.exists(dst):
        return dst
    parts = [(f"stated{i}", f"{H}/db/{u}.fasta") for i, u in enumerate(stated_ups)] + \
        [(n, f"{H}/db/{u}.fasta") for n, u in VERIFY_PANEL.items() if u not in stated_ups] + [("crap", TOOLS["crap"])]
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst + ".part", "w") as fo:
        for lab, p in parts:
            head, seq = None, []
            for line in list(open(p)) + [">"]:
                if line.startswith(">"):
                    if head:
                        s = "".join(seq)
                        fo.write(f">{LBL}{lab}|{head}\n{s}\n>DECOY_{LBL}{lab}|{head}\n{s[::-1]}\n")
                    head, seq = line[1:].strip(), []
                else:
                    seq.append(line.strip())
    os.replace(dst + ".part", dst)
    return dst


def psm_labels(idxml, fdr=0.01):
    """per target PSM at q <= fdr (rank 1): the labels of the proteins its peptide falls on"""
    txt = open(idxml).read()
    prot = dict(re.findall(r'<ProteinHit id="([^"]+)" accession="([^"]+)"', txt))
    out = []
    for m in re.finditer(r'<PeptideIdentification .*?</PeptideIdentification>', txt, re.S):
        h = re.search(r'<PeptideHit score="([^"]+)"[^>]*?protein_refs="([^"]*)"', m.group(0))
        if not h or float(h.group(1)) > fdr:
            continue
        accs = [prot.get(x, "") for x in h.group(2).split()]
        if all(a.startswith("DECOY_") for a in accs):
            continue
        out.append({a[len(LBL):].split("|")[0] for a in accs if a.startswith(LBL)})
    return out


def prose_run(cmd, log):
    with open(log, "w") as fh:
        return subprocess.run(cmd, env=dict(os.environ, **PROSE_ENV, OPENMS_TMPDIR=f"{X}/openms"),
                              stdout=fh, stderr=subprocess.STDOUT).returncode


def drawn():
    return [key(r) for r in runs()]


def verify(i):
    check_frozen([TOOLS["prose"], TOOLS["fileconverter"]])
    adj, ref = adjudicated()
    k = drawn()[i - 1]
    a, rd = adj[k], ref[k]
    o = f"{H}/annotation/verify/{fname(k)}"
    if os.path.exists(f"{o}/verify.json") or not needs_verification(k, a, rd):
        return
    PUBLISH.append(o)
    ups = [u.strip() for u in rd["proteomes"].split(";") if u.strip()]
    src = inputs()[k]
    check_inputs([src])
    cli = A.settings_cli(a)

    def go():
        sp = src
        if str(a.get("acquisition") or "").startswith("DIA"):   # DIA: DIAspeXtract's pseudo-spectra, as the arms
            sp = pseudo_spectra(src, f"{X}/pseudo/{fname(k)}.pseudo.mzML")
        sub = subset_mzml(sp, f"{X}/subset/{fname(k)}.mzML")
        db = labelled_db(ups, f"{X}/db/verify_{'+'.join(ups)}.fasta")
        idx = f"{X}/verify/{fname(k)}.idXML"
        os.makedirs(os.path.dirname(idx), exist_ok=True)
        if prose_run(A.prose_cmd(TOOLS["prose"], sub, db, idx, cli, THREADS), f"{o}/prose.log"):
            return False
        res = verify_verdict(psm_labels(idx), [f"stated{j}" for j in range(len(ups))])
        json.dump({"run": k, "proteomes": ups, "cli": cli, **res}, open(f"{o}/verify.json", "w"), indent=1)
        return True
    if run_unit(o, "verify", go) == "missing":   # final: the run is missing, never 'not verified'
        json.dump({"run": k, "proteomes": ups, "cli": cli, "error": "the verification search failed twice"},
                  open(f"{o}/verify.json", "w"), indent=1)


def needs_verification(k, a, rd):
    """§3.1 runs on the adjudicated bottom-up runs with a reference proteome and an undisputed organism"""
    return a["class"] == "bottom-up" and rd["class"] == "bottom-up" and a.get("organism_disputed") != "yes"


def pseudo_spectra(src, dst):
    """DIAspeXtract's pseudo-spectra of a DIA input (the frozen DIAspeXtract)"""
    if not os.path.exists(dst):
        from prose_auto import external as ext
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        ext.diaspextract(TOOLS["diaspextract"], src, dst + ".part.mzML", threads=THREADS, log=lambda *x: None)
        os.replace(dst + ".part.mzML", dst)
    return dst


STD = {"orbitrap": (10.0, (20.0, "ppm")), "astral": (10.0, (20.0, "ppm")), "tof": (20.0, (40.0, "ppm")),
       "iontrap": (10.0, (0.5, "Da"))}


def windows(i):
    """§3.3 by the pilot's rule (pilot.refwin_one): a wide-window mini-search of the subset with the adjudicated
    settings (DDA mzML); DIA and timsTOF runs take the instrument standard"""
    adj, ref = adjudicated()
    k = drawn()[i - 1]
    a = adj[k]
    o = f"{H}/annotation/windows/{fname(k)}.json"
    if os.path.exists(o) or a["class"] != "bottom-up" or k not in inputs():
        return
    os.makedirs(os.path.dirname(o), exist_ok=True)
    PUBLISH.append(o)
    p = inputs()[k]
    check_inputs([p])
    if a.get("instrument_class") not in STD:   # disputed (§3): no reference windows, the arms take ProSE's defaults
        res = {"prec_window": "", "prec_unit": "", "frag_window": "", "frag_unit": "", "source": "instrument class disputed"}
    elif not p.endswith(".mzML") or a.get("acquisition", "DDA") != "DDA":
        prec, frag = STD[a["instrument_class"]]
        res = {"prec_window": prec, "prec_unit": "ppm", "frag_window": frag[0], "frag_unit": frag[1],
               "source": "instrument standard (DIA / timsTOF)"}
    else:
        db = [f"{H}/db/{u.strip()}.fasta" for u in ref[k]["proteomes"].split(";") if u.strip()]
        res = refwin(p, a, db)
    json.dump({"run": k, **res}, open(o, "w"), indent=1)


def refwin(p, a, db):
    """pilot.refwin_one with the database given (the pilot read it from the auto report, which the annotation may not
    see)"""
    import numpy as np
    from prose_auto import tagsearch as ts
    from prose_auto.hypotheses import fragment_errors
    from prose_auto.detect_tolerance import fit_gauss_uniform
    from prose_auto.spectra import load_mzml
    cls = a["instrument_class"]
    hr2 = cls != "iontrap"
    mods = ts.Mods()
    mods.cys = ts.CYS_MODS.get(a["cys"], mods.cys) if a["cys"] != "none (free Cys)" else 0.0
    k_ = next((k for p_, k in (("TMTpro", "TMTpro"), ("TMT", "TMT"), ("iTRAQ8", "iTRAQ8plex"), ("iTRAQ", "iTRAQ4plex"))
               if (a.get("label") or "").startswith(p_)), None)
    if k_:
        mods.k = mods.nterm = ts.LABEL_MASS[k_]
    enz = a["protease"] if a["protease"] in ts.ENZYMES else "Trypsin"
    run = load_mzml(p, max_ms2=20000)
    scans = ts.select_scans(run, SUBSET, "random", SUBSET_SEED)
    peaks = [ts.prepare_peaks(s) for s in scans]
    name = "ref_" + re.sub(r"\W+", "_", "+".join(map(os.path.basename, db)))
    union = f"{X}/db/{name}.fasta"   # every proteome of the manifest (a mixture's, too)
    if not os.path.exists(union):
        os.makedirs(os.path.dirname(union), exist_ok=True)
        with open(union + ".part", "w") as fo:
            for f in db:
                txt = open(f).read()
                fo.write(txt if txt.endswith("\n") else txt + "\n")
        os.replace(union + ".part", union)
    idx = ts.index_for(name, union, mods, enz)
    pw, fw = 50.0, (50.0 if hr2 else 0.6)   # §3.3: 50 ppm, fragments 50 ppm (high resolution) or 0.6 Da (ion trap)
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


PANEL_SPECIES = {"human": "Homo sapiens", "mouse": "Mus musculus", "yeast": "Saccharomyces cerevisiae",
                 "ecoli": "Escherichia coli", "arabidopsis": "Arabidopsis thaliana"}


def points_to():
    """§3.1 "another organism": the runs whose verification failed and points to a panel organism ->
    H/annotation/points_to.tsv for annotate.py points-to; their verify.json is set aside (verify_first.json), so verify
    runs again once refdb.tsv is rebuilt from the updated adjudicated.tsv (one round)"""
    rows_ = []
    for f in sorted(glob.glob(f"{H}/annotation/verify/*/verify.json")):
        v = json.load(open(f))
        if not v.get("verified") and v.get("points_to") in PANEL_SPECIES:
            rows_.append({"run": v["run"], "points_to": PANEL_SPECIES[v["points_to"]],
                          "stated_share": round(v["stated_share"], 3), "psms": v["psms"]})
            os.replace(f, f.replace("verify.json", "verify_first.json"))
    with open(f"{H}/annotation/points_to.tsv", "w") as fh:
        fh.write("run\tpoints_to\tstated_share\tpsms\n" + "".join(
            f"{r['run']}\t{r['points_to']}\t{r['stated_share']}\t{r['psms']}\n" for r in rows_))
    print(len(rows_), "runs point to another organism")


REF_COLS = ["run", "class", "organism_verified", "label", "protease", "cys", "reporter_ions", "prec_window",
            "prec_unit", "frag_window", "frag_unit"]


def reference_rows(adj, ref, ver, win, conv_ok):
    """reference.tsv rows (score.py's schema): the adjudicated class (no-proteome when refdb.tsv has no database),
    organism_verified from verify.json, the windows; a run whose conversion failed has no row (missing)"""
    out = []
    for k, a in sorted(adj.items()):
        if k not in conv_ok:
            continue
        disputed = a.get("organism_disputed") == "yes"   # §3: a disputed organism or class is not verified
        cls = "bottom-up" if disputed else a["class"] if ref.get(k, {}).get("class", "bottom-up") == "bottom-up" \
            else "no-proteome"
        v = ver.get(k, {})
        if v.get("error"):   # a verification that failed twice: the run is missing (no row), never 'not verified'
            continue
        w = win.get(k, {})
        out.append({"run": k, "class": cls,
                    "organism_verified": "yes" if v.get("verified") and not disputed else "no",
                    **{c: a.get(c, "") for c in ("label", "protease", "cys", "reporter_ions")},
                    **{c: w.get(c, "") for c in ("prec_window", "prec_unit", "frag_window", "frag_unit")}})
    return out


def reference():
    adj, ref = adjudicated()
    ver = {json.load(open(f))["run"]: json.load(open(f)) for f in glob.glob(f"{H}/annotation/verify/*/verify.json")}
    win = {json.load(open(f))["run"]: json.load(open(f)) for f in glob.glob(f"{H}/annotation/windows/*.json")}
    d = set(drawn())
    ok = set(inputs())
    open_ = [k for k in sorted(d & ok) if (needs_verification(k, adj[k], ref[k]) and k not in ver)
             or (adj[k]["class"] == "bottom-up" and k not in win)]
    if open_:   # only final results count (§10.1): an unfinished verification never becomes 'not verified'
        raise SystemExit(f"{len(open_)} runs lack a final verification or windows (array A): {open_[:3]}")
    rs = [r for r in reference_rows(adj, ref, ver, win, set(inputs())) if r["run"] in d]
    with open(f"{H}/annotation/reference.tsv", "w") as fh:
        fh.write("\t".join(REF_COLS) + "\n" + "".join("\t".join(str(r[c]) for c in REF_COLS) + "\n" for r in rs))


def topup_count():
    """(refused once any arm has run) the conditional top-up's inputs (§2.1), the drawn runs only: H/annotation/topup/{reference,conversion}.tsv
    (a drawn run without a reference row and an ok conversion is listed failed when its conversion failed), and
    score.eligible_before_arms on them"""
    import score as SC
    if glob.glob(f"{H}/runs/*/*/*/*.auto.json"):
        raise SystemExit("topup-count refused: an arm has run (the top-up comes before any arm)")
    d, o = set(drawn()), f"{H}/annotation/topup"
    os.makedirs(o, exist_ok=True)
    ref = [r for r in rows(f"{H}/annotation/reference.tsv") if r["run"] in d]
    conv = [r for r in rows(CONV) if r["run"] in d]
    with open(f"{o}/reference.tsv", "w") as fh:
        fh.write("run\tclass\torganism_verified\n" + "".join(f"{r['run']}\t{r['class']}\t{r['organism_verified']}\n"
                                                           for r in ref))
    with open(f"{o}/conversion.tsv", "w") as fh:
        fh.write("run\tstatus\n" + "".join(f"{r['run']}\t{r['status']}\n" for r in conv))
    n, counts = SC.eligible_before_arms(o, RUNS_TSV)
    print(f"eligible before the arms: {n} ({counts}); {'no top-up' if n >= SC.MIN_ELIGIBLE else 'top-up batch due'}")
    return n


# ---- analyses (IBMI, §10.1 step 5a) ----------------------------------------------------------------------------------

ANALYZE_ARMS = ("auto", "spectra_only", "texts", "auto_never", "auto_always", "fast_0110")


def eligible_runs():
    """the verified bottom-up runs of reference.tsv (§2.2: only they run every arm); top-down runs: auto only"""
    out = {}
    for r in rows(f"{H}/annotation/reference.tsv"):
        if r["class"] == "bottom-up" and r["organism_verified"] == "yes":
            out[r["run"]] = ANALYZE_ARMS
        elif r["class"] == "top-down":
            out[r["run"]] = ("auto",)
    return out


def analyze_jobs():
    _, ref = adjudicated()
    inp, jobs = inputs(), []
    for k, arms in eligible_runs().items():
        orgs = [o.strip() for o in ref.get(k, {}).get("organisms", "").split(";") if o.strip()]
        for arm in arms:   # fast_0110 too: 0.1.10 --fast needs --organism (same organisms as auto, for S7)
            jobs.append({"run": k, "arm": arm, "spectra": inp[k], "organisms": orgs,
                         "out": f"{H}/runs/{k}/{arm}"})
    for r in rows(f"{H}/annotation/scenarios.tsv"):
        k = key(r)
        if r["kind"] in ("excluded", "unresolved") or r["databases"] in ("", "none") or k not in inp:
            continue
        ups = [u.strip() for u in r["proteomes"].split(";") if u.strip()]
        jobs.append({"run": k, "arm": "scenario", "kind": r["kind"], "spectra": inp[k],
                     "fasta": [f"{H}/db/{u}.fasta" for u in ups], "out": f"{H}/scenarios/{k}/{r['kind']}"})
    return jobs


def analyze_one(j):
    if glob.glob(f"{j['out']}/*.auto.json"):
        return j["out"], "exists"
    os.makedirs(j["out"], exist_ok=True)
    tools = TOOLS_0110 if j["arm"] == "fast_0110" else TOOLS
    arm = "auto" if j["arm"] == "fast_0110" else j["arm"]
    cmd = A.analyze_cmd(tools, j["spectra"], j["out"], arm, j.get("organisms", ()), j["run"].split("/")[0],
                        j.get("fasta"), THREADS)
    if no_jev():
        cmd = [x for x in cmd if x != "--jev"]
    e = env(jev_env() if "--jev" in cmd else None)
    if j["arm"] == "fast_0110":
        cmd += ["--fast"]
        e["PYTHONPATH"] = CODE_0110
    json.dump(cmd, open(f"{j['out']}/cmd_analyze.json", "w"))

    def go():
        t = time.time()
        with open(f"{j['out']}/run.log", "w") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=j["out"], env=e).returncode
        open(f"{j['out']}/seconds", "w").write(f"{time.time() - t:.1f}\n")
        return rc == 0 and bool(glob.glob(f"{j['out']}/*.auto.json"))
    st = run_unit(j["out"], "analyze", go)
    if st == "missing":
        open(f"{j['out']}/status", "w").write("missing\n")
    return j["out"], st


SEALED = "Step 4: the annotation sealed"


def sealed():
    return os.path.exists(FREEZE) and SEALED in open(FREEZE).read()


def seal():
    """§10 step 4 ends: every conversion, verification and window is final, the top-up is settled (no batch due), and
    the annotation files are hashed into FREEZE.md; analyze refuses to start before this, and a top-up batch after"""
    if sealed():
        raise SystemExit("sealed already")
    conversions_final()
    reference()   # refuses while array A is unfinished
    import score as SC
    n = topup_count()
    dry = os.environ.get("HELDOUT_DRY_SEAL") == "1" and FREEZE != f"{ACC}/FREEZE.md"   # a 2-run dry run only
    if n < SC.MIN_ELIGIBLE and len(drawn()) < 500 and not dry:
        raise SystemExit("a top-up batch is due (draw_topup.py --batch K H/annotation/topup), then rerun step 3-4")
    freeze_append(SEALED, sorted(f for f in glob.glob(f"{H}/annotation/**", recursive=True) if os.path.isfile(f)))
    print("sealed", flush=True)


def analyze(workers):
    if not sealed():
        raise SystemExit("analyze refused: the annotation and the top-up are not sealed (harness.py seal)")
    check_frozen([TOOLS[t] for t in FROZEN_TOOLS] + [TOOLS_0110["fastag"]])
    frozen_tree(CODE_0110)
    child_import_origin(env())
    child_import_origin(dict(env(), PYTHONPATH=CODE_0110), CODE_0110)
    jobs = analyze_jobs()
    check_inputs(sorted({j["spectra"] for j in jobs}))
    with ProcessPoolExecutor(workers) as ex:
        for f in as_completed([ex.submit(analyze_one, j) for j in jobs]):
            d, st = f.result()
            print("analyze", d, st, flush=True)
            if st == "retry":
                FAILED.append(f"analyze {d}")
    freeze_append("Step 5a: the analyses (IBMI), before the transfer",
                  sorted(glob.glob(f"{H}/runs/*/*/*/*.json") + glob.glob(f"{H}/scenarios/*/*/*/*.json")))


# ---- plans (§10.1 order: IBMI -> BinAC2) -----------------------------------------------------------------------------

def _report(d):
    js = sorted(glob.glob(f"{d}/*.auto.json"))
    return json.load(open(js[0])) if js else None


def spectra_of(rep, fallback):
    """the spectra the searches read: the pseudo-spectra of a DIA analysis, else the run's input"""
    i = (rep or {}).get("input") or {}
    dia = str(((rep or {}).get("decisions") or {}).get("acquisition", {}).get("value") or "").startswith("DIA")
    return i["pseudo_spectra"] if dia and i.get("pseudo_spectra") else fallback


def run_plan(k, adj_row, ref_row, win, reports, spectra, arm_spectra=None):
    """the searches of one eligible run: auto and each ablation whose ProSE command differs from auto's (§7),
    reference, sdrf, defaults, best_psm (dynamic: the task searches its subset candidates); all against the run's
    reference database (refdb proteomes + cRAP + entrapment). spectra: the comparators' input (from the adjudicated
    acquisition, never from auto's decision); arm_spectra: each prose-auto arm's own (its report's pseudo-spectra).
    A stopped auto (no_digest; also a stopped ablation) is listed in 'stopped' (0 PSMs, noncompliant in P3, §2.2/§6)
    and the other arms still run"""
    ups = [u.strip() for u in ref_row["proteomes"].split(";") if u.strip()]
    rep = reports["auto"]
    if rep is None:
        return {"run": k, "missing": "auto analysis"}
    stopped = [a for a in ["auto", *A.ABLATIONS] if reports.get(a) and "no_digest" in reports[a]]
    prec = (win["prec_window"], win["prec_unit"]) if win.get("prec_window") not in (None, "") else None
    frag = (win["frag_window"], win["frag_unit"]) if win.get("frag_window") not in (None, "") else None
    sdrf = json.loads(adj_row["sdrf_settings"]) if adj_row.get("sdrf_settings") else {}
    ref_set = {c: adj_row.get(c) for c in ("label", "protease", "cys", "enrichment")}
    searches = {**({"auto": rep["prose_cli"]} if "auto" not in stopped else {}),
                "reference": A.settings_cli(ref_set, prec, frag),
                "sdrf": A.settings_cli(sdrf, (sdrf["prec_window"], sdrf["prec_unit"]) if sdrf.get("prec_window")
                                       else None,
                                       (sdrf["frag_window"], sdrf["frag_unit"]) if sdrf.get("frag_window") else None),
                "defaults": list(A.DEFAULTS)}
    for arm in A.ABLATIONS:
        r = reports.get(arm)
        if r and "prose_cli" in r and r["prose_cli"] != rep.get("prose_cli"):
            searches[arm] = r["prose_cli"]
    settings = {"reference": {**ref_set, "prec_window": prec[0] if prec else None, "prec_unit": prec[1] if prec else None,
                              "frag_window": frag[0] if frag else None, "frag_unit": frag[1] if frag else None},
                "sdrf": {c: sdrf.get(c) for c in ("label", "protease", "cys", "prec_window", "prec_unit",
                                                  "frag_window", "frag_unit")},
                "defaults": dict(A.DEFAULT_SETTINGS)}
    return {"run": k, "spectra": spectra, "arm_spectra": {a: v for a, v in (arm_spectra or {}).items() if a in searches},
            "stopped": stopped, "proteomes": ups, "searches": searches, "settings": settings,
            "best_psm": {"start": {c: ref_set.get(c) or "" for c in ("protease", "label", "cys")},
                         "enrichment": ref_set.get("enrichment"), "prec": prec, "frag": frag},
            "reports": {arm: glob.glob(f"{H}/runs/{k}/{arm}/*.auto.json") for arm in reports if reports[arm]}}


def comparator_spectra(k, adj_row, src):
    """the comparators' spectra from the adjudicated acquisition: DIA -> DIAspeXtract's pseudo-spectra, built here
    (independent of auto's decision) and hashed with the plans; else the converted input"""
    if str(adj_row.get("acquisition") or "").startswith("DIA"):
        return pseudo_spectra(src, f"{H}/pseudo/{fname(k)}.pseudo.mzML")
    return src


def report_db_matches(rep, proteomes):
    """§3.1 binding: the FASTA sources auto searched (its report's search_database.fasta; a union is built from them)
    are the run's manifest proteomes, by SHA-256"""
    used = ((rep or {}).get("search_database") or {}).get("fasta") or []
    return all(os.path.exists(f) for f in used) and \
        sorted(sha256(f) for f in used) == sorted(sha256(f"{H}/db/{u}.fasta") for u in proteomes)


def copy_report(src_dir, o, stopped):
    os.makedirs(o, exist_ok=True)
    for f in glob.glob(f"{src_dir}/*.auto.json"):
        shutil.copy(f, o)
    if stopped:
        open(f"{o}/status", "w").write("stopped\n")


def plan():
    if not sealed():
        raise SystemExit("plan refused: the annotation is not sealed")
    adj, ref = adjudicated()
    win = {json.load(open(f))["run"]: json.load(open(f)) for f in glob.glob(f"{H}/annotation/windows/*.json")}
    inp = inputs()
    os.makedirs(f"{H}/plan", exist_ok=True)
    release = []
    for k, arms in sorted(eligible_runs().items()):
        reports = {arm: _report(f"{H}/runs/{k}/{arm}") for arm in arms}
        if arms == ("auto",):   # top-down: the auto report goes into the scoring tree (S8); no search
            rep = reports["auto"]
            copy_report(f"{H}/runs/{k}/auto", f"{H}/out/runs/{k}/auto", bool(rep and "no_digest" in rep))
            continue
        rep = reports["auto"]
        ups = [u.strip() for u in ref[k]["proteomes"].split(";") if u.strip()]
        if rep and "no_digest" not in rep and not report_db_matches(rep, ups):
            raise SystemExit(f"{k}: auto searched a database other than its manifest proteomes {ups}")
        own = {a: spectra_of(reports.get(a), inp[k]) for a in ["auto", *A.ABLATIONS] if reports.get(a)}
        p = run_plan(k, adj[k], ref[k], win[k], reports, comparator_spectra(k, adj[k], inp[k]), own)
        if p.get("missing"):   # auto's analysis failed twice: the run's auto arm is missing in the scoring tree
            os.makedirs(f"{H}/out/runs/{k}/auto", exist_ok=True)
            open(f"{H}/out/runs/{k}/auto/status", "w").write("missing\n")
        if "spectra" in p:
            p["db"] = search_db(p["proteomes"])
        json.dump(p, open(f"{H}/plan/{fname(k)}.json", "w"), indent=1)
        if "searches" in p:
            release.append(k)
    open(f"{H}/plan/release.txt", "w").write("".join(f"{k}\n" for k in release))
    scen = []
    for d in sorted(glob.glob(f"{H}/scenarios/*/*/*")):
        k, kind = "/".join(d.split("/")[-3:-1]), d.split("/")[-1]
        rep = _report(d)
        if rep and "prose_cli" in rep:
            scen.append({"run": k, "kind": kind, "cli": rep["prose_cli"], "dir": d,
                         "fasta": rep.get("search_database", {}).get("fasta"),
                         "spectra": spectra_of(rep, inp.get(k))})
    by_run = {}
    for s in scen:
        by_run.setdefault(s["run"], []).append(s)
    for k, ss in by_run.items():   # wrong settings on the correct database (single-organism runs only, §2.3)
        cor = next((s for s in ss if s["kind"] == "correct"), None)
        if cor and not any(s["kind"] == "incomplete" for s in ss):
            for kind, f in (("wrong_label", A.wrong_label), ("wrong_protease", A.wrong_protease)):
                ss.append({**cor, "kind": kind, "cli": f(cor["cli"])})
    for ss in by_run.values():
        for s_ in ss:
            search_db(sorted({os.path.basename(f).split(".")[0] for f in s_["fasta"] or []}))
    json.dump(by_run, open(f"{H}/plan/scenarios.json", "w"), indent=1)
    open(f"{H}/plan/scenario_release.txt", "w").write("".join(f"{k}\n" for k in sorted(by_run)))
    spectra = sorted({x for f in glob.glob(f"{H}/plan/*__*.json") for x in plan_spectra(json.load(open(f)))} |
                     {s_["spectra"] for ss in by_run.values() for s_ in ss if s_.get("spectra")})
    freeze_append("Step 5: the plans, search databases and the spectra they read, before the release",
                  sorted(glob.glob(f"{H}/plan/*")) + sorted(glob.glob(f"{H}/db/search/*.fasta")) +
                  [f for x in spectra if x not in set(inp.values()) for f in files_of(x)])


def plan_spectra(p):
    return [x for x in [p.get("spectra"), *(p.get("arm_spectra") or {}).values()] if x]


# ---- BinAC2 tasks (arrays B and C) -----------------------------------------------------------------------------------

def _cmd_key(cmd):
    """a search's identity: its command without the thread count (the inputs are named in it and hashed in FREEZE.md)"""
    c = list(cmd)
    if "-threads" in c:
        i = c.index("-threads")
        del c[i:i + 2]
    return c


def search(spectra, db, cli, o, name):
    """one ProSE search into o: cmd.json, log, eval.json (entrap.evaluate + corpus.proteins); the idXML stays in o.
    A result is reused only for the same command (_cmd_key); a failed search is rerun once (run_unit), then o/status
    says missing. None when there is no result (yet)"""
    os.makedirs(o, exist_ok=True)
    js, idx = f"{o}/eval.json", f"{o}/{name}.idXML"
    cmd = A.prose_cmd(TOOLS["prose"], spectra, db, idx, cli, THREADS)
    if os.path.exists(js):
        old = json.load(open(f"{o}/cmd.json")) if os.path.exists(f"{o}/cmd.json") else None
        if old is not None and _cmd_key(old) == _cmd_key(cmd):
            return json.load(open(js))
        os.replace(js, js + ".stale")   # another command's result: never reused
    if os.path.exists(f"{o}/status"):
        return None

    def go():
        sys.path.insert(0, T)
        sys.path.insert(0, f"{T}/corpus_ids")
        import corpus as C
        import entrap
        json.dump(cmd, open(f"{o}/cmd.json", "w"))
        t = time.time()
        if prose_run(cmd, f"{o}/log") or not os.path.exists(idx):
            return False
        entrap.evaluate(js + ".part", db, [idx])
        e = json.load(open(js + ".part"))
        e["protein"] = C.proteins(idx, db)
        e["seconds"] = round(time.time() - t, 1)
        json.dump(e, open(js + ".part", "w"), indent=1)
        os.replace(js + ".part", js)
        return True
    st = run_unit(o, "search", go)
    if st == "missing":
        open(f"{o}/status", "w").write("missing\n")
    return json.load(open(js)) if st == "ok" else None


def check(spectra, idxml, fasta, report, o):
    """prose-auto check into o (<stem>.check.json) and its wall time (check_seconds); its exit code counts: a failed
    check is rerun once (run_unit 'check'), then the check is missing (no .check.json; the search stays)"""
    if glob.glob(f"{o}/*.check.json"):
        return "ok"

    def go():
        t = time.time()
        rc = subprocess.run(A.check_cmd(TOOLS, spectra, idxml, fasta, report, o, THREADS), env=env(),
                            stdout=open(f"{o}/check.log", "w"), stderr=subprocess.STDOUT).returncode
        open(f"{o}/check_seconds", "w").write(f"{time.time() - t:.1f}\n")
        if rc:
            for f in glob.glob(f"{o}/*.check.json"):
                os.remove(f)
        return rc == 0 and bool(glob.glob(f"{o}/*.check.json"))
    return run_unit(o, "check", go)


def best_psm(p, o):
    """§4 best_psm: on the 5,000-MS2 subset, one factor at a time (arms.BEST_PSM_FACTORS), then the full run. A step
    selects only when every candidate has a result: a failed one is rerun first (the task exits 1), and one out of
    attempts leaves best_psm missing (never scored 0 nor dropped)"""
    os.makedirs(o, exist_ok=True)
    if os.path.exists(f"{o}/status"):
        return None
    sub = subset_mzml(p["spectra"], f"{X}/subset/{fname(p['run'])}.mzML")
    cur = dict(p["best_psm"]["start"])
    steps, under = [], None      # under: the specific protease of the label and Cys steps when unspecific won
    for factor, cands in A.BEST_PSM_FACTORS:
        scores, pending = {}, False
        for v in dict.fromkeys([cur[factor]] + list(cands)):
            s = {**cur, factor: v, "enrichment": p["best_psm"]["enrichment"]}
            if under and factor != "protease":
                s["protease"] = under
            tag = re.sub(r'[^A-Za-z0-9]+', '_', v or 'none') + (f"__{re.sub(r'[^A-Za-z0-9]+', '_', under)}" if s.get(
                "protease") == under and under else "")
            d = f"{o}/subset/{factor}/{tag}"
            ev = search(sub, p["db"], A.settings_cli(s, p["best_psm"]["prec"], p["best_psm"]["frag"]), d, "subset")
            if ev is None:
                if os.path.exists(f"{d}/status"):   # out of attempts: the ladder is unresolved, best_psm missing
                    open(f"{o}/status", "w").write("missing\n")
                    open(f"{o}/best_psm.failed", "w").write(f"subset search {factor} = {v}\n")
                    return None
                pending = True
            scores[v] = A.subset_psms(ev)
        if pending:   # no selection before every candidate has a result (the rerun is due)
            return None
        new = A.best_psm_step(cur[factor], scores)
        steps.append({"factor": factor, "scores": scores, "from": cur[factor], "to": new})
        if factor == "protease" and new == A.UNSPECIFIC:
            under = steps[-1]["label_cys_under"] = A.best_specific(cur[factor], scores)
        cur[factor] = new
    json.dump({"steps": steps, "result": cur}, open(f"{o}/best_psm.json", "w"), indent=1)
    final = {**cur, "enrichment": p["best_psm"]["enrichment"]}
    pr, fr = p["best_psm"]["prec"] or (None, None), p["best_psm"]["frag"] or (None, None)
    json.dump({**final, "prec_window": pr[0], "prec_unit": pr[1], "frag_window": fr[0], "frag_unit": fr[1]},
              open(f"{o}/settings.json", "w"))
    return search(p["spectra"], p["db"], A.settings_cli(final, p["best_psm"]["prec"], p["best_psm"]["frag"]), o, "run")


def natural_arms(evals):
    """the natural-failure cohort of one run (score.natural_cohort's rule): SDRF or defaults below half the reference"""
    ref = evals.get("reference")
    if not ref:
        return []
    return [a for a in ("sdrf", "defaults") if evals.get(a) and evals[a]["psm"]["accepted"] < 0.5 * ref["psm"]["accepted"]]


def task(i):
    """array B: every search of release run i, the auto arm's post-search check on the correct database, the natural
    failures' checks (best_psm is array B2); outputs into H/out/runs/<acc>/<file>/<arm>/ (extract.py's layout)"""
    check_frozen([TOOLS[t] for t in FROZEN_TOOLS])
    k = open(f"{H}/plan/release.txt").read().split()[i - 1]
    p = json.load(open(f"{H}/plan/{fname(k)}.json"))
    check_frozen([p["db"]] + [f"{H}/db/{u}.fasta" for u in p["proteomes"]])
    check_inputs(plan_spectra(p))
    child_import_origin(env())
    base = f"{H}/out/runs/{k}"
    PUBLISH.extend([base, f"{H}/out/natural/{k}"])
    for arm in p.get("stopped", []):   # §2.2/§6: a stopped arm has 0 PSMs (noncompliant in P3); the others still run
        copy_report(f"{H}/runs/{k}/{arm}", f"{base}/{arm}", True)
    sp = lambda arm: (p.get("arm_spectra") or {}).get(arm, p["spectra"])
    evals = {}
    for arm, cli in p["searches"].items():
        o = f"{base}/{arm}"
        evals[arm] = search(sp(arm), p["db"], cli, o, "run")
        for f in p["reports"].get(arm, []):
            shutil.copy(f, o)
        if arm in p["settings"]:
            json.dump(p["settings"][arm], open(f"{o}/settings.json", "w"))
    if evals.get("auto"):
        check(sp("auto"), f"{base}/auto/run.idXML", [p["db"]], (p["reports"].get("auto") or [None])[0],
              f"{base}/auto")
    for arm in natural_arms(evals):
        o = f"{H}/out/natural/{k}/{arm}"
        os.makedirs(o, exist_ok=True)
        check(p["spectra"], f"{base}/{arm}/run.idXML", [p["db"]], None, o)
    print("TASK_DONE", k, flush=True)


def best_psm_task(i):
    """array B2: best_psm of release run i, on a high-memory node (>= 1 TB: the unspecific subset candidate and, when
    it wins, the final full search peak at about 900 GB) and its own walltime"""
    check_frozen([TOOLS[t] for t in FROZEN_TOOLS])
    k = open(f"{H}/plan/release.txt").read().split()[i - 1]
    p = json.load(open(f"{H}/plan/{fname(k)}.json"))
    check_frozen([p["db"]] + [f"{H}/db/{u}.fasta" for u in p["proteomes"]])
    check_inputs([p["spectra"]])
    o = f"{H}/out/runs/{k}/best_psm"
    PUBLISH.append(o)
    best_psm(p, o)
    print("BEST_PSM_DONE", k, flush=True)


def scenario(i):
    """array C: every scenario of wrong-database run i: its search on the scenario database (with cRAP and
    entrapment), the post-search check with that FASTA and the scenario's report"""
    check_frozen([TOOLS[t] for t in FROZEN_TOOLS])
    child_import_origin(env())
    k = open(f"{H}/plan/scenario_release.txt").read().split()[i - 1]
    PUBLISH.append(f"{H}/out/scenarios/{k}")
    for s in json.load(open(f"{H}/plan/scenarios.json"))[k]:
        ups = sorted({os.path.basename(f).split(".")[0] for f in s["fasta"] or []})
        db = search_db(ups)
        check_frozen([db])
        check_inputs([s["spectra"]])
        o = f"{H}/out/scenarios/{k}/{s['kind']}"
        os.makedirs(o, exist_ok=True)
        for f in glob.glob(f"{s['dir']}/*.auto.json"):
            if s["kind"] not in ("wrong_label", "wrong_protease"):
                shutil.copy(f, o)
        ev = search(s["spectra"], db, s["cli"], o, "run")
        if ev is None:
            continue
        check(s["spectra"], f"{o}/run.idXML", [db], (glob.glob(f"{s['dir']}/*.auto.json") or [None])[0], o)
    print("SCENARIO_DONE", k, flush=True)


# ---- the output tree (§10 step 6) ------------------------------------------------------------------------------------

def tree():
    """H/out is already extract.py's layout; add a status file to every missing arm of a released run, and the
    timing reports (S7) of auto and fast_0110"""
    for k in open(f"{H}/plan/release.txt").read().split():
        for arm in ("auto", "reference", "sdrf", "defaults", "best_psm", *A.ABLATIONS, "fast_0110"):
            o = f"{H}/out/runs/{k}/{arm}"
            src = f"{H}/runs/{k}/{arm}"
            if arm in ("fast_0110",) and os.path.isdir(src):
                os.makedirs(o, exist_ok=True)
                for f in glob.glob(f"{src}/*.auto.json"):
                    shutil.copy(f, o)
            if os.path.isdir(o) and not os.path.exists(f"{o}/eval.json") and not os.path.exists(f"{o}/status") \
                    and arm != "fast_0110":
                open(f"{o}/status", "w").write("missing\n")


if __name__ == "__main__":
    args = sys.argv[1:]
    opt = lambda name, d=None: args[args.index(name) + 1] if name in args else d
    stage = args[0]
    frozen_prose_auto(stage)
    if stage in ("verify", "windows", "task", "best_psm", "scenario"):
        conversions_final()
    {"setup": setup, "fasta": lambda: fasta(args[1:]), "get": get, "convert": lambda: convert(int(opt("--workers", 4))),
     "sheet": sheet, "verify": lambda: verify(int(opt("--index"))), "windows": lambda: windows(int(opt("--index"))),
     "reference": reference, "topup-count": topup_count, "analyze": lambda: analyze(int(opt("--workers", 8))), "plan": plan,
     "task": lambda: task(int(opt("--index"))),
     "best_psm": lambda: best_psm_task(int(opt("--index"))), "scenario": lambda: scenario(int(opt("--index"))),
     "tree": tree, "seal": seal, "points-to": points_to}[stage]()
    if os.environ.get("HELDOUT_PUBLISH"):   # heldout.slurm copies back exactly these (relative to H)
        open(os.environ["HELDOUT_PUBLISH"], "w").write("".join(f"{os.path.relpath(x, H)}\n" for x in PUBLISH
                                                               if os.path.exists(x)))
    if FAILED:
        print(f"{len(FAILED)} units failed and are due their rerun (§10.1):", *FAILED[:10], sep="\n  ", flush=True)
        sys.exit(1)

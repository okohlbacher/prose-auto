"""FREEZE.md of PREREG 0.1.12 §10: the SHA-256 of everything frozen, written once at step 2 and appended to at steps 3-6.

    python freeze.py init --container SHA256 --jev-model NAME   step 2: the list below -> ../FREEZE.md (refuses when it
                                                                  exists)
    python freeze.py append "TITLE" PATH...                     a later step: the files (directories recursively)

A line '<sha256>  <path>' per file (harness.check_frozen reads them); a bundle also gets 'tree <sha256>  <dir>', the
SHA-256 of its sorted 'relative path, file SHA-256' lines. Step 2 freezes: the scoring, scenario, draw and power code
and their outputs and inputs; the manifests and UniProt snapshots; the harness; every FASTA (H/db, the reference panel,
cRAP); the binaries and bundles (FASTag 1.15.3 and 1.10.0, ProSE with its libraries, DIAspeXtract, both FileConverter
builds, mzpeak-convert, the pxd converter registry); prose-auto (git commit, the prose_auto tree) and the 0.1.10 code
of fast_0110; the corpus protocol (corpus.py, entrap.py); pip freeze of the environment; the Jev model, and the Jev and
text caches; the BinAC2 container image (its SHA-256, measured on BinAC2 and passed in)."""
import glob
import hashlib
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ACC = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import harness as HN   # noqa: E402

T, H = HN.T, HN.H
OUT = HN.FREEZE
ACC_FILES = ["PREREG.md", "SHA256SUMS", "extract.py", "score.py", "scenarios.py", "resolve_uniprot.py", "refdb.py",
             "power.py", "power.txt", "power_grid.txt", "power_decision.py", "power_decision.txt", "draw_sets.py",
             "draw_topup.py", "exclude.txt", "heldout.tsv", "topup.tsv", "wrongdb.tsv", "candidates.tsv",
             "topup_candidates.tsv", "scenarios.tsv", "refdb_draw.tsv", "power_inputs", "uniprot", "uniprot_heldout",
             "heldout"]
BUNDLES = {"FASTag 1.15.3": "/ceph/ibmi/abi/oliver/opt/fastag-1.15.3/FASTag",
           "FASTag 1.10.0 (fast_0110)": "/ceph/ibmi/abi/oliver/opt/fastag-1.10.0/FASTag",
           "ProSE (corpus build)": "/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr",
           "ProSE libraries": "/ceph/ibmi/abi/oliver/opt/prose-check-20261028-deps/lib",
           "FileConverter (diaweaver build, ddaPASEF)": "/ceph/ibmi/abi/oliver/spextract/diaweaver-7ed1169",
           "FileConverter (mirror recipe, Thermo)": "/ceph/ibmi/abi/oliver/opt/openms-develop-70cebcd144-vendor",
           "mzpeak-convert 0.18.0": "/ceph/ibmi/abi/oliver/opt/mzpeak-convert/0.18.0",
           "FASTag taxonomy": f"{T}/opt/FASTag-tax-ext/share-FASTag-taxonomy",
           "prose_auto 0.1.12": f"{T}/prose-auto/prose_auto",
           "prose_auto 0.1.10 (fast_0110)": f"{T}/harness/v0110c_code",
           "Jev cache": f"{H}/jev_cache", "text cache": f"{T}/text_cache_v019"}
FILES = [f"{T}/opt/bin/diaspextract", f"{T}/search/db/crap.fasta", f"{T}/entrap.py", f"{T}/corpus_ids/corpus.py",
         os.path.expanduser("~/.claude/skills/pride/tools.json")]


def walk(path):
    if os.path.isfile(path):
        return [path]
    return sorted(os.path.join(d, f) for d, _, fs in os.walk(path) for f in fs if "__pycache__" not in d)


def tree_hash(path):
    lines = "".join(f"{os.path.relpath(f, path)}\t{HN.sha256(f)}\n" for f in walk(path))
    return hashlib.sha256(lines.encode()).hexdigest()


def section(title, paths, trees=()):
    out = [f"\n## {title} ({time.strftime('%Y-%m-%d %H:%M %Z')})\n", "```"]
    for name, d in trees:
        out.append(f"tree {tree_hash(d)}  {d}  ({name})")
    for p in paths:
        for f in walk(p):
            out.append(f"{HN.sha256(f)}  {f}")
    return "\n".join(out + ["```", ""])


def init(container, jev_model):
    if os.path.exists(OUT):
        raise SystemExit(f"{OUT} exists: the freeze is written once; later steps append")
    if not jev_model or "latest" in jev_model:
        raise SystemExit("--jev-model must name a pinned Jev model")
    git = subprocess.run(["git", "-C", f"{T}/prose-auto", "rev-parse", "HEAD"], capture_output=True,
                         text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", f"{T}/prose-auto", "status", "--short"], capture_output=True,
                           text=True).stdout.strip()
    if dirty:
        raise SystemExit("the prose-auto clone has uncommitted changes: commit before the freeze")
    pip = subprocess.run([f"{T}/env/bin/pip", "freeze"], capture_output=True, text=True).stdout
    pf = f"{H}/pip_freeze.txt"
    os.makedirs(H, exist_ok=True)
    open(pf, "w").write(pip)
    exes = [f"{d}/{x}" for d, x in ((BUNDLES["FASTag 1.15.3"], "FASTag"), (BUNDLES["FASTag 1.10.0 (fast_0110)"], "FASTag"),
                                    (BUNDLES["ProSE (corpus build)"], "bin/ProSE"),
                                    (BUNDLES["FileConverter (diaweaver build, ddaPASEF)"], "bin/FileConverter"),
                                    (BUNDLES["FileConverter (mirror recipe, Thermo)"], "bin/FileConverter"),
                                    (BUNDLES["mzpeak-convert 0.18.0"], "mzpeak-convert"))]
    head = ["# FREEZE: prose-auto 0.1.12 held-out benchmark (PREREG §10)", "",
            f"- prose-auto commit: {git}", f"- Jev model: {jev_model}",
            f"- BinAC2 container image (SHA-256, measured on BinAC2): {container}",
            f"- pip freeze: {pf}", ""]
    body = section("Step 2: scoring, scenario, draw, power and harness code, manifests, snapshots",
                   [f"{ACC}/{f}" for f in ACC_FILES])
    body += section("Step 2: FASTA files (UniProt proteomes, the reference panel, cRAP)",
                    sorted(glob.glob(f"{H}/db/*.fasta")) + sorted(glob.glob(f"{H}/panel/*.fasta")) +
                    [f"{H}/db/requests.tsv", f"{T}/search/db/crap.fasta"])
    body += section("Step 2: executables, bundles, environment, caches", exes + FILES + [pf],
                    trees=[(n, d) for n, d in BUNDLES.items() if os.path.exists(d)])
    open(OUT, "w").write("\n".join(head) + body)
    print("written", OUT)


def append(title, paths):
    if not os.path.exists(OUT):
        raise SystemExit("no FREEZE.md: run init first")
    with open(OUT, "a") as fh:
        fh.write(section(title, paths))
    print("appended", title)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "init":
        init(a[a.index("--container") + 1], a[a.index("--jev-model") + 1])
    else:
        append(a[1], a[2:])

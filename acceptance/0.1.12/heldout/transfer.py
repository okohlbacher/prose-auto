"""The IBMI -> BinAC2 hand-off of PREREG 0.1.12 §10.1: the files a BinAC2 array needs, at their /ceph paths, as an
rsync --files-from list, then rsync from an IBMI node into the workspace mirror <ws>/ceph (absolute source paths).

    python transfer.py LIST            -> the list (one absolute path per line) on stdout
    python transfer.py send [DEST]     rsync -a --files-from=LIST / DEST:<ws>  (default DEST: BinAC2's login, KEY)

The list: the frozen tree the arrays run (the environment, prose-auto, the corpus protocol, ProSE and its libraries,
FASTag and its taxonomy, DIAspeXtract, cRAP), the FASTA files and search databases, the plans, the analysis reports
and the inputs every plan names (mzML in the mirror, unpacked .d), FREEZE.md (HN.FREEZE). BinAC2 then checks the hashes again
(harness.check_frozen, heldout.slurm)."""
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import harness as HN   # noqa: E402

W = "/pfs/10/work/tu_iisko01-proseauto"
DEST = "tu_iisko01@login.binac2.uni-tuebingen.de"
KEY = os.path.expanduser("~/.ssh/binac2_ibmi_ed25519")   # the IBMI-side BinAC2 key (as pxd stage)


def paths():
    T, H = HN.T, HN.H
    out = [f"{T}/env", f"{T}/prose-auto", f"{T}/corpus_ids/corpus.py", f"{T}/entrap.py", f"{T}/search/db/crap.fasta",
           "/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr", "/ceph/ibmi/abi/oliver/opt/prose-check-20261028-deps/lib",
           "/ceph/ibmi/abi/oliver/opt/fastag-1.15.3/FASTag", f"{T}/opt/FASTag-tax-ext/share-FASTag-taxonomy",
           f"{T}/opt/bin/diaspextract", "/ceph/ibmi/abi/oliver/spextract/diaweaver-7ed1169",
           f"{H}/db", f"{H}/panel", f"{H}/plan", f"{H}/runs", f"{H}/scenarios", f"{H}/annotation/conversion.tsv",
           f"{H}/annotation/adjudicated.tsv", f"{H}/annotation/refdb.tsv", f"{H}/annotation/reference.tsv",
           f"{H}/annotation/scenarios.tsv", HN.FREEZE]
    for f in glob.glob(f"{H}/plan/*.json"):
        p = json.load(open(f))
        for s in (p.values() if isinstance(p, dict) and "run" not in p else [p]):
            for x in (s if isinstance(s, list) else [s]):
                if isinstance(x, dict) and x.get("spectra"):
                    out.append(x["spectra"])
    out += [v for v in HN.inputs().values()]
    files = []
    for p in dict.fromkeys(out):
        files += HN_walk(p)
    return files


def HN_walk(p):
    if os.path.isfile(p):
        return [p]
    return sorted(os.path.join(d, f) for d, _, fs in os.walk(p) for f in fs if "__pycache__" not in d)


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print("\n".join(paths()))
    elif sys.argv[1] == "send":
        lst = f"{HN.H}/transfer.list"
        open(lst, "w").write("\n".join(p.lstrip("/") for p in paths()) + "\n")
        dest = sys.argv[2] if len(sys.argv) > 2 else DEST
        subprocess.run(["rsync", "-a", "-e", f"ssh -o BatchMode=yes -i {KEY} -o IdentitiesOnly=yes",
                        f"--files-from={lst}", "/", f"{dest}:{W}/"], check=True)

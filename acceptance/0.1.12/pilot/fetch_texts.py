"""The annotator's texts for the pilot runs: PRIDE project fields and any SDRF file of the project (PREREG §0).
python3 fetch_texts.py -> texts.json, sdrf/"""
import csv, json, os, sys, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import draw_sets as D   # noqa: E402

out = {}
os.makedirs(f"{HERE}/sdrf", exist_ok=True)
for r in csv.DictReader(open(f"{HERE}/pilot.tsv"), delimiter="\t"):
    acc = r["accession"]
    p = D.get(f"{D.API}/projects/{acc}")
    got = []
    for f in D.raw_files(acc):
        name = f.get("fileName", "")
        if "sdrf" not in name.lower():
            continue
        loc = next((l["value"] for l in f.get("publicFileLocations", []) if str(l.get("value", "")).startswith("ftp")), None)
        if loc:
            dst = f"{HERE}/sdrf/{acc}_{name}"
            urllib.request.urlretrieve(loc, dst)
            got.append(os.path.basename(dst))
    out[acc] = {k: p.get(k) for k in ("title", "projectDescription", "sampleProcessingProtocol", "dataProcessingProtocol",
                                       "instruments", "organisms", "quantificationMethods", "keywords", "experimentTypes")}
    out[acc]["sdrf"] = got
    print(acc, len(got), flush=True)
json.dump(out, open(f"{HERE}/texts.json", "w"), indent=1)

"""The reference annotation of PREREG 0.1.12 §3 (amendment 2026-11-22c): two independent agent annotators on different
model families, A (Claude) and B (Codex; Kimi as the fallback), one reconciliation round, no human adjudicator.

    python annotate.py prompts SHEETS OUT        one prompt per run from the blind sheets (harness.py sheet)
    python annotate.py run A|B OUT [--workers W] [--cli claude|codex|kimi]
    python annotate.py reconcile A|B OUT          the disputed runs: each annotator sees the other's answer and evidence
    python annotate.py final OUT                  reference_A.tsv, reference_B.tsv, reconciled_*.tsv, adjudicated.tsv,
                                                  agreement.tsv, log.json

Blindness: each annotator gets one run's sheet inside its prompt (PRIDE project texts, the SDRF rows of that file,
the raw-data indicators) and nothing else; it runs from an empty folder with every tool disabled (Claude) or a
read-only sandbox (Codex, Kimi), and never sees prose-auto's output, an arm, or the other annotator before the
reconciliation round. The rule below uses no arm value.

The rule: a field the two agree on (after normalisation: score.canon_label for labels, species names, the vocabulary
below) is the reference. A field still disputed after the reconciliation round is "not annotated":
- organisms or class disputed: the run is "organism not verified" (out of P1-P3, counted);
- any other field disputed: blank in reference.tsv, so the run leaves that S-endpoint only (S1-S4 skip a blank
  reference value); an SDRF-arm field disputed is unstated (the ProSE default, §4).
Agreement per field, raw and Cohen's kappa, is reported for the first round and after reconciliation."""
import csv
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import score as SC   # noqa: E402

VOCAB = {
    "class": ["bottom-up", "top-down", "no-proteome"],
    "label": ["label free", "TMT2plex", "TMT6plex", "TMT10/11plex", "TMTpro (16/18plex)", "TMT0 / TMT (126 only)",
              "iTRAQ4plex", "iTRAQ8plex", "SILAC (Lys8/Arg10)", "SILAC (Arg6)", "SILAC Lys4 or dimethyl (d4)",
              "SILAC (Lys4/Arg6)", "SILAC triplex (Lys4/Arg6, Lys8/Arg10)", "15N metabolic labelling",
              "Heavy water (2H2O) labelling"],
    "protease": ["Trypsin", "Trypsin/P", "Lys-C", "Arg-C", "Chymotrypsin", "Glu-C", "Asp-N", "Lys-N",
                 "unspecific cleavage"],
    "cys": ["Carbamidomethyl (C)", "none (free Cys)", "Propionamide (C)", "Nethylmaleimide (C)", "Methylthio (C)",
            "Carboxymethyl (C)"],
    "enrichment": ["none", "Phospho (STY)", "Phospho (Y)", "Acetyl (K)", "GlyGly (K)", "glyco"],
    "reporter_ions": ["yes", "no"],
    "acquisition": ["DDA", "DIA"],
    "instrument_class": ["orbitrap", "astral", "tof", "iontrap"],
}
FIELDS = ["organisms", "class"] + [f for f in VOCAB if f != "class"]
SDRF_FIELDS = ["label", "protease", "cys", "prec_window", "prec_unit", "frag_window", "frag_unit"]
CLIS = {"claude": lambda p: ["claude", "-p", p, "--output-format", "text", "--tools", ""],
        "codex": lambda p: ["codex", "exec", "-s", "read-only", "--skip-git-repo-check", "--color", "never", p],
        "kimi": lambda p: ["kimi", "-p", p]}
DEFAULT_CLI = {"A": "claude", "B": "codex"}
FALLBACK = {"codex": "kimi"}

INSTRUCTIONS = """You are annotator {who} of a preregistered proteomics benchmark. Annotate ONE mass-spectrometry run
from the material below only: the PRIDE project texts, the SDRF rows of this file (if any) and raw-data indicators
measured from the spectra (MS levels, analyser filter-string prefixes, the share of MS2/MS3 spectra carrying
isobaric reporter ions; for Bruker .d the frame types). Do not use any tool, file or web page. Decide each field from
the evidence; where the methods and the indicators disagree, the indicators win for reporter_ions, acquisition and
instrument_class.

Fields (use exactly these values):
- organisms: the species in the sample as a list of binomial names (two words each), e.g. ["Homo sapiens"];
- class: {class};
- label: {label};
- protease: {protease};
- cys: the cysteine alkylation, {cys};
- enrichment: {enrichment};
- reporter_ions: whether reporter ions are present in MS2/MS3: yes or no;
- acquisition: DDA or DIA;
- instrument_class: the MS2 analyser class: {instrument_class};
- sdrf_settings: the authors' settings AS STATED (the SDRF if one exists, else the methods): an object with any of
  label, protease, cys (same vocabularies), prec_window, prec_unit (ppm or Da), frag_window, frag_unit; leave out a
  field the authors do not state;
- evidence: one short sentence per field saying what decided it.

Answer with one JSON object and nothing else.

THE RUN
{sheet}
"""

RECONCILE = """You annotated the run below; a second, independent annotator disagrees on some fields. You see their
answer and evidence. Re-examine the same material (and nothing else: no tools, files or web pages) and give your
final answer for ALL fields in the same JSON format. Change a field only if their evidence convinces you.

Your first answer:
{mine}

The other annotator's answer:
{theirs}

Disputed fields: {fields}

THE RUN
{sheet}
"""


def prompt(sheet, who):
    return INSTRUCTIONS.format(who=who, sheet=json.dumps(sheet, indent=1),
                               **{k: " | ".join(v) for k, v in VOCAB.items()})


def parse(text):
    """the last JSON object in an agent's reply"""
    for m in reversed(list(re.finditer(r"\{", text))):
        try:
            return json.loads(text[m.start():text.rindex("}") + 1])
        except ValueError:
            continue
    return None


def species(v):
    names = v if isinstance(v, list) else str(v or "").split(";")
    return tuple(sorted({" ".join(re.sub(r"[^A-Za-z ]", " ", str(n)).split()[:2]).capitalize() for n in names
                         if str(n).strip()}))


def norm(field, v):
    """the value two annotators are compared on (and the one written to the reference)"""
    if field == "organisms":
        return species(v)
    if v is None or str(v).strip() == "":
        return None
    if field == "label":
        return SC.canon_label(v)
    if field in VOCAB:
        hit = [x for x in VOCAB[field] if x.lower() == str(v).strip().lower()]
        return hit[0] if hit else f"?{v}"   # outside the vocabulary: never equal to a valid value
    return str(v).strip()


def norm_sdrf(s):
    s = s or {}
    return {f: (norm(f, s.get(f)) if f in ("label", "protease", "cys") else
                (float(s[f]) if f in ("prec_window", "frag_window") and s.get(f) not in (None, "") else
                 (str(s.get(f)) if s.get(f) not in (None, "") else None))) for f in SDRF_FIELDS}


def disputed(a, b):
    """the fields (and sdrf_settings.<field>) on which two answers differ"""
    out = [f for f in FIELDS if norm(f, a.get(f)) != norm(f, b.get(f))]
    sa, sb = norm_sdrf(a.get("sdrf_settings")), norm_sdrf(b.get("sdrf_settings"))
    return out + [f"sdrf_settings.{f}" for f in SDRF_FIELDS if sa[f] != sb[f]]


def consensus(a, b):
    """adjudicated.tsv's row from two (reconciled) answers: agreed fields kept, disputed ones blank; a disputed
    organism or class marks the run not verified (§3 amendment 2026-11-22c)"""
    d = disputed(a, b)
    row = {f: (None if f in d else norm(f, a.get(f))) for f in FIELDS}
    sd = {f: norm_sdrf(a.get("sdrf_settings"))[f] for f in SDRF_FIELDS if f"sdrf_settings.{f}" not in d}
    row["organisms"] = "; ".join(row["organisms"]) if row["organisms"] else ""
    row["sdrf_settings"] = json.dumps({k: v for k, v in sd.items() if v is not None})
    row["organism_disputed"] = "yes" if ("organisms" in d or "class" in d) else "no"
    row["class"] = row["class"] or "bottom-up"   # disputed class: the run is not verified (above), never top-down
    row["disputed"] = ",".join(d)
    return row


def kappa(pairs):
    """(raw agreement, Cohen's kappa) of paired categorical values"""
    n = len(pairs)
    if not n:
        return None, None
    po = sum(x == y for x, y in pairs) / n
    cats = {x for p in pairs for x in p}
    pe = sum((sum(x == c for x, _ in pairs) / n) * (sum(y == c for _, y in pairs) / n) for c in cats)
    return po, (1.0 if pe == 1 else (po - pe) / (1 - pe))


def agreement(A, B):
    """per field over the runs both annotated"""
    keys = sorted(set(A) & set(B))
    out = {}
    for f in FIELDS:
        out[f] = kappa([(str(norm(f, A[k].get(f))), str(norm(f, B[k].get(f)))) for k in keys])
    for f in SDRF_FIELDS:
        out[f"sdrf_settings.{f}"] = kappa([(str(norm_sdrf(A[k].get("sdrf_settings"))[f]),
                                            str(norm_sdrf(B[k].get("sdrf_settings"))[f])) for k in keys])
    return out


# ---- the agent calls -------------------------------------------------------------------------------------------------

def call(cli, text, timeout=1800):
    """one agent call from an empty folder; returns the reply text"""
    with tempfile.TemporaryDirectory(prefix="annot_") as d:
        r = subprocess.run(CLIS[cli](text), cwd=d, capture_output=True, text=True, timeout=timeout,
                           stdin=subprocess.DEVNULL)
    return r.stdout


def answer(cli, text, raw_path):
    """an agent's parsed answer; one retry, then the fallback CLI (Codex -> Kimi), recorded in raw_path"""
    tries = []
    for c in [cli, cli] + ([FALLBACK[cli]] if cli in FALLBACK else []):
        out = call(c, text)
        tries.append({"cli": c, "reply": out})
        js = parse(out)
        if js:
            json.dump(tries, open(raw_path, "w"), indent=1)
            return {**js, "_cli": c}
    json.dump(tries, open(raw_path, "w"), indent=1)
    return None


def run(who, out, workers, cli):
    os.makedirs(f"{out}/{who}", exist_ok=True)
    todo = [p for p in sorted(glob.glob(f"{out}/prompts/*.json"))
            if not os.path.exists(f"{out}/{who}/{os.path.basename(p)}")]

    def one(p):
        job = json.load(open(p))
        js = answer(cli, prompt(job["sheet"], who), f"{out}/{who}/{os.path.basename(p)}.raw")
        json.dump({"run": job["run"], "answer": js}, open(f"{out}/{who}/{os.path.basename(p)}", "w"), indent=1)
        print(who, job["run"], "ok" if js else "NO ANSWER", flush=True)
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, todo))


def answers(out, who):
    return {j["run"]: j["answer"] for j in (json.load(open(f)) for f in glob.glob(f"{out}/{who}/*.json"))
            if j.get("answer")}


def reconcile(who, out, cli):
    other = "B" if who == "A" else "A"
    mine, theirs = answers(out, who), answers(out, other)
    os.makedirs(f"{out}/{who}_reconciled", exist_ok=True)
    for p in sorted(glob.glob(f"{out}/prompts/*.json")):
        job = json.load(open(p))
        k, dst = job["run"], f"{out}/{who}_reconciled/{os.path.basename(p)}"
        if os.path.exists(dst) or k not in mine or k not in theirs:
            continue
        d = disputed(mine[k], theirs[k])
        if not d:
            continue
        strip = lambda x: {a: b for a, b in x.items() if not a.startswith("_")}
        text = RECONCILE.format(mine=json.dumps(strip(mine[k]), indent=1), theirs=json.dumps(strip(theirs[k]), indent=1),
                                fields=", ".join(d), sheet=json.dumps(job["sheet"], indent=1))
        js = answer(cli, text, dst + ".raw")
        json.dump({"run": k, "answer": js}, open(dst, "w"), indent=1)
        print(who, "reconciled", k, flush=True)


def prompts(sheets, out):
    os.makedirs(f"{out}/prompts", exist_ok=True)
    for f in sorted(glob.glob(f"{sheets}/*.json")):
        s = json.load(open(f))
        sheet = {"run": s["run"], "pride": s["pride"], "indicators": s["indicators"],
                 "sdrf_rows": sdrf_rows(s.get("sdrf") or [], s["run"].split("/", 1)[1])}
        json.dump({"run": s["run"], "sheet": sheet}, open(f"{out}/prompts/{os.path.basename(f)}", "w"), indent=1)


def sdrf_rows(paths, file):
    """the SDRF rows of this raw file (all rows when none names it), as dicts"""
    out = []
    for p in paths:
        rs = list(csv.DictReader(open(p, encoding="utf-8", errors="replace"), delimiter="\t"))
        hit = [r for r in rs if any(os.path.basename(str(v)) == file for v in r.values())]
        out += hit or rs[:5]
    return out


def final(out):
    A0, B0 = answers(out, "A"), answers(out, "B")
    A1 = {**A0, **answers(out, "A_reconciled")}
    B1 = {**B0, **answers(out, "B_reconciled")}
    cols = ["accession", "file", "organisms", "class", "label", "protease", "cys", "enrichment", "reporter_ions",
            "acquisition", "instrument_class", "sdrf_settings", "organism_disputed", "disputed"]

    def write(name, rows_):
        with open(f"{out}/{name}", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows_:
                fh.write("\t".join("" if r.get(c) is None else str(r[c]) for c in cols) + "\n")

    def single(ans):
        return [{"accession": k.split("/", 1)[0], "file": k.split("/", 1)[1],
                 **{f: ("; ".join(norm(f, v.get(f))) if f == "organisms" else norm(f, v.get(f))) for f in FIELDS},
                 "sdrf_settings": json.dumps(v.get("sdrf_settings") or {})} for k, v in sorted(ans.items())]
    write("reference_A.tsv", single(A0))
    write("reference_B.tsv", single(B0))
    write("reconciled_A.tsv", single(A1))
    write("reconciled_B.tsv", single(B1))
    both = sorted(set(A1) & set(B1))
    adj = [{"accession": k.split("/", 1)[0], "file": k.split("/", 1)[1], **consensus(A1[k], B1[k])} for k in both]
    write("adjudicated.tsv", adj)
    ag0, ag1 = agreement(A0, B0), agreement(A1, B1)
    with open(f"{out}/agreement.tsv", "w") as fh:
        fh.write("field\traw_first\tkappa_first\traw_reconciled\tkappa_reconciled\n")
        for f in ag0:
            fh.write(f"{f}\t" + "\t".join("" if x is None else f"{x:.3f}" for x in (*ag0[f], *ag1[f])) + "\n")
    json.dump({"runs_A": len(A0), "runs_B": len(B0), "both": len(both),
               "missing_an_answer": sorted(set(A0) ^ set(B0)),
               "clis": {"A": sorted({v.get("_cli") for v in A1.values()}),
                        "B": sorted({v.get("_cli") for v in B1.values()})}}, open(f"{out}/log.json", "w"), indent=1)


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda n, d=None: a[a.index(n) + 1] if n in a else d
    if a[0] == "prompts":
        prompts(a[1], a[2])
    elif a[0] == "run":
        run(a[1], a[2], int(opt("--workers", 4)), opt("--cli", DEFAULT_CLI[a[1]]))
    elif a[0] == "reconcile":
        reconcile(a[1], a[2], opt("--cli", DEFAULT_CLI[a[1]]))
    elif a[0] == "final":
        final(a[1])

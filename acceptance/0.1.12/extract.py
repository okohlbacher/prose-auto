"""The scoring inputs from the arm outputs (PREREG 0.1.12 §10), frozen with score.py. Nothing is typed by hand.

    python extract.py OUT DIR      -> DIR/results.tsv, decisions.tsv, checks.tsv, timing.tsv

OUT is the frozen output tree of the benchmark (§10):
  OUT/runs/<accession>/<file>/<arm>/      one folder per held-out run and arm
      eval.json         the counting step of the ID corpus (corpus_ids/corpus.py): psm.accepted, psm.combined,
                        peptide.accepted, peptide.combined, protein.groups (FDP in %)
      run.auto.json     the prose-auto report of an arm that runs prose-auto (auto, auto_never, auto_always, ablations)
      settings.json     the settings of an arm that does not (reference, sdrf, defaults, best_psm): label, protease,
                        cys, prec_window, prec_unit, frag_window, frag_unit
      check.json        prose-auto check on this arm's result (auto: the post-search check on the correct database)
      status            optional, one word: missing (download, conversion or search failed)
  OUT/scenarios/<accession>/<file>/<kind>/  run.auto.json (the pre-search check), check.json (post) per scenario
  OUT/natural/<accession>/<file>/<arm>/check.json   the natural-failure checks (§6 S6)
An arm is stopped when its report has no_digest (the top-down stop) or, for a search, when status says stopped; it is
missing when it has neither counts nor a stop. Standard library only."""
import csv
import glob
import json
import os
import sys

VERDICT = {"ok": "ok", "warning": "warning", "not assessed": "not assessed", "settings likely wrong": "settings",
           "database likely wrong or incomplete": "database"}
QC_ANALYZER = {"FTMS": "ORBITRAP", "ITMS": "IT", "TOF": "TOF", "ASTMS": "ASTRAL"}
STATUS = {"auto": "auto", "review": "review", "default": "default", "given": "given", "unknown": "unknown"}
SETTING = {"label": "label", "protease": "enzyme", "cys": "cys_alkylation"}


def _json(p):
    return json.load(open(p)) if os.path.exists(p) else None


def _report(d):
    js = sorted(glob.glob(os.path.join(d, "*.auto.json")))
    return _json(js[0]) if js else None


def result(d):
    """(status, psm, pep, pg, psm_fdp, pep_fdp) of one arm folder"""
    st = open(os.path.join(d, "status")).read().strip() if os.path.exists(os.path.join(d, "status")) else None
    rep = _report(d)
    if st == "stopped" or (rep and "no_digest" in rep):
        return ("stopped",) + ("",) * 5
    e = _json(os.path.join(d, "eval.json"))
    if st == "missing" or not e:
        return ("missing",) + ("",) * 5
    return ("ok", e["psm"]["accepted"], e["peptide"]["accepted"], e["protein"]["groups"],
            round(e["psm"]["combined"], 4), round(e["peptide"]["combined"], 4))


def decisions(d):
    """[(setting, value, unit, status)] of one arm folder: from its report, else from its settings.json"""
    rep = _report(d)
    out = []
    if rep and "decisions" in rep:
        dec = rep["decisions"]
        for s, k in SETTING.items():
            if k in dec:
                out.append((s, dec[k]["value"], "", STATUS.get(dec[k]["status"], "review")))
        for s, k in (("prec_window", "precursor_tolerance"), ("frag_window", "fragment_tolerance")):
            v = (dec.get(k) or {}).get("value")
            if isinstance(v, dict) and v.get("value") is not None:
                out.append((s, v["value"], v["unit"], STATUS.get(dec[k]["status"], "review")))
        if "ms2_analyzer" in dec:
            out.append(("ms2_analyzer_header", dec["ms2_analyzer"]["value"], "", "auto"))
        cls = (rep.get("fastag_qc") or {}).get("ms2_analyzer_class") or {}
        if cls:
            top = max(cls, key=cls.get)
            out.append(("ms2_analyzer_qc", QC_ANALYZER.get(top, top), "", "auto"))
        return out
    st = _json(os.path.join(d, "settings.json"))
    for s in ("label", "protease", "cys"):
        if st and st.get(s) is not None:
            out.append((s, st[s], "", "given"))
    for s, u in (("prec_window", "prec_unit"), ("frag_window", "frag_unit")):
        if st and st.get(s) is not None:
            out.append((s, st[s], st[u], "given"))
    return out


def pre(d):
    rep = _report(d)
    v = ((rep or {}).get("database_check") or {}).get("status")
    return VERDICT.get(v) if v else None


def post(d):
    c = _json(os.path.join(d, "check.json"))
    return VERDICT.get(c["verdict"]) if c else None


def input_type(rep):
    dec = rep.get("decisions") or {}
    acq = str((dec.get("acquisition") or {}).get("value"))
    ms2 = str((dec.get("ms2_analyzer") or {}).get("value"))
    return "timsTOF" if "PASEF" in acq else "DIA" if acq.startswith("DIA") else "ion trap" if ms2 in ("IT", "LIT") \
        else "TOF" if ms2 == "TOF" else "Orbitrap DDA"


def main(out, d):
    os.makedirs(d, exist_ok=True)
    res, dec, chk, tim = [], [], [], []
    for a in sorted(glob.glob(os.path.join(out, "runs", "*", "*", "*"))):
        acc, f, arm = a.split(os.sep)[-3:]
        run = f"{acc}/{f}"
        res.append((run, arm) + result(a))
        dec += [(run, arm) + x for x in decisions(a)]
        rep = _report(a)
        if arm == "auto":
            for k, v in (("pre", pre(a)), ("post", post(a))):
                if v:
                    chk.append((run, "correct", k, v))
        if rep and "timings_s" in rep and arm in ("auto", "fast_0110") and "decisions" in rep:
            ver = "0.1.12" if arm == "auto" else "0.1.10"
            tim += [(run, ver, input_type(rep), s, x) for s, x in rep["timings_s"].items()]
    for s in sorted(glob.glob(os.path.join(out, "scenarios", "*", "*", "*"))):
        acc, f, kind = s.split(os.sep)[-3:]
        for k, v in (("pre", pre(s)), ("post", post(s))):
            if v and not (k == "pre" and kind in ("wrong_label", "wrong_protease")):
                chk.append((f"{acc}/{f}", kind, k, v))
    for s in sorted(glob.glob(os.path.join(out, "natural", "*", "*", "*"))):
        acc, f, arm = s.split(os.sep)[-3:]
        v = post(s)
        if v:
            chk.append((f"{acc}/{f}", "natural", "post", v))
    for name, head, rows in (("results", ["run", "arm", "status", "psm", "pep", "pg", "psm_fdp", "pep_fdp"], res),
                             ("decisions", ["run", "arm", "setting", "value", "unit", "status"], dec),
                             ("checks", ["run", "kind", "check", "verdict"], chk),
                             ("timing", ["run", "version", "input_type", "step", "seconds"], tim)):
        with open(os.path.join(d, f"{name}.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(head)
            w.writerows(rows)
    print(f"results {len(res)}, decisions {len(dec)}, checks {len(chk)}, timing {len(tim)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])

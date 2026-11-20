"""prose-auto check: after the search, is the database or are the settings wrong? (0.1.11)

    prose-auto check run.mzML run.idXML --fasta db.fasta --fastag FASTag [--report run.auto.json]
                     [--contaminants crap.fasta] [--fdr 0.01] [-o outdir]

The share of all MS2 identified varies with the sample (HLA ligands, single cells, DIA, ion trap). The relative
signal is the share of the spectra with a confident FASTag tag (E <= 1e-4) that the search identified: such a
spectrum carries a good peptide, so a search with the right database and settings should explain most of them. The
confident-tag spectra the search did not identify are placed on the given database and on the check panel (FASTag's
human, mouse, yeast, E. coli and the contaminants) in one FASTag pass:
- they place on the given database: its peptides are there, so the settings (label, protease, modifications,
  windows) are likely wrong;
- they place elsewhere or nowhere: the database is likely wrong or incomplete, naming a panel organism that
  explains many of them.
The thresholds are calibrated on the development set (v0111/README_0111.md).
"""

from __future__ import annotations

import json
import os
import re
import tempfile

from . import external as ext

# calibrated on the 218 release runs with correct and wrong databases and deliberately wrong settings
CHECK_MIN_CONFIDENT = 100   # spectra with a confident tag, else not assessed
CHECK_ID_OK = 0.15          # the search identifies this share of the confident-tag spectra: ok
CHECK_GIVEN_MIN = 0.20      # unidentified confident-tag spectra placing on the given database: settings, else database
CHECK_MEMBER_MIN = 0.10     # a panel member places this share of them (and the given database does not): named


def identified(path, fdr=0.01):
    """The spectra (native ids) a ProSE result identifies: idXML (PeptideIdentification spectrum_reference, first
    hit q-value <= fdr, not a decoy) or mzTab (PSM rows: spectra_ref, a q-value column, opt_global_target_decoy).
    A result without q-values (ProSE does not rescore fewer than 100 PSMs) identifies nothing. Returns (spectra,
    total PSM rows)."""
    if path.lower().endswith(".mztab"):
        return _mztab(path, fdr)
    txt = open(path, encoding="utf-8", errors="replace").read()
    out, n = set(), 0
    for m in re.finditer(r"<PeptideIdentification\b([^>]*)>(.*?)</PeptideIdentification>", txt, re.S):
        head, body = m.groups()
        n += 1
        st = re.search(r'score_type="([^"]*)"', head)
        ref = re.search(r'spectrum_reference="([^"]*)"', head)
        hit = re.search(r'<PeptideHit\b[^>]*score="([^"]+)"[^>]*>(.*?)</PeptideHit>', body, re.S)
        if not (ref and hit):
            continue
        if st and st.group(1) and "q-value" not in st.group(1).lower():
            # ProSE skips Percolator below 100 PSMs and keeps the search scores: nothing passes a q-value cutoff
            return set(), n + txt.count("<PeptideIdentification", m.end())
        decoy = re.search(r'name="target_decoy" value="decoy"', hit.group(2))
        if float(hit.group(1)) <= fdr and not decoy:
            out.add(ref.group(1).replace("&quot;", '"'))
    return out, n


def _mztab(path, fdr):
    """mzTab 1.0: the metadata declares what psm_search_engine_score[n] is; the one named a (PSM) q-value is the
    column search_engine_score[n] (0.1.12; 0.1.11 looked for "q-value" in the column names only and read 0
    identifications from OpenMS exports). Else a column named q-value."""
    head, out, n, qcol = None, set(), 0, None
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "MTD" and len(f) > 2 and f[1].startswith("psm_search_engine_score[") and qcol is None \
                and "q-value" in f[2].lower() and "peptide-level" not in f[2].lower():
            qcol = "search_engine_score[" + f[1].split("[", 1)[1]
        elif f[0] == "PSH":
            head = f
        elif f[0] == "PSM" and head:
            r = dict(zip(head, f))
            n += 1
            q = r.get(qcol) if qcol else \
                next((v for k, v in r.items() if "q-value" in k.lower() or k == "opt_global_q-value"), None)
            decoy = r.get("opt_global_target_decoy", r.get("opt_global_cv_MS:1002217_decoy_peptide", "")) in ("decoy", "1")
            try:
                ok = q is not None and float(q) <= fdr and not decoy
            except ValueError:
                ok = False
            if ok:
                out.add(re.sub(r"^ms_run\[\d+\]:", "", r.get("spectra_ref", "")))
    return out, n


def _labelled(parts, dst):
    ext.write_labelled(dst, parts)   # 0.1.12: a newline between the files


def placements(fastag, spectra, fastas, contaminants=None, frag=(20.0, "ppm"), fixed=(), threads=0, log=print,
               max_spectra=20000):
    """One FASTag pass: the confident-tag spectra and, per label (given, human, mouse, yeast, ecoli, contaminants),
    those whose confident tag places on it with an agreeing flank; and the MS2 count. Returns
    (confident, {label: spectra}, n_ms2)."""
    share = ext.fastag_share(fastag, "alkylation")
    parts = [("given", f) for f in fastas]
    parts += [(n, os.path.join(share, f"{n}.fasta.gz")) for n in ("human", "mouse", "yeast", "ecoli")] if share else []
    if contaminants:
        parts.append(("contaminants", contaminants))
    with tempfile.TemporaryDirectory(prefix="prose_auto_check_") as tmp:
        db = os.path.join(tmp, "check.fasta")
        _labelled(parts, db)
        cmd = ext._gap_cmd(fastag, spectra, db, tmp, *frag, threads, 0, fixed)
        cmd += ["-max_evalue", ext.CONFIDENT_EVALUE, "-qc_out", os.path.join(tmp, "qc.tsv")]
        cmd += ["-subsample_spectra", max_spectra] if max_spectra else []
        ext._run(cmd, log)
        conf, placed = ext.exact_placements(tmp)
        qc = ext._report(os.path.join(tmp, "qc.tsv")) if os.path.exists(os.path.join(tmp, "qc.tsv")) else {}
    # -subsample_spectra draws from all spectra (MS1 too): scale the subsample's MS2 count to the whole file
    read, total = ext._num(qc, "qc", "spectra_read"), ext._num(qc, "qc", "spectra_input")
    n_ms2 = ext._num(qc, "qc", "count_ms2") * (total / read if read else 1.0)
    return conf, {k[len(ext.DBCHECK_PREFIX):]: v for k, v in placed.items() if k}, int(round(n_ms2))


def classify(conf, ided, placed, n_ms2=0, n_ided=None):
    """The post-search verdict from the confident-tag spectra (conf), the identified spectra (ided) and the
    placements per label. Returns the report dict."""
    n = len(conf)
    hit = conf & ided
    un = conf - ided
    given = placed.get("given", set())
    crap = placed.get("contaminants", set())
    members = {m: s for m, s in placed.items() if m not in ("given", "contaminants")}
    on_given = len(un & given) / len(un) if un else 0.0
    elsewhere = {m: round(len((un & s) - given - crap) / len(un), 4) for m, s in members.items()} if un else {}
    nowhere = len(un - given - crap - set().union(*members.values())) / len(un) if un else 0.0
    out = {"ms2": n_ms2, "identified": n_ided if n_ided is not None else len(ided),
           "id_rate": round((n_ided if n_ided is not None else len(ided)) / n_ms2, 4) if n_ms2 else None,
           "confident_spectra": n, "confident_identified": round(len(hit) / n, 4) if n else None,
           "unidentified": {"spectra": len(un), "on_given": round(on_given, 4), "elsewhere": elsewhere,
                            "contaminants": round(len((un & crap) - given) / len(un), 4) if un else 0.0,
                            "nowhere": round(nowhere, 4)},
           "thresholds": {"min_confident": CHECK_MIN_CONFIDENT, "id_ok": CHECK_ID_OK, "given_min": CHECK_GIVEN_MIN,
                          "member_min": CHECK_MEMBER_MIN}}
    if n < CHECK_MIN_CONFIDENT:
        return {**out, "verdict": "not assessed", "reason": f"{n} spectra with a confident tag, "
                                                              f"{CHECK_MIN_CONFIDENT} needed"}
    if len(hit) / n >= CHECK_ID_OK:
        return {**out, "verdict": "ok", "reason": f"{len(hit) / n:.0%} of the {n} confident-tag spectra are identified"}
    top = max(elsewhere, key=elsewhere.get, default=None)
    if on_given >= CHECK_GIVEN_MIN:
        return {**out, "verdict": "settings likely wrong",
                "reason": f"only {len(hit) / n:.0%} of the confident-tag spectra are identified, and {on_given:.0%} of "
                          "the others place on the given database: check the label, protease, modifications and "
                          "windows"}
    named = top if top and elsewhere[top] >= CHECK_MEMBER_MIN else None
    return {**out, "verdict": "database likely wrong or incomplete", "organism": named,
            "reason": f"only {len(hit) / n:.0%} of the confident-tag spectra are identified, and {on_given:.0%} of the "
                      "others place on the given database"
                      + (f"; the {named} proteome places {elsewhere[named]:.0%} of them" if named else "")}


def _report_hints(path):
    """From a prose-auto report: the spectra searched (pseudo-spectra), the fragment window class and the
    isobaric label (fixed modifications FASTag needs)"""
    if not path:
        return {}
    rep = json.load(open(path))
    params = (rep.get("prose") or {}).get("params") or {}
    fixed = [m for m in params.get("Search:modifications:fixed", []) if re.match(r"(TMT|iTRAQ)", m)]
    ms2 = str(((rep.get("decisions") or {}).get("ms2_analyzer") or {}).get("value", "ORBITRAP"))
    return {"spectra": (rep.get("input") or {}).get("pseudo_spectra"), "fixed": fixed,
            "frag": (0.3, "Da") if ms2 in ("IT", "LIT") else (20.0, "ppm")}


def main(a, log=print):
    if not a.fasta:
        raise SystemExit("prose-auto check needs the --fasta the search used")
    if not a.fastag:
        raise SystemExit("prose-auto check needs --fastag")
    ext.require_fastag(a.fastag)
    hints = _report_hints(a.report)
    spectra = hints.get("spectra")
    if not (spectra and os.path.exists(spectra)):   # 0.1.12: a temporary ddaPASEF conversion is gone after the run
        if spectra:
            log(f"  {spectra} (from the report) no longer exists; tagging {a.run}")
        spectra = a.run
    ided, n_psm = identified(a.result, a.fdr)
    conf, placed, n_ms2 = placements(a.fastag, spectra, a.fasta, a.contaminants, hints.get("frag", (20.0, "ppm")),
                                     hints.get("fixed", ()), a.threads, log)
    res = classify(conf, ided, placed, n_ms2)
    res.update({"run": a.run, "result": a.result, "fasta": a.fasta, "fdr": a.fdr})
    os.makedirs(a.out, exist_ok=True)
    base = os.path.join(a.out, os.path.splitext(os.path.basename(a.run.rstrip("/")))[0])
    with open(base + ".check.json", "w") as fh:
        json.dump(res, fh, indent=1)
    u = res["unidentified"]
    print(f"{res['verdict']}: {res.get('reason', '')}\n"
          f"  identified {res['identified']} of {res['ms2']} MS2; confident-tag spectra {res['confident_spectra']}, "
          f"identified {res['confident_identified']}; the unidentified ones: on the given database {u['on_given']}, "
          + ", ".join(f"{m} {v}" for m, v in u["elsewhere"].items()) + f", nowhere {u['nowhere']}\n"
          f"written: {base}.check.json")
    return res

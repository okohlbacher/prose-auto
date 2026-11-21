"""acceptance/0.1.12: score.py's rules and scenarios.py's taxonomy rule on toy inputs (second review: the 120-run
minimum, run accounting, missing values, ratio conventions, PSM eligibility, units, P3 gates, Holm, §9 rules,
scenario proteomes and adjudication)."""

import importlib.util
import io
import math
import os
from contextlib import redirect_stdout

import pytest

D = os.path.join(os.path.dirname(__file__), "..", "acceptance", "0.1.12")
np = pytest.importorskip("numpy")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(D, f"{name}.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


SC = _load("score")
SN = _load("scenarios")
SC.B = 2000   # the tests need no 10,000 resamples


def _tsv(path, head, rows):
    path.write_text("\t".join(head) + "\n" + "".join("\t".join("" if v is None else str(v) for v in r) + "\n"
                                                      for r in rows))


REF = ["run", "class", "organism_verified", "label", "protease", "cys", "reporter_ions", "prec_window", "prec_unit",
       "frag_window", "frag_unit"]
RES = ["run", "arm", "status", "psm", "pep", "pg", "psm_fdp", "pep_fdp"]
DEC = ["run", "arm", "setting", "value", "unit", "status"]
CHK = ["run", "source", "arm", "kind", "check", "verdict", "psm", "seconds"]
own = lambda r, check, verdict: (r, "auto", "", "correct", check, verdict, "", "")


def _study(tmp, n, gain=1.1, fdp=(1.0, 1.0), stopped=0, missing=0, checks=(), extra_dec=()):
    """n drawn runs, all bottom-up; the first `missing` without a reference row, the next `stopped` with auto stopped"""
    runs = [f"PXD{i:06d}/r{i}.raw" for i in range(n)]
    _tsv(tmp / "heldout.tsv", ["accession", "file"], [r.split("/") for r in runs])
    _tsv(tmp / "reference.tsv", REF, [(r, "bottom-up", "yes", "label free", "Trypsin", "Carbamidomethyl (C)", "no",
                                       10, "ppm", 0.02, "Da") for r in runs[missing:]])
    rows = []
    for i, r in enumerate(runs[missing:]):
        a = ("stopped",) + ("",) * 5 if i < stopped else ("ok", int(1000 * gain) + i, 500, 100, *fdp)
        rows += [(r, "auto") + a, (r, "reference", "ok", 1000 + i, 500, 100, 1.0, 1.0),
                 (r, "defaults", "ok", 500 + i, 250, 50, 1.0, 1.0)]
    _tsv(tmp / "results.tsv", RES, rows)
    _tsv(tmp / "decisions.tsv", DEC, [(r, "auto", "protease", "Trypsin", "", "auto") for r in runs[missing:]]
         + list(extra_dec))
    _tsv(tmp / "checks.tsv", CHK, list(checks))
    _tsv(tmp / "scenarios.tsv", ["accession", "file", "kind", "databases"], [])
    return runs


def test_bounds_and_holm():
    lo, hi = SC.wilson(24, 25)
    assert 0.83 < lo < 0.85 and hi == pytest.approx(0.993, abs=0.002)
    assert SC.clopper_lower(95, 100) == pytest.approx(0.90, abs=0.01) and SC.clopper_lower(0, 10) == 0.0
    assert SC.binom_sf(10, 10, 0.5) == pytest.approx(0.5 ** 10)
    assert SC.holm({"a": 0.01, "b": 0.04, "c": 0.03}) == {"a": True, "b": False, "c": False}
    assert SC.holm({"a": 0.01, "b": 0.02, "c": 0.04}) == {"a": True, "b": True, "c": True}


def test_median_and_quantile_conventions():
    assert SC.median([0.8, 1.2]) == 1.0                 # arithmetic midpoint, not geometric
    assert SC.median([0.0, 2.0]) == 1.0 and SC.median([1.0, math.inf]) == math.inf
    assert SC.median([0.0, 1.0, math.inf]) == 1.0
    x = np.array([0.0] * 30 + [1.0] * 40 + [math.inf] * 30)
    assert SC.quantile(x, 0.025) == 0.0 and SC.quantile(x, 0.975) == math.inf and SC.quantile(x, 0.5) == 1.0


def test_ratio_rules():
    ok = lambda psm, pg: {"status": "ok", "psm": psm, "pg": pg}
    tab = {("r1", "auto"): ok(1000, 60), ("r1", "reference"): ok(1000, 40),
           ("r2", "auto"): {"status": "stopped"}, ("r2", "reference"): ok(800, 10),
           ("r3", "auto"): ok(300, 5), ("r3", "reference"): ok(50, 1),
           ("r4", "auto"): ok(20, 1), ("r4", "reference"): ok(30, 1),
           ("r5", "auto"): {"status": "missing"}, ("r5", "reference"): ok(900, 9),
           ("r6", "auto"): ok(500, 0), ("r6", "reference"): ok(500, 0)}
    r, cnt = SC.ratios(tab, {f"r{i}" for i in range(1, 7)}, "auto", "reference", "pg")   # eligibility by PSMs
    assert r == [1.5, 0.0, math.inf] and cnt["both zero"] == 1 and cnt["missing"] == 1 and cnt["both below"] == 1


def test_minimum_and_accounting(tmp_path):
    _study(tmp_path, 125, missing=6)        # 119 eligible
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    assert s["runs"]["missing"] == 6 and s["runs"]["bottom-up"] == 119
    assert s["pass"] is None and s["P1"]["pass"] is None and s["verdict"].startswith("no verdict: 119")
    _study(tmp_path, 126, missing=6)        # 120 eligible
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    assert s["verdict"] == "pass" and s["P1"]["median"] == pytest.approx(1.094, abs=0.002)


def test_p3_stopped_is_noncompliant(tmp_path):
    _study(tmp_path, 130, stopped=30)
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    assert s["P3"]["runs"] == 130 and s["P3"]["auto_share"] == pytest.approx(100 / 130)
    assert s["P3"]["relative_pass"] is False and s["pass"] is False


def test_p3_absolute_gate(tmp_path):
    _study(tmp_path, 130, fdp=(2.0, 2.0))   # auto noncompliant everywhere, the reference compliant
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    assert s["P3"]["absolute_pass"] is False and s["P3"]["pass"] is False


def test_strict_inputs(tmp_path):
    runs = _study(tmp_path, 130)
    rows = (tmp_path / "results.tsv").read_text().splitlines()
    (tmp_path / "results.tsv").write_text("\n".join(rows[:1] + [rows[1].rsplit("\t", 2)[0] + "\t\t"] + rows[2:]) + "\n")
    with pytest.raises(SC.InputError, match="psm_fdp"):        # a blank FDP is never a result
        SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    _study(tmp_path, 130, checks=[own(runs[0], "pre", "")])
    with pytest.raises(SC.InputError, match="verdict"):        # a blank verdict is never a detection
        SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    _study(tmp_path, 130, checks=[own(runs[0], "pre", "ok"), own(runs[0], "pre", "ok")])
    with pytest.raises(SC.InputError, match="duplicate"):
        SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    _study(tmp_path, 130, extra_dec=[("PXD999999/x.raw", "auto", "label", "TMT", "", "auto")])
    with pytest.raises(SC.InputError, match="not in the draw"):
        SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))


def test_units_and_absent_decisions(tmp_path):
    runs = [f"PXD{i:06d}/r{i}.raw" for i in range(130)]
    dec = [(r, "auto", "frag_window", 0.02, "ppm", "auto") for r in runs[:10]]          # 0.02 ppm is not 0.02 Da
    dec += [(r, "auto", "frag_window", 20, "ppm", "auto") for r in runs[10:20]]         # 20 ppm = 0.02 Da at m/z 1000
    _study(tmp_path, 130, extra_dec=dec)
    f = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))["S1-S5"]["frag_window"]
    assert f["correct"] == 10 and f["n"] == 130 and f["absent"] == 110


def test_checks_holm_and_missing(tmp_path):
    runs = [f"PXD{i:06d}/r{i}.raw" for i in range(130)]
    chk = [own(r, "pre", "ok") for r in runs[:120]] + [own(r, "pre", "not assessed") for r in runs[120:]]
    _study(tmp_path, 130, checks=chk)
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    e = s["S6"]["pre:correct"]
    assert e["assessed"] == 120 and e["not_assessed"] == 10 and e["hit"] == 120
    assert s["S6"]["post:correct"]["missing"] == 130 and s["S6"]["post:correct"]["pass"] is False
    # Holm always runs over all nine gates: the seven without scenarios keep p = 1 and do not pass
    gates = [k for k, v in s["S6"].items() if v.get("status") == "gate"]
    assert len(gates) == 9 and all(s["S6"][k]["pass"] is False for k in gates if k != "pre:correct")
    assert e["pass"] is True and e["p_value"] <= 0.05 / 9


def test_label_rule(tmp_path):
    runs = [f"PXD{i:06d}/r{i}.raw" for i in range(130)]
    dec = [(r, "auto", "label", "label free", "", "auto") for r in runs[1:]] + [(runs[0], "auto", "label", "TMT6plex", "", "auto")]
    _study(tmp_path, 130, extra_dec=dec)
    r = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))["rules"][
        "label: S1 >= always-label-free baseline, no isobaric call without reporter ions"]
    assert r["isobaric_without_reporters"] == 1 and r["pass"] is False


def _taxonomy(tmp_path):
    nodes = [(1, 1, "no rank"), (2759, 1, "superkingdom"), (33208, 2759, "kingdom"), (40674, 33208, "class"),
             (9443, 40674, "order"), (9604, 9443, "family"), (9605, 9604, "genus"), (9606, 9605, "species"),
             (9596, 9604, "genus"), (9598, 9596, "species"), (9526, 9443, "family"), (9539, 9526, "genus"),
             (9544, 9539, "species"), (9989, 40674, "order"), (10066, 9989, "family"), (10088, 10066, "genus"),
             (10090, 10088, "species"), (4751, 2759, "kingdom"), (33090, 2759, "kingdom"), (3701, 33090, "genus"),
             (3702, 3701, "species"), (6392, 33208, "family"), (116954, 6392, "genus"), (332521, 116954, "species")]
    names = {1: "root", 2759: "Eukaryota", 33208: "Metazoa", 40674: "Mammalia", 9443: "Primates", 9604: "Hominidae",
             9605: "Homo", 9606: "Homo sapiens", 9596: "Pan", 9598: "Pan troglodytes", 9526: "Cercopithecidae",
             9539: "Macaca", 9544: "Macaca mulatta", 9989: "Rodentia", 10066: "Muridae", 10088: "Mus",
             10090: "Mus musculus", 4751: "Fungi", 33090: "Viridiplantae", 3701: "Arabidopsis",
             3702: "Arabidopsis thaliana", 6392: "Lumbricidae", 116954: "Dendrobaena", 332521: "Dendrobaena veneta"}
    tax = tmp_path / "tax"
    tax.mkdir(exist_ok=True)
    (tax / "nodes.dmp").write_text("".join(f"{a}\t|\t{b}\t|\t{c}\t|\n" for a, b, c in nodes))
    (tax / "names.dmp").write_text("".join(f"{a}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for a, n in names.items()))
    _tsv(tmp_path / "cand.tsv", ["target", "level", "species", "taxid", "proteome", "proteins"],
         [(9606, "self", "Homo sapiens", 9606, "UP5640", 20400), (9606, "family", "Pan troglodytes", 9598, "UP2277", 23000),
          (9606, "order", "Macaca mulatta", 9544, "UP6718", 21000), (9606, "order", "Macaca fascicularis", 9541, "UPx", 21000),
          (10090, "self", "Mus musculus", 10090, "UP0589", 17000), (3702, "self", "Arabidopsis thaliana", 3702, "UP6548", 27000)])


def _run_scenarios(tmp_path, adjudicated=None):
    buf = io.StringIO()
    with redirect_stdout(buf):
        SN.main(str(tmp_path / "wrongdb.tsv"), str(tmp_path / "tax"), str(tmp_path / "cand.tsv"),
                str(adjudicated) if adjudicated else None)
    head, *rows = [l.split("\t") for l in buf.getvalue().splitlines()]
    return {(r[0], r[2]): dict(zip(head, r)) for r in rows}


def test_scenarios_taxonomy_rule(tmp_path):
    _taxonomy(tmp_path)
    _tsv(tmp_path / "wrongdb.tsv", ["stratum", "accession", "file", "organisms", "experiment"],
         [("x", "PXD1", "a.raw", "Homo sapiens (human)", "Shotgun proteomics"),
          ("x", "PXD2", "b.raw", "Homo sapiens (human); Mus musculus (mouse)", ""),
          ("x", "PXD3", "c.raw", "Fungi; Unknownia novus", ""),
          ("x", "PXD4", "d.d", "Homo sapiens", "Top-down proteomics"),
          ("x", "PXD5", "e.raw", "Dendrobaena veneta", "")])
    got = _run_scenarios(tmp_path)
    assert got[("PXD1", "related_close")]["databases"] == "Pan troglodytes"
    assert got[("PXD1", "related_close")]["common_rank"] == "family"
    assert got[("PXD1", "related_distant")]["databases"] == "Macaca fascicularis"   # tie: the lower taxon id
    assert got[("PXD1", "unrelated")]["proteomes"] == "UP6548" and got[("PXD1", "correct")]["proteomes"] == "UP5640"
    assert got[("PXD2", "incomplete")]["proteomes"] == "UP5640" and got[("PXD2", "related_close")]["databases"] == "none"
    assert got[("PXD2", "correct")]["proteomes"] == "UP5640; UP0589"
    assert got[("PXD4", "excluded")]["note"] == "top-down"
    assert "no reference proteome of its own: Dendrobaena veneta" in got[("PXD5", "excluded")]["note"]
    assert ("PXD3", "excluded") in got
    built = [r for k, r in got.items() if k[1] != "excluded" and r["databases"] != "none"]
    assert built and all(r["proteomes"] for r in built)


def test_scenarios_follow_the_adjudication(tmp_path):
    _taxonomy(tmp_path)
    _tsv(tmp_path / "wrongdb.tsv", ["stratum", "accession", "file", "organisms", "experiment"],
         [("x", "PXD1", "a.raw", "Homo sapiens", ""), ("x", "PXD2", "b.raw", "Homo sapiens", "")])
    _tsv(tmp_path / "adj.tsv", ["accession", "file", "organisms", "class"],
         [("PXD1", "a.raw", "Mus musculus", "bottom-up"), ("PXD2", "b.raw", "Homo sapiens", "top-down")])
    got = _run_scenarios(tmp_path, tmp_path / "adj.tsv")
    assert got[("PXD1", "correct")]["databases"] == "Mus musculus"
    assert got[("PXD2", "excluded")]["note"] == "class top-down (adjudication)"


def test_extract_tree(tmp_path):
    import json
    EX = _load("extract")
    out = tmp_path / "out"
    rep = {"decisions": {"label": {"value": "label free", "status": "auto"},
                         "enzyme": {"value": "Trypsin", "status": "auto"},
                         "cys_alkylation": {"value": "Carbamidomethyl (C)", "status": "review"},
                         "precursor_tolerance": {"value": {"value": 8, "unit": "ppm"}, "status": "auto"},
                         "fragment_tolerance": {"value": {"value": 0.02, "unit": "Da"}, "status": "default"},
                         "ms2_analyzer": {"value": "ORBITRAP", "status": "auto"},
                         "acquisition": {"value": "DDA", "status": "auto"}},
           "fastag_qc": {"ms2_analyzer_class": {"FTMS": 900, "ITMS": 3}},
           "database_check": {"status": "warning"}, "timings_s": {"total": 12.5}}
    ev = {"psm": {"accepted": 1200, "combined": 0.8}, "peptide": {"accepted": 700, "combined": 0.9},
          "protein": {"groups": 210}}
    a = out / "runs" / "PXD1" / "a.raw"
    for arm, files in (("auto", {"r.auto.json": rep, "eval.json": ev, "check.json": {"verdict": "ok"}}),
                       ("reference", {"eval.json": ev, "settings.json": {"protease": "Trypsin"}}),
                       ("defaults", {})):
        (a / arm).mkdir(parents=True)
        for n, x in files.items():
            (a / arm / n).write_text(json.dumps(x))
    td = out / "runs" / "PXD2" / "b.raw" / "auto"
    td.mkdir(parents=True)
    (td / "r.auto.json").write_text(json.dumps({"no_digest": {"class": "top-down"}}))
    sc = out / "scenarios" / "PXD1" / "a.raw" / "unrelated"
    sc.mkdir(parents=True)
    (sc / "r.auto.json").write_text(json.dumps({"database_check": {"status": "warning"}}))
    (sc / "check.json").write_text(json.dumps({"verdict": "database likely wrong or incomplete"}))
    EX.main(str(out), str(tmp_path / "in"))
    rows = lambda n: [l.split("\t") for l in (tmp_path / "in" / f"{n}.tsv").read_text().splitlines()[1:]]
    res = {(r[0], r[1]): r[2:] for r in rows("results")}
    assert res[("PXD1/a.raw", "auto")] == ["ok", "1200", "700", "210", "0.8", "0.9"]
    assert res[("PXD1/a.raw", "defaults")][0] == "missing" and res[("PXD2/b.raw", "auto")][0] == "stopped"
    dec = {(r[0], r[1], r[2]): r[3:] for r in rows("decisions")}
    assert dec[("PXD1/a.raw", "auto", "frag_window")] == ["0.02", "Da", "default"]
    assert dec[("PXD1/a.raw", "auto", "ms2_analyzer_qc")][0] == "ORBITRAP"
    assert dec[("PXD1/a.raw", "reference", "protease")] == ["Trypsin", "", "given"]
    assert sorted(tuple(r[:6]) for r in rows("checks")) == [
        ("PXD1/a.raw", "auto", "", "correct", "post", "ok"), ("PXD1/a.raw", "auto", "", "correct", "pre", "warning"),
        ("PXD1/a.raw", "scenario", "", "unrelated", "post", "database"),
        ("PXD1/a.raw", "scenario", "", "unrelated", "pre", "warning")]


def test_p3_exact_boundary(tmp_path, monkeypatch):
    # 120 runs, reference compliant everywhere, auto noncompliant on 6: the difference is exactly -0.05. The
    # bootstrap is pinned to the observed signed count, so the lower bound sits exactly on the boundary.
    pairs = [(True, True)] * 114 + [(False, True)] * 6
    assert SC.p3_relative_lower(pairs, np.random.default_rng(1), 500) <= -6
    assert -6 >= SC.P3_MARGIN_FRAC * 120                         # exact: the boundary itself passes
    monkeypatch.setattr(SC, "p3_relative_lower", lambda pr, rng=None, b=None: sum(int(a) - int(c) for a, c in pr))
    _study(tmp_path, 120, fdp=(1.0, 1.0))
    rows = (tmp_path / "results.tsv").read_text().splitlines()
    out = [rows[0]]
    for i, l in enumerate(rows[1:]):
        f = l.split("\t")
        if f[1] == "auto" and i // 3 < 6:
            f[6] = f[7] = "1.5000000001"                         # not rounded into compliance
        out.append("\t".join(f))
    (tmp_path / "results.tsv").write_text("\n".join(out) + "\n")
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    assert s["P3"]["auto_share"] == pytest.approx(114 / 120) and s["P3"]["difference"] == pytest.approx(-0.05)
    assert s["P3"]["relative_pass"] is True


def test_cys_rule_uninformative(tmp_path):
    runs = [f"PXD{i:06d}/r{i}.raw" for i in range(130)]
    dec = [(r, a, "cys", v, "", "auto") for r in runs[:6] for a, v in (("auto_never", "none (free Cys)"),
                                                                     ("auto_always", "Carbamidomethyl (C)"))]
    _study(tmp_path, 130, extra_dec=dec)
    rows = (tmp_path / "results.tsv").read_text()
    rows += "".join(f"{r}\t{a}\tstopped\t\t\t\t\t\n" for r in runs[:6] for a in ("auto_never", "auto_always"))
    (tmp_path / "results.tsv").write_text(rows)
    r = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))["rules"][
        "cys: -alkylation_out (auto_never) vs competition (auto_always)"]
    assert r["uninformative"] == 6 and r["runs"] == 0 and r["pass"] is None      # 0/0 is never evidence


def test_protease_rule_same_runs(tmp_path):
    runs = [f"PXD{i:06d}/r{i}.raw" for i in range(130)]
    dec = [(r, "auto", "protease", "Glu-C", "", "auto") for r in runs[:1]]
    _study(tmp_path, 130)
    lines = [l for l in (tmp_path / "decisions.tsv").read_text().splitlines() if not l.startswith(runs[0] + "\t")]
    (tmp_path / "decisions.tsv").write_text("\n".join(lines + ["\t".join(map(str, d)) for d in dec]) + "\n")
    r = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))["rules"][
        "protease: S2 >= stated-protease (SDRF arm) baseline"]
    assert r["n"] == 0 and r["without_sdrf_protease"] == 130 and r["pass"] is None   # no SDRF protease: no verdict


def test_label_vocabulary():
    for a, b in (("none", "label free"), ("no isobaric label", "none"), ("TMT10plex", "TMT10/11plex"),
                 ("TMT11plex", "TMT10/11plex"), ("TMTpro16plex", "TMTpro (16/18plex)"),
                 ("TMTpro18plex", "TMTpro (16/18plex)"), ("iTRAQ4plex", "itraq4plex"),
                 ("SILAC (Lys8/Arg10)", "silac (lys8/arg10)")):
        assert SC.canon_label(a) == SC.canon_label(b)
    for a, b in (("TMT6plex", "TMT10/11plex"), ("iTRAQ4plex", "iTRAQ8plex"), ("none", "TMT6plex"),
                 ("SILAC (Arg6)", "SILAC (Lys8/Arg10)"), ("TMT10/11plex", "TMTpro (16/18plex)")):
        assert SC.canon_label(a) != SC.canon_label(b)


def test_scenarios_keep_frozen_rows(tmp_path):
    _taxonomy(tmp_path)
    _tsv(tmp_path / "wrongdb.tsv", ["stratum", "accession", "file", "organisms", "experiment"],
         [("x", "PXD1", "a.raw", "Homo sapiens", ""), ("x", "PXD2", "b.raw", "Homo sapiens", "")])
    frozen = tmp_path / "frozen.tsv"
    frozen.write_text("accession\tfile\tkind\tdatabases\ttaxids\tproteomes\tcommon_rank\tnote\n"
                      "PXD1\ta.raw\tcorrect\tHomo sapiens\t9606\tUP5640\t\t\n"
                      "PXD1\ta.raw\trelated_close\tGorilla gorilla\t9595\tUPgor\tfamily\t\n")
    _tsv(tmp_path / "adj.tsv", ["accession", "file", "organisms", "class"],
         [("PXD1", "a.raw", "Homo sapiens (human)", "bottom-up"), ("PXD2", "b.raw", "Mus musculus", "bottom-up")])
    buf = io.StringIO()
    with redirect_stdout(buf):
        SN.main(str(tmp_path / "wrongdb.tsv"), str(tmp_path / "tax"), str(tmp_path / "cand.tsv"),
                str(tmp_path / "adj.tsv"), str(frozen))
    rows = [l.split("\t") for l in buf.getvalue().splitlines()[1:]]
    assert [r for r in rows if r[0] == "PXD1"] == [l.split("\t") for l in frozen.read_text().splitlines()[1:]]
    assert any(r[0] == "PXD2" and r[2] == "correct" and r[3] == "Mus musculus" for r in rows)   # changed: rebuilt


def test_extract_then_score(tmp_path):
    """the arm-output tree of a small study through extract.py into score.py: auto-arm specificity, a scenario run
    outside the draw (a multi-species wrong-database run), its correct control, natural failures of both arms,
    the header comparator and the full-precision FDP"""
    import json
    EX = _load("extract")
    out, runs = tmp_path / "out", [f"PXD{i:06d}/r{i}.raw" for i in range(121)]
    rep = lambda header: {"decisions": {
        "label": {"value": "label free", "status": "auto"}, "enzyme": {"value": "Trypsin", "status": "auto"},
        "cys_alkylation": {"value": "Carbamidomethyl (C)", "status": "auto"},
        "precursor_tolerance": {"value": {"value": 10, "unit": "ppm"}, "status": "auto"},
        "fragment_tolerance": {"value": {"value": 0.02, "unit": "Da"}, "status": "auto"},
        "ms2_analyzer": {"value": "ORBITRAP", "status": "auto",
                         "evidence": [{"source": "metadata", "method": "filter string", "value": header,
                                       "confidence": 0.4},
                                      {"source": "metadata", "method": "FASTag QC", "value": "ORBITRAP",
                                       "confidence": 0.5}]}},
        "fastag_qc": {"ms2_analyzer_class": {"FTMS": 900}}, "database_check": {"status": "ok"},
        "timings_s": {"total": 10.0, "dbcheck": 2.0}}
    ev = lambda psm, fdp: {"psm": {"accepted": psm, "combined": fdp}, "peptide": {"accepted": psm // 2, "combined": fdp},
                           "protein": {"groups": psm // 10}}

    def put(d, files):
        d.mkdir(parents=True, exist_ok=True)
        for n, x in files.items():
            (d / n).write_text(x if isinstance(x, str) else json.dumps(x))
    for i, r in enumerate(runs):
        acc, f = r.split("/")
        put(out / "runs" / acc / f / "auto", {"r.auto.json": rep("IT" if i == 0 else "ORBITRAP"),
                                              "eval.json": ev(1100, 1.5000001 if i == 1 else 1.0),
                                              "check.json": {"verdict": "ok"}})
        put(out / "runs" / acc / f / "reference", {"eval.json": ev(1000, 1.0), "settings.json": {"protease": "Trypsin"}})
        put(out / "runs" / acc / f / "defaults", {"eval.json": ev(500, 1.0)})
    for arm in ("sdrf", "defaults"):
        put(out / "natural" / "PXD000005" / "r5.raw" / arm, {"check.json": {"verdict": "settings likely wrong"}})
    for kind, verdict, psm in (("correct", "ok", 800), ("unrelated", "database likely wrong or incomplete", 80)):
        put(out / "scenarios" / "PXD9" / "mix.raw" / kind,
            {"r.auto.json": {"database_check": {"status": "ok" if kind == "correct" else "warning"},
                             "timings_s": {"dbcheck": 2.5}},
             "check.json": {"verdict": verdict}, "eval.json": ev(psm, 1.0), "check_seconds": "6.5"})
    EX.main(str(out), str(tmp_path / "in"))
    _tsv(tmp_path / "in" / "heldout.tsv", ["accession", "file"], [r.split("/") for r in runs])
    _tsv(tmp_path / "in" / "reference.tsv", REF, [(r, "bottom-up", "yes", "none", "Trypsin", "Carbamidomethyl (C)",
                                                   "no", 10, "ppm", 0.02, "Da") for r in runs])
    (tmp_path / "in" / "scenarios.tsv").write_text(
        "accession\tfile\tkind\tdatabases\nPXD9\tmix.raw\tcorrect\tHomo sapiens; Escherichia coli\n"
        "PXD9\tmix.raw\tunrelated\tArabidopsis thaliana\n")
    s = SC.score(str(tmp_path / "in"), str(tmp_path / "in" / "heldout.tsv"))
    assert s["S6"]["pre:correct"]["hit"] == 121 and s["S6"]["post:correct"]["hit"] == 121
    assert s["S6"]["post:unrelated"]["hit"] == 1 and s["S6"]["post:unrelated"]["psm_ratio_to_correct_median"] == 0.1
    assert s["S6"]["post:unrelated"]["seconds_median"] == 6.5 and s["S6"]["pre:unrelated"]["seconds_median"] == 2.5
    assert s["S6"]["scenario correct (controls)"] == {"checks": 2, "ok": 2}
    assert s["S6"]["post:natural"]["scenarios"] == 2 and s["S6"]["post:natural"]["hit"] == 2
    assert s["P3"]["auto_share"] == pytest.approx(120 / 121)       # 1.5000001 % is not compliant
    q = s["rules"]["FASTag QC: analyser class agrees with the header in >= 98 %"]
    assert q["n"] == 121 and q["agree"] == 120                       # the header evidence, not the reconciled value
    assert s["S1-S5"]["label"]["correct"] == 121                     # reference "none" = auto "label free"


def test_refdb_manifest(tmp_path):
    import sys
    sys.path.insert(0, D)
    RD = _load("refdb")
    _taxonomy(tmp_path)
    _tsv(tmp_path / "runs.tsv", ["accession", "file", "organisms"],
         [("PXD1", "a.raw", "Homo sapiens (human)"), ("PXD2", "b.raw", "Homo sapiens; Mus musculus"),
          ("PXD3", "c.raw", "Dendrobaena veneta"), ("PXD4", "d.raw", "Human gut metagenome")])
    buf = io.StringIO()
    with redirect_stdout(buf):
        RD.main(str(tmp_path / "runs.tsv"), str(tmp_path / "tax"), str(tmp_path / "cand.tsv"))
    rows = {l.split("\t")[0]: l.split("\t") for l in buf.getvalue().splitlines()[1:]}
    assert rows["PXD1"][2:6] == ["bottom-up", "Homo sapiens", "9606", "UP5640"]
    assert rows["PXD2"][5] == "UP5640; UP0589"
    assert rows["PXD3"][2] == "no-proteome" and rows["PXD3"][5] == ""
    assert rows["PXD4"][2] == "no-proteome"

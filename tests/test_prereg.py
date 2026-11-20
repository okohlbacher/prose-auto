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
    _tsv(tmp / "checks.tsv", ["run", "kind", "check", "verdict"], list(checks))
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
    _study(tmp_path, 130, checks=[(runs[0], "correct", "pre", "")])
    with pytest.raises(SC.InputError, match="verdict"):        # a blank verdict is never a detection
        SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    _study(tmp_path, 130, checks=[(runs[0], "correct", "pre", "ok"), (runs[0], "correct", "pre", "ok")])
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
    chk = [(r, "correct", "pre", "ok") for r in runs[:120]] + [(r, "correct", "pre", "not assessed") for r in runs[120:]]
    _study(tmp_path, 130, checks=chk)
    s = SC.score(str(tmp_path), str(tmp_path / "heldout.tsv"))
    e = s["S6"]["pre:correct"]
    assert e["assessed"] == 120 and e["not_assessed"] == 10 and e["hit"] == 120 and e["pass"] is True
    assert s["S6"]["post:correct"]["missing"] == 130 and s["S6"]["post:correct"]["pass"] is False


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
    assert sorted(map(tuple, rows("checks"))) == [("PXD1/a.raw", "correct", "post", "ok"),
                                                  ("PXD1/a.raw", "correct", "pre", "warning"),
                                                  ("PXD1/a.raw", "unrelated", "post", "database"),
                                                  ("PXD1/a.raw", "unrelated", "pre", "warning")]

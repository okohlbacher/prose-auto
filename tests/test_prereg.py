"""acceptance/0.1.12: score.py's rules and scenarios.py's taxonomy rule on toy inputs."""

import importlib.util
import io
import math
import os
import sys
from contextlib import redirect_stdout

import pytest

D = os.path.join(os.path.dirname(__file__), "..", "acceptance", "0.1.12")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(D, f"{name}.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


np = pytest.importorskip("numpy")
SC = _load("score")
SN = _load("scenarios")


def _tsv(path, head, rows):
    path.write_text("\t".join(head) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def test_bounds():
    lo, hi = SC.wilson(24, 25)
    assert 0.83 < lo < 0.85 and hi == pytest.approx(0.993, abs=0.002)
    assert SC.clopper_lower(95, 100) == pytest.approx(0.90, abs=0.01)
    assert SC.clopper_lower(0, 10) == 0.0


def test_ratio_rules_and_score(tmp_path):
    tab = {("r1", "auto"): {"status": "ok", "psm": "1000"}, ("r1", "reference"): {"status": "ok", "psm": "500"},
           ("r2", "auto"): {"status": "stopped", "psm": ""}, ("r2", "reference"): {"status": "ok", "psm": "800"},
           ("r3", "auto"): {"status": "ok", "psm": "300"}, ("r3", "reference"): {"status": "ok", "psm": "50"},
           ("r4", "auto"): {"status": "ok", "psm": "20"}, ("r4", "reference"): {"status": "ok", "psm": "30"},
           ("r5", "auto"): {"status": "missing", "psm": ""}, ("r5", "reference"): {"status": "ok", "psm": "900"}}
    r, cnt = SC.ratios(tab, {"r1", "r2", "r3", "r4", "r5"}, "auto", "reference")
    assert r == [2.0, 0.0, math.inf] and cnt == {"pairs": 3, "missing": 1, "both below": 1, "denominator below": 1}
    # a whole score on toy files: every eligible run gains 10 %, FDP equal, one check per kind
    runs = [f"r{i}" for i in range(40)]
    _tsv(tmp_path / "reference.tsv", ["run", "class", "organism_verified", "protease"],
         [(x, "bottom-up", "yes", "Trypsin") for x in runs] + [("td", "top-down", "yes", "")])
    _tsv(tmp_path / "results.tsv", ["run", "arm", "status", "psm", "pep", "pg", "psm_fdp", "pep_fdp"],
         [(x, a, "ok", n, n // 2, n // 10, 1.0, 1.0) for i, x in enumerate(runs)
          for a, n in (("auto", 1100 + i), ("reference", 1000 + i), ("defaults", 500 + i))]
         + [("td", "auto", "stopped", "", "", "", "", "")])
    _tsv(tmp_path / "decisions.tsv", ["run", "arm", "setting", "value", "status"],
         [(x, "auto", "protease", "Trypsin", "decided") for x in runs])
    _tsv(tmp_path / "checks.tsv", ["run", "kind", "check", "verdict"],
         [(x, "correct", "pre", "ok") for x in runs] + [(x, "unrelated", "post", "database") for x in runs])
    s = SC.score(str(tmp_path))
    assert s["runs"]["eligible"] == 40 and s["runs"]["top-down"] == 1 and s["S8 top-down"]["stopped"] == 1
    assert s["P1"]["median"] == pytest.approx(1.098, abs=0.002) and s["P1"]["pass"] and s["P2"]["pass"]
    assert s["P3"]["difference"] == 0 and s["P3"]["pass"] and s["pass"]
    assert s["S1-S5"]["protease"]["correct"] == 40
    assert s["S6"]["post:unrelated"]["exact_verdict"] == 40 and s["S6"]["pre:correct"]["pass"]


def test_scenarios_taxonomy_rule(tmp_path):
    # root 1 > Eukaryota 2759 > Metazoa 33208 > Mammalia 40674 (class) > Primates 9443 / Rodentia 9989 (orders)
    nodes = [(1, 1, "no rank"), (2759, 1, "superkingdom"), (33208, 2759, "kingdom"), (40674, 33208, "class"),
             (9443, 40674, "order"), (9604, 9443, "family"), (9605, 9604, "genus"), (9606, 9605, "species"),
             (9596, 9604, "genus"), (9598, 9596, "species"), (9526, 9443, "family"), (9539, 9526, "genus"),
             (9544, 9539, "species"), (9989, 40674, "order"), (10066, 9989, "family"), (10088, 10066, "genus"),
             (10090, 10088, "species"), (4751, 2759, "kingdom")]
    names = {1: "root", 2759: "Eukaryota", 33208: "Metazoa", 40674: "Mammalia", 9443: "Primates", 9604: "Hominidae",
             9605: "Homo", 9606: "Homo sapiens", 9596: "Pan", 9598: "Pan troglodytes", 9526: "Cercopithecidae",
             9539: "Macaca", 9544: "Macaca mulatta", 9989: "Rodentia", 10066: "Muridae", 10088: "Mus",
             10090: "Mus musculus", 4751: "Fungi"}
    tax = tmp_path / "tax"
    tax.mkdir()
    (tax / "nodes.dmp").write_text("".join(f"{a}\t|\t{b}\t|\t{c}\t|\n" for a, b, c in nodes))
    (tax / "names.dmp").write_text("".join(f"{a}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for a, n in names.items()))
    _tsv(tmp_path / "cand.tsv", ["target", "level", "species", "taxid", "proteome", "proteins"],
         [(9606, "self", "Homo sapiens", 9606, "UP5640", 20400), (9606, "family", "Pan troglodytes", 9598, "UP2277", 23000),
          (9606, "order", "Macaca mulatta", 9544, "UP6718", 21000), (9606, "order", "Macaca fascicularis", 9541, "UPx", 21000),
          (10090, "self", "Mus musculus", 10090, "UP0589", 17000)])
    _tsv(tmp_path / "wrongdb.tsv", ["stratum", "accession", "file", "organisms", "experiment"],
         [("x", "PXD1", "a.raw", "Homo sapiens (human)", "Shotgun proteomics"),
          ("x", "PXD2", "b.raw", "Homo sapiens (human); Mus musculus (mouse)", ""),
          ("x", "PXD3", "c.raw", "Fungi; Unknownia novus", ""),
          ("x", "PXD4", "d.d", "Homo sapiens", "Top-down proteomics")])
    buf = io.StringIO()
    with redirect_stdout(buf):
        SN.main(str(tmp_path / "wrongdb.tsv"), str(tax), str(tmp_path / "cand.tsv"))
    rows = [l.split("\t") for l in buf.getvalue().splitlines()[1:]]
    got = {(r[0], r[2]): (r[3], r[5]) for r in rows}
    assert got[("PXD1", "related_close")] == ("Pan troglodytes", "family")
    assert got[("PXD1", "related_distant")] == ("Macaca fascicularis", "order")   # tie: the lower taxon id
    assert got[("PXD1", "unrelated")][0] == "Arabidopsis thaliana"
    assert got[("PXD2", "incomplete")][0] == "Homo sapiens" and got[("PXD2", "related_close")][0] == "none"
    assert got[("PXD3", "excluded")] and got[("PXD4", "excluded")]
    assert [r[6] for r in rows if r[0] == "PXD4"] == ["top-down"]

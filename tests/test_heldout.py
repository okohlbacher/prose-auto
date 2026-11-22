"""The held-out harness of PREREG 0.1.12 (acceptance/0.1.12/heldout): the arm templates, the decisions the harness
makes from files (organism verification, the MS2-count check, plans, reference rows, the natural cohort), the
frozen-hash check and the Jev pin. No tool runs here; the tools are exercised by the harness dry run."""
import importlib.util
import io
import json
import os
import sys
import zipfile

import pytest

D = os.path.join(os.path.dirname(__file__), "..", "acceptance", "0.1.12", "heldout")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, D)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(D, f"{name}.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


A = _load("arms")
HN = _load("harness")


def test_settings_cli_and_wrong_settings():
    c = A.settings_cli({"label": "TMT6plex", "protease": "Glu-C", "cys": "Carbamidomethyl (C)",
                        "enrichment": "Phospho (STY)"}, (8.0, "ppm"), (0.02, "Da"))
    assert c[c.index("-Search:enzyme") + 1] == "glutamyl endopeptidase"
    fx, var = A._vals(c, "-Search:modifications:fixed"), A._vals(c, "-Search:modifications:variable")
    assert fx == ["Carbamidomethyl (C)", "TMT6plex (K)", "TMT6plex (N-term)"] and "Phospho (S)" in var
    assert c[c.index("-Search:fragment:mass_tolerance_unit") + 1] == "Da"
    s = A.settings_cli({"label": "SILAC (Arg6)", "cys": "none (free Cys)"})
    assert A._vals(s, "-Search:modifications:fixed") == [] and "Label:13C(6) (R)" in A._vals(s, "-Search:modifications:variable")
    assert "-Search:enzyme" not in s and "-Search:precursor:mass_tolerance_lower" not in s   # unstated: ProSE default
    # §2.3 wrong settings (v0111/postcal.py)
    assert A._vals(A.wrong_label(c), "-Search:modifications:fixed") == ["Carbamidomethyl (C)"]
    lf = A.settings_cli({"cys": "Carbamidomethyl (C)", "protease": "Trypsin"})
    assert A._vals(A.wrong_label(lf), "-Search:modifications:fixed")[-2:] == ["TMT6plex (K)", "TMT6plex (N-term)"]
    assert A._vals(A.wrong_protease(lf), "-Search:enzyme") == ["glutamyl endopeptidase"]
    assert A._vals(A.wrong_protease(c), "-Search:enzyme") == ["Trypsin"]


def test_prose_cmd_protocol():
    c = A.prose_cmd("ProSE", "x.mzML", "db.fasta", "o.idXML", ["-Search:precursor:max_charge", "35",
                                                                 "-Search:enzyme", "Trypsin"], 16)
    assert c[c.index("-Search:precursor:max_charge") + 1] == "20" and "-Search:peptide:missed_cleavages" in c
    assert c[-2:] == ["-threads", "16"] and all(x in c for x in A.PROTO)
    u = A.prose_cmd("ProSE", "x", "d", "o", ["-Search:enzyme", "unspecific cleavage"], 8)
    assert "-Search:peptide:missed_cleavages" not in u


def test_analyze_arms():
    t = {k: k for k in ("python", "panel", "fastag", "fastag_taxonomy", "diaspextract", "fileconverter", "crap")}
    auto = A.analyze_cmd(t, "r.mzML", "o", "auto", ["Homo sapiens"], "PXD1")
    assert "--jev" in auto and auto[auto.index("--accession") + 1] == "PXD1" and "--organism" in auto
    assert auto[auto.index("--sdrf") + 1] == "off" and auto[auto.index("--panel") + 1] == "panel"
    so = A.analyze_cmd(t, "r.mzML", "o", "spectra_only", ["Homo sapiens"], "PXD1")
    assert "--accession" not in so and "--jev" not in so
    tx = A.analyze_cmd(t, "r.mzML", "o", "texts", ["Homo sapiens"], "PXD1")
    assert "--accession" in tx and "--jev" not in tx
    nv = A.analyze_cmd(t, "r.d", "o", "auto_never", ["Homo sapiens"], "PXD1")
    assert nv[nv.index("--mini-search-cys") + 1] == "never" and "--fileconverter" in nv
    sc = A.analyze_cmd(t, "r.mzML", "o", "scenario", (), "PXD1", fasta=["a.fasta", "b.fasta"])
    assert "--organism" not in sc and sc.count("--fasta") == 2 and "--jev" in sc
    ck = A.check_cmd(t, "r.mzML", "r.idXML", ["db.fasta"], "rep.json", "o")
    assert ck[ck.index("--report") + 1] == "rep.json" and ck[ck.index("--fasta") + 1] == "db.fasta"


def test_best_psm_rule():
    assert A.best_psm_step("Trypsin", {"Trypsin": 1000, "Lys-C": 1040, "Glu-C": 900}) == "Trypsin"   # +4 %: kept
    assert A.best_psm_step("Trypsin", {"Trypsin": 1000, "Lys-C": 1050, "Glu-C": 1200}) == "Glu-C"
    assert A.best_psm_step("Trypsin", {"Trypsin": 0, "Lys-C": 10}) == "Lys-C"
    assert A.best_psm_step("Trypsin", {"Trypsin": 1000, "Lys-C": None}) == "Trypsin"
    assert A.subset_psms({"psm": {"accepted": 500, "combined": 1.2}}) == 0
    assert A.subset_psms({"psm": {"accepted": 500, "combined": 0.8}}) == 500 and A.subset_psms(None) is None
    assert [f for f, _ in A.BEST_PSM_FACTORS] == ["protease", "label", "cys"]
    assert len([p for p in A.BEST_PSM_FACTORS[0][1] if p not in ("Trypsin", "unspecific cleavage")]) == 5


def test_verify_verdict():
    v = HN.verify_verdict([{"stated0"}] * 90 + [{"stated0", "human"}] * 10 + [{"crap"}] * 50 + [{"mouse"}] * 20,
                          ["stated0"])
    assert v["psms"] == 170 and v["not_crap"] == 120 and abs(v["stated_share"] - 100 / 120) < 1e-9 and v["verified"]
    w = HN.verify_verdict([{"stated0"}] * 60 + [{"mouse"}] * 60, ["stated0"])
    assert not w["verified"] and w["points_to"] == "mouse"
    assert not HN.verify_verdict([{"stated0"}] * 99, ["stated0"])["verified"]        # fewer than 100 PSMs
    mix = [{"stated0"}] * 190 + [{"stated1"}] * 8 + [{"stated0", "stated1"}] * 2
    m = HN.verify_verdict(mix, ["stated0", "stated1"])
    assert not m["verified"] and m["component_unique_share"]["stated1"] == 8 / 200   # under 5 % uniquely


def test_conversion_status_and_ms2_count(tmp_path):
    assert HN.conversion_status((9458, 7572), (9458, 7572), "raw")[0] == "ok"
    assert HN.conversion_status((9400, 7514), (9458, 7572), "raw")[0] == "failed"   # partial conversion: missing
    assert HN.conversion_status((9458, 7570), (9458, 7572), "raw")[0] == "failed"   # MS2 compared on its own
    assert HN.conversion_status((100, 80), None, "wiff")[0] == "ok"
    assert HN.conversion_status((100, 80), None, "raw")[0] == "failed"
    assert HN.conversion_status(None, 1000, "raw")[0] == "failed" and HN.conversion_status((5, 0), 5, "raw")[0] == "failed"
    mz = tmp_path / "x.mzML"
    mz.write_text('<spectrum>\n<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>\n'
                  '<spectrum>\n<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>\n'
                  '<spectrum>\n<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="3"/>\n')
    assert HN.ms2_mzml(str(mz)) == 2 and HN.mzml_counts(str(mz)) == (3, 2)


def test_unpack(tmp_path):
    z = tmp_path / "a.d.zip"
    with zipfile.ZipFile(z, "w") as f:
        f.writestr("run.d/analysis.tdf", "x")
        f.writestr("run.d/analysis.tdf_bin", "y")
    out = HN.unpack(str(z), str(tmp_path / "in" / "a.d"))
    assert os.path.exists(os.path.join(out, "analysis.tdf"))
    bad = tmp_path / "b.d.zip"
    with zipfile.ZipFile(bad, "w") as f:
        f.writestr("../evil.d/analysis.tdf", "x")
    with pytest.raises(SystemExit, match="escapes"):
        HN.unpack(str(bad), str(tmp_path / "in" / "b.d"))


def test_run_plan_and_reference_rows():
    rep = {"prose_cli": ["-Search:enzyme", "Trypsin"], "input": {}, "decisions": {}}
    adj = {"label": "label free", "protease": "Trypsin", "cys": "Carbamidomethyl (C)", "enrichment": "",
           "sdrf_settings": json.dumps({"protease": "Lys-C"})}
    win = {"prec_window": 8.0, "prec_unit": "ppm", "frag_window": 20.0, "frag_unit": "ppm"}
    reports = {"auto": rep, "spectra_only": dict(rep), "texts": {"prose_cli": ["-Search:enzyme", "Lys-C"]},
               "auto_never": None}
    p = HN.run_plan("PXD1/a.raw", adj, {"proteomes": "UP000005640"}, win, reports, "a.mzML")
    assert set(p["searches"]) == {"auto", "reference", "sdrf", "defaults", "texts"}   # spectra_only: same command
    assert p["searches"]["sdrf"][p["searches"]["sdrf"].index("-Search:enzyme") + 1] == "Lys-C"
    assert p["best_psm"]["start"]["protease"] == "Trypsin" and p["proteomes"] == ["UP000005640"]
    st = HN.run_plan("k", adj, {"proteomes": "UP1"}, win, {"auto": {"no_digest": {}}, "texts": reports["texts"]},
                     "x.mzML", {"texts": "t.mzML"})   # a stopped auto: released, auto stopped, the comparators run
    assert st["stopped"] == ["auto"] and "auto" not in st["searches"]
    assert {"reference", "sdrf", "defaults", "texts"} <= set(st["searches"]) and st["spectra"] == "x.mzML"
    assert st["arm_spectra"] == {"texts": "t.mzML"}
    assert "missing" in HN.run_plan("k", adj, {"proteomes": "UP1"}, win, {"auto": None}, "x")
    a = {"PXD1/a.raw": {"class": "bottom-up", "label": "none"}, "PXD2/b.raw": {"class": "bottom-up"},
         "PXD3/c.raw": {"class": "bottom-up"}}
    ref = {"PXD1/a.raw": {"class": "bottom-up"}, "PXD2/b.raw": {"class": "no-proteome"}, "PXD3/c.raw": {}}
    rs = {r["run"]: r for r in HN.reference_rows(a, ref, {"PXD1/a.raw": {"verified": True}}, {"PXD1/a.raw": win},
                                                 {"PXD1/a.raw", "PXD2/b.raw"})}
    assert rs["PXD1/a.raw"]["organism_verified"] == "yes" and rs["PXD1/a.raw"]["prec_window"] == 8.0
    assert rs["PXD2/b.raw"]["class"] == "no-proteome" and "PXD3/c.raw" not in rs    # conversion failed: missing
    ev = lambda n: {"psm": {"accepted": n}}
    assert HN.natural_arms({"reference": ev(1000), "sdrf": ev(499), "defaults": ev(500)}) == ["sdrf"]


def test_manifest_helpers_and_spectra(tmp_path):
    m = tmp_path / "refdb.tsv"
    m.write_text("accession\tfile\tclass\torganisms\ttaxids\tproteomes\tnote\n"
                 "PXD1\ta.raw\tbottom-up\tHomo sapiens; Saccharomyces cerevisiae\t9606; 4932\tUP000005640; UP000002311\t\n"
                 "PXD2\tb.raw\tno-proteome\tFoo bar\t1\t\tno reference proteome of its own: Foo bar\n")
    assert HN.panel_entries(str(m)) == {"Homo_sapiens": "UP000005640", "Saccharomyces_cerevisiae": "UP000002311"}
    assert "UP000006548" in HN.proteome_ids([str(m)]) and "UP000002311" in HN.proteome_ids([str(m)])
    dia = {"input": {"pseudo_spectra": "p.mzML"}, "decisions": {"acquisition": {"value": "DIA"}}}
    assert HN.spectra_of(dia, "r.mzML") == "p.mzML" and HN.spectra_of({"decisions": {}}, "r.mzML") == "r.mzML"


def test_frozen_check_and_jev_pin(tmp_path, monkeypatch):
    f = tmp_path / "tool"
    f.write_text("binary")
    monkeypatch.setattr(HN, "FREEZE", str(tmp_path / "FREEZE.md"))
    (tmp_path / "FREEZE.md").write_text(f"```\n{HN.sha256(str(f))}  {f}\n```\n")
    HN.check_frozen([str(f)])
    f.write_text("changed")
    with pytest.raises(SystemExit, match="changed since the freeze"):
        HN.check_frozen([str(f)])
    with pytest.raises(SystemExit, match="records no Jev model"):
        HN.jev_env()


def test_jev_must_equal_the_frozen_model_endpoint_votes(tmp_path, monkeypatch):
    fz = tmp_path / "FREEZE.md"
    monkeypatch.setattr(HN, "FREEZE", str(fz))
    for v in ("JEV_MODEL", "JEV_URL", "PROSE_AUTO_JEV_VOTES"):
        monkeypatch.delenv(v, raising=False)
    fz.write_text("- Jev model: jev-latest\n- Jev endpoint: https://x/api\n- Jev votes: 3\n")
    with pytest.raises(SystemExit, match="not pinned"):
        HN.jev_env()
    fz.write_text("- Jev model: jev-2026-11-01\n- Jev endpoint: https://x/api\n- Jev votes: 3\n")
    monkeypatch.setenv("JEV_API_KEY", "k")
    assert HN.jev_env() == {"JEV_MODEL": "jev-2026-11-01", "JEV_URL": "https://x/api", "PROSE_AUTO_JEV_VOTES": "3"}
    monkeypatch.setenv("JEV_MODEL", "jev-2026-12-01")          # any other model, even a pinned one, stops
    with pytest.raises(SystemExit, match="differs from FREEZE.md"):
        HN.jev_env()
    monkeypatch.setenv("JEV_MODEL", "jev-2026-11-01")
    monkeypatch.setenv("PROSE_AUTO_JEV_VOTES", "1")
    with pytest.raises(SystemExit, match="differs from FREEZE.md"):
        HN.jev_env()
    monkeypatch.delenv("PROSE_AUTO_JEV_VOTES")
    monkeypatch.delenv("JEV_API_KEY")
    monkeypatch.delenv("TYPESAFE_API_KEY", raising=False)
    with pytest.raises(SystemExit, match="ssh stdin"):
        HN.jev_env()
    fz.write_text(fz.read_text() + f"{'a' * 64}  /x\n{'b' * 64}  /x\n")
    with pytest.raises(SystemExit, match="two digests"):
        HN.frozen_hashes()


def test_freeze_section_is_deterministic(tmp_path):
    FZ = _load("freeze")
    (tmp_path / "d").mkdir()
    (tmp_path / "d" / "b").write_text("2")
    (tmp_path / "d" / "a").write_text("1")
    assert FZ.tree_hash(str(tmp_path / "d")) == FZ.tree_hash(str(tmp_path / "d"))
    s = FZ.section("t", [str(tmp_path / "d")])
    assert s.index(str(tmp_path / "d" / "a")) < s.index(str(tmp_path / "d" / "b"))


def test_annotation_rule():
    """§3 (amendment 2026-11-22c): agreement after normalisation, disputed fields, the consensus row, kappa"""
    AN = _load("annotate")
    a = {"organisms": ["Homo sapiens (human)"], "class": "bottom-up", "label": "TMT11plex", "protease": "trypsin",
         "cys": "Carbamidomethyl (C)", "enrichment": "none", "reporter_ions": "yes", "acquisition": "DDA",
         "instrument_class": "orbitrap", "sdrf_settings": {"protease": "Trypsin", "prec_window": "10", "prec_unit": "ppm"}}
    b = {**a, "organisms": ["Homo sapiens"], "label": "TMT10/11plex", "protease": "Trypsin",
         "sdrf_settings": {"protease": "Trypsin", "prec_window": 10.0, "prec_unit": "ppm"}}
    assert AN.disputed(a, b) == []
    c = {**b, "cys": "none (free Cys)", "sdrf_settings": {"protease": "Lys-C"}}
    assert AN.disputed(a, c) == ["cys", "sdrf_settings.protease", "sdrf_settings.prec_window",
                                 "sdrf_settings.prec_unit"]
    row = AN.consensus(a, c)
    assert row["cys"] is None and row["label"] == "TMT10/11plex" and row["organism_disputed"] == "no"
    assert json.loads(row["sdrf_settings"]) == {}
    m = AN.consensus(a, {**b, "organisms": ["Mus musculus"]})
    assert m["organism_disputed"] == "yes" and m["organisms"] == "" and m["class"] == "bottom-up"
    assert AN.norm("protease", "Pepsin") == "?Pepsin" and AN.norm("label", "TMT0") == AN.norm("label", "TMT0 / TMT (126 only)")
    assert AN.kappa([("x", "x"), ("y", "y")]) == (1.0, 1.0)
    po, k = AN.kappa([("x", "x"), ("x", "y"), ("y", "y"), ("y", "x")])
    assert po == 0.5 and k == 0.0
    assert AN.parse('blah {"a": 1} more {"class": "bottom-up"}') == {"class": "bottom-up"}
    p = AN.prompt({"run": "PXD1/a.raw"}, "A")
    assert "TMTpro (16/18plex)" in p and "PXD1/a.raw" in p and "Do not use any tool" in p


def test_disputed_organism_is_not_verified():
    a = {"PXD1/a.raw": {"class": "bottom-up", "organism_disputed": "yes"}}
    ref = {"PXD1/a.raw": {"class": "no-proteome"}}   # refdb.py: no species named
    r = HN.reference_rows(a, ref, {"PXD1/a.raw": {"verified": True}}, {}, {"PXD1/a.raw"})[0]
    assert r["class"] == "bottom-up" and r["organism_verified"] == "no"


def test_frozen_prose_auto_guard(tmp_path, monkeypatch):
    """every stage imports the frozen prose_auto and, after the freeze, checks its tree hash against FREEZE.md"""
    FZ = _load("freeze")
    code = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
    monkeypatch.setattr(HN, "CODE", code)
    monkeypatch.setattr(HN, "FREEZE", str(tmp_path / "FREEZE.md"))
    HN.frozen_prose_auto("fasta")                                  # before the freeze: path check only
    with pytest.raises(SystemExit, match="only after the freeze"):
        HN.frozen_prose_auto("analyze")
    (tmp_path / "FREEZE.md").write_text(f"tree {FZ.tree_hash(code + '/prose_auto')}  {code}/prose_auto  (x)\n")
    HN.frozen_prose_auto("analyze")
    (tmp_path / "FREEZE.md").write_text(f"tree {'0' * 64}  {code}/prose_auto  (x)\n")
    with pytest.raises(SystemExit, match="differs from its FREEZE.md tree hash"):
        HN.frozen_prose_auto("task")
    monkeypatch.setattr(HN, "CODE", str(tmp_path))
    with pytest.raises(SystemExit, match="not the frozen"):
        HN.frozen_prose_auto("fasta")


def test_best_psm_is_its_own_array():
    import inspect
    assert "best_psm(" not in inspect.getsource(HN.task)          # array B: 6 h; best_psm can take days (B2)
    assert "best_psm(p, o)" in inspect.getsource(HN.best_psm_task)


def test_best_psm_unspecific_only_in_final(tmp_path, monkeypatch):
    """amendment 2026-11-22c: unspecific wins the protease step -> label and Cys are chosen with the best specific
    protease; the final full search is unspecific"""
    PSMS = {"Trypsin": 800, "Lys-C": 700, "Arg-C": 100, "Chymotrypsin": 0, "Glu-C": 0, "Asp-N": 0,
            "unspecific cleavage": 1000}
    calls = []

    def fake_search(spectra, db, cli, o, name):
        enz = A._vals(cli, "-Search:enzyme")[0]
        fixed = A._vals(cli, "-Search:modifications:fixed")
        calls.append((name, enz, o))
        n = {v: k for k, v in __import__("prose_auto.prose_params", fromlist=["x"]).ENZYME_NAMES.items()}
        base = PSMS[n.get(enz, enz)]
        return {"psm": {"accepted": base + (100 if "TMT6plex (K)" in fixed else 0), "combined": 0.5}}

    monkeypatch.setattr(HN, "search", fake_search)
    monkeypatch.setattr(HN, "subset_mzml", lambda s, out: out)
    monkeypatch.setattr(HN, "X", str(tmp_path))
    p = {"run": "PXD1/a.raw", "spectra": "a.mzML", "db": "db.fasta", "best_psm": {
        "start": {"protease": "Trypsin", "label": "label free", "cys": "Carbamidomethyl (C)"},
        "enrichment": "none", "prec": [10, "ppm"], "frag": [20, "ppm"]}}
    HN.best_psm(p, str(tmp_path / "bp"))
    j = json.load(open(tmp_path / "bp" / "best_psm.json"))
    assert j["steps"][0]["to"] == "unspecific cleavage" and j["steps"][0]["label_cys_under"] == "Trypsin"
    later = [c for c in calls if c[0] == "subset" and ("/label/" in c[2] or "/cys/" in c[2])]
    assert later and all(c[1] == "Trypsin" for c in later)
    assert j["result"] == {"protease": "unspecific cleavage", "label": "TMT6plex", "cys": "Carbamidomethyl (C)"}
    assert calls[-1][0] == "run" and calls[-1][1] == "unspecific cleavage"
    assert A.best_specific("Trypsin", PSMS) == "Trypsin"


def test_stopped_auto_is_released_and_scored_stopped(tmp_path, monkeypatch):
    """review item 1: plan() releases a stopped run; task() writes auto status=stopped (with its report) and runs the
    comparators; extract.py then reads auto as stopped (0 PSMs)"""
    H = tmp_path / "H"
    monkeypatch.setattr(HN, "H", str(H))
    k = "PXD1/a.raw"
    (H / "runs" / k / "auto").mkdir(parents=True)
    (H / "runs" / k / "auto" / "a.auto.json").write_text(json.dumps({"no_digest": {"why": "top-down"}}))
    p = {"run": k, "spectra": "x.mzML", "stopped": ["auto"], "searches": {"reference": ["-Search:enzyme", "Trypsin"]},
         "settings": {}, "reports": {}, "db": "db.fasta", "proteomes": [], "arm_spectra": {}}
    (H / "plan").mkdir()
    (H / "plan" / "release.txt").write_text(k + "\n")
    (H / "plan" / "PXD1__a.raw.json").write_text(json.dumps(p))
    for f in ("check_frozen", "check_inputs", "child_import_origin"):
        monkeypatch.setattr(HN, f, lambda *a, **kw: None)
    ran = []
    monkeypatch.setattr(HN, "search", lambda sp, db, cli, o, n: ran.append(o) or {"psm": {"accepted": 500}})
    monkeypatch.setattr(HN, "env", lambda *a: {})
    HN.task(1)
    auto = H / "out" / "runs" / k / "auto"
    assert (auto / "status").read_text().strip() == "stopped" and (auto / "a.auto.json").exists()
    assert ran == [str(H / "out" / "runs" / k / "reference")]
    EX = importlib.util.spec_from_file_location("extract", os.path.join(D, "..", "extract.py"))
    ex = importlib.util.module_from_spec(EX)
    EX.loader.exec_module(ex)
    assert ex.result(str(auto))[0] == "stopped"


def test_rerun_once_then_missing(tmp_path, monkeypatch):
    """review item 6: a failed unit is rerun once (the stage exits 1 in between), then it is missing and names its
    step; a success on the rerun counts"""
    monkeypatch.setattr(HN, "FAILED", [])
    d = str(tmp_path / "u")
    assert HN.run_unit(d, "search", lambda: False) == "retry" and HN.FAILED
    assert HN.run_unit(d, "search", lambda: False) == "missing"
    assert open(f"{d}/search.failed").read().strip() == "search"
    assert HN.run_unit(d, "search", lambda: True) == "missing"            # never a third attempt
    d2 = str(tmp_path / "v")
    HN.run_unit(d2, "check", lambda: False)
    assert HN.run_unit(d2, "check", lambda: True) == "ok"
    def boom():
        raise RuntimeError("killed")
    assert HN.run_unit(str(tmp_path / "w"), "verify", boom) == "retry"   # an exception is a failed attempt


def test_failed_search_is_rerun_then_missing_and_best_psm_waits(tmp_path, monkeypatch):
    """review items 6 and 12: search() reruns a failed ProSE once, then status=missing; best_psm selects only when
    every candidate has a result and is missing when one is out of attempts"""
    monkeypatch.setattr(HN, "FAILED", [])
    monkeypatch.setattr(HN, "prose_run", lambda cmd, log: 1)
    o = str(tmp_path / "s")
    assert HN.search("x.mzML", "db", ["-Search:enzyme", "Trypsin"], o, "run") is None
    assert not os.path.exists(f"{o}/status") and HN.FAILED
    assert HN.search("x.mzML", "db", ["-Search:enzyme", "Trypsin"], o, "run") is None
    assert open(f"{o}/status").read().strip() == "missing"
    calls = {"n": 0}

    def fake(sp, db, cli, d, name):
        calls["n"] += 1
        if "Lys_C" in d:
            os.makedirs(d, exist_ok=True)
            open(f"{d}/status", "w").write("missing\n")
            return None
        return {"psm": {"accepted": 100, "combined": 0.5}}
    monkeypatch.setattr(HN, "search", fake)
    monkeypatch.setattr(HN, "subset_mzml", lambda s, out: out)
    monkeypatch.setattr(HN, "X", str(tmp_path))
    p = {"run": "PXD1/a.raw", "spectra": "a.mzML", "db": "db", "best_psm": {
        "start": {"protease": "Trypsin", "label": "label free", "cys": "Carbamidomethyl (C)"},
        "enrichment": "none", "prec": None, "frag": None}}
    bp = str(tmp_path / "bp")
    assert HN.best_psm(p, bp) is None
    assert open(f"{bp}/status").read().strip() == "missing" and "Lys-C" in open(f"{bp}/best_psm.failed").read()
    assert not os.path.exists(f"{bp}/best_psm.json")                    # no selection from an unresolved ladder


def test_search_reuses_only_the_same_command(tmp_path, monkeypatch):
    monkeypatch.setattr(HN, "FAILED", [])
    monkeypatch.setattr(HN, "sha256", lambda p: "frozen")
    o = tmp_path / "s"
    o.mkdir()
    (o / "eval.json").write_text(json.dumps({"psm": {"accepted": 7}, "entrap_sha256": "frozen"}))
    cmd = A.prose_cmd(HN.TOOLS["prose"], "x.mzML", "db", str(o / "run.idXML"), ["-Search:enzyme", "Trypsin"], 99)
    (o / "cmd.json").write_text(json.dumps(cmd))                         # other thread count: the same search
    assert HN.search("x.mzML", "db", ["-Search:enzyme", "Trypsin"], str(o), "run")["psm"]["accepted"] == 7
    # counted by another entrap.py: counted again from the idXML, never searched again
    (o / "eval.json").write_text(json.dumps({"psm": {"accepted": 7}, "entrap_sha256": "pre22c", "seconds": 5}))
    (o / "run.idXML").write_text("")
    monkeypatch.setattr(HN, "prose_run", lambda cmd, log: pytest.fail("searched again"))
    monkeypatch.setattr(HN, "evaluate", lambda js, db, idx, sec: open(js, "w").write(
        json.dumps({"psm": {"accepted": 8}, "seconds": sec, "entrap_sha256": "frozen"})))
    e = HN.search("x.mzML", "db", ["-Search:enzyme", "Trypsin"], str(o), "run")
    assert e["psm"]["accepted"] == 8 and e["seconds"] == 5 and (o / "eval.json.entrap_old").exists()
    monkeypatch.setattr(HN, "prose_run", lambda cmd, log: 1)
    assert HN.search("x.mzML", "db", ["-Search:enzyme", "Lys-C"], str(o), "run") is None   # another command
    assert (o / "eval.json.stale").exists()


def test_stage_order_gates(tmp_path, monkeypatch):
    """review item 16: analyze and plan refuse before the seal; the seal refuses while conversions are open; the
    top-up refuses once an arm has run"""
    H = tmp_path / "H"
    monkeypatch.setattr(HN, "H", str(H))
    monkeypatch.setattr(HN, "FREEZE", str(tmp_path / "FREEZE.md"))
    monkeypatch.setattr(HN, "CONV", str(tmp_path / "conversion.tsv"))
    (tmp_path / "FREEZE.md").write_text("# FREEZE\n")
    with pytest.raises(SystemExit, match="not sealed"):
        HN.analyze(1)
    with pytest.raises(SystemExit, match="not sealed"):
        HN.plan()
    (tmp_path / "conversion.tsv").write_text("run\tstatus\tattempts\nPXD1/a.raw\tfailed\t1\n")
    with pytest.raises(SystemExit, match="due a rerun"):
        HN.seal()
    (H / "runs" / "PXD1" / "a.raw" / "auto").mkdir(parents=True)
    (H / "runs" / "PXD1" / "a.raw" / "auto" / "a.auto.json").write_text("{}")
    with pytest.raises(SystemExit, match="an arm has run"):
        HN.topup_count()
    DT = importlib.util.spec_from_file_location("draw_topup", os.path.join(D, "..", "draw_topup.py"))
    dt = importlib.util.module_from_spec(DT)
    DT.loader.exec_module(dt)
    with pytest.raises(SystemExit, match="refused"):
        dt.batch(1, str(H / "annotation" / "topup"))
    (tmp_path / "FREEZE.md").write_text("# FREEZE\n## Step 4: the annotation sealed (x)\n")
    assert HN.sealed()


def _entrap():
    f = os.path.join(HN.T, "entrap.py")
    if not os.path.exists(f):
        pytest.skip("entrap.py (corpus protocol) not here")
    spec = importlib.util.spec_from_file_location("entrap", f)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_entrapment_non_tryptic_search_by_protein_and_own_enzyme(tmp_path):
    """amendment 2026-11-22c: a Glu-C search's PSMs are classified by their source proteins (Trypsin/P peptide sets
    called all of them unclassified) and r is that of Glu-C; unspecific has r from its distinct windows"""
    E = _entrap()
    src, db = tmp_path / "src.fasta", tmp_path / "db.fasta"
    src.write_text(">sp|P1|A target one\nMSTAGHLVNPWADKLMNPQSTWVYHGLKSTAPLMNGHE\n"
                   ">sp|P2|B target two\nMGHTPLKNWYSAVDRQLMPTHGSWNYEKLPTGAMNWE\n")
    E.makedb(str(src), str(db), 0)
    ent = {h.split()[0]: s for h, s in E.read_fasta(str(db)) if h.startswith("ENTRAP_")}
    tgt = {h.split()[0]: s for h, s in E.read_fasta(str(db)) if not h.startswith(("ENTRAP_", "DECOY_"))}

    def idxml(enzyme, path):
        prots = [("PH_0", "sp|P1|A"), ("PH_1", "sp|P2|B"), ("PH_2", "ENTRAP_sp|P1|A"), ("PH_3", "DECOY_sp|P1|A")]
        hits = [("0.001", tgt["sp|P1|A"][14:30], "PH_0", "14", "target"),          # target
                ("0.002", tgt["sp|P2|B"][0:15], "PH_1", "0", "target"),            # target, protein N-term
                ("0.003", ent["ENTRAP_sp|P1|A"][3:20], "PH_2", "3", "target"),     # entrapment only
                ("0.004", ent["ENTRAP_sp|P1|A"][3:20], "PH_0 PH_2", "3 3", "target"),  # shared: target
                ("0.005", "KDAWPNVLHGATS", "PH_3", "1", "decoy")]
        body = "".join(f'<PeptideIdentification score_type="q-value">\n<PeptideHit score="{q}" sequence="{s}" charge="2" '
                       f'start="{st}" protein_refs="{r}" >\n<UserParam type="string" name="target_decoy" value="{td}"/>'
                       f'\n<UserParam type="float" name="percolator_score" value="{9 - i}"/>\n</PeptideHit>\n'
                       f'</PeptideIdentification>\n' for i, (q, s, r, st, td) in enumerate(hits))
        path.write_text(f'<IdXML>\n<SearchParameters id="SP_0" db="{db}" enzyme="{enzyme}" missed_cleavages="2" >\n'
                        '</SearchParameters>\n<ProteinIdentification>\n' +
                        "".join(f'<ProteinHit id="{i}" accession="{a}" score="0.0" sequence="" >\n</ProteinHit>\n'
                                for i, a in prots) + f'</ProteinIdentification>\n{body}</IdXML>\n')
        return str(path)

    out = tmp_path / "e.json"
    E.evaluate(str(out), str(db), [idxml("glutamyl endopeptidase", tmp_path / "gluc.idXML")])
    ev = json.loads(out.read_text())
    assert ev["enzyme"] == "glutamyl endopeptidase" and ev["classification"] == "protein"
    assert (ev["psm"]["N_T"], ev["psm"]["N_E"], ev["psm"]["N_U"]) == (3, 1, 0)
    assert ev["psm_nterm"]["N_T"] == 1
    r_glu = E.ratio(list(tgt.values()), list(ent.values()), "glutamyl endopeptidase", 2)
    assert ev["r"] == pytest.approx(r_glu[0]) and ev["r"] != pytest.approx(E.ratio(list(tgt.values()), list(ent.values()), "trypsin/p", 2)[0])
    # the old Trypsin/P sets know none of these Glu-C peptides
    t, _ = E.peptides(tgt.values())
    assert not any(E.strip(tgt["sp|P1|A"][14:30]) == p for p in t)
    E.evaluate(str(out), str(db), [idxml("unspecific cleavage", tmp_path / "unspec.idXML")])
    ev = json.loads(out.read_text())
    r_un = E.ratio(list(tgt.values()), list(ent.values()), "unspecific cleavage", 0)
    assert ev["enzyme"] == "unspecific cleavage" and ev["r"] == pytest.approx(r_un[0]) and ev["psm"]["N_E"] == 1
    with pytest.raises(SystemExit):
        E.evaluate(str(out), str(db), [idxml("pepsin", tmp_path / "x.idXML")])

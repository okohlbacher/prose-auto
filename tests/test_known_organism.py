"""0.1.11: the organism is given by default; the database check; FASTag's Cys calls and run QC; the Cys fallback."""

import argparse
import os

import pytest

from prose_auto import cli
from prose_auto import external as ext
from prose_auto.evidence import Evidence


def test_organism_required_names_the_pride_organisms(monkeypatch):
    ev = [Evidence("organism", "Glycine max", 0.6, "repository", "pride", "", {"taxid": 3847}),
          Evidence("organism", "Acyrthosiphon pisum", 0.6, "repository", "pride", "", {"taxid": 7029})]
    monkeypatch.setattr(cli.repo, "pride_evidence", lambda acc, log=print: (ev, {}))
    with pytest.raises(SystemExit) as e:
        cli.main(["analyze", "run.mzML", "--accession", "PXD039269"])
    msg = str(e.value)
    assert "PXD039269 states: Glycine max; Acyrthosiphon pisum" in msg and '--organism "Glycine max"' in msg
    with pytest.raises(SystemExit, match="needs the organism"):
        cli.main(["analyze", "run.mzML"])
    with pytest.raises(SystemExit):   # argparse: the old --fast and --infer-organism exclude each other
        cli.main(["analyze", "run.mzML", "--fast", "--infer-organism"])


def test_fasta_is_repeatable(tmp_path):
    a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
    a.write_text(">p1\nPEPTIDEK\n")
    b.write_text(">p2\nKEDITPEP")   # no final newline: the union still separates the records
    ns = argparse.Namespace(fasta=[str(a), str(b)], panel=None, out=str(tmp_path / "o"))
    cli._fasta_list(ns)
    assert ns.fastas == [str(a), str(b)] and open(ns.fasta).read() == ">p1\nPEPTIDEK\n>p2\nKEDITPEP\n"
    ns = argparse.Namespace(fasta=[str(a)], panel=None, out=str(tmp_path))
    cli._fasta_list(ns)
    assert ns.fasta == str(a)


def _tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def test_exact_placements_and_panel_labels(tmp_path):
    _tsv(tmp_path / "tags.tsv", ["spectrum", "tag", "evalue"],
         [("s1", "PEPT", 1e-6), ("s2", "LLAG", 1e-5), ("s3", "AAAA", 0.01), ("s4", "GGTK", 1e-7)])
    _tsv(tmp_path / "recon.tsv", ["spectrum", "tag", "protein", "nterm_match", "cterm_match"],
         [("s1", "PEPT", "dbcheck_human|sp|P1|X", 1, 1), ("s1", "PEPT", "sp|Q1|given", 1, 1),
          ("s2", "LLAG", "dbcheck_mouse|sp|P2|Y", 1, 0),          # one flank off: not exact
          ("s3", "AAAA", "sp|Q3|given", 1, 1),                    # not a confident tag
          ("s4", "GGTK", "dbcheck_mouse|sp|P4|Z", 1, 1)])
    conf, placed = ext.exact_placements(str(tmp_path), flanks=2)
    assert conf == {"s1", "s2", "s4"}
    assert placed == {"": {"s1"}, "dbcheck_human": {"s1"}, "dbcheck_mouse": {"s4"}}
    # the default (calibrated): one agreeing flank is enough
    assert ext.exact_placements(str(tmp_path))[1]["dbcheck_mouse"] == {"s2", "s4"}


def test_dbcheck_panel_labels_every_member(tmp_path):
    import gzip
    share = tmp_path / "bin" / "share-FASTag-alkylation"
    share.mkdir(parents=True)
    for n in ("human", "mouse", "yeast", "ecoli"):
        with gzip.open(share / f"{n}.fasta.gz", "wt") as fh:
            fh.write(f">sp|{n}1|X\nPEPTIDEK\n")
    crap = tmp_path / "crap.fasta"
    crap.write_text(">sp|K1|KRT\nKERATINK\n")
    exe = tmp_path / "bin" / "FASTag"
    exe.write_text("")
    path, names = ext.dbcheck_panel(str(exe), str(crap), str(tmp_path / "cache"))
    assert names == ["human", "mouse", "yeast", "ecoli", "contaminants"]
    heads = [l for l in open(path) if l.startswith(">")]
    assert heads[0] == ">dbcheck_human|sp|human1|X\n" and heads[-1] == ">dbcheck_contaminants|sp|K1|KRT\n"


def test_database_check_policy():
    conf = {f"s{i}" for i in range(1000)}
    s = lambda a, b: {f"s{i}" for i in range(a, b)}
    crap = s(990, 1000)
    # the right database: most confident spectra place, the panel adds nothing of its own
    ok = ext.database_check(conf, s(0, 500), {"human": s(0, 480), "mouse": s(0, 300), "contaminants": crap})
    assert ok["status"] == "ok" and ok["given_rate"] == 0.51 and ok["panel_own"]["human"] == 0
    # a wrong database (aphid for a soybean run): few placements
    low = ext.database_check(conf, s(0, 20), {"human": s(0, 20), "contaminants": crap})
    assert low["status"] == "warning" and "only 3.0%" in low["reasons"][0]
    # an incomplete one (human alone for a human/yeast/E. coli mix): yeast and E. coli have spectra of their own
    mix = ext.database_check(conf, s(0, 400), {"human": s(0, 390), "yeast": s(400, 600), "ecoli": s(600, 630)})
    assert mix["status"] == "warning" and len(mix["reasons"]) == 1 and "yeast proteome places 20.0%" in mix["reasons"][0]
    # the bar rises with the given share: 4 % of their own against a given share of 60 % is no warning, 8 % is
    assert ext.database_check(conf, s(0, 600), {"mouse": s(600, 640)})["status"] == "ok"
    assert ext.database_check(conf, s(0, 600), {"mouse": s(600, 680)})["reasons"][0].endswith("(below 6.0% expected)")
    # a member the user gave does not warn (their yeast FASTA may be another strain)
    assert ext.database_check(conf, s(0, 400), {"yeast": s(400, 600)}, organisms={"yeast"})["status"] == "ok"
    assert ext.database_check(set(list(conf)[:50]), set(), {})["status"] == "not assessed"


def _rep(rows):
    return {(s, k): v for s, k, v in rows}


def test_fastag_alkylation_and_qc_reports():
    alk = _rep([("alkylation_call", "status", "call"), ("alkylation_call", "reagent", "Propionamide"),
                ("alkylation_call", "share", "0.91")])
    alkt = _rep([("alkylation_tags_call", "status", "no_call"), ("alkylation_tags_call", "reagent", "-")])
    ev = ext.alkylation_evidence(alk, alkt)
    assert [(e.value, e.method, e.source) for e in ev] == [("Propionamide (C)", "FASTag alkylation (placements)",
                                                             "search")]
    none = _rep([("alkylation_tags_call", "status", "call"), ("alkylation_tags_call", "reagent", "unmodified"),
                 ("alkylation_tags_call", "share", "0.7")])
    assert ext.alkylation_evidence({}, none)[0].value == "none (free Cys)"
    qc = _rep([("qc", "count_ms1", "100"), ("qc", "count_ms1_FTMS", "100"), ("qc", "count_ms2", "1000"),
               ("qc", "count_ms2_ITMS", "990"), ("qc", "count_ms2_FTMS", "10"), ("qc", "resolution_ITMS_ms2", "n/a"),
               ("qc", "precursor_charge_1", "2"), ("qc", "precursor_charge_2", "600"),
               ("qc", "precursor_charge_2_share", "0.6"), ("qc", "precursor_charge_3", "390"),
               ("qc", "precursor_charge_4", "8"), ("qc", "reporter_compatible_kits", "tmt16plex,tmt18plex"),
               ("qc", "tag_yield", "0.41")])
    ev, summ = ext.qc_evidence(qc)
    got = {(e.parameter, str(e.value)) for e in ev}
    assert ("ms2_analyzer", "IT") in got and ("ms1_analyzer", "ORBITRAP") in got
    assert ("charge_range", "[2, 4]") in got and ("label", "TMTpro (16/18plex)") in got
    assert summ["charge_histogram"] == {1: 2, 2: 600, 3: 390, 4: 8} and summ["tag_yield"] == 0.41
    assert all(e.confidence <= 0.5 for e in ev)   # next to the header, never above it


def test_cys_fallback():
    c = cli._Run(argparse.Namespace(), print, "x", "x.mzML")
    assert cli._cys_fallback(c) == "FASTag and the texts make no Cys call"
    c.ev = [Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.8, "data", "FASTag alkylation (placements)")]
    assert cli._cys_fallback(c) is None                                     # FASTag calls, the texts are silent
    c.ev.append(Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.55, "text", "keywords"))
    assert cli._cys_fallback(c) is None                                     # they agree
    c.ev[1] = Evidence("cys_alkylation", "Propionamide (C)", 0.55, "text", "keywords")
    assert cli._cys_fallback(c).startswith("FASTag calls Carbamidomethyl (C), the texts say Propionamide")
    c.ev = [c.ev[1]]
    assert cli._cys_fallback(c) is None                                     # the texts alone decide
    c.ev = [Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.6, "data", "FASTag alkylation (tags)"),
            Evidence("cys_alkylation", "Nethylmaleimide (C)", 0.42, "llm", "jev:jev-latest+data")]
    assert "the texts say Nethylmaleimide" in cli._cys_fallback(c)          # Jev, read last, counts as a text
    c.ev = [Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.8, "data", "FASTag alkylation (placements)"),
            Evidence("cys_alkylation", "none (free Cys)", 0.6, "data", "FASTag alkylation (tags)")]
    assert "disagree" in cli._cys_fallback(c)


def test_post_search_check(tmp_path):
    from prose_auto import check as CK
    idx = tmp_path / "r.idXML"
    pid = lambda ref, q, td: (f'<PeptideIdentification score_type="q-value" spectrum_reference="{ref}">'
                              f'<PeptideHit score="{q}" sequence="PEPTIDEK"><UserParam type="string" '
                              f'name="target_decoy" value="{td}"/></PeptideHit></PeptideIdentification>')
    idx.write_text("<IdXML>" + pid("scan=1", 0.001, "target") + pid("scan=2", 0.2, "target")
                   + pid("scan=3", 0.001, "decoy") + "</IdXML>")
    assert CK.identified(str(idx)) == ({"scan=1"}, 3)
    raw = tmp_path / "raw.idXML"   # ProSE without rescoring (fewer than 100 PSMs): search scores, no q-values
    raw.write_text(pid("scan=1", 0.001, "target").replace("q-value", "ln(hyperscore)") * 2)
    assert CK.identified(str(raw)) == (set(), 2)
    conf = {f"s{i}" for i in range(1000)}
    s = lambda a, b: {f"s{i}" for i in range(a, b)}
    assert CK.classify(conf, s(0, 600), {"given": s(0, 800)}, 5000)["verdict"] == "ok"
    # the right database, wrong settings: the unidentified confident-tag spectra place on it
    r = CK.classify(conf, s(0, 50), {"given": s(0, 700), "human": s(0, 650)}, 5000)
    assert r["verdict"] == "settings likely wrong" and r["unidentified"]["on_given"] > 0.6
    # a wrong database: they place on a panel organism
    r = CK.classify(conf, s(0, 50), {"given": s(0, 60), "mouse": s(0, 700)}, 5000)
    assert r["verdict"] == "database likely wrong or incomplete" and r["organism"] == "mouse"
    assert CK.classify(s(0, 20), set(), {}, 100)["verdict"] == "not assessed"

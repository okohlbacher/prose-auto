"""0.1.12: regression tests for the eight defects of the 0.1.11 review (v0111/review_0111.md)."""

import argparse
import gzip
import io
import json
import os

import pytest

from prose_auto import check as CK
from prose_auto import cli
from prose_auto import external as ext
from prose_auto.evidence import Evidence

ROOT = os.environ.get("PROSE_AUTO_OPENMS_ROOT", "")
MZTAB = os.path.join(ROOT, "src", "tests", "topp", "MzTabExporter_1_output.mzTab")
HERE = os.path.dirname(__file__)


def _fastag_bundle(d):
    share = d / "share-FASTag-alkylation"
    share.mkdir(parents=True)
    for n in ("human", "mouse", "yeast", "ecoli"):
        with gzip.open(share / f"{n}.fasta.gz", "wt") as fh:
            fh.write(f">sp|{n}1|X\nPEPTIDEK\n")
    exe = d / "FASTag"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    return exe


# 1: --fastag FASTag (a bare name on PATH) from another working directory finds the bundled proteomes
def test_bare_fastag_name_on_path(tmp_path, monkeypatch):
    _fastag_bundle(tmp_path / "bin")
    elsewhere = tmp_path / "work"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    monkeypatch.setenv("PATH", f"{tmp_path / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}")
    assert ext.fastag_share("FASTag", "alkylation") == str(tmp_path / "bin" / "share-FASTag-alkylation")
    path, names = ext.dbcheck_panel("FASTag", None, str(tmp_path / "cache"))
    assert path and names == ["human", "mouse", "yeast", "ecoli"]


# 2: two --fasta with the same file name in different folders are both searched
def test_two_fastas_same_basename(tmp_path):
    (tmp_path / "A").mkdir()
    (tmp_path / "B").mkdir()
    a, b = tmp_path / "A" / "proteome.fasta", tmp_path / "B" / "proteome.fasta"
    a.write_text(">pA\nPEPTIDEK\n")
    b.write_text(">pB\nKEDITPEP\n")
    ns = argparse.Namespace(fasta=[str(a), str(b)], panel=None, out=str(tmp_path / "o"))
    cli._fasta_list(ns)
    heads = [l for l in open(ns.fasta) if l.startswith(">")]
    assert heads == [">pA\n", ">pB\n"]
    names = cli._unique_names([str(a), str(b)])
    assert list(names.values()) == [str(a), str(b)] and len(set(names)) == 2


# 3: prose-auto check separates FASTA files without a final newline
def test_check_fasta_without_final_newline(tmp_path):
    a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
    a.write_text(">a\nPEPTIDEK")   # no final newline, not the last file
    b.write_text(">b\nKEDITPEP\n")
    dst = tmp_path / "check.fasta"
    CK._labelled([("given", str(a)), ("given", str(b))], str(dst))
    assert [l for l in open(dst) if l.startswith(">")] == [">dbcheck_given|a\n", ">dbcheck_given|b\n"]


# 4: a weak text agreeing with FASTag does not hide a strong text (or Jev) contradicting it
def test_cys_fallback_mixed_statements():
    c = cli._Run(argparse.Namespace(), print, "x", "x.mzML")
    fastag = Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.86, "search", "FASTag alkylation (placements)")
    weak = Evidence("cys_alkylation", "Carbamidomethyl (C)", 0.55, "text", "keywords")
    strong = Evidence("cys_alkylation", "Nethylmaleimide (C)", 0.90, "text", "methods")
    c.ev = [fastag, weak, strong]                    # before Jev
    assert "the texts say Nethylmaleimide (C)" in cli._cys_fallback(c)
    jev = Evidence("cys_alkylation", "Nethylmaleimide (C)", 0.6, "llm", "jev:jev-latest")
    c.ev = [fastag, weak, jev]                       # after Jev: Jev's answer contradicts
    assert "Nethylmaleimide" in cli._cys_fallback(c)
    c.ev = [fastag, weak]                            # agreement alone: no competition
    assert cli._cys_fallback(c) is None


def _psm_targets(path):
    """the target spectra at q <= 0.01 of an OpenMS-exported mzTab, read independently of check._mztab"""
    head, col, out = None, None, set()
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "MTD" and f[1].startswith("psm_search_engine_score[") and "q-value" in f[2]:
            col = "search_engine_score[" + f[1].split("[")[1]
        elif f[0] == "PSH":
            head = f
        elif f[0] == "PSM":
            r = dict(zip(head, f))
            if float(r[col]) <= 0.01 and r["opt_global_cv_MS:1002217_decoy_peptide"] != "1":
                out.add(r["spectra_ref"].split(":", 1)[1])
    return out


# 5: a real OpenMS-exported mzTab: the PSM q-value is search_engine_score[n], declared in the metadata
def test_openms_mztab_qvalues():
    small = os.path.join(HERE, "data", "MzTabExporter_1_head.mzTab")
    ided, n = CK.identified(small)
    assert n == 40 and ided == _psm_targets(small) and len(ided) > 30


@pytest.mark.skipif(not os.path.exists(MZTAB), reason="OpenMS test data not found; set PROSE_AUTO_OPENMS_ROOT")
def test_openms_mztab_full():
    ided, n = CK.identified(MZTAB)
    assert n == 1620 and ided == _psm_targets(MZTAB) and len(ided) > 1500


# 6: the union and panel caches are rebuilt when a source file changes, at a new path (R3: immutable paths)
def test_cache_invalidation(tmp_path):
    a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
    a.write_text(">p1\nPEPTIDEK\n")
    b.write_text(">p2\nKEDITPEP\n")
    extra = str(tmp_path / "extra")
    u = cli._union_fasta({"a": str(a), "b": str(b)}, extra)
    b.write_text(">p3\nNEWSEQK\n")
    u2 = cli._union_fasta({"a": str(a), "b": str(b)}, extra)
    assert u2 != u and ">p3\n" in open(u2).read() and ">p2\n" in open(u).read()   # the old file is never changed
    exe = _fastag_bundle(tmp_path / "bin")
    crap = tmp_path / "crap.fasta"
    crap.write_text(">K1\nKERATINK\n")
    p, _ = ext.dbcheck_panel(str(exe), str(crap), str(tmp_path / "cache"))
    crap.write_text(">K2\nOTHERK\n")
    p2, _ = ext.dbcheck_panel(str(exe), str(crap), str(tmp_path / "cache"))
    assert p2 != p and ">dbcheck_contaminants|K2\n" in open(p2).read() and "K1" in open(p).read()


# R3: two runs with different sources under one cache name never share a file, whatever the interleaving
def test_cache_paths_by_fingerprint(tmp_path):
    a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
    a.write_text(">A\nAAAK\n")
    b.write_text(">B\nBBBK\n")
    dst = str(tmp_path / "c" / "union.fasta")
    copy = lambda src: (lambda part: open(part, "w").write(open(src).read()))
    pa = ext.cached(dst, [str(a)], copy(a))
    pb = ext.cached(dst, [str(b)], copy(b))
    assert pa != pb and open(pa).read() == ">A\nAAAK\n" and open(pb).read() == ">B\nBBBK\n"
    assert ext.cached(dst, [str(a)], lambda part: 1 / 0) == pa   # a cached fingerprint is never rewritten


# R1: a --fasta stem equal to a given organism's proteome, and one named <stem>_fasta, in both orders
@pytest.mark.parametrize("order", [0, 1])
def test_fasta_names_against_given_organism(tmp_path, monkeypatch, order):
    human = tmp_path / "Homo_sapiens.panel.fasta"
    human.write_text(">h\nHUMANK\n")
    (tmp_path / "x").mkdir()
    f1, f2 = tmp_path / "x" / "Homo_sapiens_fasta.fasta", tmp_path / "Homo_sapiens.fasta"
    f1.write_text(">f1\nFIRSTK\n")
    f2.write_text(">f2\nSECONDK\n")
    fastas = [str(f1), str(f2)] if order == 0 else [str(f2), str(f1)]
    monkeypatch.setattr(cli, "_panel", lambda a, with_fasta=False: {"Homo_sapiens": str(human)})
    a = argparse.Namespace(organism=["Homo sapiens"], fastas=fastas, panel=None, out=str(tmp_path / "o"))
    c = cli._Run(a, lambda *x: None, "x", "x.mzML")
    names, _, db = cli._organism_given(c)
    assert sorted(db["fasta"]) == sorted([str(human), str(f1), str(f2)]) and len(set(names)) == 3
    heads = {l for l in open(c.union) if l.startswith(">")}
    assert heads == {">h\n", ">f1\n", ">f2\n"}


def _mztab_with(term, tmp_path):
    src = open(os.path.join(HERE, "data", "MzTabExporter_1_head.mzTab")).read()
    old = "[MS, MS:1003115, OpenMS:Target-decoy PSM q-value, ]"
    assert old in src
    p = tmp_path / "x.mzTab"
    p.write_text(src.replace(old, term))
    return str(p)


# R2: PSM q-values by CV accession: Percolator's is read, a peptide-level one is not
def test_mztab_qvalue_cv_terms(tmp_path):
    ided, n = CK.identified(_mztab_with("[MS, MS:1001491, percolator:Q value, ]", tmp_path))
    assert n == 40 and len(ided) > 30
    ided, n = CK.identified(_mztab_with("[MS, MS:1002354, PSM-level q-value, ]", tmp_path))
    assert len(ided) > 30
    ided, n = CK.identified(_mztab_with("[MS, MS:1003116, OpenMS:Target-decoy peptide q-value, ]", tmp_path))
    assert n == 40 and ided == set()


# 7: ddaPASEF: the converted mzML is temporary; check uses the given run, not a deleted path
def test_ddapasef_converted_spectra_for_check(tmp_path, monkeypatch):
    rep = tmp_path / "r.auto.json"
    gone = tmp_path / "prose_auto_ddapasef_x" / "run.mzML"      # deleted when the analysis exited
    rep.write_text(json.dumps({"input": {"pseudo_spectra": str(gone)}, "prose": {"params": {}}, "decisions": {}}))
    seen = {}
    monkeypatch.setattr(ext, "require_fastag", lambda f: None)
    monkeypatch.setattr(CK, "identified", lambda r, fdr: (set(), 0))
    monkeypatch.setattr(CK, "placements", lambda fastag, spectra, *x, **k: seen.update(s=spectra) or (set(), {}, 0))
    a = argparse.Namespace(run=str(tmp_path / "run.d"), result="r.idXML", fasta=["db.fasta"], fastag="FASTag",
                           report=str(rep), contaminants=None, fdr=0.01, threads=1, out=str(tmp_path / "o"))
    CK.main(a, log=lambda *x: None)
    assert seen["s"] == a.run
    # the analysis does not report a temporary conversion as the searched spectra
    c = cli._Run(argparse.Namespace(), print, "x", "x.d")
    c.spectra_file, c.pseudo, c.spectra_tmp = str(gone), True, True
    assert cli._reported_spectra(c) is None
    c.spectra_tmp = False
    assert cli._reported_spectra(c) == str(gone)


# 8: the given-database and panel FASTag passes use the same fixed label
def test_db_check_passes_same_label(tmp_path, monkeypatch):
    a = argparse.Namespace(fastag="FASTag", threads=1, fast=True, check_organism="warn", contaminants=None,
                           panel=None, out=str(tmp_path), fastag_subsample=0.5)
    c = cli._Run(a, lambda *x: None, "x", "x.mzML")
    c.ev = [Evidence("label", "no isobaric label", 0.8, "data", "reporter ions")]
    used = {}

    def alk(*x, reports=None, fixed=(), **k):
        used["given"] = list(fixed)
        reports.update({"qc": {"x": 1}, "confident": {"s1"}, "placed": {"s1"}})
        return []
    monkeypatch.setattr(ext, "fastag_alkylation", alk)
    monkeypatch.setattr(ext, "qc_evidence", lambda qc: ([Evidence("label", "TMT10/11plex", 0.35, "data",
                                                                    "FASTag QC reporter kits")], {}))
    monkeypatch.setattr(ext, "dbcheck_panel", lambda *x: ("panel.fasta", ["human"]))
    monkeypatch.setattr(ext, "fastag_dbcheck", lambda *x, fixed=(), **k: used.update(panel=list(fixed)) or (set(), {}))
    cli._fastag_gaps(c, "db.fasta")
    cli._db_check(c, [])
    assert used["given"] == used["panel"] == []
    assert cli._label_mods(c).k == 0.0     # nor does the QC kit set the mini-search's fixed label


# third review: a union of compressed and plain proteomes, from a .gz organism panel and a plain --fasta
def test_union_decompresses_gz_sources(tmp_path, monkeypatch):
    panel = tmp_path / "panel"
    panel.mkdir()
    with gzip.open(panel / "Homo_sapiens.fasta.gz", "wt") as fh:
        fh.write(">sp|H1|HUMAN\nHUMANSEQK")            # compressed, no final newline
    plain = tmp_path / "my_db.fasta"
    plain.write_text(">P1\nPLAINSEQK\n")
    monkeypatch.setattr(cli, "_panel", lambda a, with_fasta=False: {"Homo_sapiens": str(panel / "Homo_sapiens.fasta.gz")})
    a = argparse.Namespace(organism=["Homo sapiens"], fastas=[str(plain)], panel=None, out=str(tmp_path / "o"))
    c = cli._Run(a, lambda *x: None, "x", "x.mzML")
    cli._organism_given(c)
    assert open(c.union).read() == ">sp|H1|HUMAN\nHUMANSEQK\n>P1\nPLAINSEQK\n"
    out = io.StringIO()
    ext.copy_fasta(out, str(panel / "Homo_sapiens.fasta.gz"), "human")
    assert out.getvalue() == ">dbcheck_human|sp|H1|HUMAN\nHUMANSEQK\n"

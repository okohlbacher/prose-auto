"""Adapters for external tools: FASTag and DIAspeXtract.

FASTag (github.com/okohlbacher/FASTag) reads sequence tags from MS2 spectra.
prose-auto uses two of its outputs:

* `-species`: a genus call from tag 7-mers against its bundled index of ~50
  reference proteomes (model organisms, contaminants, pathogens, microbiome).
  Database free, so it also works when the organism is not in the panel.
  Genus level only: close relatives (mammals) share most 7-mers, and reagent
  genera (Bos, Sus, Oryctolagus, ...) show up in almost every run.
* `-recon_fasta` + `-delta_out`: tags placed on a proteome with NO fixed
  modification; the histogram of flank mass gaps shows the Cys alkylation
  (+57.021 carbamidomethyl, +45.988 methylthio, ...) and labels as peaks.

DIAspeXtract (github.com/okohlbacher/diaspextract) turns diaPASEF (Bruker .d,
or DIA mzML) into pseudo-MS/MS spectra with one precursor each, so the DDA
detectors and the mini-search can run on DIA data too.
"""

from __future__ import annotations

import csv
import functools
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile

from .evidence import Evidence

REAGENT_GENERA = {"Bos", "Sus", "Oryctolagus", "Gallus", "Capra", "Ovis", "Equus", "Lysobacter",
                  "Achromobacter", "Staphylococcus"}

# Cys alkylation and isobaric labels as mass gaps when the search sets no fixed modification
DELTA_CANDIDATES = {
    "Carbamidomethyl (C)": 57.021464,
    "Propionamide (C)": 71.037114,
    "Nethylmaleimide (C)": 125.047679,
    "Methylthio (C)": 45.987721,
    "Carboxymethyl (C)": 58.005479,
}


def _run(cmd, log, timeout=6 * 3600):
    log("  $ " + " ".join(str(c) for c in cmd))
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        e = RuntimeError((r.stderr or r.stdout)[-1500:])
        e.returncode = r.returncode
        raise e
    return r


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _frag_args(frag_tol, frag_unit):
    return ["-fragment_tolerance", frag_tol, "-fragment_tolerance_unit", frag_unit]


def _taxonomy_args(taxonomy_dir, rank):
    """an external FASTag index (tax_k7.taxdb + pruned nodes/names.dmp in one directory) and the report rank"""
    a = ["-species_rank", rank]
    if taxonomy_dir:
        a += ["-taxdb", os.path.join(taxonomy_dir, "tax_k7.taxdb"),
              "-taxonomy_nodes", os.path.join(taxonomy_dir, "nodes.dmp"),
              "-taxonomy_names", os.path.join(taxonomy_dir, "names.dmp")]
    return a


def fastag_species(fastag, mzml, frag_tol=20.0, frag_unit="ppm", fixed=("Carbamidomethyl (C)",),
                   subsample=0.0, threads=0, log=print, taxonomy_dir=None, rank="genus"):
    """Taxon call from FASTag -species: genus with the shipped index, species with an extended one."""
    with tempfile.TemporaryDirectory(prefix="fastag_") as tmp:
        out = os.path.join(tmp, "tags.tsv")
        sp = os.path.join(tmp, "species.tsv")
        # -species_evidence best (1.10.0): a spectrum counts for the taxa of its best tag; 'all' (1.9.0) took 19 fewer
        # right species calls on the 193 validation runs and called Triticum aestivum 58 times instead of 32
        cmd = [fastag, "-in", mzml, "-out", out, "-species", "-species_evidence", "best", "-tag_length", 7,
               "-species_out", sp, "-threads", threads] + _frag_args(frag_tol, frag_unit) + _taxonomy_args(taxonomy_dir, rank)
        cmd += ["-fixed_modifications"] + list(fixed)
        if subsample:
            cmd += ["-subsample_fraction", subsample]
        _run(cmd, log)
        rows = _tsv(sp) if os.path.exists(sp) else []
    if not rows:
        return [Evidence("organism", "unknown", 0.1, "data", "FASTag species", "FASTag reported no taxon")]
    for r in rows:
        r["adjusted"] = float(r.get("adjusted") or r.get("observed") or 0)
        r["observed"] = float(r.get("observed") or 0)
    rows.sort(key=lambda r: -r["adjusted"])
    table = [{k: r.get(k) for k in ("name", "observed", "adjusted", "expected", "qvalue")} for r in rows[:15]]
    # the top non-reagent genus is the call; a reagent genus on top is reported as such
    top = rows[0]
    real = [r for r in rows if str(r["name"]).split()[0] not in REAGENT_GENERA]  # genus or species rank
    call = real[0] if real and (str(top["name"]).split()[0] in REAGENT_GENERA and real[0]["adjusted"] >= 0.5 * top["adjusted"]) else top
    nxt = next((r for r in rows if r is not call), None)
    margin = (call["adjusted"] - (nxt["adjusted"] if nxt else 0)) / max(call["adjusted"], 1)
    n = call["observed"]
    conf = 0.25 + min(0.35, 0.12 * (n ** 0.5) / 3) + min(0.3, margin)
    conf = min(0.85, conf if n >= 20 else min(conf, 0.35))
    detail = (f"FASTag {rank} call {call['name']} from {int(n)} tag-carrying spectra"
              + (f"; next {nxt['name']} ({nxt['adjusted']:.0f} vs {call['adjusted']:.0f})" if nxt else "")
              + ("; reagent genus on top skipped" if call is not top else ""))
    return [Evidence("organism", call["name"], conf, "data", "FASTag species", detail, {"taxa": table})]


FASTAG_MIN_VERSION = (1, 15, 2)   # 1.10.0: -digest_out, -label_out, -massgap_out, -species_evidence; 1.15.2: -alkylation_out,
                                  # -alkylation_tags_out, -qc_out reporter kits (1.15.3 gives 1.10.0's tags, ends and placements)
DIGEST_MODELS = os.path.join(os.path.dirname(__file__), "data", "prose_auto_013_digest_models.tsv")
_NLABEL = {"TMT6/10/11plex": ("TMT", 229.162932), "TMTpro (16/18plex)": ("TMTpro", 304.207146),
           "iTRAQ4plex": ("iTRAQ4", 144.102063)}
KR_TRYPTIC = 0.8
SILAC_MIN_TAGS = 50   # K- and R-terminal tags each before the heavy shares count (chymotrypsin: Lys8 50% of 12% K/R)
TAG_MIN = 175         # C-terminal tag ends before the termini say anything: FASTag 1.9.0 reads 1.75x (median, 111
                      # runs) the C-terminal tags of 1.8.0, so 175 keeps 1.8.0's 100; at 100 a TiO2 run read HLA-II
TAG_AAS = "GASPVTCLNDQKEMHFRYW"   # I merged into L
# model -> enzyme value; HLA-I ligands and proteasome products are class I, HLA-II ligands and unspecific cleavage
# class II (their ends barely differ from the proteome: "no protease")
TAG_MODELS = {"Trypsin": "Trypsin", "Lys-C": "Lys-C", "Arg-C": "Arg-C", "Glu-C": "Glu-C", "Chymotrypsin": "Chymotrypsin",
              "Asp-N": "Asp-N", "Lys-N": "Lys-N", "HLA-I": "unspecific cleavage", "proteasome": "unspecific cleavage",
              "HLA-II": "unspecific cleavage", "unspecific": "unspecific cleavage"}
HLA_CLASS = {"HLA-I": "I", "proteasome": "I", "HLA-II": "II", "unspecific": "II"}


@functools.lru_cache(maxsize=None)
def fastag_version(fastag):
    """(major, minor, patch) from the 'Version:' line of FASTag --help, or None"""
    r = subprocess.run([str(fastag), "--help"], capture_output=True, text=True, timeout=120)
    for line in (r.stdout + r.stderr).splitlines():
        if line.startswith("Version:"):
            try:
                return tuple(int(x) for x in line.split()[1].split(".")[:3])
            except ValueError:
                return None
    return None


def require_fastag(fastag):
    v = fastag_version(fastag)
    if not v or v < FASTAG_MIN_VERSION:
        raise RuntimeError(f"prose-auto needs FASTag >= {'.'.join(map(str, FASTAG_MIN_VERSION))} (-digest_out, "
                           f"-label_out, -massgap_out, -alkylation_out, -qc_out); {fastag} reports "
                           f"{'.'.join(map(str, v)) if v else 'no version'}")


def _report(path):
    """a FASTag long-form report (section, key, value, n, description) as {(section, key): value}"""
    return {(r["section"], r["key"]): r["value"] for r in _tsv(path)}


def _num(rep, section, key, default=0.0):
    v = rep.get((section, key), "-")
    try:
        return default if v in ("-", "", "n/a") else float(v)
    except ValueError:
        return default


def fastag_termini(fastag, mzml, frag_tol=20.0, frag_unit="ppm", subsample=0.25, threads=0, log=print):
    """Protease and labels from the termini of FASTag k=4 tags (Carbamidomethyl only, so labels are mass gaps).
    FASTag counts the tag ends (-digest_out) and the label gaps (-label_out) and scores the tag-profile models of
    data/prose_auto_013_digest_models.tsv (MEROPS, MHC Motif Atlas, proteome background) by multinomial
    log-likelihood. On 17 calibration runs: C-terminal K/R share 0.86-0.99 on tryptic HCD runs (threshold 0.80; on
    180 benchmark runs tryptic Orbitrap DDA 0.84-0.98, tryptic timsTOF 0.33-0.66), 0.27-0.42 on HLA class I and II
    (class I: hydrophobic anchors L/V/F/Y/M 0.44-0.63, class II about 0.1); N-terminal TMT/TMTpro on 5-10% of the tags
    of labelled runs (>= 0.1% otherwise); heavy Lys8 and Arg10 shares >= 0.25 only in SILAC runs (label-free <= 0.16,
    diaPASEF pseudo-spectra show an Arg10 artifact without Lys8). FASTag gives no call on MS2 scans that start above
    y1 of K (PXD066944, PXD038809: K 0.02-0.13 on tryptic runs). High-res MS2 only; 1-16 s on 25%.
    Returns evidence (enzyme, label) and the stats (C-terminal residue shares, HLA class hint)."""
    if frag_unit != "ppm":
        return [], {}
    with tempfile.TemporaryDirectory(prefix="fastag_termini_") as tmp:
        dg, lb = os.path.join(tmp, "digest.tsv"), os.path.join(tmp, "label.tsv")
        _run([fastag, "-in", mzml, "-out", os.path.join(tmp, "tags.tsv"), "-tag_length", 4, "-threads", threads,
              "-subsample_fraction", subsample, "-digest_out", dg, "-label_out", lb, "-digest_models", DIGEST_MODELS]
             + _frag_args(frag_tol, frag_unit), log)
        return termini_evidence(_report(dg), _report(lb))


def termini_evidence(d, lab):
    """The prose-auto policy on FASTag's -digest_out and -label_out reports (dicts from _report)."""
    tags, tot = int(_num(d, "count", "tags_filtered")), int(_num(d, "count", "tot"))
    if tot < TAG_MIN:
        return [], {"tags": tags, "cterm_tags": tot}
    # res: the C-terminal residue keys (labelled K and Arg10 apart), as the K/R share counts them
    res = {a: int(_num(d, "pos", f"C.P1.{a}")) for a in TAG_AAS if a not in "KR"}
    res.update({k: int(_num(d, "residue", f"C.{k}")) for k in ("K", "K+TMT", "K+TMTpro", "K+iTRAQ4", "K+8", "R", "R+10")})
    share = {a: v / tot for a, v in res.items()}
    cterm_k, n_n = _num(d, "pos", "C.P1.K"), int(_num(d, "count", "n_n"))
    stats = {"tags": tags, "cterm_tags": tot, "cterm_KR": round(_num(d, "share", "kr"), 3),
             "cterm_anchors_LVFYM": round(_num(d, "share", "anchors_LVFYM"), 3),
             "cterm_top": {a: round(v, 3) for a, v in sorted(share.items(), key=lambda x: -x[1])[:6]},
             "nterm_tags": n_n, "nterm_D": round(_num(d, "share", "nterm_D"), 3),
             "ms2_first_peak_mz": _num(d, "meta", "ms2_first_peak_mz", None)}
    ev = []
    if d.get(("call", "status")) == "call":
        model, other, margin = d[("call", "model")], d[("call", "runner_up")], _num(d, "call", "margin")
        # isobaric labels make K-ending tags easy to read (intense labelled y1): tryptic TMTpro runs show K 0.93,
        # R 0.06 (PXD037434, PXD001670), so the K/R split is no evidence there
        iso_k = res["K+TMT"] + res["K+TMTpro"] + res["K+iTRAQ4"]
        if model in ("Lys-C", "Arg-C") and TAG_MODELS.get(other) == "Trypsin" and iso_k > cterm_k / 2:
            margin = 0.0
        stats.update({"tag_call": model, "tag_runner_up": other, "tag_margin": round(margin, 3)})
        n_c = int(_num(d, "count", "n_c"))
        why = (f"tag-profile call over {n_c} C- and {n_n} N-terminal residues: {model} "
               f"(next {other}, {margin:.2f} log-likelihood per tag)")
        if model in HLA_CLASS:
            stats["hla_hint"] = f"HLA class {HLA_CLASS[model]}?"
        # the margin is the confidence: a weak call loses to a confident Jev answer (reconcile._enzyme)
        ev.append(Evidence("enzyme", TAG_MODELS[model], round(min(0.8, 0.45 + margin), 3), "data", "FASTag termini",
                           why, stats))
    else:
        stats["tag_nocall"] = d.get(("call", "reason"))
    for name, (key, m) in _NLABEL.items():
        x, bg = _num(lab, "nlabel", f"{key}.count"), _num(lab, "nlabel", f"{key}.background")
        if x >= 0.01 * tags and x >= 20 * max(bg, 1):
            ev.append(Evidence("label", name, 0.5, "data", "FASTag termini",
                               f"{x / tags:.1%} of the sequence tags start after a {m:.3f} Da N-terminal gap ({name}"
                               + ("; same mass as iTRAQ8plex, the reporter ions decide" if key == "TMTpro" else "")
                               + ")", stats))
    hk, hr = _num(lab, "heavy", "K8.share"), _num(lab, "heavy", "R10.share")
    stats.update({"heavy_K8": round(hk, 3), "heavy_R10": round(hr, 3)})
    nk, nr = _num(lab, "heavy", "K.tags"), _num(lab, "heavy", "R.tags")
    if hk >= 0.2 and hr >= 0.2 and stats["cterm_KR"] >= KR_TRYPTIC and min(nk, nr) >= SILAC_MIN_TAGS:
        # both heavy shares >= 30% occurred only in SILAC runs (5 of 164, three of them annotated label free in the
        # SDRF but SILAC by their text or file name); 20-30% stays a hint
        ev.append(Evidence("label", "SILAC (Lys8/Arg10)", 0.7 if min(hk, hr) >= 0.3 else 0.5, "data", "FASTag termini",
                           f"heavy C-terminal Lys8 {hk:.0%} and Arg10 {hr:.0%} of the tryptic C-termini", stats))
    return ev, stats


def tims_protease(stats):
    """The protease from the tag termini on timsTOF, where FASTag's tag-profile call fails: the MS2 scans start at
    m/z 200-235, above y1 of K (147.1) and R (175.1), so the C-terminal residue is read from the gap to the precursor
    and an R end often reads as R+10. The K/R share (heavy variants included) still separates: 13 tryptic timsTOF
    runs 0.33-0.75, Glu-C, Asp-N and chymotrypsin 0.006-0.025. Evidence or None (below TAG_MIN tag ends)."""
    kr, top = stats.get("cterm_KR"), stats.get("cterm_top") or {}
    if kr is None or stats.get("cterm_tags", 0) < TAG_MIN:
        return None
    if kr >= 0.30:
        val = "Trypsin"
    elif kr > 0.10:
        return None
    elif top.get("E", 0) >= 0.5:
        val = "Glu-C"
    elif stats.get("nterm_D", 0) >= 0.3:
        val = "Asp-N"
    elif sum(top.get(a, 0) for a in "FYWLM") >= 0.6:
        val = "Chymotrypsin"
    else:
        return None
    return Evidence("enzyme", val, 0.6, "data", "FASTag termini (timsTOF)",
                    f"timsTOF tag ends: K/R {kr:.2f}, C-terminal {', '.join(f'{a} {v:.2f}' for a, v in top.items())}, "
                    f"N-terminal D {stats.get('nterm_D', 0):.2f}", stats)


def fastag_shortlist(evs, n=3):
    """The top-n taxa of a FASTag species ranking (reagent genera only when nothing else is ranked), for the
    two-tier organism test."""
    taxa = [t for e in evs if e.method == "FASTag species" for t in (e.stats or {}).get("taxa", [])]
    names = [t["name"] for t in taxa if t.get("name")]
    real = [x for x in names if x.split()[0] not in REAGENT_GENERA]
    return (real or names)[:n]


CONFIDENT_EVALUE = 1e-4               # the tags the database check counts


def _gap_cmd(fastag, mzml, fasta, tmp, frag_tol, frag_unit, threads, subsample, fixed=()):
    """The tagging and placement pass of the gap stage (and of the database check, on another FASTA): k=4 tags
    extended by 6, placed with -recon_min_length 5 (the size-derived default, 6-7 on a proteome, placed 2.5-14x
    fewer spectra) and no fixed modification but a label, so the Cys alkylation shows as a mass gap."""
    cmd = [fastag, "-in", mzml, "-out", os.path.join(tmp, "tags.tsv"), "-tag_length", 4, "-extension", 6,
           "-threads", threads, "-recon_fasta", fasta, "-recon_out", os.path.join(tmp, "recon.tsv"),
           "-recon_min_length", 5, "-fixed_modifications"] + list(fixed) + _frag_args(frag_tol, frag_unit)
    return cmd + (["-subsample_fraction", subsample] if subsample else [])


def fastag_alkylation(fastag, mzml, fasta, frag_tol=20.0, frag_unit="ppm", threads=0, log=print, subsample=0.0,
                      placements=None, reports=None, fixed=()):
    """Cys alkylation and modification combinations from FASTag's mass-gap report (-massgap_out: tags placed on
    fasta without fixed modifications). placements, a list, receives (spectrum native id, peptide) of the first
    placement per spectrum: the mass errors (hypotheses.placement_accuracy). reports, a dict,
    adds -alkylation_out, -alkylation_tags_out and -qc_out to the pass and receives them ("alkylation",
    "alkylation_tags", "qc") with the database-check counts of the given database ("confident", "placed": sets of
    spectra). fixed: a label as fixed modification (an isobaric label is part of every gap otherwise, and
    -alkylation_out finds no Cys form)."""
    with tempfile.TemporaryDirectory(prefix="fastag_") as tmp:
        mg = os.path.join(tmp, "massgap.tsv")
        cmd = _gap_cmd(fastag, mzml, fasta, tmp, frag_tol, frag_unit, threads, subsample, fixed) + ["-massgap_out", mg]
        out = {"alkylation": "alk.tsv", "alkylation_tags": "alkt.tsv", "qc": "qc.tsv"} if reports is not None else {}
        for k, f in out.items():
            cmd += [f"-{k}_out", os.path.join(tmp, f)]
        _run(cmd, log)
        rc = os.path.join(tmp, "recon.tsv")
        if placements is not None and os.path.exists(rc):
            seen = set()
            for r in _tsv(rc):
                if r["spectrum"] not in seen:
                    seen.add(r["spectrum"])
                    placements.append((r["spectrum"], r["peptide"]))
        if reports is not None:
            for k, f in out.items():
                p = os.path.join(tmp, f)
                reports[k] = _report(p) if os.path.exists(p) else {}
            reports["confident"], placed = exact_placements(tmp)
            reports["placed"] = placed.get("", set())
        return massgap_evidence(_report(mg) if os.path.exists(mg) else {})


# flank masses that must agree for a placement to count: the flank that depends on the precursor mass breaks on
# isotope errors, unfixed Cys and N-terminal modifications (an exact placement for only 1 in 3 placed spectra of a
# label-free Orbitrap run), so one agreeing flank counts (calibration: v0111/README_0111.md)
DBCHECK_FLANKS = 1


def exact_placements(tmp, sep="|", flanks=DBCHECK_FLANKS):
    """The database check's counts from a FASTag tags.tsv and recon.tsv in tmp: the spectra with a confident tag
    (E <= CONFIDENT_EVALUE), and per database label the spectra whose confident tag places with at least `flanks`
    agreeing flank masses (2: both, exact). The label is the protein accession's prefix up to sep in a check panel
    (dbcheck_panel), "" in any other FASTA. FASTag lists every protein a tag places on, so one pass serves every
    member of a panel."""
    conf = {}
    for r in _tsv(os.path.join(tmp, "tags.tsv")):
        try:
            if float(r["evalue"]) <= CONFIDENT_EVALUE:
                conf.setdefault(r["spectrum"], set()).add(r["tag"])
        except (KeyError, ValueError):
            continue
    placed = {}
    rc = os.path.join(tmp, "recon.tsv")
    for r in (_tsv(rc) if os.path.exists(rc) else []):
        if (r["nterm_match"] == "1") + (r["cterm_match"] == "1") >= flanks and r["tag"] in conf.get(r["spectrum"], ()):
            prot = r["protein"]
            label = prot.split(sep, 1)[0] if prot.startswith(DBCHECK_PREFIX) else ""
            placed.setdefault(label, set()).add(r["spectrum"])
    return set(conf), placed


DBCHECK_PREFIX = "dbcheck_"   # accession prefix of the check panel members: dbcheck_human|sp|P12345|...


def fastag_share(fastag, sub):
    """FASTag's shipped data folder <sub> (alkylation, digest, taxonomy): share-FASTag-<sub> next to the executable
    (the release bundle), else <prefix>/share/FASTag/<sub>. A bare name is looked up on PATH, as subprocess runs it
    (0.1.12: realpath resolved it against the working directory)."""
    exe = shutil.which(fastag) if os.sep not in fastag else fastag
    d = os.path.dirname(os.path.realpath(exe or fastag))
    for p in (os.path.join(d, f"share-FASTag-{sub}"), os.path.join(d, "..", "share", "FASTag", sub)):
        if os.path.isdir(p):
            return os.path.normpath(p)
    return None


def sources_key(paths):
    """A fingerprint of the files a cached FASTA is built from: path and SHA-256 of the content (0.1.12)"""
    h = hashlib.sha256()
    for p in paths:
        h.update(os.path.realpath(p).encode() + b"\0")
        with open(p, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
    return h.hexdigest()


def cached(dst, paths, write):
    """The cache file for dst built from paths: dst with the first 16 hex digits of the sources' fingerprint before
    its extension, written once by write(part) and never changed afterwards. A changed source gives a new path, so a
    search that opened the old one keeps reading consistent content, and concurrent writers of one fingerprint
    publish identical files (0.1.12; 0.1.11 reused any file of the name, the first 0.1.12 kept the fingerprint in a
    side file that a concurrent run could pair with another run's content)."""
    stem, ext = os.path.splitext(dst)
    path = f"{stem}.{sources_key(paths)[:16]}{ext}"
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        part = f"{path}.{os.getpid()}.part"
        write(part)
        os.replace(part, path)
    return path


def write_labelled(dst, parts):
    """Concatenate (label, FASTA) parts into dst, every header prefixed DBCHECK_PREFIX + label + "|", every part
    ending with a newline (a FASTA without a final newline would merge its last line with the next header)"""
    with open(dst, "w") as out:
        for label, p in parts:
            last = "\n"
            with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
                for line in fh:
                    out.write(f">{DBCHECK_PREFIX}{label}|{line[1:]}" if line.startswith(">") else line)
                    last = line
            if not last.endswith("\n"):
                out.write("\n")


def dbcheck_panel(fastag, contaminants, cache):
    """The database check's panel as one labelled FASTA, written once to <cache>/dbcheck_panel.fasta: FASTag's
    shipped human, mouse, yeast and E. coli proteomes (UniProt reviewed, share/FASTag/alkylation, 17 MB) and the
    contaminants, at a path fingerprinted by their content (cached). Returns (path, member names) or (None, []) without
    the shipped proteomes."""
    share = fastag_share(fastag, "alkylation")
    members = [(n, os.path.join(share, f"{n}.fasta.gz")) for n in ("human", "mouse", "yeast", "ecoli")] if share else []
    members = [m for m in members if os.path.exists(m[1])]
    if not members:
        return None, []
    if contaminants:
        members.append(("contaminants", contaminants))
    names = [n for n, _ in members]
    dst = os.path.join(cache, "dbcheck_panel" + ("+contaminants" if contaminants else "") + ".fasta")
    return cached(dst, [p for _, p in members], lambda part: write_labelled(part, members)), names


def fastag_dbcheck(fastag, mzml, panel_fasta, frag_tol=20.0, frag_unit="ppm", threads=0, log=print, subsample=0.0,
                   fixed=()):
    """The database check's panel pass: the gap stage's tagging (same spectra, same subsample seed) placed on the
    check panel. Returns (confident spectra, {member: spectra placed exactly})."""
    with tempfile.TemporaryDirectory(prefix="fastag_dbcheck_") as tmp:
        cmd = _gap_cmd(fastag, mzml, panel_fasta, tmp, frag_tol, frag_unit, threads, subsample, fixed)
        _run(cmd + ["-max_evalue", CONFIDENT_EVALUE], log)
        conf, placed = exact_placements(tmp)
    return conf, {k[len(DBCHECK_PREFIX):]: v for k, v in placed.items() if k}


# Database check thresholds, calibrated on the development set (the 218 release runs with their correct FASTA and
# with related, unrelated and incomplete ones; in-sample): 158 of 183 correct FASTAs pass, 120/166 related, 140/167
# unrelated and 12/16 incomplete ones warn (v0111/README_0111.md). A run with fewer confident tags is not assessed.
DBCHECK_MIN_TAGS = 100
DBCHECK_MIN_RATE = 0.05    # share placed on the given FASTA below this: warn (correct FASTAs: median 0.51)
DBCHECK_OWN_MIN = 0.03     # a panel member alone places this share of the confident spectra ...
DBCHECK_OWN_REL = 0.10     # ... and at least this fraction of the given FASTA's share: warn


def database_check(confident, given, members, organisms=()):
    """The policy on the placement counts: the share of the spectra with a confident tag that place on the given
    FASTA (contaminants count as explained), and per panel member the share it places that the given FASTA does
    not ("own"). Warns when the given share is low or a member has many spectra of its own (in absolute terms and
    against the given share); with too few confident tags it abstains. organisms (taxa.key of the given organisms) leave a member out that the user gave.
    Returns a dict for the report."""
    n = len(confident)
    crap = members.get("contaminants", set())
    explained = (given | crap) & confident
    rate = len(explained) / n if n else 0.0
    rates = {m: round(len(s & confident) / n, 4) for m, s in members.items()} if n else {}
    own = {m: round(len((s & confident) - explained) / n, 4) for m, s in members.items() if m != "contaminants"} \
        if n else {}
    out = {"confident_spectra": n, "given_rate": round(rate, 4), "panel_rates": rates, "panel_own": own,
           "thresholds": {"min_tags": DBCHECK_MIN_TAGS, "min_rate": DBCHECK_MIN_RATE, "own_min": DBCHECK_OWN_MIN,
                          "own_rel": DBCHECK_OWN_REL, "flanks": DBCHECK_FLANKS}}
    if n < DBCHECK_MIN_TAGS:
        return {**out, "status": "not assessed", "reasons": [f"{n} spectra with a confident tag (E <= "
                                                             f"{CONFIDENT_EVALUE:g}), {DBCHECK_MIN_TAGS} needed"]}
    reasons = []
    if rate < DBCHECK_MIN_RATE:
        reasons.append(f"only {rate:.1%} of the {n} spectra with a confident tag place on the given database "
                       f"(at least {DBCHECK_MIN_RATE:.0%} expected)")
    bar = max(DBCHECK_OWN_MIN, DBCHECK_OWN_REL * rate)
    for m, x in sorted(own.items(), key=lambda kv: -kv[1]):
        if x >= bar and m not in organisms:
            reasons.append(f"the {m} proteome places {x:.1%} of the spectra with a confident tag that the given "
                           f"database does not (below {bar:.1%} expected)")
    return {**out, "status": "warning" if reasons else "ok", "reasons": reasons}


ALK_VALUE = {"Carbamidomethyl": "Carbamidomethyl (C)", "Propionamide": "Propionamide (C)",
             "Nethylmaleimide": "Nethylmaleimide (C)", "Methylthio": "Methylthio (C)",
             "Carboxymethyl": "Carboxymethyl (C)", "none": "none (free Cys)", "unmodified": "none (free Cys)"}


def alkylation_evidence(alk, alkt):
    """FASTag's Cys calls: -alkylation_out (placements of single-Cys peptides on the given database; experimental)
    and -alkylation_tags_out (database-free; held-out 37 of 40 calls right, 17 of 57 runs abstained). A call is
    evidence, a no-call none; the mini-search competition is the fallback (cli._cys_fallback). The placement call is a
    database result like the competition (source search, above the mass-gap peak: on the 218 release runs it agreed
    with the competition's winner on all 87 runs it called); the database-free call is data."""
    ev = []
    for rep, sec, method, conf, src in ((alk, "alkylation_call", "FASTag alkylation (placements)", 0.86, "search"),
                                        (alkt, "alkylation_tags_call", "FASTag alkylation (tags)", 0.6, "data")):
        status = (rep or {}).get((sec, "status"))
        if status != "call":
            continue
        reagent = rep.get((sec, "reagent"), "-")
        val = ALK_VALUE.get(reagent)
        if not val:
            continue
        share = _num(rep, sec, "share")
        ev.append(Evidence("cys_alkylation", val, conf, src, method,
                           f"FASTag calls {reagent}: share {share:.2f} of the binned spectra", {"share": share}))
    return ev


QC_KIT_LABEL = {"tmt6plex": "TMT6plex", "tmt10plex": "TMT10/11plex", "tmt11plex": "TMT10/11plex",
                "tmt16plex": "TMTpro (16/18plex)", "tmt18plex": "TMTpro (16/18plex)", "itraq4plex": "iTRAQ4plex",
                "itraq8plex": "iTRAQ8plex"}
QC_ANALYZER = {"FTMS": "ORBITRAP", "ITMS": "IT", "TOF": "TOF"}


def qc_evidence(qc):
    """FASTag's run QC (-qc_out, experimental) next to the header: the MS1/MS2 analyser class (from the filter
    string or the instrument configuration), the precursor charge histogram and the isobaric reporter kits. They
    are low-confidence evidence: the header and prose-auto's own detectors decide, a disagreeing analyser class
    shows as a conflict. Returns (evidence, summary for the report)."""
    if not qc:
        return [], {}
    ev, summ = [], {}
    for lvl in (1, 2):
        counts = {cls: int(_num(qc, "qc", f"count_ms{lvl}_{cls}")) for cls in QC_ANALYZER}
        tot = int(_num(qc, "qc", f"count_ms{lvl}"))
        top = max(counts, key=counts.get)
        summ[f"ms{lvl}_analyzer_class"] = {k: v for k, v in counts.items() if v}
        res = qc.get(("qc", f"resolution_{top}_ms{lvl}"), "n/a")
        if res not in ("n/a", "-", ""):
            summ[f"ms{lvl}_resolution"] = res
        if tot and counts[top] >= 0.9 * tot:
            ev.append(Evidence(f"ms{lvl}_analyzer", QC_ANALYZER[top], 0.5, "metadata", "FASTag QC",
                               f"FASTag QC: {counts[top]} of {tot} MS{lvl} spectra read out in {top}"
                               + (f", resolution {res}" if res not in ("n/a", "-", "") else "")))
    hist = {int(k.split("_")[-1]): int(float(v)) for (s, k), v in qc.items()
            if s == "qc" and k.startswith("precursor_charge_") and not k.endswith("share") and k[-1].isdigit()}
    hist = {z: n for z, n in hist.items() if z > 0 and n}
    if hist:
        summ["charge_histogram"] = hist
        tot, acc, lo, hi = sum(hist.values()), 0, None, None
        for z in sorted(hist):
            acc += hist[z]
            lo = z if lo is None and acc >= 0.005 * tot else lo
            hi = z if hi is None and acc >= 0.995 * tot else hi
        ev.append(Evidence("charge_range", [lo, hi], 0.35, "metadata", "FASTag QC",
                           f"FASTag QC: 99% of the header precursor charges lie in {lo}..{hi}", {"histogram": hist}))
    kits = [k for k in qc.get(("qc", "reporter_compatible_kits"), "-").split(",") if k in QC_KIT_LABEL]
    summ["reporter_kits"] = kits
    if kits:
        ev.append(Evidence("label", QC_KIT_LABEL[kits[0]], 0.35, "data", "FASTag QC reporter kits",
                           "FASTag QC: reporter channels compatible with " + ", ".join(kits)))
    summ["tag_yield"] = _num(qc, "qc", "tag_yield", None)
    return ev, summ


# FASTag modification name -> ProSE names (the variable modifications a dominant gap proposes)
MOD_NAMES = {
    "Oxidation": ["Oxidation (M)"], "Dioxidation": ["Dioxidation (M)"], "Deamidated": ["Deamidated (N)", "Deamidated (Q)"],
    "Gln->pyro-Glu": ["Gln->pyro-Glu (N-term Q)"], "Glu->pyro-Glu": ["Glu->pyro-Glu (N-term E)"],
    "Carbamidomethyl": ["Carbamidomethyl (C)"], "Cation:Na": ["Cation:Na (DE)"], "Cation:K": ["Cation:K (DE)"],
    "Carbamyl": ["Carbamyl (N-term)"], "Formyl": ["Formyl (N-term)"], "Acetyl": ["Acetyl (K)", "Acetyl (N-term)"],
    "Methyl": ["Methyl (K)", "Methyl (R)"], "Dimethyl": ["Dimethyl (K)", "Dimethyl (R)"], "Trimethyl": ["Trimethyl (K)"],
    "Propionyl": ["Propionyl (K)", "Propionyl (N-term)"], "Butyryl": ["Butyryl (K)"], "Crotonyl": ["Crotonyl (K)"],
    "Succinyl": ["Succinyl (K)"], "Malonyl": ["Malonyl (K)"], "GG": ["GG (K)"],
    "Phospho": ["Phospho (S)", "Phospho (T)", "Phospho (Y)"], "Sulfo": ["Sulfo (Y)"], "Nitro": ["Nitro (Y)"],
    "HexNAc": ["HexNAc (S)", "HexNAc (T)"], "Hex": ["Hex (K)"], "Cysteinyl": ["Cysteinyl (C)"],
    "Methylthio": ["Methylthio (C)"], "Propionamide": ["Propionamide (C)"], "Nethylmaleimide": ["Nethylmaleimide (C)"],
    "TMT6plex": ["TMT6plex (K)", "TMT6plex (N-term)"], "TMTpro": ["TMTpro (K)", "TMTpro (N-term)"],
    "iTRAQ4plex": ["iTRAQ4plex (K)", "iTRAQ4plex (N-term)"], "Label:13C(6)15N(2)": ["Label:13C(6)15N(2) (K)"],
    "Label:13C(6)15N(4)": ["Label:13C(6)15N(4) (R)"], "Label:2H(4)": ["Label:2H(4) (K)"],
    "Label:13C(6)": ["Label:13C(6) (R)"], "Met-loss": ["Met-loss (Protein N-term M)"],
    "Met-loss+Acetyl": ["Met-loss+Acetyl (Protein N-term M)"],
}
ENRICHMENT = {"Phospho": "Phospho (STY)", "Acetyl": "Acetyl (K)", "GG": "GlyGly (K)", "HexNAc": "Glyco (HexNAc oxonium)"}
SILAC_GAPS = {"Label:2H(4)": "Lys4", "Label:13C(6)": "Arg6", "Label:13C(6)15N(2)": "Lys8", "Label:13C(6)15N(4)": "Arg10"}


def massgap_evidence(rep, max_mods_out=5):
    """The prose-auto policy on FASTag's -massgap_out report: the Cys alkylation peak; variable_mods (modifications
    explaining >= 2% of the spectra with a mass gap), enrichment (phospho, acetyl-lysine, GlyGly, HexNAc >= 10%) and
    SILAC labels from the explained gaps (the MS Andrea idea, Buur, Winkler & Dorfer 2026, on a histogram)."""
    counts = {n: int(_num(rep, "cys", n)) for n in DELTA_CANDIDATES}
    total, zero = int(_num(rep, "cys", "spectra")), int(_num(rep, "cys", "zero_bin"))
    split = lambda v: None if v in ("none", "-", "") else v.split(" + ")
    peaks, i = [], 1
    while (("peak", f"{i}.delta")) in rep:
        alt = rep.get(("peak", f"{i}.alternatives"), "-")
        ctx = rep.get(("peak", f"{i}.context"), "-")
        peaks.append({"delta": round(_num(rep, "peak", f"{i}.delta"), 4), "spectra": int(_num(rep, "peak", f"{i}.spectra")),
                      "explained_by": split(rep.get(("peak", f"{i}.explained_by"), "none")),
                      "alternatives": [] if alt in ("-", "") else [c.split(" + ") for c in alt.split(";")],
                      **({"context": ctx} if ctx not in ("-", "") else {})})
        i += 1
    share = {k: float(v) for (sec, k), v in rep.items() if sec == "share"}
    shifted = int(_num(rep, "combo", "shifted_spectra"))
    gaps = {"placed_spectra": int(_num(rep, "combo", "placed_spectra")), "shifted_spectra": shifted,
            "peaks": peaks[:25], "share": {k: round(v, 4) for k, v in share.items()}}
    stats = {"counts": counts, "zero_bin": zero, "spectra": total}
    best = max(counts, key=counts.get)
    second = sorted(counts.values())[-2]
    if counts[best] < 10:
        ev = [Evidence("cys_alkylation", "unknown", 0.15, "data", "FASTag mass gaps",
                       f"no Cys alkylation peak in the FASTag gap histogram ({total} spectra placed)", stats)]
    else:
        ratio = (counts[best] + 1) / (second + 1)
        ev = [Evidence("cys_alkylation", best, min(0.85, 0.4 + 0.15 * ratio ** 0.5), "data", "FASTag mass gaps",
                       f"{counts[best]} placed spectra with a {DELTA_CANDIDATES[best]:+.3f} Da gap "
                       f"(next candidate {second}); searched without fixed modifications", stats)]
    if not shifted:
        return ev
    mods = [n for n, x in sorted(share.items(), key=lambda kv: -kv[1]) if x >= 0.02 and n != "Carbamidomethyl"]
    if mods:
        ev.append(Evidence("variable_mods", [x for n in mods[:max_mods_out] for x in MOD_NAMES.get(n, [n])], 0.5, "data",
                           "delta-mass combinations", "mass gaps explained as modification combinations: "
                           + ", ".join(f"{n} {share[n]:.0%}" for n in mods[:max_mods_out])
                           + f" of {shifted} spectra with a mass gap", gaps))
    # two or more SILAC labels each explaining >= 8 % of >= 500 gapped spectra (triplex PXD016465: Arg6 20 %, Lys4 18 %,
    # Lys8 8 %); label-free validation runs reach 13 % for one label at most (PXD033741, Arg10)
    seen = {v for n, v in SILAC_GAPS.items() if share.get(n, 0) >= 0.08}
    if len(seen) >= 2 and shifted >= 500:
        med, heavy = seen & {"Lys4", "Arg6"}, seen & {"Lys8", "Arg10"}
        name = "SILAC triplex (Lys4/Arg6, Lys8/Arg10)" if med and heavy else "SILAC (Lys8/Arg10)" if heavy else "SILAC (Lys4/Arg6)"
        ev.append(Evidence("label", name, 0.7, "data", "delta-mass combinations", "SILAC labels explain the mass gaps: "
                           + ", ".join(f"{n} {share[n]:.0%}" for n in SILAC_GAPS if SILAC_GAPS[n] in seen)
                           + f" of {shifted} spectra", gaps))
    for n, enr in ENRICHMENT.items():
        if share.get(n, 0) >= 0.10:
            ev.append(Evidence("enrichment", enr, min(0.7, 0.4 + share[n]), "data", "delta-mass combinations",
                               f"{n} explains {share[n]:.0%} of {shifted} spectra with a mass gap", gaps))
    return ev


# DIAspeXtract on a constant-mobility copy: mobility tracing, weighting, charge veto and twin folding off
NO_MOBILITY_OPTIONS = ["-trace:detector", "openms", "-gate:delta_im", "1.0", "-assembly:im_weight_sigma", "0",
                       "-charge:im_charge_veto", "false", "-assembly:twin_im_tolerance", "-1",
                       "-trace:im_slice_merge", "-1"]


def with_constant_mobility(src, dst, k0=1.0):
    """Copy of a DIA mzML without ion mobility in which every peak carries 1/K0 = k0, so DIAspeXtract accepts it."""
    import numpy as np
    import pyopenms as oms
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(src), exp)
    spectra = exp.getSpectra()  # iterating an MSExperiment yields copies
    for s in spectra:
        fda = oms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.full(s.size(), k0, dtype=np.float32))
        s.setFloatDataArrays([f for f in s.getFloatDataArrays() if f.getName() != "Ion Mobility"] + [fda])
        s.setDriftTimeUnit(oms.DriftTimeUnit.VSSC)
    exp.setSpectra(spectra)
    exp.updateRanges()
    oms.MzMLFile().store(str(dst), exp)
    return dst


# Bruker Frames.MsMsType: 0 MS1, 8 ddaPASEF MS/MS, 9 diaPASEF MS/MS, 10 prm-PASEF MS/MS
TIMS_MODES = {8: "ddaPASEF", 9: "diaPASEF", 10: "prm-PASEF"}


def tims_acquisition(d):
    """Acquisition mode of a Bruker timsTOF .d from analysis.tdf: 'ddaPASEF', 'diaPASEF', 'prm-PASEF', 'MS1 only',
    or None when the file is unreadable or mixes MS/MS types. The DIA window tables exist (empty) in ddaPASEF files
    too, so the frame types decide."""
    import sqlite3
    tdf = os.path.join(d, "analysis.tdf")
    if not os.path.exists(tdf):
        return None
    try:
        db = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)
        types = {t for (t,) in db.execute("SELECT DISTINCT MsMsType FROM Frames")} - {0}
    except sqlite3.Error:
        return None
    if not types:
        return "MS1 only"
    return TIMS_MODES.get(types.pop()) if len(types) == 1 else None


def tims_instrument(d):
    """Instrument name of a Bruker .d (analysis.tdf GlobalMetadata InstrumentName, e.g. 'timsTOF Pro 2'), or None"""
    import sqlite3
    try:
        db = sqlite3.connect(f"file:{os.path.join(d, 'analysis.tdf')}?mode=ro", uri=True)
        row = db.execute("SELECT Value FROM GlobalMetadata WHERE Key = 'InstrumentName'").fetchone()
    except sqlite3.Error:
        return None
    return row[0].strip() if row and row[0] and row[0].strip() else None


def fileconverter(exe, inp, out_mzml, threads=0, log=print):
    """Bruker ddaPASEF .d -> mzML with OpenMS FileConverter (built-in timsTOF reader: one MS2 per precursor)."""
    _run([exe, "-in", inp, "-out", out_mzml] + (["-threads", threads] if threads else []), log, timeout=6 * 3600)
    return out_mzml


DIASPEXTRACT_DDA_INPUT = 11   # DIAspeXtract exit code: the input is DDA (ddaPASEF), not DIA


def diaspextract(exe, inp, out_mzml, threads=0, extra=(), log=print):
    """diaPASEF .d (or DIA mzML) -> pseudo-MS/MS mzML. A RuntimeError carries DIAspeXtract's exit code as
    .returncode (DIASPEXTRACT_DDA_INPUT: the input is DDA)."""
    cmd = [exe, "-in", inp, "-out", out_mzml] + (["-threads", threads] if threads else []) + list(extra)
    _run(cmd, log, timeout=12 * 3600)
    return out_mzml

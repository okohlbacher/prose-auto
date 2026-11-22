"""The arm command templates of PREREG 0.1.12 §4 and §7, frozen with FREEZE.md (§10 step 2). Pure functions: argument
lists from settings, reports and paths; nothing here runs a tool (harness.py does).

Arms (§4): auto, reference, sdrf, defaults, best_psm; ablations (§7): spectra_only, texts, auto_never, auto_always;
S7: fast_0110 (prose-auto 0.1.10 fast mode, timing only); scenarios (§2.3, §4 S6): one prose-auto analyze with the
scenario database as --fasta (the pre-search check), its prose_cli searched against that database (with cRAP and
entrapment), then prose-auto check (post); wrong_label and wrong_protease: auto's prose_cli changed by the rules of
v0111/postcal.py, searched against the correct database, then checked. Every ProSE search uses the corpus protocol
(corpus_ids/corpus.py PROTO, 2 missed cleavages unless unspecific) against target + reversed decoys + the seed-0
entrapment twins (entrap.makedb), with -threads THREADS."""
import re

# the corpus protocol (corpus_ids/corpus.py PROTO); kept here verbatim so the frozen templates name it
PROTO = ["-Search:modifications:variable_max_per_peptide", "1", "-Search:precursor:isotope_error_min", "-2",
         "-Search:report:top_hits", "10", "-Search:report:isotope_error_convention", "theoretical_minus_observed",
         "-Search:rescore", "true", "-no_summary"]
MAX_CHARGE = 20           # ProSE's limit (as the pilot and postcal.py)
DEFAULTS = ["-Search:enzyme", "Trypsin/P", "-Search:modifications:fixed", "Carbamidomethyl (C)",
            "-Search:modifications:variable", "Oxidation (M)"]
DEFAULT_SETTINGS = {"protease": "Trypsin/P", "cys": "Carbamidomethyl (C)", "label": "label free"}
ABLATIONS = {"spectra_only": "no --accession, no description, no Jev", "texts": "--accession, no Jev",
             "auto_never": "--mini-search-cys never", "auto_always": "--mini-search-cys always"}

# best_psm (§4): one factor at a time from the reference settings, in this order; the five non-tryptic proteases are
# those of prose-auto's vocabulary other than Lys-N (Lys-C, Arg-C, Chymotrypsin, Glu-C, Asp-N)
BEST_PSM_FACTORS = (
    ("protease", ("Trypsin", "Lys-C", "Arg-C", "Chymotrypsin", "Glu-C", "Asp-N", "unspecific cleavage")),
    ("label", ("label free", "TMT6plex", "TMTpro (16/18plex)", "iTRAQ4plex", "iTRAQ8plex", "SILAC (Lys8/Arg10)")),
    ("cys", ("Carbamidomethyl (C)", "none (free Cys)", "Propionamide (C)", "Nethylmaleimide (C)", "Methylthio (C)",
             "Carboxymethyl (C)")))
BEST_PSM_GAIN = 1.05       # a candidate is kept when it beats the current value by >= 5 % in PSMs
BEST_PSM_FDP = 1.0         # ... at 1 % FDR; a candidate above 1 % entrapment FDP scores 0,
BEST_PSM_MIN_NE = 20       # but only on >= 20 entrapment PSMs (amendment 2026-11-22c item 5: below, the FDP is noise)


def _vals(c, key):
    if key not in c:
        return []
    i, out = c.index(key) + 1, []
    while i < len(c) and not c[i].startswith("-"):
        out.append(c[i])
        i += 1
    return out


def _set(c, key, vals):
    c = list(c)
    if key in c:
        i = c.index(key) + 1
        j = i
        while j < len(c) and not c[j].startswith("-"):
            j += 1
        c[i:j] = vals
    else:
        c += [key] + vals
    return c


def clamp_charge(cli):
    """ProSE accepts charges up to MAX_CHARGE"""
    c = list(cli)
    if "-Search:precursor:max_charge" in c:
        i = c.index("-Search:precursor:max_charge") + 1
        c[i] = str(min(int(c[i]), MAX_CHARGE))
    return c


def settings_cli(s, prec=None, frag=None):
    """a ProSE command from settings (label, protease, cys, enrichment) and windows ((value, unit) or None); a field
    nobody states takes the ProSE default (pilot.ref_cli)"""
    from prose_auto import prose_params as pp
    fixed = [s["cys"]] if s.get("cys") and s["cys"] != "none (free Cys)" else []
    var = ["Oxidation (M)"]
    lab = s.get("label") or "label free"
    fixed += pp.LABEL_MODS.get(lab, [])
    var += pp.SILAC_VARIABLE.get(lab, [])
    var += [m for m in pp.ENRICHMENT_MODS.get(s.get("enrichment") or "none", []) if m not in var]
    cli = []
    if prec:
        cli += ["-Search:precursor:mass_tolerance_lower", str(prec[0]), "-Search:precursor:mass_tolerance_upper",
                str(prec[0]), "-Search:precursor:mass_tolerance_unit", prec[1]]
    if frag:
        cli += ["-Search:fragment:mass_tolerance", str(frag[0]), "-Search:fragment:mass_tolerance_unit", frag[1]]
    if s.get("protease"):
        cli += ["-Search:enzyme", pp.ENZYME_NAMES.get(s["protease"], s["protease"])]
    return cli + ["-Search:modifications:fixed", *fixed, "-Search:modifications:variable", *var]


def wrong_label(cli):
    """§2.3: the isobaric label dropped, or TMT6plex added to a label-free run (v0111/postcal.py)"""
    fx = _vals(cli, "-Search:modifications:fixed")
    lab = [m for m in fx if re.match(r"(TMT|iTRAQ)", m)]
    return _set(cli, "-Search:modifications:fixed", [m for m in fx if m not in lab] if lab else
                fx + ["TMT6plex (K)", "TMT6plex (N-term)"])


def wrong_protease(cli):
    """§2.3: Trypsin replaced by Glu-C, any other protease by Trypsin (v0111/postcal.py)"""
    e = (_vals(cli, "-Search:enzyme") or ["Trypsin"])[0]
    return _set(cli, "-Search:enzyme", ["glutamyl endopeptidase"] if e.startswith("Trypsin") else ["Trypsin"])


def prose_cmd(prose_bin, spectra, db, idxml, cli, threads):
    """one ProSE search with the corpus protocol"""
    cli = clamp_charge(cli)
    enz = (_vals(cli, "-Search:enzyme") or ["Trypsin"])[0]
    proto = PROTO + ([] if enz == "unspecific cleavage" else ["-Search:peptide:missed_cleavages", "2"])
    return [prose_bin, "-in", spectra, "-database", db, "-out_idxml", idxml, *cli, *proto, "-threads", str(threads)]


def analyze_cmd(tools, spectra, out, arm, organisms=(), accession=None, fasta=None, threads=8):
    """prose-auto analyze for one arm (auto, an ablation, or a scenario pre-search with fasta set). tools: python,
    panel (the frozen reference panel), fastag, fastag_taxonomy, diaspextract, fileconverter, crap."""
    a = [tools["python"], "-m", "prose_auto.cli", "analyze", spectra, "--panel", tools["panel"],
         "--fastag", tools["fastag"], "--fastag-taxonomy", tools["fastag_taxonomy"],
         "--diaspextract", tools["diaspextract"], "--threads", str(threads), "-o", out, "--sdrf", "off",
         "--contaminants", tools["crap"]]
    if spectra.rstrip("/").endswith(".d"):
        a += ["--fileconverter", tools["fileconverter"]]
    if arm not in ("spectra_only",) and accession:
        a += ["--accession", accession]
    if arm in ("auto", "auto_never", "auto_always", "scenario"):
        a += ["--jev"]
    if arm == "auto_never":
        a += ["--mini-search-cys", "never"]
    if arm == "auto_always":
        a += ["--mini-search-cys", "always"]
    if fasta:
        for f in fasta:
            a += ["--fasta", f]
    else:
        for o in organisms:
            a += ["--organism", o]
    return a


def check_cmd(tools, spectra, idxml, fasta, report, out, threads=8):
    """prose-auto check (the post-search check, §4 S6) on one search result; it writes <stem>.check.json in out"""
    a = [tools["python"], "-m", "prose_auto.cli", "check", spectra, idxml, "--fastag", tools["fastag"],
         "--contaminants", tools["crap"], "--threads", str(threads), "-o", out]
    for f in fasta:
        a += ["--fasta", f]
    return a + (["--report", report] if report else [])


def best_psm_step(current, scores):
    """one factor of best_psm: scores maps each candidate value (the current one included) to its subset score
    (subset_psms; None when the search failed). The best candidate replaces the current value only when it beats it
    by >= 5 %; ties go to the earlier candidate"""
    base = scores.get(current) or 0
    best, val = current, base
    for v, s in scores.items():
        if s is not None and v != current and s > val and s >= BEST_PSM_GAIN * base:
            best, val = v, s
    return best


UNSPECIFIC = "unspecific cleavage"


def best_specific(current, scores):
    """amendment 2026-11-22c: when unspecific cleavage wins the protease step, the label and Cys steps search with the
    specific protease the rule would have chosen without the unspecific candidate (an unspecific subset search peaks at
    900 GB RAM and takes about 4.5 h); unspecific is searched once more, in the final full search"""
    return best_psm_step(current, {v: s for v, s in scores.items() if v != UNSPECIFIC})


def subset_psms(ev):
    """a subset search's score from its eval.json: PSMs at 1 % FDR; 0 when the entrapment FDP exceeds 1 % on at
    least BEST_PSM_MIN_NE entrapment PSMs"""
    if not ev:
        return None
    p = ev["psm"]
    return 0 if float(p["combined"]) > BEST_PSM_FDP and p["N_E"] >= BEST_PSM_MIN_NE else p["accepted"]

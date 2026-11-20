# PREREG: prose-auto 0.1.12 with a known organism, held-out benchmark (draft of 2026-11-20)

Status: **draft for sign-off by the author.** Nothing has been run on a held-out run. No held-out raw file has been
downloaded, opened or converted. No person has read the SDRF, the project record text or a result file of a held-out
project. The run lists were drawn by `draw_sets.py` from PRIDE project-level fields and file names and sizes only. The
benchmark is scored once (§10). A later change is a dated amendment at the end of this file. The changes from the
0.1.11 draft are listed under "Amendments".

## 1. Frozen code and tools

- **prose-auto:** the commit that adds this file (`git log -n1 --format=%H -- acceptance/0.1.12/PREREG.md`),
  version 0.1.12. It freezes these thresholds as committed:
  - the database check `DBCHECK_*` (`prose_auto/external.py`);
  - the post-search check `CHECK_*` (`prose_auto/check.py`);
  - the Cys fallback (`cli._cys_fallback`).
- **FASTag** 1.15.3, Linux x64 release bundle (Thermo RawFileReader included).
- **ProSE:** the build of the ID corpus (`/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr`, OpenMS 3.7.0-pre, PR
  299c9a6c99), with the corpus protocol:
  - 2 missed cleavages (none for unspecific cleavage);
  - at most 1 variable modification per peptide;
  - isotope error down to −2;
  - Percolator rescoring.
- **DIAspeXtract, FileConverter:** as in the 0.1.11 validation (`v0111/given/*/run.sh`).
- **Jev:** the model name and version the 0.1.12 runs report (`jev.model` in every report), its parameters, and the
  cache `jev_cache_v0111`. A new request is asked 3 times (`PROSE_AUTO_JEV_VOTES`) and added to a copy of the cache;
  the copy's SHA-256 is frozen after the arms (§10).
- **Binaries:** the SHA-256 of every executable and library bundle above (FASTag bundle, ProSE, DIAspeXtract,
  FileConverter, mzpeak-convert, the Python environment's `pip freeze`) is frozen in `FREEZE.md` (§10).
- **Conversion:** the ABI PRIDE mirror's default recipe (`pxd convert --to mzML`):
  - OpenMS FileConverter, vendor build, for Thermo `.raw`;
  - mzpeak-convert for Bruker `.d`;
  - FileConverter for SCIEX `.wiff`.
- **Scoring and scenario code**, in this folder, frozen by §10 before any raw file is opened:
  - `extract.py` writes the scoring inputs (`results.tsv`, `decisions.tsv`, `checks.tsv`, `timing.tsv`) from the arm
    outputs; no input is typed by hand except `reference.tsv` (§3);
  - `score.py` validates the inputs and computes every endpoint and every §9 rule;
  - `scenarios.py` builds the wrong-database scenarios from the UniProt metadata that `resolve_uniprot.py` fetched;
  - `power.py` gives the power of the criteria (`power.txt`), with `score.py`'s own functions.

## 2. Sets

### 2.1 The held-out runs

**The draw.** The sets were drawn on 2026-11-20 with seed 20261120 by `draw_sets.py`, from PRIDE's v3 search API.
They are unchanged from the 0.1.11 draft; the files are byte-identical to `acceptance/0.1.11/` of fac1a34 (`SHA256SUMS`). No fix
of 0.1.12 required a new draw.

**Excluded from the draw:** the 2,691 accessions in `exclude.txt`. These are every PXD named in the prose-auto
development tree:
- the 218-run release set and the ID corpus;
- the mass-truth, mixture, non-tryptic, HLA, intact-protein and text-benchmark sets;
- every Jev and text cache;
- FASTag's acceptance and held-out folders;
- the ABI PRIDE mirror with its legacy archive.

Whole-PRIDE screens and catalogues (candidate screens, the bigbio clone) were not counted as use.

| File | Runs | Content |
|---|---|---|
| `heldout.tsv` | 187 | one vendor file per project, 50 MB–6 GB, by quota (below), 253 GB in total |
| `wrongdb.tsv` | 40 | 30 single-organism runs of `heldout.tsv` (seeded) and 10 runs of multi-species projects, 55 GB in total |
| `candidates.tsv` | 283 | every project the draw looked at, with the file it chose, or `-` when it had none (86 projects had no vendor file of 50 MB–6 GB) |
| `draw.log` | | queries, pools and fills |

**Quotas.** Column `stratum` of `heldout.tsv` names the stratum. A project counts for the first stratum that takes it,
rarest first.

| Stratum | Runs | Selected by |
|---|---|---|
| protease Glu-C, Asp-N, chymotrypsin, Lys-C | 10 each | the enzyme named in the project texts, trypsin not named |
| TMTpro SPS-MS3 / TMTpro / TMT / iTRAQ / SILAC | 6 / 8 / 10 / 8 / 8 | quantification method or texts |
| phospho / glyco enrichment | 8 / 6 | texts |
| Astral / timsTOF DDA / timsTOF DIA / TOF (SCIEX, Waters) / ion-trap MS2 (LTQ) | 8 / 8 / 8 / 6 / 8 | instrument field, experiment type |
| Orbitrap DIA | 10 | instrument field, experiment type |
| plant / fungus / bacterium / other vertebrate / invertebrate / mouse or rat | 6 / 5 / 6 / 5 / 4 / 6 | single organism field |
| baseline: Orbitrap DDA label-free human | 13 | instrument, organism, experiment type |

The strata come from metadata and can be wrong; every endpoint uses the reference annotation (§3) instead. A stratum
with fewer than 8 runs of its class in the reference is reported as under-filled and is not topped up. This applies to
the 5 proteases and the label classes.

### 2.2 The run classes: which runs enter which endpoint

The reference annotation (§3) puts every run in one class.

- **bottom-up:** in P1–P3 and S1–S7.
- **top-down or intact protein:** out of P1–P3. Endpoint S8 instead tests whether prose-auto stopped with the
  top-down message, as it should. Five drawn runs carry "Top-down proteomics" in the experiment field; the annotation
  decides.
- **no proteome:** a metagenome, no species-level organism, or a species without a UniProt reference proteome of its
  own (prose-auto would search a relative; the reference database would be undefined). Out of P1–P3; counted.
- **organism not verified** (§3.1): out of P1–P3; counted.

These rules apply within the eligible runs:
- **A bottom-up run where prose-auto stops:** an auto failure. Its auto arm has 0 PSMs in P1 and P2, and it is
  noncompliant in P3.
- **A run whose download or conversion fails:** missing. `score.py` reads `heldout.tsv` and counts every drawn run
  without a reference row as missing; it is out of every endpoint and not replaced. A drawn run with a reference row
  but a missing arm result is left out of the endpoints that need that arm, and counted.
- **Denominators:** every endpoint reports its denominator and the counts left out, by reason.
- **The minimum:** with fewer than 120 eligible runs (bottom-up, organism verified, not missing), P1–P3 and the study
  get no verdict (`score.MIN_ELIGIBLE`, tested at 119 and 120). `power.txt` assumes 10 % missing, about 164 eligible.

### 2.3 The wrong-database scenarios (`scenarios.py`, frozen manifest)

The scenarios are built by a taxonomy rule, which replaces the hand-made map of `v0111/dbcal.py`. That map is not used
for the held-out set, and its *C. elegans* → *Drosophila* pair was not "same family or order". For each run of
`wrongdb.tsv` the rule takes:

- **Organisms:** the PRIDE organisms field is split at ";". Common names and strain suffixes are dropped; the first two
  words are kept. A taxon above species ("Fungi", "Bacteria") or a metagenome is dropped.
- **Excluded runs:** the run is excluded when:
  - the adjudicated class (below) is not bottom-up;
  - the experiment field says top-down;
  - no species is left;
  - a name is not in the frozen taxonomy;
  - or a species has no UniProt reference proteome of its own: its correct control would be a relative.

  It is reported and not replaced. The manifest excludes PXD064921 (top-down), PXD063632 (metagenome) and PXD063182
  (*Dendrobaena veneta*, no reference proteome of its own).
- **correct:** the reference proteome of each species.
- **related_close:** each species is replaced by the reference proteome of another genus of the same family.
- **related_distant:** each species is replaced by one of another family in the same order, or else in the same class.
  Each pick records the rank of the common ancestor, so the results can be stratified by distance.
- **unrelated:**
  - *Arabidopsis thaliana* for animals, fungi and other eukaryotes;
  - *Homo sapiens* for plants, bacteria and archaea.

  The group is that of the first species.
- **incomplete** (runs with two or more species): the first species alone.
- **Ties:** the candidate proteome with the most proteins wins, then the lowest taxon id.
- **No candidate:** the scenario is recorded as "none" and is not built.

**The manifest (built 2026-11-21, UniProt release 2026_04 of 04-November-2026).**
- **Inputs:** `resolve_uniprot.py wrongdb.tsv uniprot/` queried UniProt's taxonomy and proteome metadata (no sequences,
  no raw data). It wrote:
  - `uniprot/names.dmp` and `uniprot/nodes.dmp`: the lineages of the 14 organism names and of *Arabidopsis
    thaliana* (an unrelated database);
  - `uniprot/candidates.tsv`: the reference proteomes per species and level, 2,692 rows;
  - `uniprot/queries.tsv`: all 89 request URLs, with the release and the SHA-256 of each response.
- **The manifest:** `scenarios.py wrongdb.tsv uniprot uniprot/candidates.tsv > scenarios.tsv`. Both scripts are kept, so
  it can be rebuilt; a later UniProt release can change the candidates, and the frozen `scenarios.tsv` decides.
- **What it holds:**
  - 37 runs with correct, related_close, related_distant and unrelated scenarios, plus incomplete for the 8 runs with
    several species;
  - 3 runs excluded;
  - one related_close scenario that is "none" and is not built: pig (no other genus of Suidae has a reference
    proteome);
  - related_distant falls back to the class level for *S. cerevisiae*.
- **The databases:** every built row names the UniProt proteome of each of its species (column `proteomes`), the
  unrelated and incomplete rows included; `scenarios.py` stops if a built row lacks one. The FASTA files are fetched at
  the freeze (§10), and their SHA-256 is frozen.

**When the annotation changes the organism.** The manifest is built from the PRIDE organisms. After adjudication (§3),
`adjudicated.tsv` lists every wrong-database run with its adjudicated organisms and class. `resolve_uniprot.py` and
`scenarios.py` run again with it, unchanged; a run whose organisms changed gets the scenarios of its adjudicated
organisms by the same rule, a run whose class is not bottom-up is excluded. The rebuilt `scenarios.tsv` and its
proteomes are hashed in `FREEZE.md` before any arm runs; `score.py` reads that file. Runs whose organisms did not
change keep their rows.

**Wrong settings.** These scenarios apply to the single-organism runs that are not excluded, on the correct database.
The rules are those of `v0111/postcal.py`:
- **wrong_label:** the isobaric label is dropped, or TMT6plex is added to a label-free run;
- **wrong_protease:** Trypsin is replaced by Glu-C, and any other protease by Trypsin.

**Download.** After sign-off and the freeze (§10), `pxd get <PXD> --files <file>` fetches every listed file into the
mirror.

## 3. Reference annotation (before any arm runs)

**Who annotates.** Two annotators work independently. Neither sees any prose-auto, FASTag-arm or ProSE-arm output of
these runs. Each fills `reference_A.tsv` and `reference_B.tsv` from the raw data checks below and the authors' methods.

**Adjudication.** A third person settles every disagreement and records the decision and the reason in the
adjudication log. The result is `reference.tsv`, whose columns `score.py` defines:
- run, class, organism_verified;
- label, protease, cys;
- reporter_ions (yes/no);
- the windows with their units.

The adjudicated organisms and class of the wrong-database runs go into `adjudicated.tsv`. Agreement per setting
(Cohen's κ) is reported.

**Freezing (§10, step 4).** Both annotations and the adjudication log are frozen first, then `reference.tsv` and
`adjudicated.tsv`, all before the first arm runs.

1. **Organism and database: verified from the spectra.**
   - **The search:** a ProSE search of a seeded 5,000-MS2 subset with the stated settings, against a labelled
     concatenation of the stated proteome(s), the human, mouse, yeast, E. coli and *Arabidopsis* proteomes, and cRAP.
   - **Verified** when two conditions hold:
     - at least 100 target PSMs pass 1 % FDR;
     - at least 80 % of those not explained by cRAP fall on a peptide of the stated proteome(s). A shared peptide counts
       as stated.
   - **Another organism:** when the spectra point to one, the annotators record it, and it becomes the reference after
     adjudication.
   - **Mixtures:** each stated component must have at least 5 % of those PSMs uniquely.
   - **Otherwise** the run is "organism not verified" (§2.2).
   - **The database:** the frozen UniProt reference proteome of each species, plus cRAP.
2. **Categorical settings:** label, plex, protease, Cys alkylation and enrichment. Each annotator decides from:
   - the methods;
   - the raw-data indicators:
     - reporter channels in MS2/MS3 (manual check of 200 spectra);
     - SPS-MS3 from the scan headers;
     - SILAC pairs in MS1;
     - the diagnostic ions.

   The search ladder of the 0.1.11 draft is no longer the truth (see the best_psm arm, §4).
3. **Windows:** a wide-window search on the subset:
   - 50 ppm / 0.05 Da for high-resolution MS2, 0.6 Da for ion trap;
   - the precursor window is |median| + 4 robust SD;
   - the fragment window comes from a Gaussian-plus-uniform fit;
   - the same rules as the 0.1.10 mass truth.
4. **Class:** bottom-up, top-down or intact protein, or no proteome (§2.2; this includes a species without a UniProt
   reference proteome of its own).

## 4. Arms

Every arm searches the same run with the reference database (§3.1), cRAP and entrapment. The entrapment database is a
shuffled twin of every target, seed 0; the decoys are reversed. The ProSE build and protocol are those of §1.

| Arm | Settings |
|---|---|
| **auto** | prose-auto 0.1.12 default, with its `prose_cli`. Flags: `--organism` = the reference organism(s), `--accession`, `--jev`, `--sdrf off`, `--contaminants cRAP` |
| **reference** | the adjudicated annotation (§3.2) with the measured windows (§3.3) |
| **SDRF** | the authors' settings as stated: the submitted or curated SDRF if one exists, else the methods as transcribed by the annotators. A field nobody states takes the ProSE default, never auto's |
| **defaults** | ProSE defaults plus the organism: Trypsin/P, Carbamidomethyl (C) fixed, Oxidation (M) variable, ProSE's default windows and charges |
| **best_psm** (comparator, not truth) | on the subset, one factor at a time from the reference settings:<br>• protease (trypsin, the five non-tryptic proteases, unspecific);<br>• label (none, TMT, TMTpro, iTRAQ4/8, SILAC);<br>• Cys alkylation (CAM, none, propionamide, NEM, methylthio, carboxymethyl).<br>A candidate is kept when it beats the current value by ≥ 5 % in PSMs at 1 % entrapment FDP. The full run is then searched with the result |

**Wrong-database and wrong-settings searches (S6).** Every scenario (§2.3) is searched on its own database:
1. **Pre-search:** prose-auto analyze runs with the auto flags, with the scenario database given as `--fasta` instead of
   the reference organism. Its `database_check` is the pre-search verdict.
2. **The search:** ProSE searches the run with that analysis's `prose_cli`, against the scenario database plus cRAP and
   entrapment.
3. **Post-search:** `prose-auto check` runs on that result, with that FASTA and that report.

The correct-database result of the auto arm is never checked against another FASTA. For wrong settings, the scenario
settings are searched against the correct database and then checked.

## 5. Primary endpoints

**Counts:** at 1 % FDR (Percolator q-values): PSMs, unique peptides and picked-protein groups.

**Paired ratios** (`score.ratios`), per eligible run. Eligibility is by PSMs for every count (PSMs, peptides, protein
groups):
- an arm result missing: left out and counted;
- a stopped arm has 0 PSMs (0 peptides, 0 groups);
- both arms below 100 PSMs: left out and counted;
- the denominator below 100 PSMs and the numerator at or above: +∞ (the numerator is better);
- otherwise the ratio of the count; for peptides and groups a denominator of 0 gives +∞, and 0 against 0 is left out.

**The conventions** (`score.median`, `score.quantile`; `power.py` uses the same functions):
- **median of ratios:** the sample median of the ratios in [0, +∞]; for an even count the arithmetic mean of the two
  middle ones (+∞ if either is +∞). No value is clipped or replaced;
- **interval:** percentile bootstrap over runs, 10,000 resamples, `numpy.random.default_rng(20261120)`. A quantile is
  an order statistic (`inverted_cdf`, no interpolation), so 0 and +∞ stay exact;
- **two-sided 95 %:** the 2.5th and 97.5th percentiles; **one-sided 95 %:** the 5th percentile.

| # | Endpoint | Criterion |
|---|---|---|
| P1 | median per-run PSM ratio auto / reference | **non-inferior:** lower end of the two-sided 95 % interval ≥ 0.98 |
| P2 | median per-run PSM ratio auto / defaults | **superior:** lower end of the two-sided 95 % interval > 1 |
| P3 | share of runs at entrapment FDP ≤ 1.5 % (combined estimate, Wen et al. 2025) on both PSM and peptide level, over the eligible runs whose reference arm identifies (≥ 100 PSMs); an auto arm stopped or below 100 PSMs is noncompliant | **both gates:** (a) relative: the one-sided 95 % lower bound of auto's share minus the reference's share ≥ −0.05 (paired bootstrap); (b) absolute: the one-sided 95 % Clopper-Pearson lower bound of auto's share ≥ 0.75 |

**P3, the two claims.** Both are claimed and both must hold:
- **(a) relative:** auto is not materially worse than the reference settings, within 5 percentage points of the share
  of compliant runs;
- **(b) absolute:** at least 75 % of runs keep FDP ≤ 1.5 %.

**Where X = 0.75 comes from:** it is the largest bound in steps of 0.05 that the reference arm meets on development data
(share 0.836, lower bound 0.787, n = 207).

**What changed against the 0.1.11 draft:** that draft required ≥ 95 % of the runs, with no bound. On the development
corpus only 80–84 % of runs meet the limit, in every arm, so the old rule passed in 0 % of the simulations. That rule is
replaced and disclosed.

**Power** (`power.txt`, 1,000 simulations):

**Method:** each simulation draws the held-out design by quota from matching development pools:
- protease strata from the non-tryptic sets;
- Orbitrap DIA, Astral, timsTOF DDA/DIA, TOF and ion trap from their instrument pools;
- all other strata from Orbitrap DDA.

It then removes the 5 runs tagged top-down and makes each other run missing with probability 0.10 (about 164
eligible). It applies `score.py` unchanged, with 2,000 resamples.

| Endpoint | Expected pairs | Development effect | Effect shifted |
|---|---|---|---|
| P1 | 135 | 1.00 | auto PSMs × 0.98: 1.00; × 0.97: 0.53 |
| P2 | | unknown: only 16 development data sets have a defaults arm | |
| P3 (a) | 131 | 0.06 | auto as compliant as the reference: 1.00 |
| P3 (a) and (b) | 131 | 0.03 | auto as compliant as the reference: 0.77; 2 % / 5 % fewer compliant auto runs: 0.01 / 0.00 |

**Why P3's power is low on development data.** The development auto arm is prose-auto 0.1.3. In 8 of its runs it
finds fewer than 100 PSMs where the optimal settings identify, and those runs count as noncompliant. The development
share is therefore auto 0.802 against the reference's 0.836. The study passes only if the held-out auto arm is about as
compliant as the reference.

**Reported, not gating:** the peptide and protein-group ratios (auto/reference, auto/defaults) and every ratio
auto/other arm (SDRF, best_psm, the ablations of §7), with the same intervals.

**The verdict:** the study passes if P1, P2 and P3 (both gates) all hold. The three are required jointly
(intersection-union), so they are not adjusted for multiplicity.

## 6. Secondary endpoints

Against the reference annotation; exact matches unless stated. Shares come with one-sided 95 % Wilson bounds.

- **S1 label:** chemistry and plex exactly (TMT10 ≠ TMT11 ≠ TMTpro16 ≠ TMTpro18; SILAC with its channels). The family
  level (isobaric / SILAC / none) is a second line.
- **S2 protease:** the exact enzyme (Lys-C ≠ Trypsin).
- **S3 Cys alkylation.**
- **S4 mass error:** auto's precursor and fragment windows within 1.5× of the wide-window measurement (§3.3). Each
  window carries its unit; a window in Da and one in ppm are compared after conversion at m/z 1,000
  (`score.MZ_REF`; 0.02 Da = 20 ppm).
- **S1–S4, missing decisions:** an eligible run with a reference value but no auto decision counts as wrong and is
  reported as absent.
- **S5 abstention** per setting: decisions with status `default`, `unknown` or `review`.
- **S6 wrong database and wrong settings:** per check (pre, post) and per kind, over the built scenarios of the frozen
  manifest (a scenario without a check result is missing, not assessed is left out; both are counted, neither is ever a
  detection):
  - **the kinds:** correct, related_close, related_distant, unrelated, incomplete, wrong_label, wrong_protease;
  - **sensitivity:** the share flagged;
  - **specificity:** the share ok on the correct database, from the auto arm's own checks on every eligible run;
  - **exactness:** the share with the right verdict (database or settings);
  - **also reported:**
    - the identifications the wrong database loses;
    - the seconds each check costs;
    - which check catches what.
  - **natural failures:** runs whose SDRF or defaults arm finds fewer than half the reference PSMs, with the check run
    on that arm's result. Reported separately, never pooled with the constructed scenarios.
- **S7 run time** per input type (Orbitrap DDA, DIA, timsTOF, TOF, ion trap) on an idle node: prose-auto 0.1.12 total
  and per step, against 0.1.10 fast mode (arm `fast_0110`) on the same files. Descriptive (median seconds,
  `score.timing`), no criterion.
- **S8 top-down:** among the runs annotated top-down or intact protein, the share where prose-auto stopped; among the
  bottom-up runs, the number stopped.

## 7. Ablations (held-out set; identifications where the ProSE command changes)

- spectra only (no `--accession`, no description, no Jev);
- spectra plus texts (`--accession`, no Jev);
- spectra plus texts plus Jev (= auto);
- the mini-search fallback: `--mini-search-cys never`, `fallback` (= auto) and `always` (the 0.1.10 competition);
- `-alkylation_out` against the competition: the `never` and `always` arms on the runs whose Cys decision differs.

## 8. Comparators

- **Param-Medic:** tolerances and its label check, on the same mzML.
- **MSFragger 4.4.1:** mass calibration (precursor and fragment) on the same mzML and database.
- **Engine defaults:** the defaults arm (ProSE); MSFragger's default parameters as a second line.
- **best_psm** (§4).

## 9. Decision rules (fixed now)

A detector becomes or stays the default only if it passes its own criterion. Otherwise the decision abstains to the
texts or to a conservative default, and the report says so. `score.py` computes every rule (`score.rules`,
`score.checks`).

| Detector | Stays default if | Otherwise |
|---|---|---|
| FASTag `-alkylation_out` instead of the competition | arms `auto_never` vs `auto_always` (§7), on the eligible runs whose Cys decision differs: summed PSM ratio ≥ 0.99 and no run below 0.95 | the competition becomes the default again (0.1.10 behaviour) |
| mini-search fallback | arms `auto` (fallback) vs `auto_always`: as above | `always` |
| label | S1 exact ≥ the always-label-free baseline (the number of eligible runs whose reference label is none, i.e. the score of a constant "label free" decision), and no isobaric call where the reference records no reporter ions (`reporter_ions = no`) | no isobaric label without reporter ions |
| protease | S2 exact ≥ the stated-protease baseline (the SDRF arm's protease scored against the reference) | the stated protease decides |
| FASTag QC evidence | the MS2 analyser class of FASTag's QC agrees with the header's in ≥ 98 % of the eligible runs with both | removed from the evidence |
| database checks | the gates below | that check's verdicts are marked experimental |

**The check gates.** Each gate is an exact one-sided binomial test of H0 "share ≤ bound", run on the assessed
scenarios. The p-values of all nine gates are adjusted together by Holm at 0.05.

| Check | Kind | Bound | Development share | Held-out n (expected) | Power under Holm |
|---|---|---|---|---|---|
| pre | correct (specificity) | 0.70 | 0.863 (158/183) | 141 | 1.00 |
| pre | unrelated | 0.50 | 0.838 (140/167) | 32 | 0.99 |
| pre | related_close (same family) | 0.45 | 0.779 (81/104) | 31 | 0.98 |
| pre | related_distant (same order or class) | 0.45 | 0.778 (7/9) | 30 | 0.97 |
| post | correct (specificity) | 0.75 | 0.869 (159/183) | 142 | 0.97 |
| post | unrelated | 0.75 | 0.988 (163/165) | 32 | 1.00 |
| post | related_distant | 0.20 | 0.556 (5/9) | 30 | 0.99 |
| post | wrong_label | 0.65 | 0.960 (24/25) | 25 | 0.99 |
| post | wrong_protease | 0.65 | 1.000 (25/25) | 25 | 1.00 |

**How the bounds were set** (`power.py`):
- **The bound:** each is the largest bound, in steps of 0.05 from 0.10, with single-gate power ≥ 0.80 at 0.05/9
  (Bonferroni, which Holm dominates), given at least 10 expected assessed scenarios.
- **The cap:** wrong_protease is capped at wrong_label's bound, because its development 25/25 sits at the ceiling.
- **The development shares:** the related ones come from the 0.1.11 calibration, split by the rank of the common
  ancestor (`v0112/related_strata.txt`). The 9 related_distant pairs make that bound uncertain.

**Experimental (reported with interval, never a gate):** post related_close, development 0.252, which reaches no bound
≥ 0.10 with power; pre and post incomplete, 7 expected scenarios; natural failures.

**The family of decisions:** each row of the detector table is a separate deployment decision at its own criterion. No
joint confirmatory claim is made across them. Within the check gates, Holm controls the family-wise error at 0.05.

**After the score.** A default changed by these rules after the held-out score has, from then on, no independent
evidence. It needs a new held-out set before it is claimed.

## 10. Freeze and scoring once

The order is fixed. Nothing of a held-out raw file is downloaded or opened before step 3.

1. **Sign-off** of this file by the author.
2. **The freeze**, before any raw file:
   - fetch the FASTA of every reference organism and every scenario of `scenarios.tsv` (its UniProt proteome ids);
   - write `FREEZE.md` with the SHA-256 of:
     - `extract.py`, `score.py`, `scenarios.py`, `resolve_uniprot.py` and `power.py`;
     - `heldout.tsv`, `wrongdb.tsv` and `scenarios.tsv`;
     - `uniprot/` (lineages, candidates, the request log);
     - every FASTA;
     - the binaries and environment of §1 (FASTag bundle, ProSE, DIAspeXtract, FileConverter, mzpeak-convert,
       `pip freeze` of the Python environment) and the Jev model, version, parameters and cache;
     - the conversion recipe (the `pxd` recipe id and the hash of its `tools.json`);
     - the arm command templates and the output layout of `extract.py`;
   - commit it.
3. **Download and convert:** the SHA-256 of every mzML is appended to `FREEZE.md`.
4. **Annotation** (§3):
   - `reference_A.tsv`, `reference_B.tsv` and the adjudication log (each disagreement, the decision, the reason) are
     frozen;
   - then `reference.tsv` and `adjudicated.tsv` are frozen;
   - then the manifest is rebuilt from `adjudicated.tsv` (§2.3), and its FASTA files are fetched and hashed.
5. **Arms and scenario searches** (§4):
   - the SHA-256 of every arm output file is appended to `FREEZE.md`;
   - so is the Jev cache after the arms.
6. **Scoring:**
   - `python extract.py OUT DIR` writes the scoring inputs;
   - the SHA-256 of the five inputs (`reference.tsv`, `results.tsv`, `decisions.tsv`, `checks.tsv`, `scenarios.tsv`)
     and of `timing.tsv` is appended to `FREEZE.md`;
   - then `python score.py DIR` runs once.

   A bug found afterwards is fixed in a dated amendment that reports both results.

**Tests:** `tests/test_prereg.py` covers:
- `extract.py`;
- the 120-run minimum (119 against 120 runs);
- run accounting and strict inputs (blank values, unknown runs, duplicates);
- the ratio conventions with 0 and +∞;
- PSM eligibility for groups;
- units;
- stopped runs in P3;
- the absolute gate;
- Holm and the label rule;
- the scenario rules, including adjudication.

On development data (`v0112/dev_inputs.py`, in-sample), `score.py` runs end to end.

## 11. What is not claimed

- No result transfers to search engines other than ProSE.
- The strata are metadata-based; a class under-filled after annotation is reported, not claimed.
- The development results (`v0111/README_0111.md`, the 218 release runs) are in-sample and are not evidence for this
  benchmark.
- The power figures assume that the held-out runs behave like the development data, resampled by quota from pools of
  the matching instrument class. The development auto arm was prose-auto 0.1.3 and the reference was the optimal
  settings, not the adjudicated annotation; P2 has no power estimate.

## Amendments

### 2026-11-20, 0.1.12 (before sign-off; from the adversarial review of the 0.1.11 draft)

1. **S6:** a search per wrong-database scenario, with that FASTA, before both checks (§4). The 0.1.11 text could be
   read as checking the correct-database result against a wrong FASTA.
2. **Wrong-database scenarios:** a taxonomy rule (`scenarios.py`) for all 40 runs, with the manifest built now
   (`scenarios.tsv`, UniProt 2026_04). Multi-species
   runs get related, unrelated and incomplete scenarios; runs with unresolvable names are excluded by rule (§2.3).
   Related scenarios now come at two distances (close, distant). The hand-made `dbcal.py` map is retired.
3. **Top-down, stopped and missing runs:** handled by class and rule (§2.2); S8 is new.
4. **P3:** now a paired non-inferiority against the reference arm with a one-sided 95 % lower bound. Power for
   P1–P3 and the checks is in `power.txt` (§5).
5. **Checks:** a criterion with a one-sided Wilson bound per kind, set from the development power (§9). The 0.1.11
   draft pooled sensitivity at 0.70 and specificity at 0.85, with no bounds.
6. **Freeze order:** `score.py` (now written and tested), `scenarios.py`, the conversion recipe and the manifests are
   frozen and hashed before any raw file is downloaded or opened. In the 0.1.11 draft they were frozen before the first
   arm (§10).
7. **The truth:**
   - the organism is verified from the spectra (§3.1);
   - two blinded annotators decide, with adjudication (§3);
   - the best-PSM settings are a comparator arm, not the truth (§4).
8. **The sets** are unchanged (same seed, same files); `acceptance/0.1.11/` became `acceptance/0.1.12/`.

### 2026-11-21, 0.1.12 (before sign-off; from the second adversarial review)

1. **Scorer, inputs:** `score.py` reads `heldout.tsv` and accounts for every drawn run, and it validates every input:
   - columns, enumerations, numbers, duplicate keys and runs outside the draw are checked;
   - a blank value is never a result or a detection.

   The 120-run minimum is implemented as "no verdict".
2. **Scorer, endpoints:**
   - **stopped auto arms:** noncompliant in P3;
   - **eligibility:** decided by PSMs for every count, peptides and groups included;
   - **ratio conventions:** median and quantiles over 0 and +∞ are frozen, nothing is clipped, and `power.py` uses the
     same functions;
   - **units:** windows carry units;
   - **missing decisions:** counted as wrong.
3. **P3 has two gates**, relative (−5 points) and absolute (lower bound ≥ 0.75), and both are claimed (§5).
4. **§9:** every rule is computed by `score.py`:
   - the label rule against a defined baseline (a constant "label free" decision) and the reporter-ion record;
   - the protease rule against the SDRF arm;
   - the fallback and Cys rules on their arms;
   - the QC agreement.

   S7 is computed as a descriptive endpoint.
5. **Check gates:**
   - nine gates with bounds from power at 0.05/9, Holm over them, and a post-search related gate (related_distant,
     0.20);
   - related_close (pre) is now a gate;
   - post related_close and incomplete are experimental, with the reason given;
   - the family-wise policy is stated.
6. **Scenarios:**
   - every built row names its proteomes (unrelated and incomplete included);
   - a species without a reference proteome of its own excludes the run (*Dendrobaena veneta*, PXD063182);
   - after adjudication the manifest is rebuilt by the same scripts and hashed before any arm. The rebuilt manifest
     has 37 runs and 3 exclusions.
7. **Power:**
   - quota-stratified resampling and endpoint-specific eligible counts;
   - sensitivity at weaker effects;
   - P2 is labelled unknown;
   - the related power comes from the order and class strata;
   - P3's low power on the 0.1.3 development proxy is disclosed.
8. **Freeze list extended:**
   - the binaries and environment;
   - the Jev model and cache;
   - `extract.py` and the output layout;
   - both annotations and the adjudication log;
   - the arm-output hashes;
   - the hashes of the final scoring inputs.

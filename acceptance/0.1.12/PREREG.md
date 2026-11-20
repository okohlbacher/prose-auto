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
- **DIAspeXtract, FileConverter, Jev cache:** as in the 0.1.11 validation (`v0111/given/*/run.sh`). Jev answers are
  cached; a new request is asked 3 times (`PROSE_AUTO_JEV_VOTES`).
- **Conversion:** the ABI PRIDE mirror's default recipe (`pxd convert --to mzML`):
  - OpenMS FileConverter, vendor build, for Thermo `.raw`;
  - mzpeak-convert for Bruker `.d`;
  - FileConverter for SCIEX `.wiff`.
- **Scoring and scenario code:**
  - `score.py` computes every endpoint;
  - `scenarios.py` builds the wrong-database scenarios from the UniProt metadata that `resolve_uniprot.py` fetched;
  - `power.py` gives the power of the criteria (`power.txt`).

  All three are in this folder and are frozen by §10 before any raw file is opened.

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
- **no proteome:** a metagenome, or no species-level organism. Out of P1–P3; counted.
- **organism not verified** (§3.1): out of P1–P3; counted.

These rules apply within the eligible runs:
- **A bottom-up run where prose-auto stops:** an auto failure. Its auto arm has 0 PSMs in P1–P3.
- **A run whose download or conversion fails:** missing. It is out of every endpoint and counted, and it is not
  replaced.
- **Denominators:** every endpoint reports its denominator and the counts left out, by reason.
- **The minimum:** with fewer than 120 eligible runs, P1–P3 are reported without a pass or fail verdict.
  `power.txt` assumes 164 runs.

### 2.3 The wrong-database scenarios (`scenarios.py`, frozen manifest)

The scenarios are built by a taxonomy rule, which replaces the hand-made map of `v0111/dbcal.py`. That map is not used
for the held-out set, and its *C. elegans* → *Drosophila* pair was not "same family or order". For each run of
`wrongdb.tsv` the rule takes:

- **Organisms:** the PRIDE organisms field is split at ";". Common names and strain suffixes are dropped; the first two
  words are kept. A taxon above species ("Fungi", "Bacteria") or a metagenome is dropped.
- **Excluded runs:** the run is excluded when:
  - the experiment field says top-down;
  - no species is left;
  - or a name is not in the frozen taxonomy.

  It is reported and not replaced. The manifest excludes PXD064921 (top-down) and PXD063632 (metagenome).
- **correct:** the reference annotation's database.
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

**The manifest (built 2026-11-20, UniProt release 2026_04 of 04-November-2026).**
- **Inputs:** `resolve_uniprot.py wrongdb.tsv uniprot/` queried UniProt's taxonomy and proteome metadata (no sequences,
  no raw data). It wrote:
  - `uniprot/names.dmp` and `uniprot/nodes.dmp`: the lineages of the 14 organism names;
  - `uniprot/candidates.tsv`: the reference proteomes per species and level, 2,327 rows;
  - `uniprot/queries.tsv`: all 81 request URLs, with the release and the SHA-256 of each response.
- **The manifest:** `scenarios.py wrongdb.tsv uniprot uniprot/candidates.tsv > scenarios.tsv`. Both scripts are kept, so
  it can be rebuilt; a later UniProt release can change the candidates, and the frozen `scenarios.tsv` decides.
- **What it holds:**
  - 38 runs with correct, related_close, related_distant and unrelated scenarios, plus incomplete for the 8 runs with
    several species;
  - 2 runs excluded;
  - two related_close scenarios that are "none" and are not built: pig (no other genus of Suidae has a reference
    proteome) and the earthworm *Dendrobaena veneta*;
  - related_distant falls back to the class level for *Dendrobaena* and for *S. cerevisiae*.
- **The databases:** each scenario names its UniProt proteome ids. The FASTA files are fetched at the freeze (§10), and
  their SHA-256 is frozen.

**Wrong settings.** These scenarios apply to the single-organism runs that are not excluded, on the correct database.
The rules are those of `v0111/postcal.py`:
- **wrong_label:** the isobaric label is dropped, or TMT6plex is added to a label-free run;
- **wrong_protease:** Trypsin is replaced by Glu-C, and any other protease by Trypsin.

**Download.** After sign-off and the freeze (§10), `pxd get <PXD> --files <file>` fetches every listed file into the
mirror.

## 3. Reference annotation (before any arm runs)

**Who annotates.** Two annotators work independently. Neither sees any prose-auto, FASTag-arm or ProSE-arm output of
these runs. Each fills `reference_A.tsv` and `reference_B.tsv` from the raw data checks below and the authors' methods.

**Adjudication.** A third person settles every disagreement, recording the decision and the reason, in
`reference.tsv`. Agreement per setting (Cohen's κ) is reported.

**Freezing.** `reference.tsv` is frozen with its SHA-256 before the first arm runs.

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
4. **Class:** bottom-up, top-down or intact protein, or no proteome (§2.2).

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

**Paired ratios** (`score.ratios`), per eligible run:
- a stopped arm has 0 PSMs;
- both arms below 100 PSMs: left out and counted;
- the denominator below 100 and the numerator above it: +∞ (the numerator is better). The median is rank-based.

**Intervals:** percentile bootstrap over runs, 10,000 resamples, `numpy.random.default_rng(20261120)`.

| # | Endpoint | Criterion | Power (`power.txt`) |
|---|---|---|---|
| P1 | median per-run PSM ratio auto / reference | **non-inferior:** lower end of the two-sided 95 % interval ≥ 0.98 | 1.00 if the true median equals the development one (1.014) or is 2 % lower; 0.22 if 3 % lower; 0.01 if 3.5 % lower |
| P2 | median per-run PSM ratio auto / defaults | **superior:** lower end of the two-sided 95 % interval > 1 | 1.00 (development: 16 data sets only) |
| P3 | share of runs at entrapment FDP ≤ 1.5 % (combined estimate, Wen et al. 2025) on both PSM and peptide level, auto minus reference, paired over the runs with ≥ 100 PSMs in both | **non-inferior:** one-sided 95 % lower bound of the difference (5th bootstrap percentile) ≥ −0.05 | 0.81 |

**Why P3 changed.** The 0.1.11 draft required ≥ 95 % of the runs at FDP ≤ 1.5 %, with no bound. On the development
corpus only 82–85 % of runs meet the limit in every arm: auto 84.6 %, optimal settings 83.6 % (n = 214). The excess FDP
belongs to the protocol, not to the settings, so the old rule passed in 0 % of the simulations. The absolute share of
auto with its one-sided 95 % Clopper-Pearson bound is reported; it does not gate.

**Power.** The power comes from the development paired data:
- the ID corpus: prose-auto 0.1.3 against the optimal settings, 228 pairs;
- the 16-data-set defaults search.

Each simulation draws 164 eligible pairs with replacement and applies `score.py` with 2,000 resamples.

**Reported, not gating:** the peptide and protein-group ratios (auto/reference, auto/defaults, auto/SDRF) and the PSM
ratio auto/best_psm, with the same intervals.

**The verdict:** the study passes if P1, P2 and P3 all hold.

## 6. Secondary endpoints

Against the reference annotation; exact matches unless stated. Shares come with one-sided 95 % Wilson bounds.

- **S1 label:** chemistry and plex exactly (TMT10 ≠ TMT11 ≠ TMTpro16 ≠ TMTpro18; SILAC with its channels). The family
  level (isobaric / SILAC / none) is a second line.
- **S2 protease:** the exact enzyme (Lys-C ≠ Trypsin).
- **S3 Cys alkylation.**
- **S4 mass error:** auto's precursor and fragment windows within 1.5× of the wide-window measurement (§3.3).
- **S5 abstention** per setting: decisions with status `default`, `unknown` or `review`.
- **S6 wrong database and wrong settings:** per check (pre, post) and per kind:
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
  and per step, against 0.1.10 fast mode on the same files.
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
texts or to a conservative default, and the report says so.

| Detector | Stays default if | Otherwise |
|---|---|---|
| FASTag `-alkylation_out` instead of the competition | `never` vs `always`, on the runs whose Cys decision differs: summed PSM ratio ≥ 0.99 and no run below 0.95 | the competition becomes the default again (0.1.10 behaviour) |
| mini-search fallback | `fallback` vs `always`: as above | `always` |
| label | S1 exact ≥ the "always label-free" baseline, and no isobaric call without reporter ions | no isobaric label without reporter ions |
| protease | S2 ≥ the stated-protease baseline | the stated protease decides |
| pre-search database check | one-sided 95 % Wilson lower bounds: specificity ≥ 0.75, unrelated ≥ 0.65, related_distant ≥ 0.50 | the warning is marked experimental |
| post-search check | specificity ≥ 0.75, unrelated ≥ 0.85, wrong_label ≥ 0.70, wrong_protease ≥ 0.70 | the verdicts are reported as experimental |
| FASTag QC evidence | analyser class agrees with the header in ≥ 98 % of runs | removed from the evidence |

**How the check bounds were set.** Each is the largest bound (in steps of 0.05) that has ≥ 80 % power at the
development share and the expected held-out count (`power.txt`). The wrong-settings bounds are capped at 0.70, because
the development shares (24/25, 25/25) are near the ceiling. The related_distant bound comes from the development
related scenarios, of which 104 of 166 pairs were in the same genus or family (`v0112/related_strata.txt`), so it
is conservative for the distant level. Three kinds are reported only:
- related_close (no development data);
- incomplete (8 runs);
- natural failures.

**After the score.** A default changed by these rules after the held-out score has, from then on, no independent
evidence. It needs a new held-out set before it is claimed.

## 10. Freeze and scoring once

The order is fixed. Nothing of a held-out raw file is downloaded or opened before step 3.

1. **Sign-off** of this file by the author.
2. **The freeze**, before any raw file:
   - fetch the FASTA of every reference organism and every scenario of `scenarios.tsv` (its UniProt proteome ids);
   - write `FREEZE.md` with the SHA-256 of:
     - `score.py`, `scenarios.py`, `resolve_uniprot.py` and `power.py`;
     - `heldout.tsv`, `wrongdb.tsv` and `scenarios.tsv`;
     - `uniprot/` (lineages, candidates, the request log);
     - every FASTA;
     - the conversion recipe (the `pxd` recipe id and the hash of its `tools.json`);
     - the arm command templates;
   - commit it.
3. **Download and convert:** the SHA-256 of every mzML is appended to `FREEZE.md`.
4. **Annotation** (§3): `reference.tsv` is frozen.
5. **Arms and scenario searches** (§4).
6. **Scoring:** `python score.py DIR` runs once on `reference.tsv`, `results.tsv`, `decisions.tsv` and `checks.tsv`.
   A bug found afterwards is fixed in a dated amendment that reports both results.

`score.py` was tested on development data (in-sample). On the ID corpus, the defaults search and the 0.1.11 check
calibration (`v0112/dev_inputs.py`), it reproduces the development figures: post-search ok 159/183 and unrelated
163/165.

## 11. What is not claimed

- No result transfers to search engines other than ProSE.
- The strata are metadata-based; a class under-filled after annotation is reported, not claimed.
- The development results (`v0111/README_0111.md`, the 218 release runs) are in-sample and are not evidence for this
  benchmark.
- The power figures assume that the held-out runs behave like the development data. The development auto arm was
  prose-auto 0.1.3.

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

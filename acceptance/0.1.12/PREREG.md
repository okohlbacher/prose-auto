# PREREG: prose-auto 0.1.12 with a known organism, held-out benchmark (draft of 2026-11-20)

Status: **draft for sign-off by the author.** Nothing has been run on a held-out run. No held-out raw file has been
downloaded, opened or converted. No person has read the SDRF, the project record text or a result file of a held-out
project. The run lists were drawn by `draw_sets.py` from PRIDE project-level fields and file names and sizes only. The
benchmark is scored once (§10). A later change is a dated amendment at the end of this file. The changes from the
0.1.11 draft are listed under "Amendments".

## 0. Pilot (2026-11-21; independent of the held-out sets; no verdict)

The pilot is a small separate study that estimates the real power before sign-off. It sets the sample size and the
P3 gates (§5). It is not evidence for the benchmark.

**The draw.** `pilot/draw_pilot.py` drew 19 runs, one per PRIDE project, on 2026-11-21. It used seed 20261121 and the
held-out predicates and queries, one or two runs per stratum.
- **Exclusions:** every PXD in `exclude.txt`, `heldout.tsv`, `wrongdb.tsv` and `candidates.tsv`.
- **Asp-N:** no project was left, so the pilot has no Asp-N run.
- **Later draws:** `pilot/pilot.tsv` is added to the exclusions of every future draw.

**The procedure** (`pilot/pilot.py`):
1. **Download:** pxd, with PRIDE MD5 checks, into a private pxd root on node-local scratch. The group mirror's Ceph
   quota was full at the time; it was freed later the same day (§10.1).
2. **Conversion:** the mirror's default recipes.
3. **Annotation:** one annotator worked from the methods, the 4 SDRFs and the raw-data indicators (MS levels,
   analysers, reporter-ion shares, `analysis.tdf` acquisition). The annotation (`pilot/reference.tsv`) was written before
   any prose-auto output existed.
4. **Organism check:**
   - FASTag's database-free genus ranking; a stated genus must be in the top 5;
   - where FASTag cannot assess the run (DIA or Bruker input, an empty ranking, an isobaric label), the reference arm
     must reach ≥ 100 PSMs instead.
5. **Arms:** auto (prose-auto 0.1.12 release), reference (the annotation with wide-window measured windows, or the
   instrument standard for DIA and timsTOF), SDRF (the SDRF, else the methods), and defaults.
6. **Searches:** the corpus ProSE protocol with seed-0 entrapment.
7. **Scoring:** `extract.py` and `score.py --pilot`.

**Results** (`pilot/score/score.md`). Of 19 drawn runs, 12 were eligible:
- **Missing (2), no tool reads the input:**
  - PXD058482: diaPASEF; DIAspeXtract rejects two windows with the same key;
  - PXD074553: ddaPASEF; OpenTIMS does not read TDF compression 1.
- **Organism not verified (5):**
  - two runs have 0 PSMs in every arm (a targeted run with 1,217 MS2; an ion-trap run);
  - one IP is too small;
  - one E. coli expression run: FASTag ranks plant genera first;
  - one mouse run with SILAC E. coli: E. coli is not in FASTag's top 5.
- **The metadata strata were often wrong:**
  - four runs drawn as Astral, TOF, Orbitrap DIA or iTRAQ are timsTOF;
  - one drawn as ion trap is an Orbitrap DIA.
- **Annotation error:** on PXD012418 auto found 4,688 PSMs against the reference's 474. The annotation had Glu-C for a
  Lys-C + Glu-C digest.

| Endpoint | Pilot (12 runs) |
|---|---|
| P1 median PSM ratio auto / reference | 1.098 (95 % interval 1.051–1.282); 10 gain, 2 lose |
| P2 median PSM ratio auto / defaults | 1.250 (lower 1.137) |
| P3 share at FDP ≤ 1.5 % | auto 11/12, reference 10/12; difference +0.083 (lower −0.167); auto's lower bound 0.661 |
| peptides / protein groups auto / reference | 1.031 / 1.003 |
| S1–S5 | label 12/12, protease 11/12, Cys 10/11, precursor window 9/12, fragment window 9/12 |

**Power from the pilot** (`power.txt`, `power_grid.txt`; amendment 2026-11-22):
- **Method:**
  - each simulated held-out set takes the manifest's 300 runs. The 8 tagged top-down and the 17 that
    `refdb_draw.tsv` marks no-proteome are fixed out of P1–P3;
  - each of the other 275 is drawn from the 19 pilot runs as a whole: its class (12 eligible, 5 organism not verified,
    2 missing) and all its arm results together;
  - each set goes through `score.accounting` and `score.primary` unchanged.
- **No stratification:** the pilot's metadata strata were often wrong and most strata hold 0 or 1 pilot run, so no
  stratification is used. Stratum-specific power is unknown.
- **The figures are conditional on 12 eligible pilot runs.** They are a rough guide.

| Scenario at 300 drawn (X 0.75, margin −0.05; 1,000 simulations) | Eligible runs (5–95 %) | P(no verdict) | Joint power |
|---|---|---|---|
| conditional: the pilot as the truth | 174 (160–187) | 0.00 | **0.97** |
| auto PSMs ×0.97 | 174 (160–186) | 0.00 | 0.98 |
| 2 % / 5 % of compliant auto runs noncompliant | 174 | 0.00 | 0.90 / 0.71 |
| predictive: the 19 pilot runs resampled in each simulation | 175 (123–224) | 0.04 | **0.64** |
| attrition at the pilot's 95 % lower bound (eligible share 0.384), conditional | 105 (93–117) | 0.97 | 0.03 |
| the same, predictive | 106 (73–138) | 0.74 | 0.17 |

- **Project sensitivity:**
  - leaving one pilot project out gives joint power from 0.46 (without PXD024383 or PXD065759) to 1.00;
  - one pilot FDP sits next to the limit: PXD024383's reference PSM FDP is 1.509 %. Moving it to 1.500 % drops the
    conditional joint power from 0.97 to 0.50;
  - in each case P3's relative gate moves, because the paired difference rests on 2 of the 12 runs.
- **What this means:** at the pilot's own effect and attrition, 300 drawn runs suffice. Once the pilot's uncertainty is
  carried forward, the joint power is about 0.64, mostly through P3. If the held-out attrition is as high as the
  pilot's lower bound allows, the study gets no verdict in most simulations. **Power under plausible attrition is
  well below 0.8.** The draw and the gates stay as decided (§0, the author's decision); the risk is stated here, not
  removed. The stated assurance is the predictive power (§9), and the conditional top-up (§2.1) guards the 120-run
  minimum.
- **The design grid** (`power_grid.txt`, conditional, 200 simulations per cell), joint power at the pilot effect /
  ×0.98 with 2 % fewer compliant / ×0.97 with 5 % fewer:

| Drawn runs | X | Margin | Pilot effect | ×0.98, 2 % | ×0.97, 5 % | P(no verdict) |
|---|---|---|---|---|---|---|
| 187 (first draw) | 0.75 | −0.05 | 0.07 | 0.04 | 0.02 | 0.94 |
| 250 | 0.75 | −0.05 | 0.96 | 0.82 | 0.63 | 0.00 |
| **300 (adopted)** | **0.75** | **−0.05** | **0.97** | **0.93** | **0.66** | **0.00** |
| 300 | 0.75 | −0.10 | 0.99 | 0.99 | 0.98 | 0.00 |
| 300 | 0.70 | −0.10 | 1.00 | 0.99 | 0.94 | 0.00 |
| 400 | 0.75 | −0.05 | 0.99 | 0.95 | 0.78 | 0.00 |

**Adopted by the author on 2026-11-21:**
1. **300 drawn runs:** a top-up of 113 runs by `draw_topup.py` (§2.1).
2. **X = 0.75 and the margin at −0.05**, unchanged.
3. **The check bounds were recomputed at 300** by the frozen rule (§9).

**Caveats:**
- The pilot is small: 12 runs, one annotator, and a simplified organism check.
- Its loss rate (37 %) carries a wide interval, and some losses are pilot-specific:
  - the organism check is cruder than §3.1;
  - two runs fail in the format readers.
- So the power above is a rough guide, not a guarantee.

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

**The draw.** Two steps, both from PRIDE's v3 search API:
- **2026-11-20:** 187 runs, seed 20261120, by `draw_sets.py`. These sets are unchanged from the 0.1.11 draft; the files
  were byte-identical to `acceptance/0.1.11/` of fac1a34.
- **2026-11-21, the top-up to 300 (§0):** 113 runs, seed 20261122, by `draw_topup.py`:
  - the same strata, queries, predicates and file rule;
  - each stratum's 187-run quota scaled to 300 by largest remainder, minus the runs already drawn;
  - the exclusions: `exclude.txt`, the 187 runs, `wrongdb.tsv`, `candidates.tsv` and the pilot's draw and candidates;
  - where a pool runs out, its unfilled places pass one at a time to the strata in draw order. No Asp-N project was
    left (6 places), and the label-free baseline pool held 5 of 8. The 9 places went to Glu-C, chymotrypsin, Lys-C,
    TMTpro SPS-MS3, TMTpro, TMT, iTRAQ, SILAC and glyco;
  - the 113 rows are appended to `heldout.tsv`, whose first 187 rows are unchanged. `draw_topup.py` refuses to run
    unless `heldout.tsv` is still the 187-run file.

The wrong-database set and the scenarios are unchanged: `scenarios.py` rebuilds `scenarios.tsv` byte for byte.

**The conditional top-up** (amendment 2026-11-22b):
- **When it is checked:** once, after the blind annotation (§3) and the download, conversion and MS2-count check
  (§10 step 3) of all 300 runs, before any arm runs.
- **The count:** `score.eligible_before_arms` counts the eligible runs from `reference.tsv` and `conversion.tsv`
  only: bottom-up, a verified organism, a reference database (not no-proteome), and a complete conversion. It reads
  nothing of any arm or endpoint. It stops when a drawn run is neither annotated and converted nor failed.
- **The rule:** with fewer than 120 eligible runs, `draw_topup.py --batch K DIR` draws batch K of
  min(50, 500 − drawn) runs:
  - the same strata, queries, predicates and file rule;
  - the 187-run quotas scaled to the new total, minus the runs each stratum already has (at least 0); unfilled
    places pass on as before;
  - the exclusions of the top-up plus every drawn run and every candidate so far;
  - the seeds, fixed now: batch 1 20261201, batch 2 20261202, batch 3 20261203, batch 4 20261204.
- **After a batch:** its runs are annotated and converted by the same protocol, and the count is repeated.
- **When it ends:**
  - at ≥ 120 eligible runs, when the arms start;
  - at the cap of 500 drawn runs, or when every pool is exhausted, with "no verdict" for P1–P3 and the study (§2.2).
- **Safeguards:** `draw_topup.py` refuses a batch unless the count is below 120, `heldout.tsv` holds exactly the
  runs of the batches before (300 + 50 (K − 1)) and the cap is open. Each batch's runs, candidates, log and the new
  `heldout.tsv` hash go into `FREEZE.md` before its annotation.

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
| `heldout.tsv` | 300 | one vendor file per project, 50 MB–6 GB, by quota (below), 399 GB in total (187 runs of 253 GB, then the top-up of 113 runs of 146 GB) |
| `wrongdb.tsv` | 40 | 30 single-organism runs of `heldout.tsv` (seeded) and 10 runs of multi-species projects, 55 GB in total |
| `candidates.tsv` | 283 | every project the first draw looked at, with the file it chose, or `-` when it had none (86 projects had no vendor file of 50 MB–6 GB) |
| `topup.tsv`, `topup_candidates.tsv` | 113, 159 | the top-up runs (also in `heldout.tsv`) and every project it looked at |
| `draw.log`, `draw_topup.log` | | queries, pools and fills |
| `refdb_draw.tsv` | 300 | the reference-database manifest from PRIDE's organism field (`refdb.py`, UniProt snapshot `uniprot_heldout/`); rebuilt from `adjudicated.tsv` after the annotation (§3.1) |

**Quotas** (runs of the first draw; with the top-up, see `draw_topup.log`). Column `stratum` of `heldout.tsv` names the stratum. A project counts for the first stratum that takes it,
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
  top-down message, as it should. Eight drawn runs (five of the first draw, three of the top-up) carry "Top-down
  proteomics" in the experiment field; the annotation decides.
- **no proteome:** a metagenome, no species-level organism, or a species without a UniProt reference proteome of its
  own (prose-auto would search a relative; the reference database would be undefined). Out of P1–P3; counted.
- **organism not verified** (§3.1): out of P1–P3; counted.

These rules apply within the eligible runs:
- **A bottom-up run where prose-auto stops:** an auto failure. Its auto arm has 0 PSMs in P1 and P2, and it is
  noncompliant in P3.
- **A run whose download or conversion fails:** missing. `score.py` reads `heldout.tsv` and counts every drawn run
  without a reference row as missing; it is out of every endpoint and not replaced.
- **A partial conversion is a failed conversion.** At step 3 (§10) the mzML must carry every spectrum of the vendor
  file: its spectrum count at every MS level and its MS2 count must equal the vendor file's own counts, read
  directly from the vendor file by the frozen FASTag (RawFileReader for Thermo). Bruker `.d` runs are read natively
  and need MS2 spectra. SCIEX `.wiff` has no independent reader: the converter's validated index decides
  (amendment 2026-11-22c). Otherwise the run is missing. A format no frozen reader opens
  (as the pilot's TDF compression 1 and the duplicate DIA window key) is missing.
- **A run that did not execute:** a verified bottom-up run whose auto or reference arm has no result, or a missing
  one, after the reruns of §10.1 is missing, before the minimum below (`score.accounting`, "not executed"). Failed
  runs never count toward the 120.
- **Other missing arms:** a missing SDRF, defaults, best_psm or ablation result leaves the run out of the endpoints
  that need that arm, and is counted.
- **Which arms run where:** top-down runs run the auto arm only (S8). No-proteome runs and runs whose organism is not
  verified run no arm; they are counted. Every verified bottom-up run runs every arm.
- **Denominators:** every endpoint reports its denominator and the counts left out, by reason.
- **The minimum:** with fewer than 120 eligible runs (bottom-up, organism verified, not missing), P1–P3 and the study
  get no verdict (`score.MIN_ELIGIBLE`, tested at 119 and 120). `power.txt` gives the expected eligible count and
  P(no verdict) under the pilot's attrition and under higher attrition (§5).

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
`adjudicated.tsv` lists every wrong-database run with its adjudicated organisms and class. The frozen scripts then run
again, in this order:
1. `resolve_uniprot.py ... adjudicated.tsv --snapshot=uniprot/` queries UniProt only for names the snapshot never
   resolved, and keeps every snapshot row;
2. `scenarios.py ... adjudicated.tsv scenarios.tsv` copies verbatim the frozen rows of every run whose adjudicated
   organisms equal its PRIDE organisms (same species list) and whose class stays bottom-up;
3. a run whose organisms changed gets the scenarios of its adjudicated organisms by the same rule, against the
   snapshot plus the new names; a run whose class is not bottom-up is excluded.

The rebuilt `scenarios.tsv` and its proteomes are hashed in `FREEZE.md` before any arm runs; `score.py` reads that
file (tested: `test_scenarios_keep_frozen_rows`, `test_scenarios_follow_the_adjudication`).

**Wrong settings.** These scenarios apply to the single-organism runs that are not excluded, on the correct database.
The rules are those of `v0111/postcal.py`:
- **wrong_label:** the isobaric label is dropped, or TMT6plex is added to a label-free run;
- **wrong_protease:** Trypsin is replaced by Glu-C, and any other protease by Trypsin.

**Download.** After sign-off and the freeze (§10), `pxd get <PXD> --files <file>` fetches every listed file into the
mirror.

## 3. Reference annotation (before any arm runs)

**Who annotates** (amendment 2026-11-22c; `heldout/annotate.py`):
- **The annotators:** A and B are two independent agent sessions on different model families: A is Claude, B is
  Codex, with Kimi as B's fallback when Codex gives no answer twice. No human annotates, adjudicates or reviews.
- **What each one sees:** one run's annotation sheet at a time, inside its prompt (`harness.py sheet`): the PRIDE
  project texts, the SDRF rows that name that file (none when no row names it) and the raw-data indicators. It runs
  from an empty temporary folder holding only its prompt, with every tool disabled (Claude) or in a read-only sandbox
  (Codex); Kimi runs only with a tool-free agent profile, else it is refused. The annotators run on the author's
  laptop, which holds only the prompts; `annotate.py` refuses to run where the held-out outputs are reachable. It never
  sees prose-auto's output, an arm, or the other annotator's answer before the reconciliation round, which goes to
  the provider that answered round 1.
- **The output:** each fills `reference_A.tsv` or `reference_B.tsv` in a fixed vocabulary (`annotate.VOCAB`).
- **The raw-data check of §3.2:** this is the measured indicators (reporter-ion shares in 200 MS2/MS3 spectra, MS
  levels, analyser prefixes, Bruker frame types). No annotator inspects spectra by eye.

**Disagreements** (`annotate.py`). The rule uses no arm value:
- **Reconciliation:** for each run with a disputed field, A and B each see the other's answer and evidence, still
  blind to prose-auto, and may revise. There is one round.
- **A field still disputed after that round is "not annotated":**
  - organisms or class disputed: the run is "organism not verified", out of P1–P3 and counted;
  - any other field disputed: it is blank in `reference.tsv`, so the run leaves that S-endpoint only. A disputed SDRF
    field is unstated (the ProSE default, §4). A disputed instrument class gives no reference windows, so the
    reference arm takes ProSE's default windows.
- **Agreement:** raw agreement and Cohen's κ are reported per field, for the first round and after reconciliation
  (`agreement.tsv`).
- **The files:** the result is `adjudicated.tsv`, whose consensus row comes from the two reconciled answers.
  `harness.py reference` then writes `reference.tsv`, the columns `score.py` defines:
  - run, class, organism_verified;
  - label, protease, cys;
  - reporter_ions (yes/no);
  - the windows with their units.

The adjudicated organisms and class of **every** drawn run and every run of `wrongdb.tsv` go into `adjudicated.tsv`
(accession, file, organisms, class): a row for each, also when nothing changed. `refdb.py` (with `heldout.tsv`) and
`scenarios.py` stop when a run has no row. Agreement per setting (Cohen's κ) is reported.

**Freezing (§10, step 4).** Both first-round annotations, both reconciled ones and the agents' raw replies are frozen
first, then `adjudicated.tsv`, `agreement.tsv` and `reference.tsv`, all before the first arm runs.

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
   - **The database:** the UniProt reference proteome of each species, plus cRAP. `refdb.py` writes the
     reference-database manifest `refdb.tsv` for every drawn run (organisms, taxon ids, proteome ids; class
     no-proteome when a species has none of its own) from the adjudicated organisms. `resolve_uniprot.py --snapshot`
     runs first on all drawn runs. All primary arms (auto's `--organism`, reference, SDRF, defaults, best_psm) search
     that run's manifest databases. `refdb.tsv` and its FASTA files are hashed before any arm runs.
   - **One taxon once:** a species named twice in a record (two names of one taxon) is listed once, so a mixture's
     database holds each proteome once (`refdb.py`).
   - **The auto arm is bound to the manifest:** its `--organism` names come from `refdb.tsv`, and it runs with
     `--panel` set to the frozen reference panel, one FASTA per manifest proteome, named by species. prose-auto then
     takes the panel FASTA and fetches nothing (`_organism_given` matches the panel first). Before an arm's search,
     the task checks that every FASTA the analysis report lists (`fasta`, `union`) is a frozen manifest FASTA, or
     the union of exactly that run's, by SHA-256. A mismatch is a failed task.
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
shuffled twin of every target, seed 0. The decoys are ProSE's own (`-Search:decoys generate`, keep_termini, as in
the pilot and the corpus): ProSE removes the reversed decoys written into the database and regenerates them
(amendment 2026-11-22c). The ProSE build and protocol are those of §1.

| Arm | Settings |
|---|---|
| **auto** | prose-auto 0.1.12 default, with its `prose_cli`. Flags: `--organism` = the reference organism(s), `--accession`, `--jev`, `--sdrf off`, `--contaminants cRAP` |
| **reference** | the adjudicated annotation (§3.2) with the measured windows (§3.3) |
| **SDRF** | the authors' settings as stated: the submitted or curated SDRF if one exists, else the methods as transcribed by the annotators. A field nobody states takes the ProSE default, never auto's |
| **defaults** | ProSE defaults plus the organism: Trypsin/P, Carbamidomethyl (C) fixed, Oxidation (M) variable, ProSE's default windows and charges |
| **best_psm** (comparator, not truth) | on the subset, one factor at a time from the reference settings:<br>• protease (trypsin, the five non-tryptic proteases, unspecific);<br>• label (none, TMT, TMTpro, iTRAQ4/8, SILAC);<br>• Cys alkylation (CAM, none, propionamide, NEM, methylthio, carboxymethyl).<br>A candidate is kept when it beats the current value by ≥ 5 % in PSMs at 1 % entrapment FDP. If unspecific wins the protease step, label and Cys are chosen under the best specific protease (amendment 2026-11-22c). The full run is then searched with the result |

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
| P3 | share of runs at entrapment FDP ≤ 1.5 % (combined estimate, Wen et al. 2025) on both PSM and peptide level, over the eligible runs whose reference arm identifies (≥ 100 PSMs); an auto arm stopped or below 100 PSMs is noncompliant | **both gates:** (a) relative: the one-sided 95 % lower bound of auto's share minus the reference's share ≥ −0.05 (paired bootstrap of the signed per-run counts +1/0/−1, compared in integers: the boundary −0.05 itself passes); (b) absolute: the one-sided 95 % Clopper-Pearson lower bound of auto's share ≥ 0.75 |

**P3, the two claims.** Both are claimed and both must hold:
- **(a) relative:** auto is not materially worse than the reference settings, within 5 percentage points of the share
  of compliant runs;
- **(b) absolute:** at least 75 % of runs keep FDP ≤ 1.5 %.

**Where X = 0.75 comes from:** it is the largest bound in steps of 0.05 that the reference arm meets on development data
(share 0.836, lower bound 0.787, n = 207).

**What changed against the 0.1.11 draft:** that draft required ≥ 95 % of the runs, with no bound. On the development
corpus only 80–84 % of runs meet the limit, in every arm, so the old rule passed in 0 % of the simulations. That rule is
replaced and disclosed.

**Power.** The power comes from the pilot (§0: 19 drawn runs, prose-auto 0.1.12 against the annotated reference),
resampled as whole runs through `score.py`'s own decision routine (`score.primary`).
- **At 300 drawn,** X 0.75, margin −0.05:
  - conditional joint power is 0.97, and 0.71 when 5 % of the compliant auto runs become noncompliant;
  - predictive joint power is 0.64, with P(no verdict) 0.04;
  - at the pilot's lower attrition bound, P(no verdict) is 0.74–0.97;
  - leaving one pilot project out gives 0.46–1.00.
- **The first draw alone** (187) would have had P(no verdict) 0.94.

The development corpus (prose-auto 0.1.3 against the optimal settings) gave a P3 power of 0.03. It is replaced by the
pilot as the power basis and kept only as a disclosure.

**Reported, not gating:** the peptide and protein-group ratios (auto/reference, auto/defaults) and every ratio
auto/other arm (SDRF, best_psm, the ablations of §7), with the same intervals.

**The verdict:** the study passes if P1, P2 and P3 (both gates) all hold. The three are required jointly
(intersection-union), so they are not adjusted for multiplicity.

## 6. Secondary endpoints

Against the reference annotation; exact matches unless stated. Shares come with one-sided 95 % Wilson bounds.

- **S1 label:** the category, at the granularity prose-auto reports (`score.LABEL_CANON`, its schema vocabulary):
  - none: "label free", "none", "no isobaric label" and "no SILAC" are one category;
  - TMT2plex, TMT6plex, TMT10/11plex (TMT10 and TMT11 together), TMTpro (16/18plex) (TMTpro16 and 18 together) and
    TMT0;
  - iTRAQ4plex and iTRAQ8plex;
  - each SILAC channel set;
  - 15N and 2H2O.

  The annotators use the same names. TMT10 and TMT11, or TMTpro16 and 18, are not told apart (prose-auto does not
  report it). The family level (isobaric / SILAC / none; 15N and 2H2O stay themselves) is a second line
  (`score.label_family`). Each canonical name is also accepted in its short form (TMT0, 15N, 2H2O).
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
    - the identifications the wrong database loses: per kind, the median ratio of the scenario search's PSMs to the
      correct scenario's. It comes from `scenario_counts.tsv`, every built scenario's search whatever its checks
      gave. Scenarios without a count and controls without a count or with 0 PSMs are counted;
    - the seconds each check costs: per kind, the median (pre: the report's `dbcheck` step; post: the harness's
      `check_seconds`);
    - the correct-scenario controls (checks on the correct database of each wrong-database run), per check against
      the manifest's built correct scenarios, apart from the auto arm's specificity;
    - which check catches what.
  - **natural failures:** the eligible runs whose SDRF or defaults arm finds strictly fewer than half the reference
    PSMs, with the check run on that arm's result. The cohort comes from the search counts (`score.natural_cohort`). A
    cohort arm without a check is missing; checks outside the cohort are counted and ignored. Reported separately,
    never pooled with the constructed scenarios.
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
| FASTag `-alkylation_out` instead of the competition | arms `auto_never` vs `auto_always` (§7), on the eligible runs whose Cys decision differs: summed PSM ratio ≥ 0.99 and no run below 0.95. A run with 0 PSMs in both arms is uninformative (left out, counted, never evidence); x/0 with x > 0 is +∞. Fewer than 5 informative runs: no verdict, the default stays and the report says why | the competition becomes the default again (0.1.10 behaviour) |
| mini-search fallback | arms `auto` (fallback) vs `auto_always`: as above | `always` |
| label | S1 exact ≥ the always-label-free baseline (the number of eligible runs whose reference label is none, i.e. the score of a constant "label free" decision), and no isobaric call where the reference records no reporter ions (`reporter_ions = no`) | no isobaric label without reporter ions |
| protease | on the same runs, those with a reference protease and an SDRF-arm protease: auto's exact hits ≥ the SDRF arm's (the stated-protease baseline); runs without an SDRF protease are counted; none left: no verdict | the stated protease decides |
| FASTag QC evidence | FASTag QC's own MS2 analyser call (its ≥ 90 % rule) agrees with the header's own evidence (the scan header or `analysis.tdf`, read from the report's evidence before reconciliation, never the reconciled decision) in ≥ 98 % of the eligible runs with header evidence. A run where FASTag QC made no call, or did not run, stays in the denominator as a disagreement; calls, no-calls and runs without header evidence are reported | removed from the evidence |
| database checks | the gates below | that check's verdicts are marked experimental |

**The stated assurance for P1–P3** (`power_decision.txt`, `power_decision.py`, amendment 2026-11-22b). The
predictive power, with the 19 pilot runs resampled in each simulation, is the power this study claims:
- **The adopted design** (300 drawn, X 0.75, margin −0.05):
  - joint power 0.64, P(no verdict) 0.05;
  - 0.60 and 0.52 with 2 % and 5 % fewer compliant auto runs;
  - a minimum of 0.45 with one pilot project left out (PXD065759);
  - at the pilot's lower attrition bound, 0.18, with P(no verdict) 0.74, before the conditional top-up of §2.1.
- **The limit is the pilot's uncertainty about the P3 effect, not the size.** P1 and P2 reach 0.95–1.00 from 300
  drawn runs, while P3's paired compliance difference rests on 12 pilot runs and two projects.

| Drawn | Margin | Joint (P1 / P2 / P3) | P(no verdict) | Eligible 5–95 % | 2 % / 5 % fewer compliant | One project left out, minimum |
|---|---|---|---|---|---|---|
| **300** | **−0.05** | **0.64 (0.95 / 0.92 / 0.66)** | **0.05** | **119–226** | **0.60 / 0.52** | **0.45** |
| 300 | −0.10 | 0.75 (0.96 / 0.94 / 0.76) | 0.04 | 122–226 | 0.68 / 0.65 | 0.57 |
| 400 | −0.05 | 0.69 (1.00 / 0.96 / 0.72) | 0.00 | 156–302 | 0.65 / 0.53 | 0.51 |
| 400 | −0.10 | 0.80 (1.00 / 0.98 / 0.81) | 0.00 | 164–297 | 0.73 / 0.69 | 0.64 |
| 500 | −0.05 | 0.74 (1.00 / 0.98 / 0.76) | 0.00 | 200–376 | 0.69 / 0.59 | 0.49 |
| 500 | −0.10 | 0.82 (1.00 / 0.98 / 0.83) | 0.00 | 211–385 | 0.76 / 0.69 | 0.64 |

1,000 simulations per 300 cell and 500 per 400 and 500 cell, 200 per left-out project. The author kept 300 drawn
runs and the −0.05 margin, with the conditional top-up.

**The check gates.** Each gate is an exact one-sided binomial test of H0 "share ≤ bound", run on the assessed
scenarios. The p-values of all nine gates are adjusted together by Holm at 0.05. The family is always all nine: a
gate without an assessed scenario keeps p = 1 and does not pass.

| Check | Kind | Bound | Development share | Held-out n (expected, conservative) | Exact power at 0.05/9 (≤ power under Holm) |
|---|---|---|---|---|---|
| pre | correct (specificity) | 0.65 | 0.863 (158/183) | 63 | 0.918 |
| pre | unrelated | 0.50 | 0.838 (140/167) | 32 | 0.938 |
| pre | related_close (same family) | 0.45 | 0.779 (81/104) | 31 | 0.872 |
| pre | related_distant (same order or class) | 0.50 | 0.838 (rank-weighted; order 6/7, class 1/2) | 30 | 0.899 |
| post | correct (specificity) | 0.65 | 0.869 (159/183) | 63 | 0.937 |
| post | unrelated | 0.75 | 0.988 (163/165) | 32 | 0.943 |
| post | related_distant | 0.25 | 0.568 (rank-weighted; order 4/7, class 1/2) | 30 | 0.824 |
| post | wrong_label | 0.65 | 0.960 (24/25) | 25 | 0.924 |
| post | wrong_protease | 0.65 | 1.000 (25/25) | 25 | 1.000 |

**How the bounds were set** (`power.py`, amendment 2026-11-22):
- **The bound:** each is the largest bound, in steps of 0.05 from 0.10, with exact single-gate power ≥ 0.80 at 0.05/9
  (Bonferroni), given at least 10 expected assessed scenarios. The power is an exact binomial sum, not a simulation.
  Holm rejects every hypothesis Bonferroni rejects, so each exact power is a lower bound of that gate's power under
  Holm, whatever the dependence between the checks. The simulated "power under Holm" of earlier drafts assumed
  independent gates and is dropped.
- **Specificity at high attrition:** the expected count is 73 eligible runs, the 5th percentile in the predictive
  scenario at the pilot's lower attrition bound (§0). Times the development share assessed, that gives 63. At that
  size the rule gives 0.65 for both checks; at the conditional count (about 174 eligible) it would give 0.75. The
  development rows are the 0.1.11 correct-database scenario checks (source scenario). Development has no auto-arm
  checks. The held-out gate uses the auto arm's own checks (§6).
- **related_distant:** the development share is weighted by the manifest's rank mix (35 order, 2 class; the first
  species' common ancestor) with the rank-specific development shares (`power_inputs/related_ranks.tsv`, from
  `v0112/related_strata.txt`). The class rank has 2 development pairs and 2 held-out scenarios. At that rank alone
  (0.5) the exact power would be 0.003 (pre, at 0.50) and 0.572 (post, at 0.25). The 9 development pairs make these
  bounds uncertain.
- **The cap:** wrong_protease is capped at wrong_label's bound, because its development 25/25 sits at the ceiling.
- **Inputs:** `power_inputs/` holds the calibration rows `power.py` read (`dev_checks.tsv`, `related_ranks.tsv`), frozen
  with it.

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
     - `extract.py`, `score.py`, `scenarios.py`, `resolve_uniprot.py`, `refdb.py` and `power.py`, with `power.txt`,
       `power_grid.txt` and `power_inputs/` (the calibration rows `power.py` read);
     - `heldout.tsv` (300 runs), `topup.tsv`, `wrongdb.tsv`, `scenarios.tsv` and `refdb_draw.tsv`, with
       `draw_sets.py` and `draw_topup.py` (with the conditional top-up's seeds);
     - `power_decision.py` and `power_decision.txt`;
     - `uniprot/` and `uniprot_heldout/` (lineages, candidates, the request log);
     - every FASTA;
     - the binaries and environment of §1 (FASTag bundle, ProSE, DIAspeXtract, FileConverter, mzpeak-convert,
       `pip freeze` of the Python environment) and the Jev model, version, parameters and cache; also the frozen
       reference panel (§3.1) and the BinAC2 container image (§10.1);
     - the conversion recipe (the `pxd` recipe id and the hash of its `tools.json`);
     - the arm command templates and the output layout of `extract.py` (its docstring; `checks.tsv` keys each check
       by run, source (auto, scenario, natural), arm, kind and check; FDP values keep full precision);
   - commit it.
3. **Download and convert:** the SHA-256 of every mzML is appended to `FREEZE.md`, with each run's MS2 count
   against its vendor index (§2.2: a partial conversion is missing).
4. **Annotation** (§3), then the conditional top-up (§2.1), whose batches are annotated the same way:
   - `reference_A.tsv`, `reference_B.tsv`, `reconciled_A.tsv`, `reconciled_B.tsv` and every raw agent reply are
     frozen;
   - then `reference.tsv` and `adjudicated.tsv` are frozen;
   - then both manifests are rebuilt from `adjudicated.tsv`, which has a row for every drawn and every
     wrong-database run: the scenarios (`scenarios.py`, §2.3) and the reference databases (`refdb.py` with
     `heldout.tsv`, §3.1, as `refdb.tsv`);
   - their FASTA files and the frozen reference panel built from `refdb.tsv` are fetched and hashed, and `refdb.tsv`
     and the rebuilt `scenarios.tsv` are hashed, all before any arm runs.
5. **Arms and scenario searches** (§4):
   - the SHA-256 of every arm output file is appended to `FREEZE.md`;
   - so is the Jev cache after the arms.
6. **Scoring:**
   - `python extract.py OUT DIR` writes the scoring inputs;
   - the SHA-256 of the six inputs (`reference.tsv`, `results.tsv`, `decisions.tsv`, `checks.tsv`,
     `scenario_counts.tsv`, `scenarios.tsv`) and of `timing.tsv` is appended to `FREEZE.md`;
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
- the fourth-review regressions: failed runs and the minimum, the natural cohort, scenario counts, the QC no-call,
  label aliases, strict JSON, refdb coverage and taxa, power through the scorer, and the output guard.

On development data (`v0112/dev_inputs.py`, in-sample), `score.py` runs end to end.

### 10.1 Execution plan (where each step runs; nothing started)

| Step (§10) | Where | Load |
|---|---|---|
| 3. Download and conversion | IBMI: `pxd get` and `pxd convert` into the shared PRIDE mirror `/ceph/ibmi/abi/data/PXD` (the mirror's recipe), on a claimed IBMI node | 300 vendor files, 399 GB, and their mzML |
| 4. Annotation sheets (`harness.py sheet`) | IBMI | one sheet per run, from PRIDE and the converted input |
| 4. Annotators A and B, reconciliation (`annotate.py`) | the author's laptop, where the agent CLIs (Claude, Codex, Kimi) run; the sheets are copied there and the answers back | two agent calls per run, a third and fourth for each disputed run |
| 4. Verification and window searches (`harness.py verify`, `windows`) | BinAC2, array A | 300 tasks, one subset search and one mini-search each |
| 5a. Every prose-auto `analyze`: the auto arm, the ablations of §7, the S6 pre-search analyses | IBMI, on a claimed node, 8 at a time with 8 threads | about 2,000 analyses at a median 145 s (pilot): about 10 h. They call Jev, whose key reaches only IBMI and only over ssh stdin, so they do not run on BinAC2 |
| 5b. Every ProSE search: auto, reference, SDRF, defaults, best_psm (subset and full), the ablation searches where the command changes, the S6 scenario searches | BinAC2, arrays B (one task per verified bottom-up held-out run, §2.2, 6 h), B2 (its best_psm, on high-memory nodes: 1,000 GB, 2 days) and C (one task per scenario run of `wrongdb.tsv`) | about 30 searches per run and 213 scenario searches |
| 5c. Post-search checks (`prose-auto check`, which runs FASTag) | BinAC2, in the task of their search | — |
| 6. `extract.py`, `score.py` | IBMI | minutes |

**The mirror's quota** (`/ceph/ibmi/abi`, one CephFS quota realm of 80 TiB):
- **Then and now:** full on 2026-11-21 during the pilot (22 GB free, §0). After a cleanup on 2026-11-21 it had
  3.88 TB free (2026-11-22). The 300 vendor files and their mzML need about 1.2 TB.
- **The 50 GB margin:** `pxd get` refuses a download unless the realm keeps 50 GB free beyond twice its size.
- **Fallback when it refuses:** a private pxd root on the claimed node's local scratch, as in the pilot, with the same
  PRIDE MD5 checks and the same recipe. The files then go to BinAC2 from there.

**The order of a run** (per run; nothing of a run starts before the steps before it are done and hashed):
1. **IBMI:**
   - conversion and its MS2-count check (§2.2);
   - every analysis of 5a, which writes per run and arm `cmd.json` and the report `*.auto.json` into
     `/ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/out/<accession>/<file>/<arm>/`;
   - their SHA-256 is appended to `FREEZE.md`.
2. **Transfer:** `rsync` from that absolute source path to `binac2j:/pfs/10/project/bw26g016/prose-auto-0112/in/`,
   with the mzML (`pxd stage --to binac2j:…`) and the run's manifest FASTA files. BinAC2 checks the SHA-256 of each
   received file against `FREEZE.md`.
3. **Release:** a run's array-B task is released only when all its files have arrived and been checked. Slurm holds
   it until then; the release list is a file of task ids.
4. **Back to Ceph:** results return by `rsync` from an absolute source path and are hashed (§10 step 5).

**The hash check on both sites:** before any step, a task compares the SHA-256 of every executable and bundle it
runs with `FREEZE.md`, and stops on a mismatch, which is a failed task:
- **IBMI:** FileConverter and mzpeak-convert for conversion. For analysis: prose-auto (the git tree hash and `pip
  freeze`), FASTag, ProSE, DIAspeXtract, the Jev model name and version, and the Jev cache.
- **BinAC2:** the container image, ProSE, DIAspeXtract, the FASTag bundle (for `prose-auto check`), prose-auto, the
  `cmd.json` and the FASTA files.

**Failures, the same rule on both sites:** a failed step (download, conversion, analysis, search or check, on IBMI or
on BinAC2) is rerun once. A Slurm task is rerun with twice the memory and time. A run whose step fails again is
missing (§2.2), and its `status` file names the step.

**BinAC2 (RV bw26g016):**
- **Paths:**
  - inputs: `/pfs/10/project/bw26g016/prose-auto-0112/in`; outputs: `.../out`. The project directory is permanent;
    each task's outputs are copied there at its end;
  - intermediates go to the node-local `$TMPDIR`;
  - the container and the frozen tree live in the workspace `proseauto`, with a copy in
    `/pfs/10/project/bw26g016/prose-auto-0112/tree`.
- **The workspace:**
  - owner `tu_iisko01`; extended on 2026-11-22 to 2026-12-21, with 4 extensions of 30 days left (`ws_extend proseauto 30`);
  - it is extended whenever fewer than 7 days remain while runs are pending;
  - if it lapses, the tasks run from the project copy of the tree, which has the same hashes;
  - the workspace `prose-heldout` (expires 2026-12-14) is not used.
- **Environment:** the Ubuntu 24.04 container of the October 2026 development run (below), with the frozen IBMI tree
  bound at its own paths.
- **Slurm:** partition `compute`, arrays with at most 40 tasks running. Per task:
  - `--cpus-per-task=32`, that is 16 physical cores, since hyperthreading is on;
  - ProSE `-threads 16`, with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` all 16;
  - `--mem=128gb`, and `--time=06:00:00`;
  - the rerun asks for 256 GB and 12 h, on the same partition;
  - the memory figures: the pilot's largest search peaked at 18 GB. In the October 2026 development run on BinAC2
    (`bench2`, 127 runs, job 6309342), 34 tasks ran out of memory at 96 GB, with the conversion inside the task.
    All 127 projects are in `exclude.txt`, and none is in a 0.1.12 set.
- **Nodes:** the standard `compute` nodes have 512 GB, of which 503 GB are usable, so 128 GB and 256 GB tasks fit on
  them. The 2 TB nodes are reached on the same partition by the memory request; there is no separate high-memory
  partition.

**Budget and ETA:**

| Item | Estimate |
|---|---|
| Pilot searches (8 threads) | median 5 s, 90th percentile 21 s, longest 91 s |
| Arrays B and C | about 50 node-hours; the cap is 150 |
| Wall time with 40 tasks running | about 6 h, plus the queue (compute was 199 of 204 nodes allocated on 2026-11-21) |
| IBMI analyses | about 10 h |
| Conversion | about 4 h |
| Download | depends on Globus |

The first 10 tasks of array B recalibrate the time limit before the rest are released.

## 11. What is not claimed

- No result transfers to search engines other than ProSE.
- The strata are metadata-based; a class under-filled after annotation is reported, not claimed.
- The development results (`v0111/README_0111.md`, the 218 release runs) are in-sample and are not evidence for this
  benchmark.
- The power figures (§5) assume that the held-out runs behave like the 19 pilot runs, resampled as whole runs. The
  pilot is too small for stratum-specific power, and the figures under higher attrition and with the pilot resampled
  show how far they can fall.

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

### 2026-11-21b, 0.1.12 (before sign-off; from the third review and the pilot)

1. **Scorer and extractor:**
   - checks are keyed by source (the auto arm, a scenario, a natural failure) and arm, so correct-scenario controls and
     both natural-failure arms coexist, and scenario runs outside the draw are validated against the manifest;
   - FDP keeps full precision;
   - the QC header comparator is the header's own evidence;
   - an extract → score integration test was added.
2. **P3 relative:** compared in signed integer counts, so −0.05 itself passes; `power.py` shares the function.
3. **Holm** always runs over all nine gates.
4. **The Cys rules:** 0/0 is uninformative, with at least 5 informative runs. The protease rule uses the runs with both
   a reference and an SDRF protease.
5. **`power.py`:** validates with `score.validate_checks`, derives the bounds exactly, resamples eligible runs only, and
   takes the loss rate from its data.
6. **Adjudication:** keeps the frozen scenario rows of unchanged organisms (`--snapshot`, frozen manifest).
   `refdb.py` writes the reference-database manifest of every drawn run.
7. **S1** uses prose-auto's label vocabulary. **S6** reports identification losses and check times.
8. **§0 Pilot** was added. It proposes a top-up to 300 drawn runs and keeps the P3 gates (X 0.75, −0.05). Not yet
   adopted: the author decides.

### 2026-11-21c, 0.1.12 (before sign-off; the author's decision on the pilot)

1. **Sample size:** the author adopted the top-up to 300 drawn runs and kept X = 0.75 and the margin at −0.05.
   `draw_topup.py` drew 113 runs (seed 20261122) and appended them to `heldout.tsv` (§2.1). Asp-N stays at 10 runs:
   no project was left.
2. **Power at 300** (superseded by 2026-11-22, which adds the sensitivity and attrition ranges): joint power 0.98 at the pilot effect, about 184
   eligible runs. The grid now counts every drawn run, the top-down share included. Before, its 187 row counted 187
   runs not tagged top-down.
3. **Check bounds at 300** by the frozen rule (§9): pre-search specificity 0.70 → 0.75; the other eight are unchanged (superseded by 2026-11-22).
4. **Manifests:**
   - `scenarios.tsv` is unchanged; it rebuilds byte for byte from the unchanged `wrongdb.tsv`;
   - `refdb_draw.tsv` gives the reference databases of all 300 runs from PRIDE's organism field;
   - its UniProt snapshot `uniprot_heldout/` is the `--snapshot` of the rebuild after adjudication;
   - on the PRIDE field, 282 runs have a database and 18 are no-proteome: 12 species without a reference proteome of
     their own, and 6 names UniProt's taxonomy does not resolve (synonyms such as *Neosartorya fumigata*, and
     non-species names). The annotation decides (§3.1);
   - `refdb.py` now puts its own folder on the import path, as `resolve_uniprot.py` does, so it runs under
     `python -I`.
5. **Execution plan** (§10.1): searches and checks on BinAC2, and every analysis that calls Jev on IBMI, with the node
   budget and the ETA.

### 2026-11-22, 0.1.12 (before sign-off; from the fourth adversarial review)

The 300-run draw, its seed and both P3 gates are unchanged.

1. **Scorer** (`score.py`, each with a test in `tests/test_prereg.py`):
   - a verified bottom-up run whose auto or reference arm has no result, or a missing one, is missing before the
     120-run minimum ("not executed"; `test_failed_runs_never_reach_the_minimum`);
   - the natural-failure cohort comes from the search counts, strictly below half the reference PSMs; a cohort arm
     without a check is missing, and checks outside the cohort are counted and ignored;
   - S6 losses come from the new `scenario_counts.tsv`, independent of the checks. The correct-scenario controls are
     reported per check against the manifest;
   - the QC rule counts every eligible run with header evidence, and a FASTag no-call as a disagreement;
   - S1 accepts the canonical short names (TMT0, 15N, 2H2O) and reports the family line;
   - `score.json` is strict JSON (+∞ as "inf");
   - P1–P3 and the verdict are one routine (`score.primary`), which `power.py` calls.
2. **Extractor** (`extract.py`):
   - it reads the CLI's `<stem>.check.json` (exactly one per folder);
   - it writes `scenario_counts.tsv`;
   - the QC value is FASTag QC's own call, or "no call".
3. **Manifests:**
   - `adjudicated.tsv` covers every drawn and wrong-database run; `refdb.py` (with `heldout.tsv`) and `scenarios.py`
     stop on a run without a row;
   - `refdb.py` lists a taxon once, and `refdb_draw.tsv` is regenerated: 5 rows lose a duplicate species, for example
     UP000002311 twice in PXD051453 and PXD086884;
   - §10 step 4 rebuilds and hashes `refdb.tsv` and the frozen reference panel;
   - the auto arm runs on that panel, and every task checks its FASTA files.
4. **The top-up's replay** (`draw_topup.py --replay`, 2026-11-22):
   - the draw reads PRIDE's live search results and file lists, which were not archived;
   - the replay against today's PRIDE gives the same first 19 runs. It diverges from the TMTpro SPS-MS3 stratum on,
     whose pool grew from 289 to 291 projects overnight, and 89 of the 113 projects agree;
   - **only the append is verifiable:** `topup.tsv` equals rows 188–300 of `heldout.tsv`, and the first 187 rows are
     the 187-run file. The replay's output (`topup_replay.tsv`, `draw_topup_replay.log`) is kept as the record.
5. **Power** (`power.py` rewritten; §0, §5, §9):
   - whole pilot runs are resampled, with their attrition, and the manifest's top-down and no-proteome runs are fixed;
   - there are no metadata pools, and stratum-specific power is unknown;
   - the scorer's own routine is used;
   - conditional, predictive, attrition and project-level sensitivity are reported, with P(no verdict);
   - the grid uses the exact 187- and 300-run manifests and states its 200 simulations per cell;
   - **the result:** conditional joint power 0.97; predictive 0.64; at the pilot's lower attrition bound,
     P(no verdict) 0.74–0.97; leaving one project out, 0.46–1.00. **Power under plausible attrition is well below
     0.8.** It is stated, and the draw is not changed.
6. **Check bounds** (§9):
   - specificity at the high-attrition count: 0.75 → 0.65 for both checks;
   - related_distant at the manifest's rank mix: pre 0.45 → 0.50, post 0.20 → 0.25;
   - the simulated Holm power is replaced by exact Bonferroni powers, a lower bound under any dependence;
   - the calibration rows are frozen in `power_inputs/`;
   - `power.py` records its command and refuses to overwrite its outputs without `--overwrite`.
7. **Run classes and conversion** (§2.2):
   - top-down runs run the auto arm only, and no-proteome and organism-not-verified runs run no arm;
   - a partial conversion is missing: the MS2 count must match the vendor index (§10 step 3).
8. **Execution plan** (§10.1):
   - the IBMI → BinAC2 hand-off of each run's analyses, with its order and release;
   - a hash check of every executable on both sites (FASTag on BinAC2 included);
   - one rerun for every failed step on either site;
   - the node RAM for the 256 GB rerun;
   - the workspace owner, extension (`proseauto` to 2026-12-21) and fallback;
   - the quota then and now (3.88 TB free on 2026-11-22), with a scratch fallback.
9. **README_0112.md:** the development-proxy power section is marked superseded.

### 2026-11-22b, 0.1.12 (before sign-off; the author's decision on the predictive power)

1. **Kept:** 300 drawn runs, X 0.75, margin −0.05, the seed and every gate. This was the author's choice among the
   predictive cells for 300, 400 and 500 drawn runs and the −0.05 and −0.10 margins.
2. **The conditional top-up** (§2.1):
   - it is checked after the blind annotation and conversion of all 300 runs, before any arm;
   - with fewer than 120 eligible runs, batches of 50 are drawn with seeds 20261201–20261204, fixed now;
   - it ends at ≥ 120, or with "no verdict" at the 500-run cap or when the pools are exhausted;
   - the count reads annotation and conversion only (`score.eligible_before_arms`), and `draw_topup.py --batch`
     accepts a batch only by that count;
   - tested in `test_conditional_topup_rule`.
3. **The stated assurance** (§9):
   - the predictive power and its decision table (`power_decision.txt`, `power_decision.py`), now in this folder and
     in `SHA256SUMS`;
   - the limit is the pilot's uncertainty about the P3 effect, not the size.

### 2026-11-22c, 0.1.12 (before sign-off of the freeze; the annotators and the harness)

1. **The annotators** (§3; the author's decision: "Two agents has to be enough"):
   - A is Claude and B is Codex, with Kimi as B's fallback; both are blind; one reconciliation round; no human
     adjudicator or review;
   - a field still disputed is not annotated: a disputed organism or class makes the run not verified, and any other
     disputed field leaves only its S-endpoint;
   - raw agreement and κ are reported per field;
   - the raw-data check is the measured indicators.
2. **The harness** (`heldout/`: `arms.py` templates, `harness.py` stages, `annotate.py`, `freeze.py`,
   `heldout.slurm`, `transfer.py`; `tests/test_heldout.py`), frozen at step 2. The facts it fixes:
   - **Bruker runs** are unpacked (zip CRCs checked) and read natively by prose-auto, ProSE and FASTag, as in the
     pilot; the mirror's mzpeak-convert recipe writes mzPeak, not mzML (§1 corrected);
   - **SCIEX `.wiff` runs** are fetched with their `.wiff.scan` companion, which the draw did not list, and converted
     to mzML by mzpeak-convert (§1 named FileConverter, which reads Thermo only);
   - **the conversion check** (§2.2) compares two counts of the mzML with the same two counts of the vendor file
     itself: the spectra at every MS level and the MS2 spectra. Both must agree, as §2.2 asked for MS2. The vendor
     counts come from the frozen FASTag reading the `.raw` (Thermo RawFileReader) or `.d` directly, in one full pass
     (`-census` input_spectra, `-qc_out` count_ms2), never from the converted file. A partial conversion loses spectra
     and fails;
     No frozen tool reads `.wiff` independently, so the `.wiff` runs (8 drawn, 1 wrong-database) have only the
     mzML-side count: the converter's validated index, with no vendor count to compare;
   - **verification of DIA runs** (§3.1) searches the subset of DIAspeXtract's pseudo-spectra. The reference windows
     of DIA and timsTOF runs are the instrument standard (the pilot's rule; §3.3 named no exception);
   - **the five non-tryptic proteases** of best_psm are Lys-C, Arg-C, Chymotrypsin, Glu-C and Asp-N; a candidate
     whose entrapment FDP exceeds 1 % scores 0;
   - **Jev:** `JEV_MODEL` must name a pinned model, never `jev-latest`, and it is recorded in `FREEZE.md`. The harness
     refuses Jev otherwise;
   - **`FREEZE.md`** also holds the corpus protocol (`corpus.py`, `entrap.py`), the 0.1.10 code and FASTag 1.10.0
     of `fast_0110`, both FileConverter builds, the pxd converter registry and the BinAC2 container image. The
     analyses (step 5a), the plans and the search databases are appended before the BinAC2 release;
   - **best_psm: unspecific cleavage only in the final search** (the author's decision, option 2: "Unspecific in
     final only"). Unspecific stays a protease candidate. When it wins the protease step, the label and Cys steps
     search with the specific protease the rule would have chosen without it (`arms.best_specific`), and unspecific is
     searched once more, in the final full search. Label and Cys are then chosen under a different protease than the
     final search uses; `best_psm.json` records it (`label_cys_under`). The reason: in the dry run (PXD000561, 5,000
     MS2, human + entrapment, 8 threads) one unspecific subset search peaked at 900 GB RAM and took about 4.5 h,
     and a run where unspecific wins needed about 12 of them;
   - **best_psm runs in its own array (B2) on high-memory nodes** (`--mem=1000gb`, 2 days per task): the unspecific
     subset candidate of every run, and the final search where it wins, need about 900 GB;
   - **`fast_0110`** gets the same `--organism` as auto: 0.1.10 fast mode refuses to run without it.
3. **Fixes from the short review of the harness** (Codex and Kimi, 2026-11-22):
   - **a stopped auto** (no_digest) is released: auto is scored stopped (0 PSMs, noncompliant in P3, §2.2/§6) and the
     comparators still run; the auto reports of top-down runs enter the scoring tree (S8);
   - **the comparators' spectra** come from the adjudicated acquisition (DIA: DIAspeXtract's pseudo-spectra, built
     by the harness and hashed with the plans), never from auto's decision; each prose-auto arm searches its own;
   - **Jev** must equal the model id, endpoint and vote count in `FREEZE.md` (`freeze.py init --jev-model ID
     --jev-url URL --jev-votes N`; ID is the concrete model jev-latest resolves to at the freeze);
   - **the ion-trap wide window** of §3.3 is 0.6 Da (the harness had the pilot's 1.0 Da); the DDA window search
     uses every manifest proteome of a mixture;
   - **reruns (§10.1)** are counted per unit (search, check, analysis, verification, conversion) across reruns of a
     stage or Slurm task; a unit due its rerun makes the task exit 1, a unit failing twice is missing and names its
     step; a check's exit code counts; conversions are tried at most twice;
   - **best_psm** selects a step only when every candidate has a result; a candidate out of attempts leaves best_psm
     missing. A search result is reused only for the same command;
   - **verification:** only final results count. A verification that fails twice makes the run missing, never "not
     verified"; `reference.tsv` is not written while one is open. A run whose spectra point to a panel organism goes
     to both annotators once more (`harness.py points-to`, `annotate.py points-to`, §3.1 "another organism"): agreed,
     that organism is the reference and the run is verified again; otherwise the organism is disputed;
   - **freeze coverage:** the converted inputs are hashed into `FREEZE.md` at step 3 and the pseudo-spectra with the
     plans, and every task checks them before use; a path listed twice with two digests stops; the converter
     registry, FASTag, the 0.1.10 tree and the prose_auto a child process imports are checked; the auto report's
     FASTA sources must be the run's manifest proteomes;
   - **outputs:** each Slurm task copies back exactly the paths it wrote (array A's verification and windows too);
   - **stage order:** `harness.py seal` ends step 4 (all conversions, verifications and windows final, no top-up
     due) and hashes the annotation; analyze and plan refuse before it, and a top-up batch after it or after any arm;
   - **score.py** requires `scenario_counts.tsv` when scenarios were built, with unique (run, kind) rows covering
     every built scenario; a failed conversion still gets an annotation sheet (from the texts).


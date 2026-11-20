# PREREG: prose-auto 0.1.11 with a known organism, held-out benchmark (frozen 2026-11-20)

Status: **draft for sign-off by the author.** Nothing has been run on a held-out run. No held-out raw file has been
downloaded, opened or converted; no SDRF, project record text or result file of a held-out project has been read by
a person. The run lists were drawn by `draw_sets.py` from PRIDE project-level fields and file names and sizes only.
Scored once; later changes are dated amendments at the end of this file.

## 1. Frozen code and tools

- **prose-auto:** the commit that adds this file (`git log -n1 --format=%H -- acceptance/0.1.11/PREREG.md`),
  version 0.1.11. Thresholds as committed: database check `DBCHECK_*` and post-search check `CHECK_*`
  (`prose_auto/external.py`, `prose_auto/check.py`), the Cys fallback (`cli._cys_fallback`).
- **FASTag** 1.15.3, Linux x64 release bundle (Thermo RawFileReader included).
- **ProSE:** the build of the ID corpus (`/ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr`, OpenMS 3.7.0-pre, PR
  299c9a6c99), with the corpus protocol: 2 missed cleavages (none for unspecific cleavage), at most 1 variable
  modification per peptide, isotope error down to −2, Percolator rescoring.
- **DIAspeXtract, FileConverter, Jev cache:** as in the 0.1.11 validation (`v0111/given/*/run.sh`). Jev answers are
  cached; a new request is asked 3 times (`PROSE_AUTO_JEV_VOTES`).
- **Conversion:** the ABI PRIDE mirror's default recipe (`pxd convert --to mzML`: OpenMS FileConverter vendor build
  for Thermo `.raw`, mzpeak-convert for Bruker `.d`); SCIEX `.wiff` through FileConverter. The recipe hash and the
  SHA-256 of every mzML are archived before any arm runs.

## 2. Sets

All drawn on 2026-11-20 with seed 20261120 by `draw_sets.py` (this folder), from PRIDE's v3 search API, excluding
the 2,691 accessions in `exclude.txt`: every PXD named in the prose-auto development tree (the 218-run release
set, the ID corpus, the mass-truth, mixture, non-tryptic, HLA, intact-protein and text-benchmark sets, every Jev
and text cache), FASTag's acceptance and held-out folders, and the ABI PRIDE mirror with its legacy archive.
Whole-PRIDE screens and catalogues (candidate screens, the bigbio clone) were not counted as use.

| File | Runs | Content |
|---|---|---|
| `heldout.tsv` | 187 | one vendor file per project, 50 MB–6 GB, by quota (below), 253 GB |
| `wrongdb.tsv` | 40 | 30 single-organism runs of `heldout.tsv` (seeded) and 10 runs of multi-species projects, 55 GB |
| `candidates.tsv` | 283 | every project the draw looked at, with the file it chose or `-` (86 had no vendor file of 50 MB–6 GB) |
| `draw.log` | | queries, pools and fills |
| `SHA256SUMS` | | the checksums of the files above, the freeze |

**Quotas** (`heldout.tsv` column `stratum`; a project counts for the first stratum that takes it, rarest first):

| Stratum | Runs | Selected by |
|---|---|---|
| protease Glu-C, Asp-N, chymotrypsin, Lys-C | 10 each | the enzyme named in the project texts, trypsin not named |
| TMTpro SPS-MS3 / TMTpro / TMT / iTRAQ / SILAC | 6 / 8 / 10 / 8 / 8 | quantification method or texts |
| phospho / glyco enrichment | 8 / 6 | texts |
| Astral / timsTOF DDA / timsTOF DIA / TOF (SCIEX, Waters) / ion-trap MS2 (LTQ) | 8 / 8 / 8 / 6 / 8 | instrument field, experiment type |
| Orbitrap DIA | 10 | instrument field, experiment type |
| plant / fungus / bacterium / other vertebrate / invertebrate / mouse or rat | 6 / 5 / 6 / 5 / 4 / 6 | single organism field |
| baseline: Orbitrap DDA label-free human | 13 | instrument, organism, experiment type |

The strata come from metadata and can be wrong. Every endpoint uses the **reference annotation** (§3), not the
stratum. A stratum that the reference shows to have fewer than 8 runs of its class (the 5 proteases, the label
classes) is reported as under-filled; it is not topped up after the freeze.

**Download:** after sign-off, `pxd get <PXD> --files <file>` for every listed file into the mirror. A file that
fails to download or convert is replaced by nothing; it is reported as missing.

## 3. Reference annotation (before any arm runs)

Made from the raw data and the authors' methods, by one annotator who has not seen any prose-auto, FASTag or arm
output of these runs. Recorded in `reference.tsv` (one row per run) and frozen with its SHA-256 before the first arm
runs.

1. **Organism and database:** the species the project states; the database is its UniProt reference proteome
   (frozen release, one per species; mixtures concatenated) plus cRAP. For the wrong-database set the related,
   unrelated and incomplete databases are chosen by the rules of `v0111/dbcal.py` (same family or order; another
   kingdom; a mixture's first component only) before any search.
2. **Settings ladder:** on a seeded 5,000-MS2 subset, ProSE searches each candidate of: protease (the stated one,
   trypsin, the five non-tryptic proteases, unspecific), label (none, TMT, TMTpro, iTRAQ4/8, SILAC by the reporter
   region and MS1 pairs), Cys alkylation (CAM, none, propionamide, NEM, methylthio, carboxymethyl) and the stated
   enrichment, one factor at a time from the stated settings; the candidate with the most PSMs at 1 % entrapment FDP
   wins when it beats the stated value by ≥ 5 %, else the stated value stays. Label plex from the reporter channels
   observed in MS2/MS3 (manual check of 200 spectra), SPS-MS3 from the scan headers.
3. **Windows:** a wide-window search (50 ppm / 0.05 Da high-resolution MS2, 0.6 Da ion trap) on the subset;
   precursor |median| + 4 robust SD, fragment from a Gaussian-plus-uniform fit, the same rules as the 0.1.10 mass
   truth.
4. **Disagreements** between the ladder and the methods are resolved by the annotator and recorded with the reason.

## 4. Arms

All arms search the same run with the same database (reference organism, cRAP, entrapment: a shuffled twin of every
target, seed 0; decoys reversed), the ProSE build and protocol of §1.

| Arm | Settings |
|---|---|
| **auto** | prose-auto 0.1.11 default: `--organism` = the reference organism(s), `--accession`, `--jev`, `--sdrf off`, `--contaminants cRAP`; its `prose_cli` |
| **reference** | §3: the annotated settings with the measured windows |
| **SDRF** | the authors' settings as stated: the submitted or curated SDRF if one exists, else the methods transcribed by the annotator; a field nobody states takes the ProSE default, never auto's |
| **defaults** | ProSE defaults plus the organism: Trypsin/P, Carbamidomethyl (C) fixed, Oxidation (M) variable, ProSE's default windows and charges |

## 5. Primary endpoints

Counts at 1 % FDR (Percolator q-values): PSMs, unique peptides, picked-protein groups. Ratios are paired per run;
runs where the denominator arm has fewer than 100 PSMs are left out of that ratio and counted. Confidence intervals:
percentile bootstrap over runs, 10,000 resamples, `numpy.random.default_rng(20261120)`.

| # | Endpoint | Criterion |
|---|---|---|
| P1 | median per-run PSM ratio auto / reference | **non-inferior:** lower 95 % bound ≥ 0.98 |
| P2 | median per-run PSM ratio auto / defaults | **superior:** lower 95 % bound > 1 |
| P3 | per-run entrapment FDP (combined estimate, Wen et al. 2025) of auto, PSM and peptide level | ≤ 1.5 % in ≥ 95 % of the runs with ≥ 100 PSMs |

Peptide and protein-group ratios (auto/reference, auto/defaults, auto/SDRF) are reported with the same intervals and
the same bounds; they do not gate. The study **passes** if P1, P2 and P3 all hold.

## 6. Secondary endpoints

Against the reference annotation; exact matches unless stated.

- **S1 label:** chemistry and plex exactly (TMT10 ≠ TMT11 ≠ TMTpro16 ≠ TMTpro18; SILAC with its channels);
  family level (isobaric / SILAC / none) as a second line.
- **S2 protease:** the exact enzyme (Lys-C ≠ Trypsin).
- **S3 Cys alkylation.**
- **S4 mass error:** auto's precursor and fragment tolerance within 1.5× of the wide-window measurement (§3.3).
- **S5 abstention** per setting: decisions with status `default`, `unknown` or `review`.
- **S6 wrong database and wrong settings, pre- against post-search** (on `wrongdb.tsv`; each run with the correct,
  related, unrelated and, for the 10 multi-species runs, incomplete database):
  - pre-search: the `database_check` warning of the auto report;
  - post-search: `prose-auto check` on the auto arm's ProSE result, verdict `database likely wrong or incomplete`;
  - wrong settings: `prose-auto check` on deliberately wrong settings of the 30 single-organism runs with the
    correct database (isobaric label dropped, or TMT6plex added to a label-free run; Trypsin replaced by Glu-C, any
    other protease by Trypsin; the rules of `v0111/postcal.py`), verdict `settings likely wrong`; and on every
    held-out run whose SDRF or defaults arm finds fewer than half the reference PSMs;
  - reported: sensitivity per wrong-database kind, specificity on the correct database, which check catches what,
    the identifications lost by the wrong database (the loss a warning prevents), and the seconds each check costs.
- **S7 run time** per input type (Orbitrap DDA, DIA, timsTOF, TOF, ion trap) on an idle node: prose-auto 0.1.11
  total and per step, against 0.1.10 fast mode on the same files.

## 7. Ablations (held-out set; identifications where the ProSE command changes)

- spectra only (no `--accession`, no description, no Jev);
- spectra plus texts (`--accession`, no Jev);
- spectra plus texts plus Jev (= auto);
- mini-search fallback: `--mini-search-cys never`, `fallback` (= auto), `always` (the 0.1.10 competition);
- `-alkylation_out` against the competition: the `never` and `always` arms on the runs whose Cys decision differs.

## 8. Comparators

- **Param-Medic** (tolerances; its label check) on the same mzML.
- **MSFragger 4.4.1** mass calibration (precursor and fragment) on the same mzML and database.
- **Engine defaults:** the defaults arm (ProSE); MSFragger's default parameters as a second line.

## 9. Decision rules (fixed now)

A detector becomes or stays the default only if it passes its own criterion; otherwise the decision abstains to the
texts or to a conservative default, and the report says so.

| Detector | Stays default if | Otherwise |
|---|---|---|
| FASTag `-alkylation_out` instead of the competition | `never` vs `always` on the runs whose Cys decision differs: summed PSM ratio ≥ 0.99 and no run below 0.95 | the competition becomes the default again (0.1.10 behaviour) |
| mini-search fallback | `fallback` vs `always`: as above | `always` |
| label | S1 exact ≥ the "always label-free" baseline and no isobaric call without reporter ions | no isobaric label without reporter ions |
| protease | S2 ≥ the stated-protease baseline | the stated protease decides |
| pre-search database check | sensitivity ≥ 0.70 over the wrong databases at specificity ≥ 0.85 | warning text only marked experimental |
| post-search check | the same bounds, separately for database and settings | verdicts reported as experimental |
| FASTag QC evidence | analyser class agrees with the header in ≥ 98 % of runs | removed from the evidence |

## 10. Scoring once

`score.py` (written before the first arm runs, its SHA-256 recorded in `FREEZE.md` with `reference.tsv` and the mzML
checksums) computes every endpoint above from the arm outputs. It is run once. A bug found afterwards is fixed in a
dated amendment that reports both results.

## 11. What is not claimed

- No result transfers to search engines other than ProSE.
- The strata are metadata-based; a class under-filled after annotation is reported, not claimed.
- The development results (`v0111/README_0111.md`, the 218 release runs) are in-sample and are not evidence for
  this benchmark.

## Amendments

(none)

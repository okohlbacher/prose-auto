# Changelog

## 0.1.10 (2026-11-05)

Mass errors from FASTag
- Tolerances and offsets come from the FASTag tag placements of the gap stage (`-recon_out`), in both modes:
  - the precursor error recomputed at full precision;
  - the b/y fragment errors of the placed peptides;
  - |offset| + 4 SD (precursor) and a Gauss-plus-uniform fit (fragment), as the wide-window reference.
- The mini-search PSM errors are no longer used.
- Mass truth (120 DDA cases, within 1.5×): precursor 59/59, high-res fragment 41/41, ion-trap fragment 19/20 (0.1.9:
  59/59, 40/41, 19/20).
- ProSE IDs on the runs whose command changed: PSMs +0.77 % (full, 140 runs) and +0.14 % (fast, 154 runs). On runs
  whose only change is the mass errors: +0.13 % and +0.15 %.
- Chance placements on a wrong given proteome spread over the window. They give no mass error: the precursor SD must
  stay below a quarter of the window, the fragment σ below 0.4 of it.

Mini-search
- **What it still does:** the organism competition (full mode) and the Cys alkylation competition (both modes). On
  timsTOF it runs the protease competition only when the tag ends give no call (2 of 13 tryptic timsTOF runs).
- **Its PSM mass errors:** no longer used.
- **The Cys competition stays:**
  - Without it, the SDRF Cys score rose by 2 runs.
  - But the 11 runs whose decision changed lost 9.5 % of their PSMs, the three that then matched the SDRF included.
- **Runtime:** most of the mini-search time is the Cys competition, so the changed steps take only 2–3 % less
  (timsTOF full mode: −19 %).

FASTag 1.10.0
- Requires FASTag >= 1.10.0. The species call uses `-species_evidence best`: on the 193 validation runs it gets 131
  species right, against 112 for `all` (1.9.0), and the organism decisions stay at 187/193.

Runs without a digest
- The recommendation depends on the class:
  - "use FLASHTnT (OpenMS)" for top-down, middle-down and native runs with MS2;
  - "use intact-mass deconvolution, e.g. FLASHDeconv (OpenMS)" for intact-mass runs and runs without MS2, because
    FLASHTnT needs MS2.

## 0.1.9 (2026-11-05)

Protease on timsTOF from the spectra
- The MS2 scans start at m/z 200–235, above y1 of K and R. FASTag's tag-profile call fails there, and an R end reads
  as R+10.
- Two signals now decide, with Jev and the texts as before (Jev still picks the member of the K/R family):
  - **The tag ends:** K/R ≥ 0.30 means Trypsin. K/R ≤ 0.10 means Glu-C (E ≥ 0.5), Asp-N (N-terminal D ≥ 0.3) or
    chymotrypsin (FYWLM ≥ 0.6).
  - **The mini-search protease competition**, on timsTOF only.
- Validation: 17/17 timsTOF runs, 13 tryptic and 4 non-tryptic. The tags call 14 and the competition 15. No HLA
  timsTOF run was available.
- The heavy-Arg share of the tag ends (0.29–0.91 on label-free timsTOF) reaches Jev with a timsTOF caveat. No label
  decision depended on it.

Mixture unions: 3 organisms at most (was 6)
- The HYE runs and PXD012636 keep their unions. PXD017035 (8 species) drops its fourth member, E. coli: PSMs −7.4 %,
  protein groups −12.1 %.
- A cap of 2 would have cost the five affected runs 8.0 % of their PSMs and 14.4 % of their protein groups.
- `--fast` with every `--organism` still searches the full list.

Jev 3-vote
- Jev takes no temperature or seed. A new request is now asked 3 times (`PROSE_AUTO_JEV_VOTES`). The majority choice
  wins per question, with the probabilities averaged; if all three answers differ, the vote abstains. The vote and
  the single answers are cached under a new key.
- On the 730 requests of the validation, two single answers differ in 1.15 % of the choices, and a single answer
  differs from the vote in 0.51 %.
- SDRF accuracy is unchanged. The cost is 3× the requests (about 0.03 USD per run) and +0.36 s per run.

Texts and downloads
- `repository.fetch` retries 5xx three times. `PROSE_AUTO_TEXT_CACHE` keeps fetched texts, so a PRIDE or Europe PMC
  outage does not change a rerun. Europe PMC answered 503 for hours on 2026-11-05.
- Proteome downloads wait as long as UniProt's Retry-After asks on HTTP 429.

## 0.1.8 (2026-11-05)

Fast mode
- `--fast` with `--organism` (repeatable) and/or `--fasta`: the given organisms are taken as true. Skipped:
  - the ProSE signature-panel organism search;
  - the panel mini-search organism competition;
  - the FASTag species call;
  - tier 2;
  - mixture detection and union building;
  - Jev's organism question.

  Every other detector runs against the given proteome, or against the union of several. The report shows the
  organism as "given (fast mode)" and leaves it out of the agreement score. See README "Fast and full mode" and
  docs/decision_flow.md.
- Database: the `--panel` file, else the UniProt reference proteome by name. Without one, the nearest relative's
  proteome within 3 parent taxa is used, named in the report and on the console. If nothing is found, fast mode
  stops with an error and never infers the organism.
- `--organism` without `--fast` is a user statement (confidence 0.9) for full mode.
- Fast against full on 217 validation runs:
  - the organism-dependent steps take 75 % less time;
  - the whole run takes 66 % less on DDA, 56 % less on DIA, and 13 % less on timsTOF, where reading the .d files
    dominates;
  - with the same database, PSMs change by −0.12 % and protein groups by −0.09 %;
  - on non-organism parameters, fast mode is right exactly where full mode is (curated SDRF).
  - The remaining decision differences come from the database: mass accuracy from the PSMs of the given proteome,
    not of full mode's search winner (median offset change 0.05 ppm).

Proteome cache
- UniProt's name search lists viruses of a host before the host itself. For 12 host species, 0.1.0–0.1.7 searched a
  virus proteome: Tomato torrado virus for tomato, Triticum mosaic virus for wheat, and papillomaviruses and other
  viruses filed under the host name.
  - A species query now takes the species itself or its strains. A renamed genus with the same epithet also counts
    (Lactobacillus → Lactiplantibacillus plantarum).
  - It never takes a virus or phage named after the species.
  - Effect on 0.1.6/0.1.7: the decided organism's database was a virus in 6 validation runs, which affected their
    mass accuracy, Cys competition and mass gaps. Corrected values are in the v018 report.
- `<cache>/resolved.tsv` maps each query and taxon id to its cached file, so a cached proteome needs no request.
- Downloads:
  - a cut chunked stream (IncompleteRead) is retried;
  - when the REST stream refuses (HTTP 429/5xx), the proteome comes from UniProt's FTP mirror, checked against the
    taxon id in its headers.

## 0.1.7 (2026-11-05)

FASTag 1.9.0
- Requires FASTag >= 1.9.0. 1.9.0 keeps reversed (b-ion) placements and changes its default tags; options and
  report formats are unchanged.
- TAG_MIN 100 -> 175 C-terminal tags. 1.9.0 reads a median 1.75x the C-terminal tags of 1.8.0, so 175 keeps
  1.8.0's operating point. At 100, PXD031090 (TiO2) called HLA-II against trypsin.
- Against 0.1.6 (218 runs, Jev cache on), 63 runs change a decision:
  - 58 modification lists from the new mass-gap peaks (4 of them together with a new Jev answer);
  - 3 organism labels of the same species;
  - one tolerance, one Cys and one enzyme decision (the last two from new Jev answers).
  No organism or union decision changes.
- ProSE IDs on the 71 runs whose command changed: PSMs -0.23 %, protein groups -0.08 %.
- FASTag time -40 to -43 %, whole run -16 to -29 %.

## 0.1.6 (2026-11-04)

Complete mixture unions
- Every organism that the user, the record, Jev (>= 0.3) or the file name names is mini-searched, together with
  the panel mini-search's own winner. The file name counts as HYE tokens, amount lists such as `05Y_50H_45E`,
  or words such as `Ecoli`. An organism outside the panel is fetched from UniProt by name, or by its NCBI taxon id
  from the PRIDE record (walking up to 3 parents), and cached under `<panel>/extra`.
- Each candidate's spectra are tested against the union so far, so a third or later species can join. The
  candidate with the most unexplained spectra joins first. To join it needs at least that many spectra the union
  does not explain:
  - 1 when the file name names it;
  - max(10, 0.5 %) when an indicator or the winner names it;
  - max(40, 2 %) when only the texts, or a record beyond its first three organisms, name it.
- The cap rises from 4 to 6 organisms. On PXD017035 (8 species), PSMs level off at 5–6 members while the combined
  PSM FDP stays at 0.85–1.09 %.
- Mixset:
  - PXD017035 goes from 2 species and 639 PSMs to 4 species and 28,391 PSMs; the reference arm has 24,446.
  - PXD039576 adds E. coli.
  - Over the 19 union runs: 406,802 PSMs, 104 % of the true-organism search.
  - Pure files are unchanged.
- Benchmark: PXD066672 adds its 5 % yeast. No other organism decision changed.

Jev cache
- Jev accepts no temperature or seed, and identical requests differ by up to 0.05 in probability. Answers are now
  cached per study (`--jev-cache DIR` or `PROSE_AUTO_JEV_CACHE`, `DIR/<accession>/<sha256>.json`).
- The cache key is the request: the texts, the questions, and the measured profile with its numbers rounded to
  2 significant digits. The API key is never part of it.
- A second pass over valset and borg_ho (180 runs) asked Jev nothing and changed no decision.

## 0.1.5 (2026-11-03)

Mixture union search
- When organism_mix is suspected, ProSE, the alkylation probe and the FASTag mass-gap pass search the union of the
  suspected proteomes (`<panel>/extra/union/A+B.fasta`). An organism joins if an indicator of strength >= 0.3 names
  it and the panel mini-search gives it at least max(10, 0.5 % of the spectra) spectra the winner does not explain;
  winner first, then by those spectra, at most 4. Reagent genera join only when the user or the record names them.
- `.auto.json` carries `organisms` and `search_database`; `.auto.md` names the union; the SDRF organism column
  lists them separated by ";".
- Mixset (38 runs): 19/19 mixtures get a union, no pure file does. On the union runs PSMs +24 %, peptides +27 %,
  protein groups +47 % over 0.1.3 (97 %, 94 %, 97 % of the true-organism search), combined PSM FDP <= 1.5 %.
  Benchmark (180): 3 unions (PXD012636, PXD053406, PXD066672), PSMs +32 %.

Charge range
- The charge range counts only precursors up to 9,000 Da (ProSE's peptide:max_mass) and is clamped to 1..20,
  ProSE's limit. PXD062749 (heavy non-peptide ions at 8+ to 60+) went from 2..56 to 2..20; five other runs narrow.

## 0.1.4 (2026-11-03)

Runs without a digest
- Top-down, middle-down, native and intact-mass runs are detected right after reading the spectra
  (`topdown.py`). Such a run gets no ProSE parameters, no Jev questions and no mini-search. prose-auto writes the
  verdict, the reason with its numbers and "use FLASHTnT (OpenMS)" to `.auto.json` and `.auto.md`, and exits with
  status 3.
- Signals: charged MS2 precursors (median mass and charge, share at 8+); without precursor charges FLASHDeconv
  (`--flashdeconv`) on 40 MS1 and 40 MS2 spectra, and a charge-state series among the MS1 peaks scored against
  chance. Native: a charge series above m/z 2000, or precursors of >= 10 kDa at most at 0.078 sqrt(M).
- On 123 public runs (13 studies) against 268 bottom-up runs: intact mass 13/13, top-down 59/60, native 30/38,
  middle-down 8/12 aborted; no bottom-up run aborted. Two benchmark runs are top-down and are now aborted
  (PXD033100 intact EF-Tu, PXD024405 protamines). Thresholds set on these runs.

FASTag 1.8.0
- Requires FASTag >= 1.8.0. The termini and mass-gap passes read FASTag's own reports (`-digest_out`,
  `-label_out`, `-massgap_out`) with the 0.1.3 model table (`data/prose_auto_013_digest_models.tsv`); the tag-TSV
  parsing, the model scoring, the mass-gap histogram code, `delta_combos.py` and the MS2 scan-start gate (now
  FASTag's `ms2_first_peak` rule) are gone. The policy (confidences, TMTpro margin override, label, SILAC,
  alkylation and modification rules, HLA mapping, PASEF gate) stays in prose-auto.
- On 213 rerun runs every FASTag value equals 0.1.3's (terminal tags within 17, margins within 0.006); eight
  decisions moved because Jev answered differently. The two FASTag passes take 0.70x the time (59.5 s -> 41.5 s
  on 10 runs).

## 0.1.3 (2026-11-02)

Proteases and HLA ligands
- Protease from the published matrices: the FASTag tag termini (C-terminal residue and, from 100 tags, the
  N-terminal one; window +-0.008 Da) are scored by multinomial log-likelihood against MEROPS P1/P1' frequencies
  (trypsin, Lys-C, Lys-N, chymotrypsin, Glu-C, proteasome), PeptideCutter rules for Arg-C and Asp-N, the MHC
  Motif Atlas ligand ends (HLA-I, HLA-II) and the human proteome background ("unspecific") (5 % non-specific
  cleavage mixed into every model; `data/protease_models.tsv`). The margin per tag is the confidence; a call needs
  100 C-terminal tags and a margin >= 0.15. HLA-I/proteasome give class I ligands, HLA-II/unspecific class II.
- No tag call on timsTOF (ddaPASEF, diaPASEF: C-terminal R is nearly invisible), when the MS2 scans start above
  the y1 ion of K (147.11), or between trypsin and Lys-C/Arg-C when isobaric-labelled K dominates (TMTpro).
- Enzyme decision: user or repository, else the stronger of the tag call and Jev "+data" at p >= 0.8, else the
  keywords, else Trypsin. Inside the trypsin family (Trypsin, Trypsin/P, Lys-C, Arg-C) a Jev answer picks the
  member. The mini-search protease competition is removed.
- Keywords: Lys-N, lysyl and glutamyl endopeptidase, V8 protease; "unspecific" needs immunopeptidomics or HLA
  ligands, a bare "HLA" no longer counts; "no enzyme" in a control phrasing and texts naming three or more
  proteases give no call; a protease in the file name (`AspN`, `GluC`, `LysC`, `LysN`, `ArgC`, `chymo`) is evidence.
- SILAC needs >= 50 tags on both K and R; a low K/R share alone no longer votes for HLA ligands.
- Non-tryptic Orbitrap DDA (16 runs): enzyme 9 -> 15/16, label 14 -> 15/16, Cys alkylation 7 -> 8/10.
  DIA/timsTOF (11 non-tryptic runs, 12 tryptic controls): 5 -> 8/11, controls 12/12. HLA ligands 14/14.

Organisms
- Organism names resolve to NCBI taxon ids (`taxa.py`, `data/organism_taxa.tsv`, 165 taxa from the UniProt
  taxonomy); "yeast" is *S. cerevisiae* (4932), not *Candida albicans*.
- Tier 2 also searches every stated organism (user, repository, text from 0.5). A record's pick that was searched
  and explains under max(10, 1 %) of the spectra yields to another named species the search finds. A database
  under 10 % of the largest searched one's peptides tests nothing; an organism searched only for being stated
  does not keep a record's pick on the data shortlist.
- `organism_mix`: the "A + B" search indicator needs >= 40 unique spectra and 0.12 x the winner's PSMs, without
  PSMs of trypsin, keratins and serum albumin. A mixture only the records and texts describe is reported as
  study-level (`scope`), and not suspected when every stated second species was searched and is absent.
- Label agreement: "no SILAC" and "no isobaric label" agree unless they negate the decided label class.
- Mixture set (38 runs): sensitivity 19/19, specificity 15 -> 18/19; single-species files right 18 -> 19/19;
  suspected mixtures among the 180 validation runs 19 -> 9. Agreement at 0.2 on all data: precision
  0.16 -> 0.25, recall 0.71 -> 0.65.

Search windows
- ProSE windows per analyzer, as multiples of the measured error with ppm floors: Orbitrap precursor 4x, fragment
  2x (10/20 ppm); Astral 8x, 2x (10/10 ppm); TOF 1x, 1x (10/15 ppm); ion trap unchanged.

Fixes
- Target-decoy FDR of the mini-search: tied scores pass or fail together; the 1 % cut no longer falls inside a group
  of equal scores.
- A tier-2 organism winner also redoes the FASTag mass-gap evidence (Cys alkylation, modification combinations)
  against its database; before, the gaps against the losing database were kept.
- The FASTag species call runs once per run: after a failed ProSE organism search the panel stages reuse it.
- `--threads` also caps the OpenMP, OpenBLAS and MKL thread pools (unless set in the environment).

180 validation and held-out runs: enzyme 174/175 and every other setting as in 0.1.2; organism 174 -> 177/180.

## 0.1.2 (2026-11-02)

- Jev reads last: after the mini-search, the ProSE organism search, FASTag and the mass gaps, with the texts and
  every indicator (counts and reference ranges, no verdicts, no confidences). Its second request judges each final
  setting against the texts and the indicators (supports / texts contradict / indicators contradict / not
  mentioned); tolerances and modifications against the texts only.
- New output `agreement`: per setting, the confidence-weighted share of independent detectors that disagree
  (inconsistency, level, dissenting detectors, Jev's verdict), and an overall score. An inconsistency >= 0.2
  flags 53 of the 73 wrong decisions and 16 % of the correct ones.
- New output `organism_mix`: a suspected mixture of organisms from at least two kinds of indicator (a second
  database with spectra of its own, several genera stated, Jev's p(mixed species)).
- No single indicator decides an enrichment: a strong diagnostic-ion signal or two independent indicators are
  needed; on ion-trap MS2 the neutral-loss test no longer counts as tested. Enrichment 160 -> 165/180.
- HLA ligand call from several indicators (tag termini where valid, hydrophobic C-termini, precursor masses,
  singly charged share); tag termini alone no longer claim an unspecific digest.
- Lock-mass share counts scans with a hit, not hits.

## 0.1.1 (2026-11-02)

- pyopenms 3.5 (the last release for Python 3.10): idXML is read into a `PeptideIdentificationList`; a plain list
  failed there and broke the ProSE organism search.
- Tests: the shell-script stand-ins of DIAspeXtract and FileConverter are skipped on Windows.

## 0.1.0 (2026-11-02)

First release of the prototype.

- `auto` organism mode: the spectra check a stated organism (panel mini-search, FASTag, tier 2). Without a
  stated organism, one ProSE signature-panel search decides. 174/180 correct on the validation and held-out studies.
- HLA ligands: the unspecific-cleavage claim of the ProSE ladder needs a fully-tryptic PSM share below 0.5
  (tryptic runs 77-98 %, HLA ligand files 4-5 %).
- Persistent peptide indexes, a random 2,000-spectrum mini-search sample and query-time Cys and protease
  variants: a median run takes about 1-2 min.
- Readers (LLM, Jev) see all metadata and the measured data features. `--jev-priors` adds parameter
  distributions to the Jev prompt (no measurable effect; off by default).
- Repository fetches retry once on DNS, timeout and 5xx failures.

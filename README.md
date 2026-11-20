# prose-auto: an Auto Mode prototype for ProSE

[![tests](https://github.com/okohlbacher/prose-auto/actions/workflows/tests.yml/badge.svg)](https://github.com/okohlbacher/prose-auto/actions/workflows/tests.yml)


`prose-auto` proposes ProSE search settings from an mzML file and the organism (or the FASTA) the user
gives. It uses the spectra and the file metadata, and optionally a free-text description, a methods
section, and the PRIDE/MassIVE/SDRF records of a dataset. It checks that the database fits the run, before
the search and, with `prose-auto check`, after it. Every proposed setting comes with a confidence, its
source and the evidence behind it, and a setting is flagged for review
when the sources disagree.

This is a temporary exploration tool. It shows which detectors work and what
the GUI needs. It is not meant to be the final implementation: the detectors
that work belong in ProSE (C++) or in a prose-gui sidecar, as described in the
Auto Mode report.

## Install and run

```bash
pip install pyopenms numpy pyyaml          # or: pip install -e .
python -m prose_auto.cli panel proteomes/  # once: UniProt Swiss-Prot of 12 model organisms
python -m prose_auto.cli analyze run.mzML --organism "Homo sapiens" --panel proteomes/ \
    --fastag FASTag --contaminants crap.fasta \
    --description "HeLa, TMTpro 16plex, IAA, Exploris 480" \
    --accession PXD012345 --llm claude -o out/
python -m prose_auto.cli check run.mzML run.idXML --fasta db.fasta --fastag FASTag --report out/run.auto.json
```

Since 0.1.11 the organism is an input: `--organism` (repeatable) or `--fasta` (repeatable) is required. Without
them prose-auto stops; with `--accession` it names the organisms the PRIDE record states. `--infer-organism`
infers the organism instead (the 0.1.10 full mode).

The outputs, written per input to `out/`:

| File | Content |
|---|---|
| `<run>.auto.json` | everything: decisions, evidence, ProSE parameters, SDRF columns, LLM review (the GUI contract) |
| `<run>.auto.md` | readable report: settings table, conflicts, ProSE command line, evidence |
| `<run>.prose.ini` | ProSE INI fragment (`ProSE:1:Search:...`) |
| `<run>.sdrf_technical.tsv` | SDRF technical columns (instrument, label, tolerances, modifications) |

`--json` prints the report JSON to stdout and the log to stderr, so prose-gui
can run it as a sidecar.

## External tools (optional)

| Flag | Tool | What it adds |
|---|---|---|
| `--fastag <exe>` | [FASTag](https://github.com/okohlbacher/FASTag) >= 1.15.2 | on the given FASTA: the protease from the tag ends, the mass errors of the tag placements, the Cys alkylation (mass gaps, `-alkylation_out` on the placements, the database-free `-alkylation_tags_out`), the run QC (`-qc_out`: analyser class, charges, reporter kits), and the database check; with `--infer-organism` also the species call |
| `--contaminants <fasta>` | | cRAP: part of the database-check panel (and of the ProSE organism search) |
| `--fastag-subsample <f>` | FASTag | analyse only a fraction of the MS2 spectra |
| `--diaspextract <exe>` | [DIAspeXtract](https://github.com/okohlbacher/diaspextract) | converts a diaPASEF `.d` folder (required for `.d` input) or a DIA mzML with ion mobility into pseudo-MS/MS spectra (Orbitrap DIA without ion mobility skips it and stays on the DIA path; DDA-PASEF `.d` is not supported and stops with an error), so the mini-search and FASTag also run on DIA data |
| `--threads <n>` | both | threads for the external tools |

FASTag needs high-resolution MS2 to be useful: on the ion-trap BSA example it
places only a handful of tags. Without FASTag the mass errors come from the database-free estimates only.

## Runs without a digest

Top-down, middle-down, native and intact-mass runs are recognised right after reading the spectra
(`topdown.py`): from the MS2 precursor charges and masses, else FLASHDeconv on a few MS1 and MS2 spectra
(`--flashdeconv`) and charge-state series in MS1. prose-auto then sets no ProSE parameters, writes the verdict with
its reason and a recommendation, and exits with status 3. Top-down, middle-down and native runs with MS2 get "use
FLASHTnT (OpenMS)". Intact-mass runs and runs without MS2 get "use intact-mass deconvolution, e.g. FLASHDeconv
(OpenMS)", because FLASHTnT needs MS2.

## Pipeline

| Step | Module | What it decides | Basis |
|---|---|---|---|
| 1 metadata | `detect_instrument.py` | instrument model and class, MS1/MS2 analyzer (filter string `FTMS`/`ITMS`, mzML CV), activation, NCE, ion mobility | mzML CV terms, Thermo filter strings |
| 2 acquisition | `detect_instrument.py` | DDA / DIA / diaPASEF / PRM, DIA window scheme, charge range | isolation-window width, number of distinct window centres, repeats per cycle, charge annotation |
| 3 mass precision | `detect_tolerance.py` | fragment and precursor precision without a database | Param-Medic pairs (May et al. 2017), extended to DIA window repeats and MS1 scan-to-scan; Gaussian + uniform EM |
| 3 mass offset | `detect_tolerance.py` | systematic offset without a database | lock-mass polysiloxane ions (MS1), immonium and y1 ions (HCD MS2) |
| 4 chemistry | `detect_chemistry.py` | TMT/TMTpro/iTRAQ plex, SILAC, phospho/glyco/acetyl enrichment, Cys alkylation (weak) | reporter and diagnostic ions against control m/z values (Param-Medic 2019 idea); MS1 heavy/light pairs with an intensity check |
| 5 mini-search | `tagsearch.py`, `hypotheses.py` | organism (full mode), Cys alkylation; protease on timsTOF without a tag-rule call | binomial b/y scoring with target/reversed-decoy FDR on the top-TIC subset; competing hypotheses compared by PSMs at 1% FDR |
| 5 mass errors | `hypotheses.placement_accuracy` (FASTag `-recon_out`) | calibrated offsets and tolerances | the precursor and b/y fragment errors of the tag placements of the gap stage: \|offset\| + 4 SD (precursor), Gauss-plus-uniform fit (fragment) |
| 5 protease | `external.py` (FASTag `-digest_out`) | protease, HLA ligand class | FASTag scores the tag termini against MEROPS, PeptideCutter, MHC Motif Atlas and proteome models (`data/prose_auto_013_digest_models.tsv`); on timsTOF the K/R and residue shares of the tag ends (`tims_protease`) and the mini-search protease competition instead |
| 5 records and text | `repository.py`, `text_extract.py` | what the submitters and the user say | PRIDE API v3, PROXI, bigbio/sdrf-annotated-datasets, Europe PMC methods; keyword rules with verbatim quotes |
| 6 LLM | `llm.py` | the same from free text, literature priors, a final review | any CLI or OpenAI-compatible endpoint; quotes must occur in the text |
| 6b Jev (last) | `jev.py`, `schema.py` | every head from the texts and all indicators; afterwards a consistency verdict per setting | typed choice questions over the texts plus the measured and search profiles (`--jev`) |
| 7 reconcile | `reconcile.py` | one decision per parameter, confidence, conflicts, status; indicator agreement; suspected organism mix | measured beats stated for instrument-determined settings; a clear spectral result beats text for sample prep; an enrichment needs a strong signal or two indicators |
| (opt) ProSE | `prose_runner.py` | confirmation with ProSE's own calibration pass and an open search | `-summary_out` YAML (`precursor_error`, `fragment_error`), `-out_mod_analysis_dir` |

## LLM backends

`--llm` takes `none` (the default; the rule-based extractor still runs),
`claude` (`claude -p`), `codex` (`codex exec -`), `cli:<any command reading
stdin>`, or `openai:<model>` (`OPENAI_BASE_URL`, `OPENAI_API_KEY`). It
answers three prompts: extract settings from text, typical tolerances for an
instrument/acquisition (`--llm-priors`), and a review of the final settings
(`--llm-review`).

## Reading the texts: one schema, many heads

`schema.py` defines the ProSE settings that texts can state as *heads*
(organism, sample type, protease, missed cleavages, Cys alkylation, label,
enrichment, acquisition, fragmentation, MS2 analyzer, instrument, precursor
and fragment tolerance), each with a closed list of values, plus yes/no flags
for variable modifications and for search types ProSE does not cover
(open, cross-link, glycan, top-down). All texts of a study (the user's
description, SDRF, PRIDE title/description/protocols, methods) go into one
labelled state and are read in **one request**:

* `--llm`: the generative LLM gets the heads, values and consistency rules in
  one prompt, answers every head with a verbatim quote (checked), and lists
  where the texts contradict each other or a rule.
* `--jev`: [Jev](https://docs.typesafe.ai) answers every head as a typed
  `choice` question (with "not stated") and every flag as a `noul` question.
  The head answers are decoded jointly: the most probable combination that
  breaks no hard rule wins, so a near-tie in one head is settled by the
  others. Jev runs after all searches and sees the texts plus every indicator
  as profiles with reference ranges. After reconciling, it judges each final
  setting against the texts and the indicators (supports / texts contradict /
  indicators contradict / not mentioned). Two requests per run, each asked 3 times, about $0.03. A TypeSafe console key
  (`TYPESAFE_API_KEY` or `JEV_API_KEY`) uses `api.typesafe.ai/v1/systemone`; a
  `jv_live_` key uses the jevtypesafeai.com proxy. About $0.001 per request.
  Jev takes no temperature or seed, and two identical requests differ by up to 0.05 in probability, and in
  1.2 % of the answers by their choice. So a new request is asked 3 times (`PROSE_AUTO_JEV_VOTES`). Per question
  the majority choice wins, with the probabilities averaged. If all three answers differ, the vote abstains
  (not stated). The vote is cached (`--jev-cache DIR` or `PROSE_AUTO_JEV_CACHE`;
  `DIR/<accession>/<sha256>.json`, with the single answers). The key is the request: the texts, the questions
  and the measured profile with its numbers rounded to 2 significant digits (the reference ranges stay as
  written), so a rerun asks nothing and decides the same. Fetched texts can be kept too
  (`PROSE_AUTO_TEXT_CACHE`), so an outage of PRIDE or Europe PMC does not change a rerun.

**Consistency rules** (`schema.RULES`) say which values exclude each other
(*hard*: TMT/iTRAQ with ion-trap CID MS2, PASEF on a non-timsTOF, ETD on an
instrument without ETD, trypsin on an immunopeptidome, ppm fragment windows
on ion-trap MS2, Orbitrap vs TOF analyzer), which are unusual (*soft*: TMT
with DIA, unspecific search of a digest, K-GG without trypsin, 0.3 Da fragment
windows on high-res MS2, metabolic 15N/heavy-water labels, mixed species,
out-of-scope search types) and which imply others (*imply*: phospho
enrichment adds Phospho (STY), PASEF implies TOF and ion mobility). They are
applied to the text answers, to the final settings (`consistency` in the
report) and to the tuning loop's proposals.

## Tuning loop (`python -m prose_auto.tune`)

Starts from a report's ProSE settings and improves them with ProSE on a
top-TIC subset: each round reads ProSE's diagnostics (PSMs, calibration error
of precursor and fragment, open-search mass shifts), builds a shortlist of
changes (widen/narrow tolerances, add a modification seen as a mass shift,
drop surplus modifications, Trypsin/P, more missed cleavages, calibration),
drops any that break a hard rule, lets Jev order them (`--jev`; rule order
otherwise), searches the top few and keeps a change that gains at least 1%.
`--entrap entrap.py` scores by PSMs at 1% entrapment FDP instead of raw PSMs.
Every step goes to `<run>.tune.trace.jsonl` for later calibration of the
chooser.

```bash
python -m prose_auto.tune out/run.auto.json --prose ProSE --database td_entrap.fasta \
    --entrap entrap.py --jev --threads 16
```

## Tests

```bash
git clone --depth 1 --filter=blob:none --sparse https://github.com/OpenMS/OpenMS openms
git -C openms sparse-checkout set share/OpenMS/examples/BSA share/OpenMS/examples/TOPPAS/data src/tests/topp
python -m pytest tests -q                 # or set PROSE_AUTO_OPENMS_ROOT to an OpenMS checkout
```

The tests use OpenMS test data: BSA1 (LTQ Orbitrap XL, CID/ion trap), the
DIAuditor and OpenSWATH DIA inputs, and the IsobaricAnalyzer iTRAQ input. A
fake LLM covers the quote guard.

## Known organism (the default) and `--infer-organism`

**The default** takes the organism as known:

```bash
python -m prose_auto.cli analyze run.mzML --organism "Homo sapiens" --panel proteomes/ --fastag FASTag -o out/
python -m prose_auto.cli analyze mix.mzML --organism "Homo sapiens" --organism "Escherichia coli" --panel proteomes/
python -m prose_auto.cli analyze run.mzML --fasta my_db.fasta --fasta crap.fasta
```

- **What is given:** `--organism` (repeatable) and/or `--fasta` (repeatable). Without either, prose-auto stops;
  with `--accession` the message names the organisms of the PRIDE record (`--organism "Glycine max" ...`).
- **The search database:**
  - each `--organism` resolves to its panel proteome, else to its UniProt reference proteome (cached under
    `<panel>/extra`);
  - a species without a UniProt reference proteome is searched against its nearest relative's (the first reference
    proteome under its NCBI taxon or up to 3 parent taxa). The organism stays as given; the report names the
    database used. Nothing within 3 parents, or a failed download, stops the run;
  - every `--fasta` is added as it is; several databases are concatenated once into `<panel>/extra/union/`.
- **Skipped:** everything that only serves the organism: the ProSE organism search, the organism competition of the
  mini-search, the FASTag species call, tier 2, the mixture test and union, Jev's organism question. Every other
  organism statement is dropped.
- **One FASTag pass on the given database** (k=4 tags, the isobaric label the reporter ions show as fixed
  modification) gives the mass errors of the placements, the mass gaps (Cys, variable modifications, enrichment),
  FASTag's Cys calls and the run QC.
- **The mini-search is a fallback** (`--mini-search-cys fallback`):
  - the Cys competition runs only when FASTag's Cys call (`-alkylation_out` on the placements, else the
    database-free `-alkylation_tags_out`) and the texts disagree, or when both abstain; Jev's answer counts as a
    text and is checked again after Jev reads;
  - on timsTOF the protease competition runs when the tag ends give no call (fewer than 175);
  - `always` restores the 0.1.10 competition, `never` uses FASTag's call alone (ablations).
- **Report:** `mode: "given"`; the organism decision has status `given`; `database_check`, `fastag_qc` and
  `mini_search` (whether it ran, and the Cys competition) are reported.

**`--infer-organism`** is the 0.1.10 full mode: the organism and the database come from the spectra and the texts
(ProSE or panel organism search, FASTag species, tier 2, mixture union); `--organism` and `--fasta` are then
statements and candidates. It keeps the 0.1.10 FASTag pass (no label fixed, no FASTag Cys call, no QC).

### Database check (before the search)

A wrong or incomplete FASTA is the main failure when the organism is an input. A tolerance estimate does not catch
it: on PXD039269, given as aphid (truth soybean), the wrong proteome still gave 8 ppm.

- **Signal:** the spectra with a confident FASTag tag (E ≤ 1e-4) of the gap stage; the share whose tag places on
  the given database with an agreeing flank mass (one flank: the precursor-dependent flank breaks on isotope
  errors, unfixed Cys and N-terminal modifications), against a fixed panel placed in one more FASTag pass on the
  same spectra: FASTag's shipped human, mouse, yeast and E. coli proteomes (UniProt reviewed, 17 MB) and the
  `--contaminants`.
- **Warning** (`--check-organism warn`, the default) when fewer than 5 % of the confident-tag spectra place on the
  given database, or a panel organism alone places at least 3 % of them and at least a tenth of the given share.
  Fewer than 100 confident-tag spectra: not assessed. `strict` infers the organism after a warning (needs
  `--panel`); `off` skips the check.
- **Development calibration** (the 218 release runs, in-sample; see the 0.1.11 report): correct FASTA 158 of 183
  pass; warnings on 120 of 166 related, 140 of 167 unrelated and 12 of 16 incomplete FASTAs, and on 2 of the 3 known
  wrong SDRF organisms. Several warnings on "correct" FASTAs look real (pig and zebrafish runs that place more on the
  human proteome; human runs with a third of their tags on mouse). Cost: one FASTag pass, about 2.4 s per DDA run.

### Post-search check (`prose-auto check`)

```bash
python -m prose_auto.cli check run.mzML run.idXML --fasta db.fasta --fastag FASTag --contaminants crap.fasta \
    [--report out/run.auto.json] [--fdr 0.01] -o out/
```

- **Signal:** the share of the confident-tag spectra the search identified (q ≤ 0.01; raw ID rates vary too much
  with the sample type), and for the unidentified ones the share whose tags place on the given database against
  the panel.
- **Verdicts:** `ok` when at least 15 % of the confident-tag spectra are identified; otherwise `settings likely
  wrong` when at least 20 % of the unidentified ones place on the given database (check label, protease,
  modifications, windows), else `database likely wrong or incomplete`, naming a panel organism that places at least
  10 % of them. ProSE results without q-values (no rescoring below 100 PSMs) identify nothing.
- **Development calibration** (in-sample): correct FASTA 159 of 183 ok; deliberately wrong settings (label dropped
  or added, wrong protease) flagged 49 of 50, 35 of them as settings; unrelated FASTAs flagged 163 of 165, related
  70 of 165, incomplete 2 of 16. The pre-search check is better for related and incomplete FASTAs; only the
  post-search check sees wrong settings.

The order of the decisions, and which detectors feed each parameter, are in
[docs/decision_flow.md](docs/decision_flow.md).

## Organism inference (`--infer-organism`) and validation

With `--infer-organism`, `--organism-source auto` (the default there):
- **An organism is stated and `--panel` is given:** the spectra check the stated organism. The panel mini-search, FASTag and tier 2 run with the stated organisms as evidence, and the data overrule the record only from 5 % of the mini-search spectra.
- **With `--fasta` or without a panel:** the stated organism is used directly.
- **Nothing stated:** one ProSE search against the signature panel (`--organism-panel`), else the panel mini-search.

On 180 public studies (62 validation, 118 held-out, one run each; truth from the curated SDRFs) the default reaches:

| Parameter | Correct |
|---|---|
| organism | 177/180 |
| instrument | 163/180 |
| acquisition | 168/180 |
| label | 160/178 |
| enrichment | 165/180 |
| enzyme | 174/175 |
| Cys alkylation | 122/128 |

The mass tolerances are measured, not copied. A wide-window search against the true proteome agrees with the measured error in 117 of 120 cases. The curated SDRF value is the authors' search-engine window, not the instrument's accuracy, so it is not used as truth.

Every report also carries `agreement` (per setting: the inconsistency of the independent indicators, the
dissenting detectors and Jev's verdict; overall score and level) and `organism_mix` (suspected, score,
indicators, scope); a suspected mixture is searched against the union of the proteomes (`organisms`,
`search_database`). Every organism the user, the record, Jev (>= 0.3) or the file name (`HYE`, `05Y_50H_45E`,
`Ecoli`) names, and the panel mini-search's own winner, is mini-searched against its proteome (panel, or
UniProt by name or NCBI taxon id, cached under `<panel>/extra`). Greedily, the candidate with the most spectra
the union so far does not explain joins while it has at least max(10, 0.5 % of the spectra) of them (1 when the
file name names it; max(40, 2 %) when only the texts or a record beyond its first three organisms name it), up
to 3 organisms. An inconsistency of 0.2 or more flags 65 % of the wrong decisions (precision 0.25).

Other sets (0.1.3): non-tryptic Orbitrap DDA enzyme 15/16, DIA/timsTOF non-tryptic 8/11 with 12/12 tryptic
controls, HLA ligand files 14/14 (unspecific, human); organism mixtures (38 runs) sensitivity 19/19, specificity
18/19.

## Known limits (prototype)

* The mini-search is pure Python; with persistent peptide indexes it takes 10-30 s per run. Since 0.1.11 it is a
  fallback in the default mode (Cys when FASTag and the texts disagree or both abstain, timsTOF protease); with
  `--infer-organism` it decides the organism and the Cys alkylation as in 0.1.10.
* Vendor input: FASTag 1.15.3 reads Thermo `.raw` and Bruker `.d` directly, but prose-auto's own detectors (MS1
  SILAC pairs and scan-to-scan precision, MS3 reporter ions, the mini-search) read the spectra through pyopenms, and
  FASTag passes on only the tagged MS2 spectra (`-out_spectra`). So `.raw` still needs an mzML; ddaPASEF `.d` is still
  converted by FileConverter.
* Repeated-ion and MS1 scan-to-scan estimates give precision without the offset and are about 2x too small
  alone; the FASTag placement errors or the immonium/y1 probe should decide. On the 120 mass-truth cases the
  placement errors are within 1.5x of a wide-window reference in 58/58 (precursor), 41/41 (high-res fragment) and
  19/20 (ion-trap fragment) cases; the mini-search PSM errors were in 55/56, 39/39 and 18/19.
* Proteases on timsTOF: the MS2 scans start at m/z 200-235, above y1 of K and R, so FASTag's tag-profile call
  fails and an R end reads as R+10 (heavy Arg). The total K/R share of the tag ends (13 tryptic runs 0.33-0.75;
  Glu-C, Asp-N, chymotrypsin 0.006-0.025) and the residue shares decide instead; with fewer than 175 tag ends,
  the mini-search protease competition (17/17 timsTOF runs). There are no HLA timsTOF runs in the validation sets.
  N-side cutters (Lys-N, Asp-N) need 100 N-terminal tags; LysargiNase has no model and reads as "unspecific".
* The non-tryptic sets are small (16 Orbitrap DDA, 11 DIA/timsTOF runs).
* `organism_mix` from the records and texts alone is study-level: a pure file of a mixed study is flagged unless
  the stated second species was searched and is absent.
  A flagged run searches a union only when a further species has spectra of its own in the mini-search; a
  species at a few percent with no such spectra (PXD020394, SARS-CoV-2) is missed.
* Unions stop at 3 organisms (0.1.6-0.1.8: 6). The fourth and later members of PXD017035 (8 species) are no longer
  searched; `--fast` with every `--organism` searches them all.

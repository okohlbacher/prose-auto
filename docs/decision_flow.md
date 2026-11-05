# Decision flow (prose-auto 0.1.9)

prose-auto has two modes:
- **Full mode (the default)** infers the organism and the search database from the spectra and the texts.
- **Fast mode (`--fast` with `--organism` and/or `--fasta`)** takes the given organisms as true. It runs every other
  detector against their proteome, or against the union of their proteomes when several are given.

Both modes decide the parameters in the same order:
1. **What the file is:** top-down check, then instrument, analyzers, fragmentation, acquisition and charge range.
2. **What needs no database:** mass accuracy, chemistry and the tag termini.
3. **What is stated:** repository records, the user's description, methods sections, the file name.
4. **The organism and the search database.** Only this step differs between the modes.
5. **What needs the database:** the mini-search (mass accuracy, Cys alkylation) and the FASTag mass gaps.
6. **The readers and reconciling:** Jev reads last, then reconcile decides each parameter.

In full mode, a suspected organism mixture adds two more steps after reconciling: the union search, then the
database-dependent steps run again against the union.

## Order of the stages

Grey dashed boxes run in full mode only. Fast mode replaces them with the "organism given" box.

```mermaid
flowchart TD
    A[mzML / .d input] --> TD{top-down,<br/>native, intact?}
    TD -- yes --> STOP[abort: no ProSE settings]
    TD -- no --> M[1 metadata: instrument, MS1/MS2 analyzer,<br/>fragmentation, acquisition, charge range<br/>DIA: DIAspeXtract pseudo-spectra]
    M --> TOL[2a mass accuracy without a database:<br/>repeated ions, lock mass, immonium/y1]
    TOL --> CH[2b chemistry: reporter ions, SILAC pairs,<br/>diagnostic ions, digest hint, FASTag termini]
    CH --> TX[3 texts: repository record, user description,<br/>methods, file name, LLM reader]
    TX --> MODE{--fast?}
    MODE -- no --> SRC[4 organism source: auto / data / prose / metadata]
    SRC --> PO[ProSE organism search<br/>signature panel, HLA-first]:::full
    SRC --> PM[panel mini-search<br/>organism competition]:::full
    PO --> FS[FASTag species call]:::full
    PM --> FS
    FS --> T2[tier 2: shortlist, UniProt proteomes,<br/>second mini-search]:::full
    MODE -- yes --> GIV[4 organism given:<br/>proteome per --organism, + --fasta,<br/>union of several]:::fast
    T2 --> DB[(search database)]
    GIV --> DB
    DB --> MS[5a mini-search on the database:<br/>PSM mass errors, Cys alkylation competition]
    MS --> GAP[5b FASTag mass gaps on the database:<br/>Cys alkylation, variable mods, enrichment]
    GAP --> JEV[6a Jev reads last: texts + measured profiles<br/>full: + organism profile and organism question]
    JEV --> REC[6b reconcile: one decision per parameter]
    REC --> MIX{organism mix<br/>suspected?}:::full
    MIX -- yes --> UN[union: candidates mini-searched,<br/>own spectra, UniProt]:::full
    UN --> REDO[Cys competition and FASTag gaps<br/>against the union, reconcile again]:::full
    MIX -- no --> OUT
    REDO --> OUT
    REC -- fast --> OUT[ProSE parameters, SDRF columns,<br/>Jev cross-check, indicator agreement]
    classDef full fill:#eeeeee,stroke:#888888,stroke-dasharray: 5 5
    classDef fast fill:#ffffff,stroke:#000000,stroke-width:2px
```

## Which detectors feed which decision

An edge means that the detector gives evidence for that parameter. Grey dashed detectors run in full mode only.
Detectors marked "db" need the search database, so in fast mode they run against the given proteome.

```mermaid
flowchart LR
    subgraph stated[stated]
        USR[user: --organism, --fasta,<br/>description]
        REP[repository record / SDRF]
        TXT[texts, file name]
        JEV[Jev, last]
    end
    subgraph data[spectra, no database]
        META[mzML metadata]
        PM2[repeated ions, lock mass]
        DIAG[reporter / diagnostic ions,<br/>SILAC pairs]
        TERM[FASTag termini]
    end
    subgraph db[spectra against the database]
        MINI[mini-search, db]
        GAPS[FASTag mass gaps, db]
    end
    subgraph orgdet[organism detectors]
        PROSE[ProSE organism search]:::full
        PANEL[panel mini-search competition]:::full
        FSP[FASTag species]:::full
        TIER2[tier 2]:::full
        MIXU[mixture test and union]:::full
    end
    ORG[organism / database]
    TOLD[tolerances, offsets]
    CHG[charge range]
    ENZ[enzyme]
    CYS[Cys alkylation]
    LAB[label]
    ENR[enrichment]
    VAR[variable mods]
    INS[instrument, acquisition,<br/>fragmentation, analyzers]

    USR --> ORG
    REP --> ORG
    TXT --> ORG
    JEV -. full only .-> ORG
    PROSE --> ORG
    PANEL --> ORG
    FSP --> ORG
    TIER2 --> ORG
    MIXU --> ORG
    ORG --> MINI
    ORG --> GAPS

    META --> INS
    META --> CHG
    REP --> INS
    TXT --> INS
    JEV --> INS

    PM2 --> TOLD
    MINI --> TOLD
    JEV --> TOLD
    TXT --> TOLD

    REP --> ENZ
    USR --> ENZ
    TERM --> ENZ
    JEV --> ENZ
    TXT --> ENZ
    PROSE -. HLA-first .-> ENZ

    DIAG --> CYS
    MINI --> CYS
    GAPS --> CYS
    TXT --> CYS
    JEV --> CYS

    DIAG --> LAB
    TERM --> LAB
    TXT --> LAB
    JEV --> LAB

    DIAG --> ENR
    GAPS --> ENR
    TXT --> ENR
    JEV --> ENR

    GAPS --> VAR
    TXT --> VAR
    JEV --> VAR
    ENR --> VAR

    classDef full fill:#eeeeee,stroke:#888888,stroke-dasharray: 5 5
```

## Per parameter

| Parameter | Evidence, in decision order | Fast mode |
|---|---|---|
| organism, search database | user, record, texts, then the organism source: the ProSE organism search when nothing is stated, else the panel mini-search with the FASTag species call and tier 2; Jev's organism answer; after reconciling, the mixture test and the union (3 organisms at most) | given; no detector runs; several organisms are one union |
| instrument, analyzers, fragmentation, acquisition | mzML metadata (Bruker `analysis.tdf`) first, then the record, the texts and Jev | unchanged |
| charge range | precursor charges up to 9 kDa (DIA: the pseudo-spectra), clamped to 1–20 | unchanged |
| precursor / fragment tolerance and offset | database-free measurement (repeated ions, lock mass, immonium/y1), then the mini-search PSM errors on the database; texts and Jev only as coarse windows | mini-search on the given proteome |
| enzyme | user and record first; then the FASTag tag-termini call (margin ≥ 0.15; on timsTOF the K/R and residue shares of the tag ends and the mini-search protease competition) or Jev (p ≥ 0.8), the more confident one; then keywords. In full mode, the ProSE organism search (HLA-first) also decides when nothing else does | no ProSE organism search |
| Cys alkylation | diagnostic ions (weak), the mini-search competition (CAM / none / others) on the database, the FASTag mass gaps on the database, texts, Jev | competition and gaps on the given proteome |
| label | reporter ions, SILAC pairs, FASTag N-terminal gaps, texts, Jev | unchanged |
| enrichment | diagnostic ions, FASTag mass-gap combinations on the database, texts, Jev | gaps on the given proteome |
| variable mods | FASTag mass-gap combinations on the database, texts, Jev, implied by the enrichment | gaps on the given proteome |

The mini-search decides the protease only on timsTOF, where the tag-profile call fails (MS2 scans start above y1 of K
and R). Elsewhere its protease competition lost to the tag termini on tryptic and HLA runs and is not called.

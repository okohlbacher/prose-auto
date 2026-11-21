# Held-out score (PREREG 0.1.12)

Runs: {'drawn': 19, 'bottom-up': 12, 'top-down': 0, 'no-proteome': 0, 'organism not verified': 5, 'missing': 2}

**Verdict: pilot: no verdict**

| Endpoint | n | estimate | bound | pass |
|---|---|---|---|---|
| P1 median PSM ratio | 12 | 1.098 | 1.051 | None |
| P2 median PSM ratio | 12 | 1.250 | 1.137 | None |
| P3 relative (auto - reference) | 12 | 0.083 | -0.167 | None |
| P3 absolute (auto share) | 12 | 0.917 | 0.661 | None |

| Check:kind | status | assessed | share | Wilson lower | bound | Holm pass |
|---|---|---|---|---|---|---|
| post:correct | gate | 0 | n/a | n/a | 0.750 | False |
| post:related_distant | gate | 0 | n/a | n/a | 0.200 | False |
| post:unrelated | gate | 0 | n/a | n/a | 0.750 | False |
| post:wrong_label | gate | 0 | n/a | n/a | 0.650 | False |
| post:wrong_protease | gate | 0 | n/a | n/a | 0.650 | False |
| pre:correct | gate | 12 | 0.750 | 0.513 | 0.700 | False |
| pre:related_close | gate | 0 | n/a | n/a | 0.450 | False |
| pre:related_distant | gate | 0 | n/a | n/a | 0.450 | False |
| pre:unrelated | gate | 0 | n/a | n/a | 0.500 | False |

| §9 rule | pass |
|---|---|
| cys: -alkylation_out (auto_never) vs competition (auto_always) | None |
| mini-search fallback (auto) vs always (auto_always) | None |
| label: S1 >= always-label-free baseline, no isobaric call without reporter ions | True |
| protease: S2 >= stated-protease (SDRF arm) baseline | False |
| FASTag QC: analyser class agrees with the header in >= 98 % | True |

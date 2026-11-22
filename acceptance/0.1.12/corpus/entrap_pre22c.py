"""Entrapment FDP for ProSE idXML (doubled search space, as eeval.py of the PR's earlier evaluation).

  entrap.py makedb SRC_TD_FASTA OUT_FASTA SEED   targets + shuffled twin per target + all decoys (eeval.makedb)
  entrap.py eval OUT_JSON DB_FASTA IDXML ...     counts of one arm x shuffle (all files pooled)
  entrap.py selfcheck

Entrapment twin: residues inside each Trypsin/P segment shuffled, C-terminal K/R and the first residue fixed, seeded per
protein (7919 i + 17 + 1000003 SEED, as eeval.py). A target PSM is entrapment if its peptide (I = L, unmodified) is an
entrapment peptide and not a target peptide. Wen et al. (2025): combined = N_E (1 + 1/r) / (N_T + N_E),
lower = N_E / (N_T + N_E), r = entrapment / target peptides of the searchable space (Trypsin/P, <= 2 missed cleavages,
length 7-40, Met-clipped variants). Protein N-terminal stratum: peptides starting at residue 0 or 1 (after M) of a protein.
"""
import json, random, re, sys


def read_fasta(path):
    out, hdr, seq = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(seq)))
            hdr, seq = line[1:].rstrip("\n"), []
        else:
            seq.append(line.strip())
    out.append((hdr, "".join(seq)))
    return out


def shuffled(seq, seed):  # = eeval.shuffled
    rng = random.Random(seed)
    res, start = list(seq), 1
    for i in range(1, len(seq) + 1):
        if i == len(seq) or seq[i - 1] in "KR":
            end = i - 1 if seq[i - 1] in "KR" else i
            mid = res[start:end]
            rng.shuffle(mid)
            res[start:end] = mid
            start = i
    return "".join(res)


def makedb(src, out, seed):  # = eeval.makedb
    entries = read_fasta(src)
    targets = [(h, s) for h, s in entries if not h.startswith("DECOY_")]
    decoys = [(h, s) for h, s in entries if h.startswith("DECOY_")]
    entrap = [(f"ENTRAP_{h.split()[0]} entrapment of {h.split()[0]}", shuffled(s, 7919 * i + 17 + 1000003 * seed))
              for i, (h, s) in enumerate(targets)]
    with open(out, "w") as f:
        for h, s in targets + entrap + decoys + [(f"DECOY_{h}", s[::-1]) for h, s in entrap]:
            f.write(f">{h}\n")
            for j in range(0, len(s), 60):
                f.write(s[j:j + 60] + "\n")


def peptides(seqs):
    """(all, protein N-terminal) Trypsin/P peptides, <= 2 missed cleavages, length 7-40, Met-clipped variants; I = L"""
    allp, nterm = set(), set()
    for s in seqs:
        s = s.replace("I", "L")
        sites = [0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR"] + [len(s)]
        for a in range(len(sites) - 1):
            for b in range(a + 1, min(a + 4, len(sites))):
                p = s[sites[a]:sites[b]]
                if 7 <= len(p) <= 40:
                    allp.add(p)
                    if a == 0:
                        nterm.add(p)
                if a == 0 and s.startswith("M") and 7 <= len(p) - 1 <= 40:
                    allp.add(p[1:]); nterm.add(p[1:])
    return allp, nterm


def strip(seq):
    return re.sub(r"[^A-Z]", "", re.sub(r"\([^)]*\)|\[[^\]]*\]", "", seq)).replace("I", "L")


PID = re.compile(r"<PeptideIdentification\b.*?</PeptideIdentification>", re.S)
HIT = re.compile(r'<PeptideHit score="([^"]+)" sequence="([^"]+)"(.*?)</PeptideHit>', re.S)
UP = lambda body, k: (re.search(rf'name="{k}" value="([^"]*)"', body) or [None, None])[1]


def psms(path):
    """rank-1 hit per spectrum: (q, percolator_score, canonical peptide, is_decoy)"""
    txt = open(path).read()
    for m in PID.finditer(txt):
        h = HIT.search(m.group(0))
        if h:
            q, seq, body = h.groups()
            yield float(q), float(UP(body, "percolator_score") or "nan"), strip(seq), "decoy" in (UP(body, "target_decoy") or "")


def fdp(nt, ne, r):
    n = nt + ne
    return (100 * ne * (1 + 1 / r) / n if n else 0.0), (100 * ne / n if n else 0.0)


def equal_fdp(rows, r, target=1.0):
    """rows: (q, cls) of target-side PSMs; N_T at the largest q threshold whose monotonised combined FDP <= target"""
    rows = sorted(rows)
    nt = ne = 0; curve = []
    for i, (q, c) in enumerate(rows):
        nt += c == "T"; ne += c == "E"
        if i + 1 == len(rows) or rows[i + 1][0] != q:
            curve.append((nt, fdp(nt, ne, r)[0]))
    # q_emp(t) = min over t' >= t of FDP(t'); the largest t with q_emp(t) <= target
    run_min, best = float("inf"), 0
    qemp = []
    for nt_, f in reversed(curve):
        run_min = min(run_min, f); qemp.append((nt_, run_min))
    for nt_, qe in qemp:
        if qe <= target:
            best = max(best, nt_)
    return best


def evaluate(out, db, files):
    entries = read_fasta(db)
    t, tn = peptides(s for h, s in entries if not h.startswith(("DECOY_", "ENTRAP_")))
    e, en = peptides(s for h, s in entries if h.startswith("ENTRAP_"))
    e -= t; en -= t
    r, rn = len(e) / len(t), len(en) / len(tn)
    cls = lambda p: "T" if p in t else ("E" if p in e else "U")
    rows, pep_best = [], {}
    for f in files:
        for q, sc, p, dec in psms(f):
            if not dec:
                rows.append((q, cls(p), p in tn or p in en))
            key = ("D" if dec else "T", p)
            if sc == sc and sc > pep_best.get(key, -1e300):
                pep_best[key] = sc
    res = dict(db=db, files=len(files), r=r, r_nterm=rn)
    acc = [x for x in rows if x[0] <= 0.01]
    nt, ne = sum(c == "T" for _, c, _ in acc), sum(c == "E" for _, c, _ in acc)
    res["psm"] = dict(N_T=nt, N_E=ne, N_U=sum(c == "U" for _, c, _ in acc), accepted=len(acc),
                      combined=fdp(nt, ne, r)[0], lower=fdp(nt, ne, r)[1],
                      N_T_at_equal_fdp_1pct=equal_fdp([(q, c) for q, c, _ in rows if c != "U"], r))
    acc_n = [x for x in acc if x[2]]
    ntn, nen = sum(c == "T" for _, c, _ in acc_n), sum(c == "E" for _, c, _ in acc_n)
    res["psm_nterm"] = dict(N_T=ntn, N_E=nen, accepted=len(acc_n), share_of_accepted=100 * len(acc_n) / max(1, len(acc)),
                            combined=fdp(ntn, nen, rn)[0], lower=fdp(ntn, nen, rn)[1])
    # peptide level: best percolator score per distinct peptide (target and decoy side), target-decoy q, monotonised
    peps = sorted(((sc, side, p) for (side, p), sc in pep_best.items()), reverse=True)
    nd = ntg = 0; qs = []
    for sc, side, p in peps:
        nd += side == "D"; ntg += side == "T"
        qs.append(nd / max(1, ntg))
    for i in range(len(qs) - 2, -1, -1):
        qs[i] = min(qs[i], qs[i + 1])
    acc_p = [p for (sc, side, p), q in zip(peps, qs) if side == "T" and q <= 0.01]
    ntp, nep = sum(cls(p) == "T" for p in acc_p), sum(cls(p) == "E" for p in acc_p)
    res["peptide"] = dict(N_T=ntp, N_E=nep, accepted=len(acc_p), combined=fdp(ntp, nep, r)[0], lower=fdp(ntp, nep, r)[1])
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "makedb":
        makedb(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    elif cmd == "eval":
        evaluate(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "selfcheck":
        s = "MAGKLLPRSTVWQKEEK"; e = shuffled(s, 3)
        assert sorted(e) == sorted(s) and e[0] == "M" and [i for i, c in enumerate(e) if c in "KR"] == [i for i, c in enumerate(s) if c in "KR"]
        assert strip("PEPT(Oxidation)IDE") == "PEPTLDE" and strip(".(Acetyl)MAGK") == "MAGK"
        a, n = peptides(["MAAAAAAAKBBBBBBBR"])
        assert "MAAAAAAAK" in n and "AAAAAAAK" in n and "BBBBBBBR" in a and "BBBBBBBR" not in n
        rows = [(0.001, "T")] * 99 + [(0.002, "E")] + [(0.02, "E")] * 5 + [(0.03, "T")] * 10
        assert equal_fdp(rows, 1.0, 2.5) == 99, equal_fdp(rows, 1.0, 2.5)
        print("selfcheck ok")

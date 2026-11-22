"""Entrapment FDP for ProSE idXML (doubled search space, as eeval.py of the PR's earlier evaluation).

  entrap.py makedb SRC_TD_FASTA OUT_FASTA SEED   targets + shuffled twin per target + all decoys (eeval.makedb)
  entrap.py eval OUT_JSON DB_FASTA IDXML ...     counts of one arm x shuffle (all files pooled)
  entrap.py selfcheck

Entrapment twin: residues inside each Trypsin/P segment shuffled, C-terminal K/R and the first residue fixed, seeded per
protein (7919 i + 17 + 1000003 SEED, as eeval.py). A target PSM is entrapment when all its source proteins are
entrapment twins (ENTRAP_), target when any is a target. Wen et al. (2025): combined = N_E (1 + 1/r) / (N_T + N_E),
lower = N_E / (N_T + N_E), r = entrapment / target peptides of the searchable space of the search's own enzyme and
missed cleavages (the idXML's SearchParameters; length 7-40, Met-clipped variants; unspecific: every distinct window). Protein
N-terminal stratum: PSMs starting at residue 0, or 1 after M, of a source protein. Amendment 2026-11-22c; until then
both the classification and r used Trypsin/P peptide sets for every search (entrap_pre22c.py).
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


# the cleavage rules of the frozen ProSE build (pr/share/OpenMS/CHEMISTRY/Enzymes.xml, RegEx), by idXML enzyme name
ENZYMES = {"trypsin": r"(?<=[KRX])(?!P)", "trypsin/p": r"(?<=[KRX])", "lys-c": r"(?<=[KX])(?!P)",
           "lys-c/p": r"(?<=[KX])", "lys-n": r"(?=[KX])", "arg-c": r"(?<=[RX])(?!P)", "arg-c/p": r"(?<=[RX])",
           "asp-n": r"(?=[DBX])", "chymotrypsin": r"(?<=[FYWLJX])(?!P)", "chymotrypsin/p": r"(?<=[FYWLJX])",
           "glutamyl endopeptidase": r"(?<=[DBEZX])", "unspecific cleavage": None}
LMIN, LMAX = 7, 40


def peptides(seqs, enzyme="trypsin/p", mc=2):
    """(all, protein N-terminal) peptides of the enzyme, <= mc missed cleavages, length 7-40, Met-clipped variants;
    I = L after digestion (I is no cleavage site, L is one of chymotrypsin's)"""
    rx = re.compile(ENZYMES[enzyme])
    allp, nterm = set(), set()
    for s in seqs:
        sites = [0] + [m.start() for m in rx.finditer(s) if 0 < m.start() < len(s)] + [len(s)]
        L = s.replace("I", "L")
        for a in range(len(sites) - 1):
            for b in range(a + 1, min(a + mc + 2, len(sites))):
                p = L[sites[a]:sites[b]]
                if LMIN <= len(p) <= LMAX:
                    allp.add(p)
                    if a == 0:
                        nterm.add(p)
                if a == 0 and s.startswith("M") and LMIN <= len(p) - 1 <= LMAX:
                    allp.add(p[1:]); nterm.add(p[1:])
    return allp, nterm


def _hashes(seqs, k_range):
    """64-bit hashes of every distinct window (I = L) per length k: (all, protein N-terminal: start 0, or 1 after M)"""
    import numpy as np
    seqs = [x.replace("I", "L") for x in seqs]
    a = np.frombuffer(("\0".join(seqs) + "\0").encode(), np.uint8).astype(np.uint64)
    n = len(a)
    sep = np.concatenate(([0], np.cumsum(a == 0)))
    nt = np.zeros(n, bool)
    o = 0
    for x in seqs:
        nt[o] = True
        nt[o + 1] = x.startswith("M") and len(x) > 1
        o += len(x) + 1
    h = np.zeros(n, np.uint64)
    B = np.uint64(0x9E3779B97F4A7C15)
    with np.errstate(over="ignore"):
        for k in range(1, max(k_range) + 1):
            m = max(0, n - k + 1)
            h = h[:m] * B + a[k - 1:k - 1 + m]
            if k in k_range:
                ok = sep[k:k + m] - sep[:m] == 0
                yield k, np.unique(h[ok]), np.unique(h[ok & nt[:m]])


def windows(targets, entraps):
    """unspecific cleavage: (r, r_nterm) from distinct windows of length 7-40 (64-bit hashes, one length at a time),
    entrapment windows equal to a target window not counted"""
    import numpy as np
    ks = range(LMIN, LMAX + 1)
    t = tn = e = en = 0
    for (k, ht, htn), (_, he, hen) in zip(_hashes(targets, ks), _hashes(entraps, ks)):
        t += len(ht); tn += len(htn)
        e += int((~np.isin(he, ht, assume_unique=True)).sum())
        en += int((~np.isin(hen, ht, assume_unique=True)).sum())
    return e / t, en / tn


def ratio(targets, entraps, enzyme, mc):
    """r and r_nterm of the search's own enzyme: entrapment / target peptides of the searchable space, entrapment
    peptides equal to a target peptide not counted"""
    if ENZYMES[enzyme] is None:
        return windows(targets, entraps)
    t, tn = peptides(targets, enzyme, mc)
    e, en = peptides(entraps, enzyme, mc)
    e -= t; en -= t
    return len(e) / len(t), len(en) / len(tn)


def strip(seq):
    return re.sub(r"[^A-Z]", "", re.sub(r"\([^)]*\)|\[[^\]]*\]", "", seq)).replace("I", "L")


PID = re.compile(r"<PeptideIdentification\b.*?</PeptideIdentification>", re.S)
HIT = re.compile(r'<PeptideHit score="([^"]+)" sequence="([^"]+)"(.*?)</PeptideHit>', re.S)
PROT = re.compile(r'<ProteinHit id="([^"]+)" accession="([^"]+)"')
SP = re.compile(r"<SearchParameters\b[^>]*>")
UP = lambda body, k: (re.search(rf'name="{k}" value="([^"]*)"', body) or [None, None])[1]
ATTR = lambda tag, k: (re.search(rf'\b{k}="([^"]*)"', tag) or [None, None])[1]


def search_enzyme(txt):
    """(enzyme, missed cleavages) of the search, from the idXML's SearchParameters; an unknown enzyme stops"""
    tag = SP.search(txt).group(0)
    enz = (ATTR(tag, "enzyme") or "").lower()
    if enz not in ENZYMES:
        raise SystemExit(f"entrap: unknown enzyme {enz!r} in SearchParameters")
    mc = ATTR(tag, "missed_cleavages")
    if mc is None and ENZYMES[enz] is not None:
        raise SystemExit("entrap: no missed_cleavages in SearchParameters")
    return enz, int(mc or 0)


def psms(path, seqs=None):
    """rank-1 hit per spectrum: (q, percolator_score, canonical peptide, is_decoy, class, protein N-terminal).
    class by the source proteins (Wen et al. 2025): T when any protein is a target, E when all are entrapment
    twins (ENTRAP_), U without proteins. N-terminal: starts at residue 0, or 1 of a protein starting with M."""
    txt = open(path).read()
    prot = dict(PROT.findall(txt))
    seqs = seqs or {}
    for m in PID.finditer(txt):
        h = HIT.search(m.group(0))
        if h:
            q, seq, body = h.groups()
            tag = body.split(">", 1)[0]
            accs = [prot.get(r, "") for r in (ATTR(tag, "protein_refs") or "").split()]
            starts = [int(x) for x in (ATTR(tag, "start") or "").split()]
            cls = ("U" if not accs else "T" if any(not a.startswith(("ENTRAP_", "DECOY_")) for a in accs)
                   else "E" if all(a.startswith("ENTRAP_") for a in accs) else "U")
            nterm = any(s == 0 or (s == 1 and seqs.get(a, "").startswith("M")) for a, s in zip(accs, starts))
            yield (float(q), float(UP(body, "percolator_score") or "nan"), strip(seq),
                   "decoy" in (UP(body, "target_decoy") or ""), cls, nterm)


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
    """amendment 2026-11-22c: PSMs classified by their source proteins and r for the search's own enzyme (the idXML's
    SearchParameters), unspecific included; before, both used Trypsin/P peptide sets whatever the search's enzyme"""
    entries = read_fasta(db)
    seqs = {h.split()[0]: s for h, s in entries}
    enz = {search_enzyme(open(f).read()) for f in files}
    if len(enz) != 1:
        raise SystemExit(f"entrap: pooled files of different searches {sorted(enz)}")
    (enzyme, mc), = enz
    r, rn = ratio([s for h, s in entries if not h.startswith(("DECOY_", "ENTRAP_"))],
                  [s for h, s in entries if h.startswith("ENTRAP_")], enzyme, mc)
    rows, pep_best, pep_cls = [], {}, {}
    for f in files:
        for q, sc, p, dec, c, nt_ in psms(f, seqs):
            if not dec:
                rows.append((q, c, nt_))
                if c != "U" and pep_cls.get(p) != "T":
                    pep_cls[p] = c
            key = ("D" if dec else "T", p)
            if sc == sc and sc > pep_best.get(key, -1e300):
                pep_best[key] = sc
    res = dict(db=db, files=len(files), enzyme=enzyme, missed_cleavages=mc, classification="protein", r=r,
               r_nterm=rn)
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
    ntp, nep = sum(pep_cls.get(p) == "T" for p in acc_p), sum(pep_cls.get(p) == "E" for p in acc_p)
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
        a, _ = peptides(["AAAAAAAKPGGGGGGGRCCCCCCCK"], "lys-c", 0)   # Lys-C: after K, not before P; R no site
        assert a == {"AAAAAAAKPGGGGGGGRCCCCCCCK"}, a
        a, _ = peptides(["AAAAAAIGGGGGGGLCCCCCCCF"], "chymotrypsin", 0)   # I no site, L one; I = L afterwards
        assert a == {"AAAAAALGGGGGGGL", "CCCCCCCF"}, a
        T, En = ["MACDEFGHKLACDEFGHKLW", "ACDEFGHKLMN"], ["MHGFEDCAKLWDEFGHKLAC", "NMLKHGFEDCA"]
        brute = lambda xs: ({x[i:i + k].replace("I", "L") for x in xs for k in range(7, 41) for i in range(len(x) - k + 1)},
                            {x[o:o + k] for x in xs for o in ((0, 1) if x.startswith("M") else (0,)) for k in range(7, 41)
                             if o + k <= len(x)})
        (bt, btn), (be, ben) = brute(T), brute(En)
        assert windows(T, En) == (len(be - bt) / len(bt), len(ben - bt) / len(btn)), (windows(T, En), len(bt), len(be))
        assert windows(T, T) == (0.0, 0.0)
        rows = [(0.001, "T")] * 99 + [(0.002, "E")] + [(0.02, "E")] * 5 + [(0.03, "T")] * 10
        assert equal_fdp(rows, 1.0, 2.5) == 99, equal_fdp(rows, 1.0, 2.5)
        print("selfcheck ok")
